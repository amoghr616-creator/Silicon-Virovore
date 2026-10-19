"""Run an isolated Vina-only 1YCR/p53 positive-control redocking benchmark.

This is not the HybriDock/RAPiDock pipeline and its Vina affinity is not an
experimental binding free energy. The output directory must be new or empty.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser, is_aa
from Bio.SeqUtils import seq1
from rdkit import Chem


ROOT = Path(__file__).resolve().parents[1]
BENCHMARK_ROOT = ROOT / "prototype2" / "benchmark_1ycr"
REFERENCE_PDB = BENCHMARK_ROOT / "1YCR.pdb"
REFERENCE_SEQUENCE = "SQETFSDLWKLLPEN"
OBSERVED_PEPTIDE_SEQUENCE = "ETFSDLWKLLPEN"
VINA_DEFAULT = Path(
    "/opt/homebrew/Caskroom/miniforge/base/envs/sv-docking/bin/vina"
)
MEEKO_PREP_DEFAULT = Path(
    "/Library/Frameworks/Python.framework/Versions/3.14/bin/mk_prepare_receptor.py"
)
CONFIG = {
    "backend": "AutoDock Vina 1.2.7 only",
    "backend_comparison_scope": "separate from HybriDock/RAPiDock; not score-equivalent",
    "scoring_function": "vina",
    "vina_score_direction": "more negative is more favorable; kcal/mol score, not measured affinity",
    "exhaustiveness": 64,
    "num_modes": 20,
    "energy_range_kcal_mol": 5.0,
    "cpu": 4,
    "seed": 20261016,
    "box_center_source": "centroid of the experimentally observed peptide chain B heavy atoms",
    "box_size_A": [30.0, 30.0, 30.0],
    "box_margin_note": "native peptide coordinate span is 21.03 x 17.51 x 15.49 A; 30 A cubic box leaves at least 4.48 A per side along the largest span",
    "native_contact_cutoff_A": 5.0,
    "interface_residue_cutoff_A": 10.0,
    "severe_clash_cutoff_A": 2.0,
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _parse_pdb_sequences(path: Path) -> tuple[str, str, list[dict], list[str], list[str]]:
    text = path.read_text()
    seqres_names: list[str] = []
    for line in text.splitlines():
        if line.startswith("SEQRES") and line[11:12] == "B":
            seqres_names.extend(line[19:70].split())
    deposited = "".join(seq1(name) for name in seqres_names)

    structure = PDBParser(QUIET=True).get_structure("1ycr", str(path))
    model = next(structure.get_models())
    observed_residues = [
        residue
        for residue in model["B"]
        if is_aa(residue, standard=True)
    ]
    observed = "".join(seq1(residue.resname) for residue in observed_residues)

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -2.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(deposited, observed)[0]

    mapping: list[dict] = []
    for deposited_slice, observed_slice in zip(*alignment.aligned):
        if deposited_slice[1] - deposited_slice[0] != observed_slice[1] - observed_slice[0]:
            raise ValueError("alignment block does not provide one-to-one residue mapping")
        for deposited_index, observed_index in zip(
            range(deposited_slice[0], deposited_slice[1]),
            range(observed_slice[0], observed_slice[1]),
        ):
            residue = observed_residues[observed_index]
            mapping.append({
                "observed_sequence_index_1based": observed_index + 1,
                "deposited_sequence_index_1based": deposited_index + 1,
                "deposited_residue": deposited[deposited_index],
                "observed_residue": observed[observed_index],
                "author_residue_number": int(residue.id[1]),
                "insertion_code": str(residue.id[2]).strip(),
            })

    if deposited != REFERENCE_SEQUENCE:
        raise ValueError(
            f"unexpected deposited chain B sequence: {deposited!r}"
        )
    if observed != OBSERVED_PEPTIDE_SEQUENCE:
        raise ValueError(
            f"unexpected observed chain B sequence: {observed!r}"
        )
    if len(mapping) != len(observed_residues):
        raise ValueError("not every observed peptide residue aligned")

    receptor_lines = [
        line for line in text.splitlines()
        if line.startswith("ATOM  ") and line[21:22] == "A"
    ]
    ligand_lines = [
        line for line in text.splitlines()
        if line.startswith("ATOM  ")
        and line[21:22] == "B"
        and line[16:17] in (" ", "A")
    ]
    if not receptor_lines or not ligand_lines:
        raise ValueError("1YCR must provide coordinate records for receptor A and peptide B")

    chains = sorted({line[21:22] for line in receptor_lines})
    if chains != ["A"]:
        raise ValueError(f"clean receptor extraction retained unexpected chains: {chains}")
    return deposited, observed, mapping, receptor_lines, ligand_lines


def _parse_index_map(pdbqt_text: str) -> dict[int, int]:
    """Return source 1-based RDKit atom index -> output PDBQT serial."""
    mapping: dict[int, int] = {}
    for line in pdbqt_text.splitlines():
        if not line.startswith("REMARK INDEX MAP"):
            continue
        values = [int(value) for value in line.split()[3:]]
        if len(values) % 2:
            raise ValueError("malformed Meeko REMARK INDEX MAP line")
        for input_index, output_serial in zip(values[::2], values[1::2]):
            if input_index in mapping:
                raise ValueError(f"duplicate source atom index in map: {input_index}")
            mapping[input_index] = output_serial
    if not mapping:
        raise ValueError("Meeko did not emit an explicit atom index map")
    return mapping


def _parse_pdbqt_atoms(text: str, *, first_model_only: bool = True) -> dict[int, dict]:
    atoms: dict[int, dict] = {}
    model_seen = False
    model_number = None
    for line in text.splitlines():
        if line.startswith("MODEL"):
            model_seen = True
            model_number = int(line[10:14].strip() or "1")
            continue
        if line.startswith("ENDMDL"):
            if first_model_only and model_seen:
                break
            continue
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if first_model_only and model_seen and model_number != 1:
            continue
        serial = int(line[6:11])
        if serial in atoms:
            raise ValueError(f"duplicate PDBQT atom serial: {serial}")
        atoms[serial] = {
            "line": line,
            "element": line[76:78].strip().upper(),
            "xyz": np.array([
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            ], dtype=float),
        }
    return atoms


def _prepare_ligand(ligand_lines: list[str], output_path: Path) -> tuple[dict[int, int], dict[int, dict]]:
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    block = "\n".join(ligand_lines + ["TER", "END"]) + "\n"
    molecule = Chem.MolFromPDBBlock(
        block,
        sanitize=True,
        removeHs=False,
        proximityBonding=True,
    )
    if molecule is None:
        raise ValueError("RDKit failed to construct the observed peptide graph")
    if molecule.GetNumAtoms() < len(ligand_lines):
        raise ValueError("RDKit omitted input peptide atoms")
    for atom_index, line in enumerate(ligand_lines):
        if molecule.GetAtomWithIdx(atom_index).GetSymbol().upper() != line[76:78].strip().upper():
            raise ValueError(f"RDKit atom order disagrees with source PDB at atom {atom_index + 1}")

    molecule = Chem.AddHs(molecule, addCoords=True)
    setups = MoleculePreparation().prepare(molecule)
    if len(setups) != 1:
        raise ValueError(f"expected one Meeko ligand setup, received {len(setups)}")
    pdbqt_text, success, error = PDBQTWriterLegacy.write_string(
        setups[0],
        add_index_map=True,
    )
    if not success:
        raise ValueError(f"Meeko ligand preparation failed: {error}")
    output_path.write_text(pdbqt_text)
    input_to_serial = _parse_index_map(pdbqt_text)
    prepared_atoms = _parse_pdbqt_atoms(pdbqt_text, first_model_only=False)
    required_input_indices = set(range(1, len(ligand_lines) + 1))
    if not required_input_indices.issubset(input_to_serial):
        missing = sorted(required_input_indices - set(input_to_serial))
        raise ValueError(f"Meeko atom map misses observed peptide heavy atoms: {missing[:10]}")
    for input_index, source_line in enumerate(ligand_lines, start=1):
        serial = input_to_serial[input_index]
        if serial not in prepared_atoms:
            raise ValueError(f"mapped ligand PDBQT serial is absent: {serial}")
        original_xyz = np.array([
            float(source_line[30:38]),
            float(source_line[38:46]),
            float(source_line[46:54]),
        ])
        if not np.allclose(prepared_atoms[serial]["xyz"], original_xyz, atol=0.01):
            raise ValueError(
                f"prepared coordinates do not map to source atom {input_index}; "
                "refusing to score unverified atom correspondence"
            )
    return input_to_serial, prepared_atoms


def _format_atom_line(source_line: str, serial: int, xyz: np.ndarray) -> str:
    return (
        source_line[:6]
        + f"{serial:5d}"
        + source_line[11:30]
        + f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
        + source_line[54:]
    )


def _assemble_pose(
    receptor_lines: list[str],
    ligand_lines: list[str],
    input_to_serial: dict[int, int],
    docked_atoms: dict[int, dict],
) -> str:
    missing = sorted(set(input_to_serial.values()) - set(docked_atoms))
    if missing:
        raise ValueError(f"docked model is missing mapped ligand atom serials: {missing[:10]}")
    assembled = []
    next_serial = 1
    for line in receptor_lines:
        assembled.append(_format_atom_line(
            line,
            next_serial,
            np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
        ))
        next_serial += 1
    assembled.append("TER")
    for input_index, line in enumerate(ligand_lines, start=1):
        pose_atom = docked_atoms[input_to_serial[input_index]]
        assembled.append(_format_atom_line(line, next_serial, pose_atom["xyz"]))
        next_serial += 1
    assembled.extend(["TER", "END"])
    return "\n".join(assembled) + "\n"


def _heavy_atoms(residue) -> list:
    return [atom for atom in residue if atom.element.upper() not in {"H", "D"}]


def _minimum_distance(left: list, right: list) -> float:
    left_xyz = np.array([atom.coord for atom in left], dtype=float)
    right_xyz = np.array([atom.coord for atom in right], dtype=float)
    distances = np.linalg.norm(left_xyz[:, None, :] - right_xyz[None, :, :], axis=2)
    return float(distances.min())


def _contacts(receptor_chain, peptide_chain, cutoff: float) -> set[tuple[int, int]]:
    contact_pairs = set()
    for receptor_residue in receptor_chain:
        receptor_atoms = _heavy_atoms(receptor_residue)
        for peptide_index, peptide_residue in enumerate(peptide_chain, start=1):
            peptide_atoms = _heavy_atoms(peptide_residue)
            if receptor_atoms and peptide_atoms and _minimum_distance(receptor_atoms, peptide_atoms) <= cutoff:
                contact_pairs.add((int(receptor_residue.id[1]), peptide_index))
    return contact_pairs


def _kabsch_transform(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mobile_center = mobile.mean(axis=0)
    target_center = target.mean(axis=0)
    mobile_centered = mobile - mobile_center
    target_centered = target - target_center
    covariance = mobile_centered.T @ target_centered
    left, _singular, right_t = np.linalg.svd(covariance)
    rotation = left @ np.diag([1.0, 1.0, np.linalg.det(left @ right_t)]) @ right_t
    translation = target_center - mobile_center @ rotation
    return rotation, translation


def _analyze_pose(reference_pdb: Path, pose_pdb: Path) -> dict:
    parser = PDBParser(QUIET=True)
    reference = next(parser.get_structure("reference", str(reference_pdb)).get_models())
    pose = next(parser.get_structure("pose", str(pose_pdb)).get_models())
    receptor_ref = reference["A"]
    peptide_ref = reference["B"]
    receptor_pose = pose["A"]
    peptide_pose = pose["B"]

    reference_sequence = "".join(seq1(residue.resname) for residue in peptide_ref if is_aa(residue, standard=True))
    pose_sequence = "".join(seq1(residue.resname) for residue in peptide_pose if is_aa(residue, standard=True))
    if reference_sequence != OBSERVED_PEPTIDE_SEQUENCE or pose_sequence != reference_sequence:
        raise ValueError(
            f"peptide sequence mismatch: ref={reference_sequence}, pose={pose_sequence}"
        )

    ref_residues = [residue for residue in peptide_ref if is_aa(residue, standard=True)]
    pose_residues = [residue for residue in peptide_pose if is_aa(residue, standard=True)]
    backbone_names = ("N", "CA", "C", "O")
    missing_backbone = []
    for index, (ref_residue, pose_residue) in enumerate(zip(ref_residues, pose_residues), start=1):
        missing = [name for name in backbone_names if name not in pose_residue]
        if missing:
            missing_backbone.append({"peptide_index": index, "residue": pose_residue.resname, "atoms": missing})
        if ref_residue.resname != pose_residue.resname:
            raise ValueError(f"residue mismatch at peptide index {index}")
    if missing_backbone:
        raise ValueError(f"pose missing peptide backbone atoms: {missing_backbone}")

    ref_ca = np.array([residue["CA"].coord for residue in ref_residues], dtype=float)
    pose_ca = np.array([residue["CA"].coord for residue in pose_residues], dtype=float)
    direct_ca_rmsd = float(np.sqrt(np.mean(np.sum((pose_ca - ref_ca) ** 2, axis=1))))

    native_contacts = _contacts(receptor_ref, peptide_ref, CONFIG["native_contact_cutoff_A"])
    predicted_contacts = _contacts(receptor_pose, peptide_pose, CONFIG["native_contact_cutoff_A"])
    recovered_contacts = native_contacts & predicted_contacts

    native_interface_receptor = []
    for receptor_residue in receptor_ref:
        nearest = min(
            _minimum_distance(_heavy_atoms(receptor_residue), _heavy_atoms(peptide_residue))
            for peptide_residue in peptide_ref
        )
        if nearest <= CONFIG["interface_residue_cutoff_A"] and "CA" in receptor_residue:
            native_interface_receptor.append(receptor_residue)
    if not native_interface_receptor:
        raise ValueError("native interface has no receptor C-alpha atoms")
    receptor_ref_ca = np.array([residue["CA"].coord for residue in native_interface_receptor])
    receptor_pose_ca = np.array([receptor_pose[residue.id]["CA"].coord for residue in native_interface_receptor])
    rotation, translation = _kabsch_transform(receptor_pose_ca, receptor_ref_ca)
    native_interface_peptide_indices = sorted({peptide_index for _receptor_index, peptide_index in native_contacts})
    interface_ref = np.array([
        residue[name].coord
        for peptide_index in native_interface_peptide_indices
        for residue in [ref_residues[peptide_index - 1]]
        for name in backbone_names
    ])
    interface_pose = np.array([
        residue[name].coord
        for peptide_index in native_interface_peptide_indices
        for residue in [pose_residues[peptide_index - 1]]
        for name in backbone_names
    ])
    interface_pose_aligned = interface_pose @ rotation + translation
    interface_rmsd = float(np.sqrt(np.mean(np.sum((interface_pose_aligned - interface_ref) ** 2, axis=1))))

    receptor_atoms = [atom for residue in receptor_pose for atom in _heavy_atoms(residue)]
    peptide_atoms = [atom for residue in peptide_pose for atom in _heavy_atoms(residue)]
    receptor_xyz = np.array([atom.coord for atom in receptor_atoms])
    peptide_xyz = np.array([atom.coord for atom in peptide_atoms])
    distances = np.linalg.norm(receptor_xyz[:, None, :] - peptide_xyz[None, :, :], axis=2)
    minimum_distance = float(distances.min())
    severe_clash_count = int(np.sum(distances < CONFIG["severe_clash_cutoff_A"]))

    arginine_metrics = {}
    for residue_number in (203, 205):
        residue = receptor_pose[residue_number]
        residue_atoms = _heavy_atoms(residue)
        residue_xyz = np.array([atom.coord for atom in residue_atoms])
        residue_distances = np.linalg.norm(residue_xyz[:, None, :] - peptide_xyz[None, :, :], axis=2)
        arginine_metrics[str(residue_number)] = {
            "minimum_heavy_atom_distance_A": float(residue_distances.min()),
            "peptide_heavy_atom_counts_within_4A": int(np.sum(residue_distances <= 4.0)),
            "peptide_heavy_atom_counts_within_5A": int(np.sum(residue_distances <= 5.0)),
        }

    geometry_valid = (
        pose_sequence == OBSERVED_PEPTIDE_SEQUENCE
        and not missing_backbone
        and np.isfinite(pose_xyz).all()
        and severe_clash_count == 0
    )
    return {
        "chain_mapping": {
            "receptor_chain": "A",
            "peptide_chain": "B",
            "reference_peptide_sequence": reference_sequence,
            "pose_peptide_sequence": pose_sequence,
            "peptide_residue_count": len(pose_residues),
            "sequence_mapping_valid": True,
        },
        "pose_ca_rmsd_direct_A": direct_ca_rmsd,
        "interface_rmsd_A": interface_rmsd,
        "interface_rmsd_definition": "RMSD of N/CA/C/O atoms for peptide residues with native heavy-atom contacts within 5 A, after receptor-interface C-alpha Kabsch alignment; receptor interface residues are within 10 A of native peptide.",
        "native_contact_cutoff_A": CONFIG["native_contact_cutoff_A"],
        "native_residue_pair_contacts": len(native_contacts),
        "predicted_residue_pair_contacts": len(predicted_contacts),
        "recovered_native_residue_pair_contacts": len(recovered_contacts),
        "native_contact_recovery": len(recovered_contacts) / len(native_contacts) if native_contacts else None,
        "minimum_receptor_peptide_heavy_atom_distance_A": minimum_distance,
        "severe_clash_cutoff_A": CONFIG["severe_clash_cutoff_A"],
        "severe_clash_count": severe_clash_count,
        "arg203_arg205_contacts": arginine_metrics,
        "complete_peptide_backbone": not missing_backbone,
        "missing_backbone_atoms": missing_backbone,
        "finite_coordinates": bool(np.isfinite(peptide_xyz).all()),
        "geometry_valid": geometry_valid,
    }


def _run(args: argparse.Namespace) -> dict:
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"refusing to overwrite non-empty output: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    if not REFERENCE_PDB.exists():
        raise FileNotFoundError(REFERENCE_PDB)
    vina = args.vina.resolve()
    receptor_prep = args.meeko_receptor_prep.resolve()
    if not vina.is_file() or not receptor_prep.is_file():
        raise FileNotFoundError("Vina binary or Meeko receptor-preparation script missing")

    deposited_sequence, observed_sequence, sequence_mapping, receptor_lines, ligand_lines = _parse_pdb_sequences(REFERENCE_PDB)
    if observed_sequence != args.peptide:
        raise ValueError(f"configured peptide {args.peptide} != observed chain B {observed_sequence}")

    source_pdb = output_dir / "reference_complex.pdb"
    clean_receptor_pdb = output_dir / "receptor_chainA_only.pdb"
    source_pdb.write_bytes(REFERENCE_PDB.read_bytes())
    clean_receptor_pdb.write_text("\n".join(receptor_lines + ["TER", "END"]) + "\n")

    ligand_pdbqt = output_dir / "ligand_meeko.pdbqt"
    input_to_serial, prepared_atoms = _prepare_ligand(ligand_lines, ligand_pdbqt)
    receptor_pdbqt = output_dir / "receptor_meeko.pdbqt"
    receptor_command = [
        str(args.python),
        str(receptor_prep),
        "--read_pdb", str(clean_receptor_pdb),
        "-o", str(output_dir / "receptor_meeko"),
        "-p", str(receptor_pdbqt),
    ]
    receptor_result = subprocess.run(receptor_command, capture_output=True, text=True, check=False)
    (output_dir / "receptor_prepare.stdout.txt").write_text(receptor_result.stdout)
    (output_dir / "receptor_prepare.stderr.txt").write_text(receptor_result.stderr)
    if receptor_result.returncode != 0 or not receptor_pdbqt.is_file():
        raise RuntimeError(
            f"Meeko receptor preparation failed ({receptor_result.returncode}); see saved stdout/stderr"
        )

    reference_structure = PDBParser(QUIET=True).get_structure("reference", str(source_pdb))
    reference_model = next(reference_structure.get_models())
    native_peptide_heavy = [
        atom.coord
        for residue in reference_model["B"]
        if is_aa(residue, standard=True)
        for atom in _heavy_atoms(residue)
    ]
    native_coords = np.array(native_peptide_heavy, dtype=float)
    box_center = native_coords.mean(axis=0)
    config = {
        **CONFIG,
        "peptide_sequence": args.peptide,
        "deposited_chain_B_sequence": deposited_sequence,
        "observed_chain_B_sequence": observed_sequence,
        "observed_to_deposited_residue_mapping": sequence_mapping,
        "missing_coordinate_residues": [
            {
                "deposited_sequence_index_1based": index + 1,
                "deposited_residue": deposited_sequence[index],
                "author_residue_number": 15 + index,
            }
            for index in range(len(deposited_sequence))
            if index not in {entry["deposited_sequence_index_1based"] - 1 for entry in sequence_mapping}
        ],
        "peptide_heavy_atom_count": len(native_coords),
        "box_center_A": box_center.tolist(),
        "receptor_chain": "A only; extracted from 1YCR, peptide chain B excluded",
        "receptor_preparation_command": receptor_command,
        "vina_executable": str(vina),
        "vina_version": subprocess.run([str(vina), "--version"], capture_output=True, text=True, check=True).stdout.strip(),
        "python_executable": str(args.python),
        "meeko_version": __import__("importlib.metadata", fromlist=["version"]).version("meeko"),
        "reference_sha256": _sha256(source_pdb),
        "receptor_pdbqt_sha256": _sha256(receptor_pdbqt),
        "ligand_pdbqt_sha256": _sha256(ligand_pdbqt),
        "pose_atom_map_source": "Meeko REMARK INDEX MAP, validated against experimental ligand heavy-atom coordinates before docking",
        "scientific_scope": "computational positive-control pose recovery only",
    }
    box_size = CONFIG["box_size_A"]
    pose_path = output_dir / "vina_raw_out.pdbqt"
    command = [
        str(vina),
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", f"{box_center[0]:.3f}",
        "--center_y", f"{box_center[1]:.3f}",
        "--center_z", f"{box_center[2]:.3f}",
        "--size_x", str(box_size[0]),
        "--size_y", str(box_size[1]),
        "--size_z", str(box_size[2]),
        "--exhaustiveness", str(CONFIG["exhaustiveness"]),
        "--num_modes", str(CONFIG["num_modes"]),
        "--energy_range", str(CONFIG["energy_range_kcal_mol"]),
        "--cpu", str(CONFIG["cpu"]),
        "--seed", str(CONFIG["seed"]),
        "--out", str(pose_path),
    ]
    config["vina_command"] = command
    (output_dir / "inputs.json").write_text(json.dumps(config, indent=2) + "\n")
    (output_dir / "vina_command.txt").write_text(" ".join(command) + "\n")
    started = time.perf_counter()
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    runtime_seconds = time.perf_counter() - started
    (output_dir / "vina.stdout.txt").write_text(result.stdout)
    (output_dir / "vina.stderr.txt").write_text(result.stderr)
    if result.returncode != 0 or not pose_path.is_file():
        failure = {
            "status": "docking_failed",
            "returncode": result.returncode,
            "runtime_seconds": runtime_seconds,
            "stdout_path": str(output_dir / "vina.stdout.txt"),
            "stderr_path": str(output_dir / "vina.stderr.txt"),
            "raw_pose_exists": pose_path.exists(),
        }
        (output_dir / "validation.json").write_text(json.dumps(failure, indent=2) + "\n")
        return failure

    raw_pose_text = pose_path.read_text()
    docked_atoms = _parse_pdbqt_atoms(raw_pose_text, first_model_only=True)
    if set(input_to_serial.values()) != set(docked_atoms):
        raise ValueError(
            "Vina pose atom serials do not match the Meeko source atom map; "
            "pose conversion is intentionally blocked"
        )
    for input_index, source_line in enumerate(ligand_lines, start=1):
        serial = input_to_serial[input_index]
        input_type = prepared_atoms[serial]["line"][77:79].strip()
        output_type = docked_atoms[serial]["line"][77:79].strip()
        if input_type != output_type:
            raise ValueError(f"ligand PDBQT atom type changed at mapped input atom {input_index}")

    pose_pdb = output_dir / "docked_complex_mapped.pdb"
    pose_pdb.write_text(_assemble_pose(receptor_lines, ligand_lines, input_to_serial, docked_atoms))
    validation = _analyze_pose(source_pdb, pose_pdb)
    vina_score = None
    for line in raw_pose_text.splitlines():
        if line.startswith("REMARK VINA RESULT:"):
            vina_score = float(line.split()[3])
            break
    if vina_score is None:
        validation["vina_score_kcal_mol"] = None
        validation["vina_score_parse_status"] = "missing_remark"
    else:
        validation["vina_score_kcal_mol"] = vina_score
        validation["vina_score_parse_status"] = "parsed_from_raw_pdbqt_remark"
    validation.update({
        "status": "completed",
        "backend": "AutoDock Vina-only; separate from Hybridock/RAPiDock",
        "runtime_seconds": runtime_seconds,
        "raw_pose_pdbqt": str(pose_path),
        "mapped_pose_pdb": str(pose_pdb),
        "raw_pose_sha256": _sha256(pose_path),
        "config": config,
    })
    (output_dir / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    return validation


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=BENCHMARK_ROOT / "vina_positive_control_20261016")
    parser.add_argument("--peptide", default=OBSERVED_PEPTIDE_SEQUENCE)
    parser.add_argument("--python", type=Path, default=Path(sys.executable))
    parser.add_argument("--vina", type=Path, default=VINA_DEFAULT)
    parser.add_argument("--meeko-receptor-prep", type=Path, default=MEEKO_PREP_DEFAULT)
    args = parser.parse_args()
    try:
        result = _run(args)
    except FileExistsError:
        raise
    except Exception as exc:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        failure = {
            "status": "workflow_failed",
            "error_type": type(exc).__name__,
            "error": str(exc),
            "raw_docking_started": (args.output_dir / "vina_raw_out.pdbqt").exists(),
            "backend": "Vina-only; separate from Hybridock/RAPiDock",
        }
        (args.output_dir / "validation.json").write_text(json.dumps(failure, indent=2) + "\n")
        raise
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()