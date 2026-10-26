"""Independent target-aware structural scoring for supplied complex PDB files.

This module scores a *pose already represented in one coordinate frame*. It does
not generate poses, predict affinity, or modify the legacy native c_score.
Scoring version: target_aware_structural_score_v1.1_clash_gated.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

SCORE_VERSION = "target_aware_structural_score_v1.1_clash_gated"
CONTACT_CUTOFFS_A = (3.5, 4.0, 5.0, 8.0)
SEVERE_CLASH_CUTOFF_A = 2.0
FAR_SEPARATION_WARNING_A = 25.0

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "HID": "H", "HIE": "H", "HIP": "H",
}
TWO_LETTER_ELEMENTS = {"CL", "BR", "SE", "ZN", "FE", "MG", "MN", "CU", "CO", "NI", "NA", "K", "CA"}


class StructuralInputError(ValueError):
    """Invalid or incompatible structural input; no score should be emitted."""


@dataclass(frozen=True)
class HotspotSpec:
    chain: str
    residue_number: int
    residue_name: str = "ARG"
    insertion_code: str = ""

    def label(self) -> str:
        ins = self.insertion_code or ""
        return f"{self.chain}:{self.residue_number}{ins}:{self.residue_name.upper()}"


@dataclass(frozen=True)
class StructuralScoreConfig:
    contact_cutoffs_A: tuple[float, ...] = CONTACT_CUTOFFS_A
    severe_clash_cutoff_A: float = SEVERE_CLASH_CUTOFF_A
    far_separation_warning_A: float = FAR_SEPARATION_WARNING_A
    hotspot_residues: tuple[HotspotSpec, ...] = field(default_factory=tuple)
    # Formula constants are fixed in v1; change them only under a new version.
    residue_pair_saturation_count: int = 20
    interface_residue_saturation_count: int = 12
    clash_saturation_pair_count: int = 5
    weights: tuple[float, float, float] = (0.45, 0.30, 0.25)
    clash_penalty_weight: float = 0.35
    geometry_optimum_A: float = 3.5
    geometry_sigma_A: float = 1.5


@dataclass(frozen=True)
class Atom:
    chain: str
    residue_number: int
    insertion_code: str
    residue_name: str
    atom_name: str
    element: str
    x: float
    y: float
    z: float

    @property
    def residue_key(self) -> tuple[str, int, str, str]:
        return (self.chain, self.residue_number, self.insertion_code, self.residue_name)

    @property
    def xyz(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _element_from_line(line: str) -> str:
    explicit = line[76:78].strip().upper() if len(line) >= 78 else ""
    if explicit:
        return explicit
    raw = line[12:16]
    stripped = raw.strip().upper()
    if not stripped:
        return ""
    # PDB atom-name alignment distinguishes e.g. " CA " (carbon alpha) from
    # right-justified two-letter elements. If alignment is absent, prefer a
    # conservative one-letter interpretation.
    if raw and raw[0] != " " and len(stripped) >= 2 and stripped[:2] in TWO_LETTER_ELEMENTS:
        return stripped[:2]
    return stripped[0]


def _parse_pdb_atoms(path: str | Path, chains: set[str]) -> list[Atom]:
    pdb_path = Path(path)
    if not pdb_path.is_file():
        raise FileNotFoundError(f"Structure file not found: {pdb_path}")
    try:
        lines = pdb_path.read_text(encoding="utf-8", errors="strict").splitlines()
    except (OSError, UnicodeError) as exc:
        raise StructuralInputError(f"Cannot read PDB text: {exc}") from exc

    has_models = any(line.startswith("MODEL ") for line in lines)
    in_first_model = not has_models
    atoms: list[Atom] = []
    for line_no, line in enumerate(lines, start=1):
        if has_models and line.startswith("MODEL "):
            if not in_first_model:
                in_first_model = True
            continue
        if has_models and line.startswith("ENDMDL") and in_first_model:
            break
        if has_models and not in_first_model:
            continue
        record = line[:6].strip().upper()
        if record not in {"ATOM", "HETATM"}:
            continue
        if len(line) < 22:
            # Only error if the chain could have been identified as requested.
            raise StructuralInputError(f"Malformed atom record at line {line_no}: fewer than 22 columns")
        chain = line[21:22].strip()
        if chain not in chains:
            continue
        if len(line) < 54:
            raise StructuralInputError(f"Malformed atom coordinates at line {line_no}: fewer than 54 columns")
        altloc = line[16:17]
        if altloc not in {" ", "A", ""}:
            continue
        residue_name = line[17:20].strip().upper()
        # Restrict this v1 protein/peptide scorer to recognized amino-acid
        # residues. It does not silently count waters or arbitrary ligands.
        if residue_name not in AA3_TO_1:
            continue
        atom_name = line[12:16].strip().upper()
        element = _element_from_line(line)
        if not atom_name or not element:
            raise StructuralInputError(f"Missing atom name/element at line {line_no}")
        if element in {"H", "D"}:
            continue
        try:
            residue_number = int(line[22:26].strip())
            insertion_code = line[26:27].strip()
            x, y, z = (float(line[a:b].strip()) for a, b in ((30, 38), (38, 46), (46, 54)))
        except (ValueError, IndexError) as exc:
            raise StructuralInputError(f"Malformed residue/coordinates at line {line_no}") from exc
        if not all(math.isfinite(v) for v in (x, y, z)):
            raise StructuralInputError(f"Non-finite coordinate at line {line_no}")
        atoms.append(Atom(chain, residue_number, insertion_code, residue_name,
                          atom_name, element, x, y, z))
    return atoms


def _group_residues(atoms: Sequence[Atom]) -> dict[tuple[str, int, str, str], list[Atom]]:
    groups: dict[tuple[str, int, str, str], list[Atom]] = {}
    for atom in atoms:
        groups.setdefault(atom.residue_key, []).append(atom)
    return groups


def _chain_sequence(atoms: Sequence[Atom], chain: str) -> str:
    seen: set[tuple[int, str, str]] = set()
    letters: list[str] = []
    for atom in atoms:
        if atom.chain != chain:
            continue
        key = (atom.residue_number, atom.insertion_code, atom.residue_name)
        if key in seen:
            continue
        seen.add(key)
        letters.append(AA3_TO_1[atom.residue_name])
    return "".join(letters)


def inspect_complex_structure(complex_structure: str | Path) -> dict[str, Any]:
    """List chains/residue numbering so callers can verify hotspot IDs first."""
    path = Path(complex_structure)
    # Parse all plausible protein chains by discovering chain IDs from columns.
    if not path.is_file():
        raise FileNotFoundError(f"Structure file not found: {path}")
    text = path.read_text(encoding="utf-8", errors="strict")
    chains: set[str] = set()
    for line in text.splitlines():
        if line[:6].strip().upper() in {"ATOM", "HETATM"} and len(line) >= 22:
            if line[21:22].strip():
                chains.add(line[21:22].strip())
    atoms = _parse_pdb_atoms(path, chains)
    by_chain: dict[str, list[dict[str, Any]]] = {}
    for chain in sorted(chains):
        residues = []
        seen = set()
        for atom in atoms:
            if atom.chain != chain:
                continue
            key = (atom.residue_number, atom.insertion_code, atom.residue_name)
            if key in seen:
                continue
            seen.add(key)
            residues.append({"residue_number": atom.residue_number,
                             "insertion_code": atom.insertion_code,
                             "residue_name": atom.residue_name,
                             "one_letter": AA3_TO_1[atom.residue_name]})
        by_chain[chain] = residues
    return {
        "path": str(path), "sha256": sha256_file(path),
        "chains": {chain: {"residue_count": len(residues),
                            "sequence": "".join(r["one_letter"] for r in residues),
                            "residues": residues}
                   for chain, residues in by_chain.items()},
        "note": "Use these author residue identifiers to verify chain and hotspot mapping before scoring.",
    }


def _atom_chain_checksum(atoms: Sequence[Atom], chain: str) -> str:
    canonical = [
        [a.chain, a.residue_number, a.insertion_code, a.residue_name,
         a.atom_name, a.element, round(a.x, 5), round(a.y, 5), round(a.z, 5)]
        for a in atoms if a.chain == chain
    ]
    payload = json.dumps(canonical, separators=(",", ":"), ensure_ascii=True).encode()
    return hashlib.sha256(payload).hexdigest()


def make_structural_cache_key(peptide_sequence: str, receptor_checksum: str,
                              pose_checksum: str, config: StructuralScoreConfig | None = None) -> str:
    """Cache key contract: sequence + receptor + pose + score version + config."""
    cfg = config or StructuralScoreConfig()
    payload = {
        "peptide_sequence": peptide_sequence,
        "receptor_checksum": receptor_checksum,
        "pose_checksum": pose_checksum,
        "score_version": SCORE_VERSION,
        "config": asdict(cfg),
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=lambda x: asdict(x)).encode()
    return hashlib.sha256(raw).hexdigest()


def _failed_result(path: str | Path, target_structure_id: str | None,
                   receptor_chain: str, peptide_chain: str, status: str,
                   reason: str, provenance: str) -> dict[str, Any]:
    return {
        "c_score": None,
        "c_score_status": "not_computed_by_structural_scorer",
        "target_aware_structural_score": None,
        "structural_score_status": status,
        "score_valid": False,
        "target_structure_id": target_structure_id,
        "peptide_sequence": None,
        "score_version": SCORE_VERSION,
        "receptor_chain": receptor_chain,
        "peptide_chain": peptide_chain,
        "structure_provenance": {"kind": provenance, "path": str(path),
                                 "coordinate_frame": "single_complex_pdb_not_transformed",
                                 "transformations_applied": []},
        "structural_metrics": {},
        "eligibility_reasons": [reason],
    }


def score_peptide_receptor_complex(
    receptor_structure: str | Path,
    peptide_structure: str | Path | None = None,
    config: StructuralScoreConfig | None = None,
    *,
    receptor_chain: str,
    peptide_chain: str,
    target_structure_id: str | None = None,
    structure_provenance: str = "user_supplied_complex_pdb",
) -> dict[str, Any]:
    """Score chains within one complex PDB, ensuring a shared coordinate frame.

    `receptor_structure` must be a single PDB containing both chains in the same
    coordinate frame. `peptide_structure` is reserved for API clarity and must
    be None or the same path; separate files are rejected because relative pose
    cannot be inferred safely without an explicit documented transform.
    """
    path = Path(receptor_structure)
    if peptide_structure is not None and Path(peptide_structure).resolve() != path.resolve():
        return _failed_result(path, target_structure_id, receptor_chain, peptide_chain,
                              "invalid_input", "receptor_and_peptide_must_share_one_complex_pdb_or_explicit_alignment", structure_provenance)
    cfg = config or StructuralScoreConfig()
    if receptor_chain == peptide_chain or not receptor_chain or not peptide_chain:
        return _failed_result(path, target_structure_id, receptor_chain, peptide_chain,
                              "invalid_input", "receptor_and_peptide_chain_ids_must_be_distinct_and_nonempty", structure_provenance)
    try:
        atoms = _parse_pdb_atoms(path, {receptor_chain, peptide_chain})
        receptor_atoms = [a for a in atoms if a.chain == receptor_chain]
        peptide_atoms = [a for a in atoms if a.chain == peptide_chain]
        if not receptor_atoms:
            raise StructuralInputError(f"no_valid_heavy_atoms_in_receptor_chain:{receptor_chain}")
        if not peptide_atoms:
            raise StructuralInputError(f"no_valid_heavy_atoms_in_peptide_chain:{peptide_chain}")
        receptor_residues = _group_residues(receptor_atoms)
        peptide_residues = _group_residues(peptide_atoms)
        if not receptor_residues or not peptide_residues:
            raise StructuralInputError("both chains must contain recognized amino-acid residues")
        for hotspot in cfg.hotspot_residues:
            if hotspot.chain != receptor_chain:
                raise StructuralInputError(f"hotspot_chain_mismatch:{hotspot.label()} expected receptor chain {receptor_chain}")
            matches = [k for k in receptor_residues if k[0] == hotspot.chain and k[1] == hotspot.residue_number
                       and k[2] == hotspot.insertion_code and k[3] == hotspot.residue_name.upper()]
            if not matches:
                raise StructuralInputError(f"hotspot_not_found_or_renumbered:{hotspot.label()}")
    except FileNotFoundError as exc:
        return _failed_result(path, target_structure_id, receptor_chain, peptide_chain,
                              "preparation_failure", str(exc), structure_provenance)
    except (StructuralInputError, UnicodeError) as exc:
        return _failed_result(path, target_structure_id, receptor_chain, peptide_chain,
                              "invalid_input", str(exc), structure_provenance)
    except Exception as exc:  # Keep scoring failures explicit in the output schema.
        return _failed_result(path, target_structure_id, receptor_chain, peptide_chain,
                              "scoring_failure", f"unexpected_parse_error:{type(exc).__name__}:{exc}", structure_provenance)

    try:
        rec_xyz = np.asarray([a.xyz for a in receptor_atoms], dtype=float)
        pep_xyz = np.asarray([a.xyz for a in peptide_atoms], dtype=float)
        cutoffs = tuple(sorted(set(float(c) for c in cfg.contact_cutoffs_A)))
        if not cutoffs or any((not math.isfinite(c) or c <= 0) for c in cutoffs):
            raise StructuralInputError("contact cutoffs must be finite positive distances")
        if not set(CONTACT_CUTOFFS_A).issubset(set(cutoffs)):
            raise StructuralInputError(f"v1 requires all prespecified cutoffs: {CONTACT_CUTOFFS_A}")
        if (not math.isfinite(cfg.severe_clash_cutoff_A) or cfg.severe_clash_cutoff_A <= 0
                or not math.isfinite(cfg.far_separation_warning_A) or cfg.far_separation_warning_A <= 0
                or not math.isfinite(cfg.geometry_sigma_A) or cfg.geometry_sigma_A <= 0
                or cfg.residue_pair_saturation_count < 1 or cfg.interface_residue_saturation_count < 1
                or cfg.clash_saturation_pair_count < 1
                or len(cfg.weights) != 3 or any(not math.isfinite(w) or w < 0 for w in cfg.weights)
                or not math.isfinite(cfg.clash_penalty_weight) or cfg.clash_penalty_weight < 0):
            raise StructuralInputError("invalid scoring parameters: scales must be positive and weights/penalties finite nonnegative")
        atom_counts = {c: 0 for c in cutoffs}
        residue_pairs: dict[float, set[tuple[tuple[str, int, str, str], tuple[str, int, str, str]]]] = {c: set() for c in cutoffs}
        receptor_interface: dict[float, set[tuple[str, int, str, str]]] = {c: set() for c in cutoffs}
        peptide_interface: dict[float, set[tuple[str, int, str, str]]] = {c: set() for c in cutoffs}
        min_distance = float("inf")
        clash_count = 0
        rkeys = [a.residue_key for a in receptor_atoms]
        pkeys = [a.residue_key for a in peptide_atoms]
        # Chunk receptor atoms to bound peak memory on large receptor structures.
        max_cutoff = max(max(cutoffs), cfg.severe_clash_cutoff_A)
        for start in range(0, len(rec_xyz), 256):
            stop = min(start + 256, len(rec_xyz))
            delta = rec_xyz[start:stop, None, :] - pep_xyz[None, :, :]
            dist2 = np.einsum("ijk,ijk->ij", delta, delta)
            if dist2.size:
                local_min = float(np.sqrt(float(np.min(dist2))))
                min_distance = min(min_distance, local_min)
            clash_count += int(np.count_nonzero(dist2 < cfg.severe_clash_cutoff_A ** 2))
            for cutoff in cutoffs:
                mask = dist2 <= cutoff ** 2
                atom_counts[cutoff] += int(np.count_nonzero(mask))
                rows, cols = np.where(mask)
                for rr, cc in zip(rows.tolist(), cols.tolist()):
                    ri, pi = start + rr, cc
                    residue_pairs[cutoff].add((rkeys[ri], pkeys[pi]))
                    receptor_interface[cutoff].add(rkeys[ri])
                    peptide_interface[cutoff].add(pkeys[pi])
        if not math.isfinite(min_distance):
            raise StructuralInputError("no_finite_interchain_distance")

        # Hotspots are validated by chain, author residue number, insertion code,
        # and residue name before these metrics are calculated.
        hotspot_metrics = []
        for hotspot in cfg.hotspot_residues:
            hkey = (hotspot.chain, hotspot.residue_number, hotspot.insertion_code, hotspot.residue_name.upper())
            hatoms = [a for a in receptor_atoms if a.residue_key == hkey]
            hxyz = np.asarray([a.xyz for a in hatoms], dtype=float)
            hd = hxyz[:, None, :] - pep_xyz[None, :, :]
            hd2 = np.einsum("ijk,ijk->ij", hd, hd)
            hmask = hd2 <= 5.0 ** 2
            _, pidx = np.where(hmask)
            hotspot_metrics.append({
                "hotspot": hotspot.label(),
                "atom_pairs_within_5A": int(np.count_nonzero(hmask)),
                "distinct_peptide_residues_within_5A": len({pkeys[int(i)] for i in pidx.tolist()}),
                "minimum_distance_A": float(np.sqrt(float(np.min(hd2)))),
            })

        cutoff_key = lambda c: f"{c:g}A"
        atom_contacts = {cutoff_key(c): atom_counts[c] for c in cutoffs}
        residue_pair_counts = {cutoff_key(c): len(residue_pairs[c]) for c in cutoffs}
        unique_residue_pair_5 = len(residue_pairs.get(5.0, set())) if 5.0 in residue_pairs else len(residue_pairs[min(cutoffs, key=lambda x: abs(x - 5.0))])
        iface_rec_5 = len(receptor_interface.get(5.0, set())) if 5.0 in receptor_interface else len(receptor_interface[min(cutoffs, key=lambda x: abs(x - 5.0))])
        iface_pep_5 = len(peptide_interface.get(5.0, set())) if 5.0 in peptide_interface else len(peptide_interface[min(cutoffs, key=lambda x: abs(x - 5.0))])
        # v1 score: normalized residue-pair diversity + interface coverage + a
        # smooth near-contact geometry term - explicit severe-clash penalty.
        # No SASA, electrostatics, affinity, or learned weights are asserted.
        contact_component = min(unique_residue_pair_5 / max(cfg.residue_pair_saturation_count, 1), 1.0)
        interface_component = min((iface_rec_5 + iface_pep_5) / max(cfg.interface_residue_saturation_count, 1), 1.0)
        geometry_component = math.exp(-0.5 * ((min_distance - cfg.geometry_optimum_A) / cfg.geometry_sigma_A) ** 2)
        clash_component = min(clash_count / max(cfg.clash_saturation_pair_count, 1), 1.0)
        w_contact, w_interface, w_geometry = cfg.weights
        raw_score = (w_contact * contact_component + w_interface * interface_component
                     + w_geometry * geometry_component - cfg.clash_penalty_weight * clash_component)
        score = float(min(1.0, max(0.0, raw_score)))
        reasons = []
        if clash_count:
            reasons.append("severe_interchain_overlaps_lt_2A;structural_score_rejected_and_excluded_from_ranking")
        contacts_5 = atom_contacts.get("5A", 0)
        if contacts_5 == 0:
            reasons.append("no_interchain_heavy_atom_contacts_within_5A")
        if min_distance > cfg.far_separation_warning_A:
            reasons.append("peptide_far_from_receptor_gt_25A_coordinate_frame_or_decoy_check_required")
        input_sha = sha256_file(path)
        receptor_hash = _atom_chain_checksum(atoms, receptor_chain)
        pose_hash = _atom_chain_checksum(atoms, peptide_chain)
        peptide_sequence = _chain_sequence(peptide_atoms, peptide_chain)
        cache_key = make_structural_cache_key(peptide_sequence, receptor_hash, pose_hash, cfg)
        return {
            "c_score": None,
            "c_score_status": "not_computed_by_structural_scorer",
            # The composite remains visible only as a diagnostic when hard-gated.
            # Official structural score is null for any severe inter-chain overlap.
            "target_aware_structural_score": score if clash_count == 0 else None,
            "structural_score_status": "success" if clash_count == 0 else "invalid_complex",
            "score_valid": clash_count == 0,
            "target_structure_id": target_structure_id,
            "peptide_sequence": peptide_sequence,
            "score_version": SCORE_VERSION,
            "receptor_chain": receptor_chain,
            "peptide_chain": peptide_chain,
            "structure_provenance": {
                "kind": structure_provenance, "path": str(path), "sha256": input_sha,
                "coordinate_frame": "single_complex_pdb_not_transformed",
                "transformations_applied": [],
                "receptor_chain_coordinate_sha256": receptor_hash,
                "peptide_pose_coordinate_sha256": pose_hash,
                "cache_key": cache_key,
                "model_policy": "first_MODEL_only_if_multiple_models",
                "altloc_policy": "blank_or_A_only",
                "protonation_assumptions": "none; supplied heavy-atom coordinates used as-is",
            },
            "structural_metrics": {
                "composite_score_diagnostic_before_clash_gate": score,
                "units": "angstroms_for_distances",
                "interchain_heavy_atom_contact_pairs": atom_contacts,
                "distinct_interchain_residue_pairs": residue_pair_counts,
                "minimum_interchain_heavy_atom_distance_A": min_distance,
                "severe_clash_cutoff_A": cfg.severe_clash_cutoff_A,
                "severe_clash_heavy_atom_pairs": clash_count,
                "interface_receptor_residue_count_within_5A": iface_rec_5,
                "interface_peptide_residue_count_within_5A": iface_pep_5,
                "hotspot_metrics": hotspot_metrics,
                "buried_surface_area_A2": None,
                "buried_surface_area_status": "not_calculated_no_SASA_method_used",
                "components": {
                    "contact_residue_pair_diversity_0_1": contact_component,
                    "interface_residue_coverage_0_1": interface_component,
                    "near_contact_geometry_0_1": geometry_component,
                    "severe_clash_penalty_0_1": clash_component,
                },
                "simple_baseline_objectives_0_1": {
                    "contact_only": contact_component,
                    "geometry_only": geometry_component,
                    "clash_only_no_clash_is_best": 1.0 - clash_component,
                },
                "score_formula": (
                    f"clip({w_contact:g}*min(residue_pairs_5A/{cfg.residue_pair_saturation_count},1) + "
                    f"{w_interface:g}*min((receptor_interface_residues_5A+peptide_interface_residues_5A)/{cfg.interface_residue_saturation_count},1) + "
                    f"{w_geometry:g}*exp(-0.5*((min_distance_A-{cfg.geometry_optimum_A:g})/{cfg.geometry_sigma_A:g})^2) - "
                    f"{cfg.clash_penalty_weight:g}*min(severe_clash_atom_pairs_lt_2A/{cfg.clash_saturation_pair_count},1),0,1)"
                ),
                "score_weights": {"residue_pair_diversity": w_contact, "interface_residue_coverage": w_interface,
                                  "near_contact_geometry": w_geometry, "severe_clash_penalty": cfg.clash_penalty_weight},
            },
            "eligibility_reasons": reasons,
        }
    except StructuralInputError as exc:
        return _failed_result(path, target_structure_id, receptor_chain, peptide_chain,
                              "invalid_input", str(exc), structure_provenance)
    except Exception as exc:
        return _failed_result(path, target_structure_id, receptor_chain, peptide_chain,
                              "scoring_failure", f"unexpected_scoring_error:{type(exc).__name__}:{exc}", structure_provenance)


def aligned_peptide_ca_rmsd(
    reference_complex_pdb: str | Path,
    moving_complex_pdb: str | Path,
    *, receptor_chain: str = "A",
    peptide_chain: str = "B",
    moving_receptor_chain: str | None = None,
    moving_peptide_chain: str | None = None,
) -> dict[str, Any]:
    """Align receptors, then calculate peptide C-alpha RMSD.

    Receptor residues match by author ID, insertion code and residue name.
    Peptide residues match in sequence order only if their complete observed
    sequences are identical. This allows numbering offsets without pairing
    different peptide sequences.
    """
    moving_receptor_chain = moving_receptor_chain or receptor_chain
    moving_peptide_chain = moving_peptide_chain or peptide_chain

    ref_atoms = _parse_pdb_atoms(
        reference_complex_pdb, {receptor_chain, peptide_chain}
    )
    mov_atoms = _parse_pdb_atoms(
        moving_complex_pdb, {moving_receptor_chain, moving_peptide_chain}
    )

    def ca_map(atoms: Sequence[Atom], chain: str):
        out = {}
        for atom in atoms:
            if atom.chain == chain and atom.atom_name == "CA":
                key = (atom.residue_number, atom.insertion_code,
                       atom.residue_name)
                out[key] = np.asarray(atom.xyz, dtype=float)
        return out

    def ordered_peptide_ca(atoms: Sequence[Atom], chain: str):
        records = {}
        for atom in atoms:
            if atom.chain == chain and atom.atom_name == "CA":
                key = (atom.residue_number, atom.insertion_code,
                       atom.residue_name)
                records[key] = np.asarray(atom.xyz, dtype=float)

        ordered = sorted(
            records.items(),
            key=lambda item: (item[0][0], item[0][1], item[0][2]),
        )
        try:
            sequence = "".join(AA3_TO_1[key[2]] for key, _ in ordered)
        except KeyError as exc:
            raise StructuralInputError(
                f"unsupported_residue_for_peptide_RMSD:{exc}"
            ) from exc

        coords = (
            np.vstack([xyz for _, xyz in ordered])
            if ordered else np.empty((0, 3), dtype=float)
        )
        return sequence, coords, [key for key, _ in ordered]

    ref_rec = ca_map(ref_atoms, receptor_chain)
    mov_rec = ca_map(mov_atoms, moving_receptor_chain)
    rec_keys = sorted(set(ref_rec) & set(mov_rec))

    if len(rec_keys) < 3:
        raise StructuralInputError(
            "receptor_alignment_requires_at_least_3_corresponding_CA_atoms;"
            f"found:{len(rec_keys)}"
        )

    ref_seq, ref_pep, ref_keys = ordered_peptide_ca(
        ref_atoms, peptide_chain
    )
    mov_seq, mov_pep, mov_keys = ordered_peptide_ca(
        mov_atoms, moving_peptide_chain
    )

    if len(ref_keys) < 3 or len(mov_keys) < 3:
        raise StructuralInputError(
            "peptide_RMSD_requires_at_least_3_CA_atoms_in_each_chain;"
            f"reference:{len(ref_keys)};moving:{len(mov_keys)}"
        )

    if ref_seq != mov_seq:
        raise StructuralInputError(
            "peptide_sequence_mismatch_for_RMSD;"
            f"reference:{ref_seq};moving:{mov_seq}"
        )

    moving_rec = np.vstack([mov_rec[key] for key in rec_keys])
    reference_rec = np.vstack([ref_rec[key] for key in rec_keys])
    moving_center = moving_rec.mean(axis=0)
    reference_center = reference_rec.mean(axis=0)

    covariance = (
        (moving_rec - moving_center).T
        @ (reference_rec - reference_center)
    )
    u, _, vt = np.linalg.svd(covariance)
    rotation = u @ vt

    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt

    aligned_peptide = (
        (mov_pep - moving_center) @ rotation + reference_center
    )
    delta = aligned_peptide - ref_pep
    rmsd = float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))

    return {
        "reference_complex_sha256": sha256_file(reference_complex_pdb),
        "moving_complex_sha256": sha256_file(moving_complex_pdb),
        "reference_receptor_chain": receptor_chain,
        "reference_peptide_chain": peptide_chain,
        "moving_receptor_chain": moving_receptor_chain,
        "moving_peptide_chain": moving_peptide_chain,
        "receptor_CA_correspondences": len(rec_keys),
        "peptide_CA_correspondences": len(ref_keys),
        "peptide_reference_sequence": ref_seq,
        "peptide_moving_sequence": mov_seq,
        "peptide_correspondence_method":
            "exact_sequence_order_ignore_author_residue_number",
        "receptor_aligned_peptide_CA_RMSD_A": rmsd,
        "alignment_method":
            "Kabsch fit on common receptor CA atoms; transform applied to peptide",
        "residue_correspondence":
            "receptor: author ID/insertion/name; peptide: ordered exact sequence",
    }

def translate_peptide_chain_pdb(source_pdb: str | Path, output_pdb: str | Path,
                                peptide_chain: str, dx: float, dy: float, dz: float) -> dict[str, Any]:
    """Create a geometry-control decoy by translating only the selected peptide chain."""
    src, dst = Path(source_pdb), Path(output_pdb)
    if not src.is_file():
        raise FileNotFoundError(src)
    if not all(math.isfinite(v) for v in (dx, dy, dz)):
        raise ValueError("translation values must be finite")
    out_lines = []
    changed = 0
    has_models = any(line.startswith("MODEL ") for line in src.read_text().splitlines())
    in_first_model = not has_models
    for line in src.read_text().splitlines(keepends=True):
        bare = line.rstrip("\r\n")
        if has_models and bare.startswith("MODEL "):
            if not in_first_model:
                in_first_model = True
            out_lines.append(line)
            continue
        if has_models and bare.startswith("ENDMDL") and in_first_model:
            out_lines.append(line)
            break
        if has_models and not in_first_model:
            continue
        if bare[:6].strip().upper() in {"ATOM", "HETATM"} and len(bare) >= 54 and bare[21:22].strip() == peptide_chain:
            try:
                xyz = [float(bare[a:b].strip()) for a, b in ((30, 38), (38, 46), (46, 54))]
            except ValueError as exc:
                raise StructuralInputError(f"malformed peptide coordinates while creating decoy: {bare}") from exc
            moved = [xyz[0] + dx, xyz[1] + dy, xyz[2] + dz]
            if not all(math.isfinite(v) for v in moved):
                raise StructuralInputError("non-finite translated coordinate")
            bare = bare[:30] + "".join(f"{v:8.3f}" for v in moved) + bare[54:]
            line = bare + ("\n" if line.endswith("\n") else "")
            changed += 1
        out_lines.append(line)
    if changed == 0:
        raise StructuralInputError(f"no atom records found for peptide chain {peptide_chain}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text("".join(out_lines), encoding="utf-8")
    return {"source": str(src), "output": str(dst), "peptide_chain": peptide_chain,
            "translation_A": {"x": dx, "y": dy, "z": dz}, "translated_atom_records": changed,
            "control_type": "computational_translated_pose_decoy_not_experimental_nonbinder"}


def parse_hotspot_spec(value: str) -> HotspotSpec:
    """Parse CHAIN:RESNUM[:RESNAME], e.g. A:203:ARG or A:203A:ARG."""
    import re
    match = re.fullmatch(r"([^:]):(-?\d+)([A-Za-z]?):([A-Za-z]{3})", value.strip())
    if not match:
        raise ValueError(f"invalid hotspot {value!r}; expected CHAIN:RESNUM[:INSERTION]:RESNAME, e.g. A:203:ARG")
    chain, number, insertion, resname = match.groups()
    return HotspotSpec(chain, int(number), resname.upper(), insertion.upper())
