#!/usr/bin/env python3
"""
Standalone HERV-K SU / peptide docking verification workflow.

Purpose
-------
Re-run one 21-residue peptide against a user-defined HERV-K SU docking box
with AutoDock Vina at controlled, higher search effort, while preserving the
input predicted peptide conformation as a rigid body by default.

The script deliberately keeps three things separate:
  1. raw AutoDock Vina affinity from stdout,
  2. the Vina PDBQT REMARK affinity cross-check,
  3. any score stored by a larger pipeline.

Inputs
------
--receptor-pdb   HERV-K SU PDB (e.g. herv_k_su.pdb)
--peptide-pdb    predicted 3-D peptide PDB for the candidate
--vina           local Vina executable; auto-discovers `vina` if omitted
--obabel         Open Babel executable; auto-discovers `obabel`/`babel`
--box-center     explicit center_x center_y center_z in Angstrom
--box-size       explicit size_x size_y size_z in Angstrom

Optional
--------
--anchor-residue CHAIN:RESNAME:RESNUM ...
    Computes a box center from the heavy atoms of the named receptor residues.
    This is a convenience for generating/checking the box center; the box
    size remains explicit.
--sequence       sequence annotation for the report; does not modify docking
--num-modes      number of Vina poses to retain (default 9)
--exhaustiveness Vina exhaustiveness (default 32)
--cpu             Vina CPU count; 0 means Vina chooses available CPUs
--seed            Vina random seed
--rigid           make the prepared peptide completely rigid (default)
--keep-flex       keep Meeko's default ligand torsion tree instead

Notes
-----
- A rigid 21-mer is a useful diagnostic first pass because unrestricted peptide
  flexibility can create a very large search problem.
- pLDDT is not itself evidence that a peptide is alpha-helical. If an alpha
  helix is the model you intend to test, verify the secondary structure from
  the actual predicted coordinates (e.g. DSSP/PyMOL) before interpreting it.
- The target-side HERV-K CD98HC interface is not fully resolved. Use a box
  supported by the structure/literature and the actual PDB coordinates, rather
  than assuming a literature residue number maps directly between structures.

Outputs
-------
<outdir>/
  receptor.pdbqt
  peptide_prepared.pdbqt
  peptide_rigid.pdbqt       (if --rigid)
  vina_out.pdbqt
  vina.log
  vina.stdout.txt
  vina.stderr.txt
  vina_command.txt
  vina_pose1.pdbqt
  vina_pose1.pdb
  docked_complex_mode1.pdb
  docking_box.pdb
  docking_summary.json
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Iterable

try:
    from Bio.PDB import PDBParser
except ImportError:
    PDBParser = None


AFFINITY_REMARK = re.compile(
    r"REMARK\s+VINA RESULT:\s*(-?\d+(?:\.\d+)?)"
)
AFFINITY_TABLE = re.compile(
    r"^\s*1\s+(-?\d+(?:\.\d+)?)\s+",
    re.MULTILINE,
)
MODE_TABLE = re.compile(
    r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s*$",
    re.MULTILINE,
)


def run_checked(
    cmd: list[str],
    *,
    cwd: Path | None = None,
    stdout_path: Path | None = None,
    stderr_path: Path | None = None,
) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        cmd,
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )
    if stdout_path is not None:
        stdout_path.write_text(result.stdout, encoding="utf-8")
    if stderr_path is not None:
        stderr_path.write_text(result.stderr, encoding="utf-8")
    if result.returncode != 0:
        detail = result.stderr.strip() or result.stdout.strip()
        raise RuntimeError(
            f"Command failed with exit code {result.returncode}: {' '.join(cmd)}\n"
            f"{detail[-1500:]}"
        )
    return result


def find_executable(explicit: str | None, names: Iterable[str]) -> str:
    if explicit:
        path = Path(explicit).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(f"Executable not found: {path}")
        if not os.access(path, os.X_OK):
            raise PermissionError(f"Executable is not executable: {path}")
        return str(path)
    for name in names:
        found = shutil.which(name)
        if found:
            return found
    raise FileNotFoundError(
        f"Could not find executable. Tried: {', '.join(names)}"
    )


def parse_triplet(values: list[str], label: str) -> tuple[float, float, float]:
    try:
        nums = tuple(float(x) for x in values)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Invalid {label}: {values}") from exc
    if len(nums) != 3 or any(not math.isfinite(x) for x in nums):
        raise argparse.ArgumentTypeError(f"{label} must contain 3 finite numbers")
    return nums  # type: ignore[return-value]


def parse_anchor(value: str) -> tuple[str, str, int]:
    parts = value.split(":")
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(
            "Anchor must be CHAIN:RESNAME:RESNUM, e.g. A:ARG:203"
        )
    chain, resname, resnum_s = parts
    try:
        resnum = int(resnum_s)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Invalid residue number: {value}") from exc
    return chain, resname.upper(), resnum


def anchor_center(receptor_pdb: Path, anchors: list[tuple[str, str, int]]) -> tuple[float, float, float]:
    if PDBParser is None:
        raise RuntimeError(
            "Biopython is required for --anchor-residue. Install biopython or "
            "use explicit --box-center values."
        )
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", str(receptor_pdb))
    points: list[tuple[float, float, float]] = []
    requested = set(anchors)

    for model in structure:
        for chain in model:
            for residue in chain:
                hetflag, resseq, _icode = residue.id
                if hetflag.strip():
                    continue
                key = (chain.id, residue.resname.upper(), int(resseq))
                if key not in requested:
                    continue
                for atom in residue:
                    if atom.element and atom.element.upper() == "H":
                        continue
                    coord = atom.coord
                    points.append((float(coord[0]), float(coord[1]), float(coord[2])))
        # first model is sufficient for a receptor box
        break

    missing = requested - {
        (chain.id, residue.resname.upper(), int(residue.id[1]))
        for model in structure
        for chain in model
        for residue in chain
        if not residue.id[0].strip()
    }
    if missing:
        missing_text = ", ".join(f"{c}:{r}:{n}" for c, r, n in sorted(missing))
        raise ValueError(f"Anchor residue(s) not found in receptor PDB: {missing_text}")
    if not points:
        raise ValueError("Anchor residues contained no heavy-atom coordinates")

    return tuple(sum(p[i] for p in points) / len(points) for i in range(3))  # type: ignore[return-value]


def strip_end_records(lines: list[str]) -> list[str]:
    return [line for line in lines if not line.startswith(("END", "ENDMDL"))]


def make_rigid_pdbqt(source: Path, destination: Path) -> int:
    """Flatten a Meeko PDBQT torsion tree so all ligand atoms are rigid."""
    atom_lines: list[str] = []
    remarks: list[str] = []
    for raw in source.read_text(encoding="utf-8").splitlines():
        line = raw.rstrip("\n")
        if line.startswith(("ATOM  ", "HETATM")):
            atom_lines.append(line)
        elif line.startswith("REMARK"):
            remarks.append(line)

    if not atom_lines:
        raise ValueError(f"No ATOM/HETATM records found in {source}")

    out = [
        *remarks,
        "REMARK  RIGIDIFIED_BY_VERIFY_HERVK_VINA",
        "ROOT",
        *atom_lines,
        "ENDROOT",
        "TORSDOF 0",
        "END",
    ]
    destination.write_text("\n".join(out) + "\n", encoding="utf-8")
    return len(atom_lines)


def first_model_pdbqt(source: Path, destination: Path) -> None:
    """Extract MODEL 1 from a multi-pose PDBQT output."""
    lines = source.read_text(encoding="utf-8").splitlines()
    has_model = any(line.startswith("MODEL") for line in lines)
    if not has_model:
        destination.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")
        return

    in_first = False
    collected: list[str] = []
    for line in lines:
        if line.startswith("MODEL"):
            if in_first:
                break
            in_first = True
            collected.append(line)
            continue
        if in_first:
            collected.append(line)
            if line.startswith("ENDMDL"):
                break

    if not collected:
        raise ValueError(f"Could not extract first model from {source}")
    destination.write_text("\n".join(collected) + "\n", encoding="utf-8")


def rewrite_ligand_chain(pdb_path: Path, chain_id: str = "X") -> str:
    """Return ligand PDB text with a separate chain ID for visual inspection."""
    out: list[str] = []
    for line in pdb_path.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM  ", "HETATM", "TER")) and len(line) >= 22:
            chars = list(line)
            chars[21] = chain_id
            line = "".join(chars)
        if line.startswith("END"):
            continue
        out.append(line)
    return "\n".join(out) + "\n"


def make_complex_pdb(receptor_pdb: Path, ligand_pdb: Path, destination: Path) -> None:
    receptor_lines = []
    for line in receptor_pdb.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM  ", "HETATM", "ANISOU", "TER", "CONECT")):
            receptor_lines.append(line)
    ligand_text = rewrite_ligand_chain(ligand_pdb)
    ligand_lines = [
        line for line in ligand_text.splitlines()
        if line.startswith(("ATOM  ", "HETATM", "ANISOU", "TER", "CONECT"))
    ]
    destination.write_text(
        "\n".join(receptor_lines + ["TER"] + ligand_lines + ["END"]) + "\n",
        encoding="utf-8",
    )


def make_box_pdb(center: tuple[float, float, float], size: tuple[float, float, float], destination: Path) -> None:
    cx, cy, cz = center
    sx, sy, sz = size
    xmin, xmax = cx - sx / 2, cx + sx / 2
    ymin, ymax = cy - sy / 2, cy + sy / 2
    zmin, zmax = cz - sz / 2, cz + sz / 2
    corners = [
        (xmin, ymin, zmin), (xmax, ymin, zmin),
        (xmax, ymax, zmin), (xmin, ymax, zmin),
        (xmin, ymin, zmax), (xmax, ymin, zmax),
        (xmax, ymax, zmax), (xmin, ymax, zmax),
    ]
    edges = [(1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5),(1,5),(2,6),(3,7),(4,8)]
    lines = []
    for idx, (x, y, z) in enumerate(corners, 1):
        lines.append(
            f"HETATM{idx:5d}  V   BOX X   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          C"
        )
    for a, b in edges:
        lines.append(f"CONECT{a:5d}{b:5d}")
    lines.append("END")
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_affinities(stdout: str, pose_pdbqt: str) -> tuple[float, float | None, list[dict[str, float]]]:
    table_match = AFFINITY_TABLE.search(stdout)
    remark_match = AFFINITY_REMARK.search(pose_pdbqt)
    stdout_aff = float(table_match.group(1)) if table_match else None
    remark_aff = float(remark_match.group(1)) if remark_match else None

    if stdout_aff is None and remark_aff is None:
        raise ValueError("Could not find a Vina affinity in stdout or output PDBQT")

    if stdout_aff is not None and remark_aff is not None:
        if abs(stdout_aff - remark_aff) > 1e-3:
            raise ValueError(
                f"Vina affinity mismatch: stdout={stdout_aff}, PDBQT REMARK={remark_aff}"
            )

    primary = stdout_aff if stdout_aff is not None else remark_aff
    assert primary is not None

    modes = [
        {
            "mode": float(mode),
            "affinity_kcal_mol": float(affinity),
            "rmsd_lb": float(rmsd_lb),
            "rmsd_ub": float(rmsd_ub),
        }
        for mode, affinity, rmsd_lb, rmsd_ub in MODE_TABLE.findall(stdout)
    ]
    return primary, remark_aff, modes


def prepare_receptor(vina_prep: str, receptor_pdb: Path, output_base: Path, center: tuple[float, float, float], size: tuple[float, float, float]) -> Path:
    """Use Meeko receptor preparation. Returns rigid receptor PDBQT path."""
    cmd = [
        vina_prep,
        "-i", str(receptor_pdb),
        "-o", str(output_base),
        "-p", "-v",
        "--box_size", *(f"{x:.3f}" for x in size),
        "--box_center", *(f"{x:.3f}" for x in center),
        "-a",
    ]
    run_checked(cmd)
    candidates = [
        output_base.with_name(output_base.name + ".pdbqt"),
        output_base.with_name(output_base.name + "_rigid.pdbqt"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        f"Meeko receptor preparation completed, but no receptor PDBQT was found near {output_base}"
    )


def prepare_ligand(meeko_script: str, obabel: str, peptide_pdb: Path, output_dir: Path) -> Path:
    sdf = output_dir / "peptide_input.sdf"
    pdbqt = output_dir / "peptide_prepared.pdbqt"

    # Preserve the supplied 3-D coordinates while ensuring explicit hydrogens.
    run_checked([
        obabel,
        str(peptide_pdb),
        "-O", str(sdf),
        "-h",
    ])

    run_checked([
        meeko_script,
        "-i", str(sdf),
        "-o", str(pdbqt),
        "--charge_model", "gasteiger",
        "--remove_smiles",
    ])
    return pdbqt


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify HERV-K SU peptide docking with AutoDock Vina")
    parser.add_argument("--receptor-pdb", default="herv_k_su.pdb", type=Path)
    parser.add_argument("--peptide-pdb", required=True, type=Path)
    parser.add_argument("--sequence", default=None)
    parser.add_argument("--vina", default=None)
    parser.add_argument("--obabel", default=None)
    parser.add_argument("--meeko-ligand", default=None, help="Path to mk_prepare_ligand.py")
    parser.add_argument("--meeko-receptor", default=None, help="Path to mk_prepare_receptor.py")
    parser.add_argument("--box-center", nargs=3, metavar=("X", "Y", "Z"), default=None)
    parser.add_argument("--box-size", nargs=3, metavar=("SX", "SY", "SZ"), required=True)
    parser.add_argument("--anchor-residue", action="append", type=parse_anchor, default=[])
    parser.add_argument("--outdir", type=Path, default=Path("results/vina_verify"))
    parser.add_argument("--exhaustiveness", type=int, default=32)
    parser.add_argument("--num-modes", type=int, default=9)
    parser.add_argument("--cpu", type=int, default=0)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--energy-range", type=float, default=3.0)
    parser.add_argument("--rigid", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--max-evals", type=int, default=None)
    args = parser.parse_args()

    receptor_pdb = args.receptor_pdb.expanduser().resolve()
    peptide_pdb = args.peptide_pdb.expanduser().resolve()
    if not receptor_pdb.exists():
        raise FileNotFoundError(f"Receptor PDB not found: {receptor_pdb}")
    if not peptide_pdb.exists():
        raise FileNotFoundError(f"Peptide PDB not found: {peptide_pdb}")
    if args.exhaustiveness < 1:
        raise ValueError("--exhaustiveness must be >= 1")
    if args.num_modes < 1:
        raise ValueError("--num-modes must be >= 1")
    if args.cpu < 0:
        raise ValueError("--cpu must be >= 0")

    vina = find_executable(args.vina, ["vina"])
    obabel = find_executable(args.obabel, ["obabel", "babel"])
    meeko_ligand = find_executable(args.meeko_ligand, ["mk_prepare_ligand.py"])
    meeko_receptor = find_executable(args.meeko_receptor, ["mk_prepare_receptor.py"])

    size = parse_triplet(args.box_size, "--box-size")
    if any(x <= 0 for x in size):
        raise ValueError("All box dimensions must be > 0")

    if args.box_center is not None:
        center = parse_triplet(args.box_center, "--box-center")
    elif args.anchor_residue:
        center = anchor_center(receptor_pdb, args.anchor_residue)
    else:
        raise ValueError("Provide either --box-center X Y Z or at least one --anchor-residue")

    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Save exact input metadata.
    (outdir / "inputs.json").write_text(
        json.dumps(
            {
                "receptor_pdb": str(receptor_pdb),
                "peptide_pdb": str(peptide_pdb),
                "sequence": args.sequence,
                "box_center_angstrom": center,
                "box_size_angstrom": size,
                "anchor_residues": [f"{c}:{r}:{n}" for c, r, n in args.anchor_residue],
                "exhaustiveness": args.exhaustiveness,
                "num_modes": args.num_modes,
                "cpu": args.cpu,
                "seed": args.seed,
                "energy_range": args.energy_range,
                "rigid": args.rigid,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    make_box_pdb(center, size, outdir / "docking_box.pdb")

    receptor_base = outdir / "receptor_prepared"
    receptor_pdbqt = prepare_receptor(
        meeko_receptor,
        receptor_pdb,
        receptor_base,
        center,
        size,
    )
    # Stable user-facing name.
    stable_receptor = outdir / "receptor.pdbqt"
    shutil.copy2(receptor_pdbqt, stable_receptor)

    prepared_ligand = prepare_ligand(meeko_ligand, obabel, peptide_pdb, outdir)
    ligand_for_docking = prepared_ligand
    if args.rigid:
        rigid_ligand = outdir / "peptide_rigid.pdbqt"
        n_atoms = make_rigid_pdbqt(prepared_ligand, rigid_ligand)
        ligand_for_docking = rigid_ligand
        print(f"Prepared rigid ligand with {n_atoms} atoms: {rigid_ligand}")
    else:
        print(f"Keeping Meeko ligand torsion tree: {prepared_ligand}")

    pose_pdbqt = outdir / "vina_out.pdbqt"
    vina_stdout = outdir / "vina.stdout.txt"
    vina_stderr = outdir / "vina.stderr.txt"
    vina_log = outdir / "vina.log"
    command = [
        vina,
        "--receptor", str(stable_receptor),
        "--ligand", str(ligand_for_docking),
        "--center_x", f"{center[0]:.3f}",
        "--center_y", f"{center[1]:.3f}",
        "--center_z", f"{center[2]:.3f}",
        "--size_x", f"{size[0]:.3f}",
        "--size_y", f"{size[1]:.3f}",
        "--size_z", f"{size[2]:.3f}",
        "--exhaustiveness", str(args.exhaustiveness),
        "--num_modes", str(args.num_modes),
        "--energy_range", f"{args.energy_range:.3f}",
        "--cpu", str(args.cpu),
        "--out", str(pose_pdbqt),
        "--log", str(vina_log),
    ]
    if args.seed is not None:
        command.extend(["--seed", str(args.seed)])
    if args.max_evals is not None:
        command.extend(["--max_evals", str(args.max_evals)])

    (outdir / "vina_command.txt").write_text(
        " ".join(command) + "\n",
        encoding="utf-8",
    )

    result = run_checked(
        command,
        stdout_path=vina_stdout,
        stderr_path=vina_stderr,
    )
    if not pose_pdbqt.exists() or pose_pdbqt.stat().st_size == 0:
        raise FileNotFoundError("Vina exited successfully but produced no output PDBQT")

    pose_text = pose_pdbqt.read_text(encoding="utf-8")
    primary, remark, modes = parse_affinities(result.stdout, pose_text)

    pose1_pdbqt = outdir / "vina_pose1.pdbqt"
    pose1_pdb = outdir / "vina_pose1.pdb"
    complex_pdb = outdir / "docked_complex_mode1.pdb"
    first_model_pdbqt(pose_pdbqt, pose1_pdbqt)
    run_checked([
        obabel,
        "-ipdbqt", str(pose1_pdbqt),
        "-opdb",
        "-O", str(pose1_pdb),
    ])
    make_complex_pdb(receptor_pdb, pose1_pdb, complex_pdb)

    summary = {
        "sequence": args.sequence,
        "receptor": str(receptor_pdb),
        "peptide_structure": str(peptide_pdb),
        "vina_affinity_kcal_mol": primary,
        "vina_affinity_from_pdbqt_remark_kcal_mol": remark,
        "vina_stdout_affinity_kcal_mol": primary if AFFINITY_TABLE.search(result.stdout) else None,
        "modes": modes,
        "exhaustiveness": args.exhaustiveness,
        "num_modes": args.num_modes,
        "cpu": args.cpu,
        "seed": args.seed,
        "box_center_angstrom": center,
        "box_size_angstrom": size,
        "rigid_ligand": args.rigid,
        "receptor_pdbqt": str(stable_receptor),
        "ligand_pdbqt": str(ligand_for_docking),
        "pose_pdbqt": str(pose_pdbqt),
        "complex_pdb": str(complex_pdb),
    }
    (outdir / "docking_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    print("\n=== Vina verification complete ===")
    print(f"Affinity (kcal/mol): {primary:.3f}")
    if remark is not None:
        print(f"PDBQT REMARK cross-check: {remark:.3f} kcal/mol")
    print(f"Exhaustiveness: {args.exhaustiveness}")
    print(f"Box center: {center}")
    print(f"Box size (A): {size}")
    print(f"Docked complex: {complex_pdb}")
    print(f"Raw Vina stdout: {vina_stdout}")
    print(f"Raw Vina log: {vina_log}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
