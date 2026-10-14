"""Write a non-destructive, standard-library 9SYA pocket crop."""

from pathlib import Path
import argparse
import math

try:
    from .target_audit import AA3_TO_1, audit_residue_mapping, parse_pdb, residues_from_atoms
except ImportError:  # direct script execution
    from target_audit import AA3_TO_1, audit_residue_mapping, parse_pdb, residues_from_atoms


ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "data" / "receptor" / "9SYA.pdb"
INPUT = ROOT / "prototype2" / "target_aware_audit" / "9SYA_chainA_repaired.pdb"
OUTPUT = ROOT / "prototype2" / "target_aware_audit" / "9SYA_pocket15.pdb"
SITE = (-8.405, -5.823, 1.523)
RADIUS = 15.0


def make_pocket(input_path: Path = INPUT, output_path: Path = OUTPUT) -> Path:
    mapping = audit_residue_mapping(input_path, reference_path=REFERENCE, chain="A")
    if not mapping["mapping_valid"]:
        raise ValueError(f"Refusing to crop receptor with unresolved target mapping: {mapping['failure_reasons']}")
    atoms, lines = parse_pdb(input_path)
    residues = residues_from_atoms(atoms, "A")
    retained = set()
    for residue in residues:
        if residue.resname not in AA3_TO_1:
            continue
        if any(math.dist((atom.x, atom.y, atom.z), SITE) <= RADIUS for atom in residue.atoms):
            retained.add(residue.key)
    kept_line_numbers = {
        atom.line_number
        for atom in atoms
        if atom.residue_key in retained
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        for line_number, line in enumerate(lines, 1):
            if line.startswith(("ATOM  ", "HETATM")):
                if line_number not in kept_line_numbers:
                    continue
            handle.write(line)
    print(f"Saved: {output_path}")
    print(f"Pocket radius: {RADIUS:.1f} A")
    print(f"Residues retained: {len(retained)}")
    print("Residues:", [(key[1], next(r.resname for r in residues if r.key == key)) for key in sorted(retained)])
    return output_path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    make_pocket(args.input, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
