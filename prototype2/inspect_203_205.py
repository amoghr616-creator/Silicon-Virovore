"""Run the fail-closed ARG203/ARG205 mapping audit."""

from pathlib import Path
import argparse
import json

try:
    from .target_audit import compare_receptor_mappings
except ImportError:  # direct script execution
    from target_audit import compare_receptor_mappings


ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "data" / "receptor" / "9SYA.pdb"
DEFAULT_RECEPTORS = [
    ROOT / "data" / "receptor" / "9SYA.pdb",
    ROOT / "prototype2" / "9SYA_chainA.pdb",
    ROOT / "prototype2" / "9SYA_chainA_final.pdb",
    ROOT / "prototype2" / "9SYA_chainA_meeko_clean.pdb",
    ROOT / "prototype2" / "9SYA_chainA_docking.pdb",
    ROOT / "prototype2" / "target_aware_audit" / "9SYA_chainA_repaired.pdb",
    ROOT / "prototype2" / "p2_rapidock" / "P2-01-clean" / "receptor_for_rapidock_full.pdb",
]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receptor", action="append", type=Path)
    parser.add_argument("--output", type=Path, default=ROOT / "prototype2" / "target_aware_audit" / "residue_mapping_audit.json")
    args = parser.parse_args()
    receptors = args.receptor or DEFAULT_RECEPTORS
    result = compare_receptor_mappings(receptors, reference_path=REFERENCE, chain="A")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "all_representations_valid": result["all_representations_valid"],
        "output": str(args.output),
    }, indent=2))
    for item in result["representations"]:
        print(f"{Path(item['structure_path']).name}: mapping_valid={item['mapping_valid']}")
        for target in item["target_residues"]:
            print(
                f"  author ARG{target['requested_author_residue_number']} -> "
                f"{target.get('residue_name', target.get('reference_residue_name'))} "
                f"{target.get('pdb_author_residue_number', 'missing')}"
            )
    return 0 if result["all_representations_valid"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
