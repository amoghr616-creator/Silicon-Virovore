"""Check a relaxed P2 pose with the same independent geometry criteria."""

from pathlib import Path
import argparse
import json

try:
    from .target_audit import evaluate_pose
except ImportError:  # direct script execution
    from target_audit import evaluate_pose


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "prototype2" / "p2_rapidock" / "P2-01-clean"
REFERENCE = ROOT / "data" / "receptor" / "9SYA.pdb"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receptor", type=Path, default=BASE / "receptor_for_rapidock_full.pdb")
    parser.add_argument("--pose", type=Path, default=BASE / "best_pose_vina_relaxed.pdb")
    parser.add_argument("--sequence", default="MKQAVNALLVFFAGSSDAIRR")
    parser.add_argument("--reference", type=Path, default=REFERENCE)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = evaluate_pose(args.receptor, args.pose, args.sequence, reference_path=args.reference)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print("=== P2 RELAXED POSE INDEPENDENT GEOMETRY ===")
    print("Peptide heavy atoms:", result["peptide_atom_count_heavy"])
    print("Target contact atom count:", result["target_contact_atom_count"])
    print("Minimum receptor-peptide heavy-atom distance:", result["minimum_receptor_peptide_heavy_atom_distance_A"], "A")
    print("Severe clash count:", result["severe_clash_count"])
    print("Geometry status:", result["pose_geometry_status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
