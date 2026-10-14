"""Compare ARG203/ARG205 coordinates across receptor representations."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

try:
    from .target_audit import coordinate_frame_audit
except ImportError:  # direct script execution
    from target_audit import coordinate_frame_audit


ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "data" / "receptor" / "9SYA.pdb"
DEFAULT_RECEPTORS = [
    ROOT / "data" / "receptor" / "9SYA.pdb",
    ROOT / "prototype2" / "9SYA_chainA_final.pdb",
    ROOT / "prototype2" / "target_aware_audit" / "9SYA_chainA_repaired.pdb",
    ROOT / "prototype2" / "p2_rapidock" / "P2-01-clean" / "receptor_for_rapidock_full.pdb",
    ROOT / "prototype2" / "p2_rapidock" / "P2-01-clean" / "receptor_for_rapidock.pdb",
]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, default=REFERENCE)
    parser.add_argument("--receptor", action="append", type=Path)
    parser.add_argument("--site", nargs=3, type=float, default=(-8.405, -5.823, 1.523))
    parser.add_argument("--output", type=Path, default=ROOT / "prototype2" / "target_aware_audit" / "coordinate_frame_audit.json")
    args = parser.parse_args()
    receptors = [path for path in (args.receptor or DEFAULT_RECEPTORS) if path.exists()]
    result = coordinate_frame_audit(
        receptors,
        reference_path=args.reference,
        site_center=args.site,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "all_mappings_valid": result["all_mappings_valid"],
        "all_coordinate_frames_comparable": result["all_coordinate_frames_comparable"],
        "output": str(args.output),
    }, indent=2))
    return 0 if result["all_mappings_valid"] and result["all_coordinate_frames_comparable"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
