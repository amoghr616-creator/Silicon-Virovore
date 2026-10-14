"""Prepare a separate, auditable 9SYA receptor copy.

The previous implementation deleted GLN A:109. This wrapper performs a
conservative reference-based repair instead: missing heavy atoms are copied
from the matching deposited conformer, while residues, existing coordinates,
author numbering, and original files are preserved.
"""

from pathlib import Path
import json

try:
    from .target_audit import repair_from_reference
except ImportError:  # direct ``python prototype2/clean_9sya_receptor.py``
    from target_audit import repair_from_reference


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "prototype2" / "9SYA_chainA_final.pdb"
REFERENCE = ROOT / "data" / "receptor" / "9SYA.pdb"
OUTPUT_DIR = ROOT / "prototype2" / "target_aware_audit"
OUTPUT = OUTPUT_DIR / "9SYA_chainA_repaired.pdb"
REPORT = OUTPUT_DIR / "receptor_preparation_report.json"


def prepare_receptor(
    input_path: Path = INPUT,
    reference_path: Path = REFERENCE,
    output_path: Path = OUTPUT,
    report_path: Path = REPORT,
) -> Path:
    result = repair_from_reference(input_path, reference_path, output_path, chain="A")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"Prepared copy: {output_path}")
    print(f"Preparation report: {report_path}")
    print(f"Status: {result['preparation_status']} | added heavy atoms: {len(result['added_atoms'])}")
    return output_path


if __name__ == "__main__":
    prepare_receptor()
