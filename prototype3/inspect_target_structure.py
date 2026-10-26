#!/usr/bin/env python3
"""Inspect author chain/residue numbering before applying hotspot identifiers."""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.target_aware_scoring import inspect_complex_structure


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", required=True, help="Local PDB containing the intended receptor")
    parser.add_argument("--chain", help="Optional chain ID to show in detail")
    args = parser.parse_args()
    report = inspect_complex_structure(args.pdb)
    if args.chain:
        if args.chain not in report["chains"]:
            parser.error(f"chain {args.chain!r} not found; available={list(report['chains'])}")
        report["chains"] = {args.chain: report["chains"][args.chain]}
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
