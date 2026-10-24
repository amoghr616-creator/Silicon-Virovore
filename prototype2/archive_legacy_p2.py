#!/usr/bin/env python3
"""Reversibly archive prior Prototype 2 experiment outputs, preserving final evidence."""
from __future__ import annotations
import argparse, json, shutil
from datetime import datetime
from pathlib import Path

KEEP_DIRS = {
    "p2_six_method_smoke_20261021",
    "p2_uniform_random_reset_vs_v2_100_20261023",
    "structural_validation",
}
ARCHIVE_PREFIXES = ("p2_", "substitution_", "patch_backup_", "validation_checkpoint_")
ARCHIVE_EXACT = {"benchmark_1f4v", "debug_receptor"}
TARGET_AUDIT_OLD_RESULTS = {
    "P2-01_geometry.json",
    "P2-01_relaxed_geometry.json",
    "p2_rapidock_target_rank_report.json",
    "p2_target_rank_report.json",
    "p2_v2_target_rank_report.json",
}
BENCHMARK_1YCR_OLD_OUTPUTS = {"docking", "rapidock_20", "vina_positive_control_20261016"}
DIRECT_OLD_RESULT_PREFIXES = ("p2_", "adaptive_gen10_pool", "baseline_gen10_pool", "pilot_runs_preflight", "readiness_", "target_aware_audit_preflight")

def plan(root: Path) -> list[Path]:
    proto = root / "prototype2"
    if not proto.is_dir():
        raise SystemExit(f"ERROR: no prototype2 directory under {root}")
    items: list[Path] = []
    for p in sorted(proto.iterdir()):
        if p.is_dir():
            if p.name in KEEP_DIRS:
                continue
            if p.name.startswith(ARCHIVE_PREFIXES) or p.name in ARCHIVE_EXACT:
                items.append(p)
        elif p.is_file():
            # Archive generated result/data files, but preserve .py source scripts and target PDB assets.
            if p.suffix.lower() != ".py" and (p.name.startswith(DIRECT_OLD_RESULT_PREFIXES) or p.name == "P2_PILOT_RUNBOOK.md"):
                items.append(p)
    # Only archive previous generated results under the 1YCR benchmark; keep 1YCR.pdb and code.
    ycr = proto / "benchmark_1ycr"
    if ycr.is_dir():
        for p in sorted(ycr.iterdir()):
            if p.is_dir() and p.name in BENCHMARK_1YCR_OLD_OUTPUTS:
                items.append(p)
    audit = proto / "target_aware_audit"
    if audit.is_dir():
        for p in sorted(audit.iterdir()):
            if p.is_file() and p.name in TARGET_AUDIT_OLD_RESULTS:
                items.append(p)
    # Deduplicate while retaining order.
    return list(dict.fromkeys(items))

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=Path.home() / "Desktop" / "Silicon-Virovore")
    ap.add_argument("--apply", action="store_true", help="perform moves; without this, only preview")
    args = ap.parse_args()
    root = args.root.expanduser().resolve()
    items = plan(root)
    print(f"Repo: {root}")
    print(f"Items to archive: {len(items)}")
    for p in items:
        print(("MOVE " if args.apply else "WOULD MOVE ") + str(p.relative_to(root)))
    print("PRESERVE prototype2/p2_six_method_smoke_20261021")
    print("PRESERVE prototype2/p2_uniform_random_reset_vs_v2_100_20261023 (final 100-seed evidence)")
    print("PRESERVE prototype2/target_aware_audit receptor/mapping/preparation assets")
    print("PRESERVE prototype2/benchmark_1ycr/1YCR.pdb and benchmark source scripts")
    if not args.apply:
        print("DRY RUN ONLY. Re-run with --apply to move these items to a sibling archive.")
        return
    archive = root.parent / f"{root.name}_P2_legacy_archive_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    archive.mkdir(parents=True, exist_ok=False)
    moved = []
    for p in items:
        dest = archive / p.relative_to(root)
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists():
            raise SystemExit(f"ABORT: archive destination already exists: {dest}")
        shutil.move(str(p), str(dest))
        moved.append({"source_relative_to_repo": str(p.relative_to(root)), "archive_relative_path": str(dest.relative_to(archive))})
    manifest = {"repo": str(root), "archive": str(archive), "moves": moved, "preserved": sorted(KEEP_DIRS), "operation": "reversible move; no deletion"}
    (archive / "archive_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"\nARCHIVED {len(moved)} items to:\n{archive}\nManifest: {archive / 'archive_manifest.json'}")

if __name__ == "__main__":
    main()
