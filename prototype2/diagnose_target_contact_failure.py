"""Audit archived Prototype 2 poses before target-aware optimization.

This script never launches a docking backend.  It evaluates every parseable
pose file in the supplied result roots against the independently mapped
ARG203/ARG205 site and writes raw measurements plus failure classifications.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Iterable

try:
    from .target_audit import audit_residue_mapping, coordinate_frame_audit, diagnose_pose
except ImportError:  # direct script execution
    from target_audit import audit_residue_mapping, coordinate_frame_audit, diagnose_pose


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ROOTS = (
    ROOT / "prototype2" / "p2_rapidock",
    ROOT / "prototype2" / "p2_heldout",
    ROOT / "prototype2" / "p2_v2_rapidock",
)
DEFAULT_REFERENCE = ROOT / "data" / "receptor" / "9SYA.pdb"
FAILURE_CATEGORIES = (
    "wrong_mapping",
    "coordinate_mismatch",
    "search_region_miss",
    "malformed_pose",
    "sequence_mismatch",
    "ambiguous_chain",
    "severe_clash",
    "missing_output",
    "geometrically_valid_but_misses_target_contacts",
)


def _metadata(path: Path) -> dict:
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}


def _pose_files(candidate_dir: Path) -> list[tuple[str, Path]]:
    """Collect each existing output representation once."""

    candidates: list[tuple[str, Path]] = []
    for directory_name in ("poses", "poses_minimized", "poses_scored"):
        directory = candidate_dir / directory_name
        if directory.exists():
            candidates.extend((directory_name, path) for path in sorted(directory.glob("*.pdb")))
    for name in ("best_pose.pdb", "best_pose_vina_relaxed.pdb"):
        path = candidate_dir / name
        if path.exists():
            candidates.append((name.removesuffix(".pdb"), path))
    seen: set[Path] = set()
    result = []
    for kind, path in candidates:
        resolved = path.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        result.append((kind, path))
    return result


def _ranked_scores(candidate_dir: Path) -> dict[str, float]:
    ranked = candidate_dir / "ranked_poses.csv"
    if not ranked.exists():
        return {}
    scores: dict[str, float] = {}
    try:
        with ranked.open(newline="") as handle:
            for row in csv.DictReader(handle):
                try:
                    scores[row["pose_filename"]] = float(row["delta_g"])
                except (KeyError, TypeError, ValueError):
                    continue
    except OSError:
        return {}
    return scores


def _candidate_dirs(roots: Iterable[Path]) -> list[tuple[Path, Path]]:
    return [
        (root, candidate)
        for root in roots
        if root.exists()
        for candidate in sorted(root.iterdir())
        if candidate.is_dir()
    ]


def diagnose_candidate(
    root: Path,
    candidate_dir: Path,
    *,
    reference: Path,
) -> dict[str, object]:
    metadata = _metadata(candidate_dir / "run_metadata.json")
    cli_args = metadata.get("cli_args", {})
    sequence = cli_args.get("peptide_sequence") or metadata.get("peptide_sequence") or metadata.get("sequence")
    receptor = candidate_dir / "receptor_for_rapidock_full.pdb"
    if not receptor.exists():
        receptor = candidate_dir / "receptor_for_rapidock.pdb"
    if not receptor.exists():
        receptor = candidate_dir / "receptor.pdb"
    site_center = cli_args.get("site_coords")
    box_size = cli_args.get("box_size")
    base = {
        "result_root": str(root),
        "candidate": candidate_dir.name,
        "candidate_dir": str(candidate_dir),
        "sequence": sequence,
        "receptor_path": str(receptor),
        "docking_status": metadata.get("status"),
        "site_center_A": site_center,
        "declared_box_size_A": box_size,
        "metadata_provenance": {
            "rapidock_commit_sha": metadata.get("rapidock_commit_sha"),
            "pipeline_git_sha": metadata.get("git_sha"),
            "effective_box_size_A": metadata.get("effective_box_size"),
            "long_checkpoint_threshold": cli_args.get("long_checkpoint_threshold"),
            "blind": cli_args.get("blind"),
            "pocket_search": cli_args.get("pocket_search"),
        },
    }
    pose_files = _pose_files(candidate_dir)
    base["pose_file_count"] = len(pose_files)
    if not sequence or not receptor.exists() or not pose_files:
        base["failure_category"] = "missing_output"
        base["poses"] = []
        return base

    mapping = audit_residue_mapping(receptor, reference_path=reference, chain="A")
    try:
        frame = coordinate_frame_audit(
            [receptor], reference_path=reference, chain="A", site_center=site_center,
        )["representations"][0]
    except (OSError, ValueError) as exc:
        frame = {"coordinate_frame_comparable": False, "coordinate_frame_error": str(exc)}
    base["receptor_mapping"] = mapping
    base["coordinate_frame"] = frame
    scores = _ranked_scores(candidate_dir)
    pose_rows = []
    for pose_kind, pose_path in pose_files:
        row: dict[str, object] = {
            "pose_kind": pose_kind,
            "pose_path": str(pose_path),
            "docking_score_kcal_mol": scores.get(pose_path.name),
        }
        try:
            diagnostic = diagnose_pose(
                receptor,
                pose_path,
                sequence,
                reference_path=reference,
                site_center=site_center,
                box_size_A=box_size,
                mapping_result=mapping,
                coordinate_frame_result=frame,
            )
            # These are constant for every pose from one candidate/receptor.
            # Keep the complete records at candidate level and retain only
            # the per-pose comparability result here.
            diagnostic.pop("target_mapping", None)
            diagnostic.pop("coordinate_frame", None)
            diagnostic["coordinate_frame_comparable"] = frame.get("coordinate_frame_comparable", False)
            row["diagnostic"] = diagnostic
        except (OSError, ValueError, RuntimeError) as exc:
            row["diagnostic"] = {
                "failure_categories": ["malformed_pose"],
                "primary_failure_category": "malformed_pose",
                "diagnostic_error": str(exc),
            }
        pose_rows.append(row)
    base["poses"] = pose_rows
    categories = Counter(
        category
        for row in pose_rows
        for category in row["diagnostic"].get("failure_categories", [])
    )
    base["failure_category_counts"] = dict(sorted(categories.items()))
    base["parseable_pose_count"] = sum(
        1 for row in pose_rows if "diagnostic_error" not in row["diagnostic"]
    )
    return base


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", action="append", type=Path, dest="roots")
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    roots = tuple(args.roots or DEFAULT_ROOTS)
    candidates = _candidate_dirs(roots)
    candidate_rows = [diagnose_candidate(root, candidate, reference=args.reference) for root, candidate in candidates]
    all_pose_rows = [row for candidate in candidate_rows for row in candidate["poses"]]
    category_counts = Counter(
        category
        for row in all_pose_rows
        for category in row["diagnostic"].get("failure_categories", [])
    )
    primary_counts = Counter(
        row["diagnostic"].get("primary_failure_category")
        for row in all_pose_rows
        if row["diagnostic"].get("primary_failure_category")
    )
    report = {
        "purpose": "diagnose_target_contact_failure_before_target_aware_optimization",
        "docking_launched": False,
        "favorable_docking_score_is_binding_proof": False,
        "reference_path": str(args.reference),
        "results_roots": [str(root) for root in roots],
        "target_site": {
            "chain": "A",
            "author_residues": [203, 205],
            "contact_cutoffs_A": [3.0, 4.0, 5.0, 8.0],
            "severe_clash_cutoff_A": 2.0,
            "structural_contact_cutoff_A": 4.5,
        },
        "candidate_count": len(candidate_rows),
        "pose_count": len(all_pose_rows),
        "parseable_pose_count": sum(candidate["parseable_pose_count"] for candidate in candidate_rows),
        "failure_categories": list(FAILURE_CATEGORIES),
        "failure_category_counts": {
            category: category_counts.get(category, 0)
            for category in FAILURE_CATEGORIES
        },
        "primary_failure_category_counts": dict(sorted(primary_counts.items())),
        "candidates": candidate_rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "candidate_count": report["candidate_count"],
        "pose_count": report["pose_count"],
        "parseable_pose_count": report["parseable_pose_count"],
        "failure_category_counts": report["failure_category_counts"],
        "output": str(args.output),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
