"""Evidence-gated ranking for Prototype 2 docking outputs.

This is an external analysis tool, not a replacement for the main sequence
prior in ``src/ranking.py``. It refuses to produce a target-aware ranking if
the receptor mapping, peptide coordinates, or independent geometry checks are
not valid. A Vina/HybriDock score is retained as a descriptive field only;
it is never treated as proof of binding.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

try:
    from .target_audit import TargetMappingError, PoseValidationError, diagnose_pose
except ImportError:  # direct script execution
    from target_audit import TargetMappingError, PoseValidationError, diagnose_pose


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ROOT = ROOT / "prototype2" / "p2_rapidock"
DEFAULT_REFERENCE = ROOT / "data" / "receptor" / "9SYA.pdb"
DEFAULT_SEQUENCE_METRICS = ROOT / "prototype2" / "p2_sequence_metrics.tsv"


def _metadata(path: Path) -> dict:
    try:
        return json.loads(path.read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        return {}


def _best_score(path: Path, pose_name: str) -> float | None:
    ranked = path / "ranked_poses.csv"
    if not ranked.exists():
        return None
    with ranked.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        if row.get("pose_filename") == pose_name:
            try:
                return float(row["delta_g"])
            except (KeyError, TypeError, ValueError):
                return None
    return None


def _sequence_baseline(path: Path) -> dict[str, float]:
    if not path.exists():
        return {}
    with path.open() as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        result = {}
        for row in rows:
            try:
                result[row["candidate"]] = float(row["old_best_score"])
            except (KeyError, TypeError, ValueError):
                continue
        return result


def audit_candidate(candidate_dir: Path, reference: Path, sequence: str | None) -> dict:
    metadata = _metadata(candidate_dir / "run_metadata.json")
    cli_args = metadata.get("cli_args", {})
    sequence = sequence or cli_args.get("peptide_sequence")
    # Prefer the explicitly emitted Vina-relaxed pose when present.  Keep the
    # selected representation in the report; never mix raw and relaxed scores
    # without saying which coordinates were audited.
    pose = candidate_dir / "best_pose_vina_relaxed.pdb"
    pose_variant = "best_pose_vina_relaxed"
    if not pose.exists():
        pose = candidate_dir / "best_pose.pdb"
        pose_variant = "best_pose"
    receptor = candidate_dir / "receptor_for_rapidock_full.pdb"
    if not receptor.exists():
        receptor = candidate_dir / "receptor_for_rapidock.pdb"
    row = {
        "candidate": candidate_dir.name,
        "sequence": sequence,
        "pose_path": str(pose),
        "pose_variant": pose_variant,
        "receptor_path": str(receptor),
        "docking_score_kcal_mol": _best_score(candidate_dir, "best_pose.pdb"),
        "docking_status": metadata.get("status"),
        "mapping_valid": False,
        "geometry_valid": False,
        "target_aware_evidence_eligible": False,
    }
    if not sequence:
        row["failure_reason"] = "peptide sequence is missing from run metadata and arguments"
        return row
    if not pose.exists() or not receptor.exists():
        row["failure_reason"] = "separate receptor and best_pose coordinate files are required"
        return row
    try:
        result = diagnose_pose(
            receptor,
            pose,
            sequence,
            reference_path=reference,
            site_center=cli_args.get("site_coords"),
            box_size_A=cli_args.get("box_size"),
        )
    except (TargetMappingError, PoseValidationError) as exc:
        row["failure_reason"] = str(exc)
        return row
    row.update({
        "mapping_valid": result.get("mapping_valid", False),
        "geometry_valid": result.get("pose_geometry_status") == "valid",
        "pose_coordinates_valid": result.get("coordinates_valid", False),
        "peptide_sequence_valid": result.get("sequence_valid", False),
        "peptide_chain": result.get("peptide_chain"),
        "target_contact_atom_count": result.get("target_contact_atom_count", 0),
        "minimum_receptor_peptide_heavy_atom_distance_A": result.get("minimum_receptor_peptide_heavy_atom_distance_A"),
        "severe_clash_count": result.get("severe_clash_count"),
        "target_metrics": result.get("target_metrics"),
        "search_volume": result.get("search_volume"),
        "failure_categories": result.get("failure_categories", []),
        "coordinate_frame": result.get("coordinate_frame"),
        "target_contact_evidence_status": result.get("target_contact_evidence_status"),
        "target_aware_evidence_eligible": (
            result.get("mapping_valid", False)
            and result.get("coordinates_valid", False)
            and result.get("sequence_valid", False)
            and row["docking_status"] == "complete"
            and result.get("pose_geometry_status") == "valid"
            and result.get("target_contact_atom_count", 0) > 0
            and "search_region_miss" not in result.get("failure_categories", [])
            and "coordinate_mismatch" not in result.get("failure_categories", [])
            and not any(
                category in result.get("failure_categories", [])
                for category in ("malformed_pose", "sequence_mismatch", "ambiguous_chain")
            )
        ),
    })
    return row


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--sequence", help="Use one sequence for a directory without run metadata")
    parser.add_argument("--sequence-metrics", type=Path, default=DEFAULT_SEQUENCE_METRICS)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--allow-target-aware", action="store_true", help="Explicitly enable evidence-gated target-aware ranking")
    args = parser.parse_args()
    candidates = [path for path in sorted(args.results_root.iterdir()) if path.is_dir()]
    rows = [audit_candidate(path, args.reference, args.sequence) for path in candidates]
    baseline_scores = _sequence_baseline(args.sequence_metrics)
    for row in rows:
        row["sequence_only_baseline_score"] = baseline_scores.get(row["candidate"])
    eligible = [row for row in rows if row["target_aware_evidence_eligible"]]
    report = {
        "target_aware_ranking_requested": args.allow_target_aware,
        "target_aware_ranking_performed": False,
        "target_aware_ranking_status": "disabled_without_explicit_flag",
        "ranking_basis": "independent target contacts + clash-free geometry; docking score descriptive only",
        "favorable_docking_score_is_binding_proof": False,
        "candidate_count": len(rows),
        "eligible_candidate_count": len(eligible),
        "candidates": rows,
        "ranking": [],
        "sequence_only_baseline": sorted(
            [row for row in rows if row["sequence_only_baseline_score"] is not None],
            key=lambda row: (-row["sequence_only_baseline_score"], row["candidate"]),
        ),
        "sequence_only_baseline_source": str(args.sequence_metrics),
    }
    if args.allow_target_aware:
        failed_mapping = [row for row in rows if not row["mapping_valid"]]
        if failed_mapping:
            report["target_aware_ranking_status"] = "stopped_mapping_or_pose_failure"
        elif not eligible:
            report["target_aware_ranking_status"] = "stopped_no_eligible_structural_evidence"
        else:
            # Contact count is a geometry descriptor, not a binding affinity.
            report["target_aware_ranking_performed"] = True
            report["target_aware_ranking_status"] = "complete_evidence_gated"
            report["ranking"] = sorted(
                eligible,
                key=lambda row: (
                    -row["target_contact_atom_count"],
                    row["minimum_receptor_peptide_heavy_atom_distance_A"],
                    row["docking_score_kcal_mol"] if row["docking_score_kcal_mol"] is not None else float("inf"),
                    row["candidate"],
                ),
            )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": report["target_aware_ranking_status"],
        "eligible_candidate_count": report["eligible_candidate_count"],
        "output": str(args.output),
    }, indent=2))
    return 0 if report["target_aware_ranking_status"] in {"disabled_without_explicit_flag", "complete_evidence_gated"} else 2


if __name__ == "__main__":
    raise SystemExit(main())
