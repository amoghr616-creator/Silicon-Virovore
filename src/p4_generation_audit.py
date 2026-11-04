"""Consolidated, evidence-gated P4 generation report."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

SCHEMA_VERSION = "p4-generation-assessment-v1"


def finite_number(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def file_hash(path: Path) -> str | None:
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def epistasis_provenance() -> dict[str, Any]:
    configured = os.environ.get("SILICON_VIROVORE_P4_EPISTASIS_REPORT")
    if not configured:
        return {
            "status": "separate_benchmark_not_attached",
            "candidate_level_epistasis_used": False,
        }

    path = Path(configured).expanduser()
    if not path.is_file():
        return {
            "status": "configured_report_missing",
            "path": str(path),
            "candidate_level_epistasis_used": False,
        }

    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return {
            "status": "invalid_external_report",
            "path": str(path),
            "error": f"{type(exc).__name__}: {exc}",
            "candidate_level_epistasis_used": False,
        }

    return {
        "status": "external_benchmark_attached",
        "path": str(path),
        "sha256": file_hash(path),
        "report_type": payload.get("report_type")
            if isinstance(payload, dict) else None,
        "candidate_level_epistasis_used": False,
        "interpretation": (
            "Benchmark provenance only; not candidate-specific "
            "HERV-K binding evidence."
        ),
    }


def write_p4_generation_assessment(
    candidates: list[Any],
    target_report: dict[str, Any],
    output_directory: str | Path,
) -> dict[str, Any]:
    if not isinstance(target_report, dict):
        raise TypeError("target_report must be a dictionary")

    generation = int(target_report["generation"])
    rows = []

    for index, candidate in enumerate(candidates):
        metadata = getattr(candidate, "metadata", None)
        if not isinstance(metadata, dict):
            raise TypeError(f"Candidate {index} has invalid metadata")

        target = metadata.get("p4_target_evidence") or {}
        reliability = metadata.get("p4_reliability") or {}
        epistasis = metadata.get("p4_epistasis_evidence") or {}
        ood = metadata.get("p4_ood_decision") or {}
        gate = metadata.get("p4_integrated_gate") or {}
        epistasis_operator = metadata.get("p4_epistasis_operator") or {}

        target_status = target.get("status", "missing")
        calibrated = reliability.get("calibrated") is True
        feedback = reliability.get("feedback_active") is True

        if target_status != "success":
            decision = "ABSTAIN_TARGET_EVIDENCE_INVALID"
        elif target_report.get("pareto_reordered_candidates") is True:
            decision = "TARGET_PARETO_FRONT_ASSIGNED"
        else:
            decision = "TARGET_ASSESSED_LEGACY_ORDER_RETAINED"

        if gate.get("decision"):
            decision = gate["decision"]

        rows.append({
            "candidate_index": index,
            "candidate_id": metadata.get(
                "candidate_id", f"generation-{generation}-candidate-{index}"
            ),
            "rank": getattr(candidate, "rank", None),
            "target_status": target_status,
            "target_eligible": target_status == "success",
            "target_geometry_score": finite_number(
                target.get("target_geometry_score_mean")
            ),
            "hotspot_proximity_score": finite_number(
                target.get("hotspot_proximity_score")
            ),
            "pose_source_groups": target.get("pose_source_groups"),
            "reliability_status": (
                "externally_calibrated" if calibrated
                else "diagnostic_only"
                if reliability.get("status") == "diagnostic_only"
                else "uncalibrated_or_unavailable"
            ),
            "reliability_calibrated": calibrated,
            "reliability_feedback_active": feedback,
            "epistasis_status": epistasis.get(
                "status", "no_candidate_specific_evidence"
            ),
            "ood_status": ood.get("status", "not_evaluated"),
            "candidate_decision": decision,
            "integrated_gate_decision": gate.get(
                "decision", "NOT_EVALUATED"
            ),
            "integrated_gate_reason": gate.get("reason"),
            "target_gate_valid": bool(
                gate.get("target_optimization", {}).get("valid")
            ),
            "reliability_gate_valid": bool(
                gate.get("reliability", {}).get("valid")
            ),
            "epistasis_gate_valid": bool(
                gate.get("epistasis", {}).get("valid")
            ),
            "epistasis_operator_active": bool(
                epistasis_operator.get("active")
            ),
            "epistasis_operator_mode": epistasis_operator.get("mode"),
        })

    count = len(rows)
    target_count = sum(row["target_eligible"] for row in rows)
    calibrated_count = sum(row["reliability_calibrated"] for row in rows)
    feedback_count = sum(row["reliability_feedback_active"] for row in rows)

    evaluator_path = Path(__file__).with_name("p4_pipeline_integration.py")

    report = {
        "schema_version": SCHEMA_VERSION,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "generation": generation,
        "run_id": os.environ.get("SILICON_VIROVORE_RUN_ID"),
        "candidate_count": count,
        "target_evidence_count": target_count,
        "target_evidence_coverage": target_count / count if count else 0.0,
        "target_manifest": target_report.get("manifest"),
        "target_selection": {
            "pareto_requested": target_report.get("pareto_requested", False),
            "pareto_reordered_candidates": target_report.get(
                "pareto_reordered_candidates", False
            ),
            "reason": target_report.get("reason"),
            "legacy_fitness_scores_modified": target_report.get(
                "legacy_fitness_scores_modified", False
            ),
        },
        "reliability": {
            "diagnostic_candidate_count": target_report.get(
                "reliability_diagnostic_candidate_count", 0
            ),
            "calibrated_candidate_count": calibrated_count,
            "feedback_active_candidate_count": feedback_count,
            "calibrated": count > 0 and calibrated_count == count,
            "feedback_active": feedback_count > 0,
            "interpretation": (
                "Repeatability diagnostics are not structural-accuracy "
                "calibration."
            ),
        },
        "integrated_gate": {
            "status": "candidate_level_adjudication_executed",
            "eligible_count": sum(
                row["integrated_gate_decision"]
                == "ELIGIBLE_FOR_MULTI_OBJECTIVE_REVIEW"
                for row in rows
            ),
            "abstain_count": sum(
                row["integrated_gate_decision"] == "ABSTAIN"
                for row in rows
            ),
            "target_valid_count": sum(row["target_gate_valid"] for row in rows),
            "reliability_valid_count": sum(
                row["reliability_gate_valid"] for row in rows
            ),
            "epistasis_valid_count": sum(
                row["epistasis_gate_valid"] for row in rows
            ),
            "scalar_fitness_generated": False,
            "abstention_blocks_multi_objective_review_not_legacy_search": True,
        },
        "epistasis_live_loop": {
            "active_candidate_count": sum(
                row["epistasis_operator_active"] for row in rows
            ),
            "modes": sorted({
                row["epistasis_operator_mode"]
                for row in rows
                if row["epistasis_operator_mode"]
            }),
            "candidate_level_delta_epsilon_validated_count": sum(
                row["epistasis_gate_valid"] for row in rows
            ),
            "interpretation": (
                "Pairwise epistasis/entropy modes guide mutation proposals "
                "and generation feedback. They are not candidate-level "
                "experimental epistasis validation."
            ),
        },
        "epistasis_benchmark": epistasis_provenance(),
        "ood_assessment": {
            "status": (
                "candidate_metadata_available"
                if any(row["ood_status"] != "not_evaluated" for row in rows)
                else "not_evaluated"
            ),
            "candidate_level_decisions_recorded": sum(
                row["ood_status"] != "not_evaluated" for row in rows
            ),
        },
        "evsi_and_pruning": {
            "status": "ABSTAINED_NO_EXPERIMENTAL_EVIDENCE",
            "pruning_allowed": False,
            "reason": (
                "Generation-level computational scores do not satisfy "
                "the independent experimental evidence requirements "
                "of candidate_prune_decision."
            ),
        },
        "optimizer_feedback": {
            "reliability_feedback_active": feedback_count > 0,
            "experimental_pruning_performed": False,
            "legacy_fitness_scores_modified": target_report.get(
                "legacy_fitness_scores_modified", False
            ),
        },
        "source_provenance": {
            "target_report_path": target_report.get("report_path"),
            "target_evaluator_sha256": file_hash(evaluator_path),
            "assessment_module_sha256": file_hash(Path(__file__)),
        },
        "candidates": rows,
    }

    output = Path(output_directory) / "p4_integrated_assessment"
    output.mkdir(parents=True, exist_ok=True)

    json_path = output / f"generation_{generation:04d}.json"
    csv_path = output / f"generation_{generation:04d}.csv"

    json_tmp = json_path.with_suffix(".json.tmp")
    json_tmp.write_text(
        json.dumps(json_safe(report), indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    json_tmp.replace(json_path)

    fields = list(rows[0]) if rows else [
        "candidate_index", "candidate_id", "rank", "target_status",
        "target_eligible", "target_geometry_score", "hotspot_proximity_score",
        "pose_source_groups",
        "reliability_status", "reliability_calibrated",
        "reliability_feedback_active", "epistasis_status", "ood_status",
        "candidate_decision",
    ]

    csv_tmp = csv_path.with_suffix(".csv.tmp")
    with csv_tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()

        csv_rows = []
        for row in rows:
            csv_row = dict(row)
            csv_row["pose_source_groups"] = json.dumps(
                json_safe(row.get("pose_source_groups")),
                sort_keys=True,
                separators=(",", ":"),
                allow_nan=False,
            )
            csv_rows.append(csv_row)

        writer.writerows(csv_rows)
    csv_tmp.replace(csv_path)

    report["report_path"] = str(json_path)
    report["candidate_csv_path"] = str(csv_path)
    return report
