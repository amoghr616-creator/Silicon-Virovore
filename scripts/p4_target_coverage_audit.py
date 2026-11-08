#!/usr/bin/env python3
"""Audit P4 target-evidence coverage from project-generated run reports.

No structures are synthesized and no missing target scores are imputed. The
script consumes only p4_target_evidence/generation_*.json reports produced by
run_pipeline.py and writes candidate-level CSV plus a summary JSON.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from statistics import mean
from typing import Any


def _read_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return data


def audit_run(run_dir: Path, out_dir: Path) -> dict[str, Any]:
    evidence_dir = run_dir / "p4_target_evidence"
    reports = sorted(evidence_dir.glob("generation_*.json"))
    if not reports:
        raise FileNotFoundError(
            f"No P4 target-evidence generation reports found under {evidence_dir}"
        )

    candidate_rows: list[dict[str, Any]] = []
    generation_rows: list[dict[str, Any]] = []
    all_shortlist = all_evidence = all_repeatable = 0
    all_candidates = 0

    for report_path in reports:
        report = _read_json(report_path)
        generation = int(report.get("generation", 0))
        records = report.get("candidates", [])
        shortlist = [
            r for r in records
            if (r.get("structure_shortlist") or {}).get("selected") is True
        ]
        valid = [r for r in shortlist if r.get("status") == "success"]
        repeatable = [r for r in shortlist if r.get("target_selection_eligible") is True]
        geometries = [
            float(r["target_geometry_score_mean"])
            for r in valid
            if r.get("target_geometry_score_mean") is not None
        ]

        for record in records:
            shortlist_info = record.get("structure_shortlist") or {}
            evidence = record.get("complex_generation") or {}
            reliability = record.get("reliability") or {}
            candidate_rows.append({
                "generation": generation,
                "candidate_id": record.get("candidate_id"),
                "sequence": record.get("sequence"),
                "shortlisted": bool(shortlist_info.get("selected")),
                "shortlist_rank": shortlist_info.get("rank"),
                "complex_generation_status": evidence.get("status"),
                "registered_complex_count": record.get("registered_complex_count", 0),
                "target_evidence_status": record.get("status"),
                "target_geometry_score_mean": record.get("target_geometry_score_mean"),
                "hotspot_proximity_score": record.get("hotspot_proximity_score"),
                "independent_run_count": reliability.get(
                    "n_independent_runs", reliability.get("independent_run_count", 0)
                ),
                "repeatability_status": reliability.get("status"),
                "selection_eligible": bool(record.get("target_selection_eligible")),
                "selection_reason": record.get("target_selection_reason"),
                "errors": json.dumps(record.get("errors", []), sort_keys=True),
            })

        n_short = len(shortlist)
        all_shortlist += n_short
        all_evidence += len(valid)
        all_repeatable += len(repeatable)
        all_candidates += len(records)
        generation_rows.append({
            "generation": generation,
            "population_candidates": len(records),
            "shortlist_candidates": n_short,
            "valid_target_evidence": len(valid),
            "repeatability_eligible": len(repeatable),
            "target_evidence_coverage_of_shortlist": len(valid) / n_short if n_short else 0.0,
            "repeatability_coverage_of_shortlist": len(repeatable) / n_short if n_short else 0.0,
            "mean_shortlist_geometry_score": mean(geometries) if geometries else None,
            "best_shortlist_geometry_score": max(geometries) if geometries else None,
            "pareto_reordered_candidates": bool(report.get("pareto_reordered_candidates")),
            "exploratory_override_used": bool(report.get("exploratory_override_used")),
            "selection_validation_status": report.get("selection_validation_status"),
            "report_path": str(report_path),
        })

    out_dir.mkdir(parents=True, exist_ok=True)
    candidate_csv = out_dir / "p4_target_coverage_candidates.csv"
    generation_csv = out_dir / "p4_target_coverage_generations.csv"

    def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
        if not rows:
            path.write_text("", encoding="utf-8")
            return
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    write_csv(candidate_csv, candidate_rows)
    write_csv(generation_csv, generation_rows)

    summary = {
        "run_dir": str(run_dir.resolve()),
        "generation_count": len(reports),
        "population_candidate_records": all_candidates,
        "shortlist_candidate_records": all_shortlist,
        "valid_target_evidence_records": all_evidence,
        "repeatability_eligible_records": all_repeatable,
        "target_evidence_coverage_of_shortlists": all_evidence / all_shortlist if all_shortlist else 0.0,
        "repeatability_coverage_of_shortlists": all_repeatable / all_shortlist if all_shortlist else 0.0,
        "generation_metrics": generation_rows,
        "candidate_csv": str(candidate_csv),
        "generation_csv": str(generation_csv),
        "affinity_validation": False,
        "interpretation": (
            "Coverage and structural repeatability accounting only. Geometry scores are computational "
            "measurements and are not experimental binding affinity or efficacy."
        ),
    }
    summary_path = out_dir / "p4_target_coverage_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    summary["summary_path"] = str(summary_path)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--out-dir", type=Path)
    args = parser.parse_args()
    run_dir = args.run_dir.expanduser().resolve()
    out_dir = args.out_dir.expanduser().resolve() if args.out_dir else run_dir / "p4_target_coverage_audit"
    summary = audit_run(run_dir, out_dir)
    print(json.dumps({k: v for k, v in summary.items() if k != "generation_metrics"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
