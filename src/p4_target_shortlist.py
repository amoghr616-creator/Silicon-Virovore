"""Budget-aware P4 target-evidence shortlists.

This module keeps two operations distinct:

* pre-structure selection chooses a deterministic, diverse subset for the
  expensive structure/HADDOCK budget using sequence fitness and novelty;
* post-scoring reporting ranks only candidates with valid target evidence.

Neither operation changes the population's primary parent-selection order.
Structural geometry/proximity are heuristic diagnostics, not affinity values.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any


def _fitness(candidate: Any) -> float:
    value = getattr(candidate, "c_score", None)
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float("-inf")
    return number if math.isfinite(number) else float("-inf")


def _edit_distance(left: str, right: str) -> int:
    """Levenshtein distance; sequences are short fixed-length peptides."""
    if len(left) < len(right):
        left, right = right, left
    previous = list(range(len(right) + 1))
    for i, a in enumerate(left, start=1):
        current = [i]
        for j, b in enumerate(right, start=1):
            current.append(min(
                current[-1] + 1,
                previous[j] + 1,
                previous[j - 1] + (a != b),
            ))
        previous = current
    return previous[-1]


def select_p4_structure_shortlist(
    candidates: list[Any],
    limit: int,
) -> list[Any]:
    """Select a deterministic sequence-fitness/diversity shortlist.

    The first slot is the highest finite sequence-fitness candidate. Further
    slots maximize minimum edit distance from already selected sequences,
    with sequence fitness and stable input order as tie-breakers.
    """
    if not isinstance(limit, int) or limit < 1:
        raise ValueError("P4 shortlist limit must be an integer >= 1")

    indexed = []
    for index, candidate in enumerate(candidates):
        sequence = str(getattr(candidate, "sequence", "")).strip().upper()
        if not sequence:
            candidate.metadata["p4_structure_shortlist"] = {
                "selected": False,
                "status": "invalid_sequence",
            }
            continue
        indexed.append({
            "index": index,
            "candidate": candidate,
            "sequence": sequence,
            "fitness": _fitness(candidate),
        })

    # The seed candidate is best by fitness, with deterministic tie-breaking.
    indexed.sort(key=lambda row: (
        -row["fitness"], row["sequence"], row["index"]
    ))
    selected = indexed[:1]
    remaining = indexed[1:]

    while remaining and len(selected) < limit:
        remaining.sort(key=lambda row: (
            -min(
                _edit_distance(row["sequence"], chosen["sequence"])
                for chosen in selected
            ),
            -row["fitness"],
            row["sequence"],
            row["index"],
        ))
        selected.append(remaining.pop(0))

    selected_by_index = {row["index"]: rank for rank, row in enumerate(selected, 1)}
    selected_rows = {row["index"]: row for row in selected}

    for index, candidate in enumerate(candidates):
        rank = selected_by_index.get(index)
        if rank is None:
            candidate.metadata["p4_structure_shortlist"] = {
                "selected": False,
                "status": "not_selected_budget",
                "limit": limit,
            }
            continue

        row = selected_rows[index]
        rationale = (
            "highest_sequence_fitness"
            if rank == 1
            else "maximum_minimum_edit_distance_then_fitness"
        )
        min_distance = None if rank == 1 else min(
            _edit_distance(row["sequence"], prior["sequence"])
            for prior in selected[: rank - 1]
        )
        candidate.metadata["p4_structure_shortlist"] = {
            "selected": True,
            "status": "selected",
            "rank": rank,
            "limit": limit,
            "selection_method": "fitness_seeded_greedy_max_min_sequence_diversity",
            "rationale": rationale,
            "sequence_fitness": (
                row["fitness"] if math.isfinite(row["fitness"]) else None
            ),
            "minimum_edit_distance_to_prior_selected": min_distance,
        }

    return [row["candidate"] for row in selected]


def write_p4_target_shortlist(
    candidates: list[Any],
    *,
    generation: int,
    output_directory: str | Path,
) -> dict[str, Any]:
    """Write a diagnostic ranking of candidates with valid target evidence.

    This is deliberately separate from population ranking and does not alter
    candidate rank, overall_score, c_score, or parent selection.
    """
    scored = []
    statuses = {}

    for candidate in candidates:
        meta = candidate.metadata if isinstance(candidate.metadata, dict) else {}
        evidence = meta.get("p4_target_evidence", {})
        status = str(evidence.get("status", "missing"))
        statuses[status] = statuses.get(status, 0) + 1
        if status != "success":
            continue

        geometry = evidence.get("target_geometry_score_mean")
        proximity = evidence.get("hotspot_proximity_score")
        try:
            geometry = float(geometry)
            proximity = float(proximity)
        except (TypeError, ValueError):
            continue
        if not math.isfinite(geometry) or not math.isfinite(proximity):
            continue

        scored.append({
            "candidate_id": meta.get("candidate_id"),
            "sequence": candidate.sequence,
            "target_geometry_score_mean": geometry,
            "hotspot_proximity_score": proximity,
            "sequence_fitness_c_score": _fitness(candidate),
            "reliability_status": meta.get("p4_reliability", {}).get(
                "status", "unavailable"
            ),
            "reliability_calibrated": False,
            "target_evidence_status": status,
        })

    # Keep metrics visible and use a transparent lexicographic order. No
    # arbitrary weighted composite is presented as an affinity estimate.
    scored.sort(key=lambda row: (
        -row["target_geometry_score_mean"],
        -row["hotspot_proximity_score"],
        -row["sequence_fitness_c_score"],
        row["sequence"],
    ))
    for rank, row in enumerate(scored, start=1):
        row["diagnostic_rank"] = rank

    report = {
        "generation": generation,
        "status": "available" if scored else "no_valid_target_evidence",
        "candidate_count": len(candidates),
        "valid_target_evidence_count": len(scored),
        "target_evidence_coverage": (
            len(scored) / len(candidates) if candidates else 0.0
        ),
        "status_counts": statuses,
        "ranking_method": (
            "lexicographic_target_geometry_then_hotspot_proximity_then_sequence_fitness"
        ),
        "ranking_scope": "valid_target_evidence_candidates_only",
        "population_reordered": False,
        "parent_selection_changed": False,
        "reliability_calibrated": False,
        "binding_affinity_claim": False,
        "interpretation": (
            "Diagnostic shortlist only. Target geometry and hotspot proximity "
            "are heuristic structural measurements, not experimentally "
            "validated binding affinity. Unscored candidates are not mixed "
            "into this ranking, and the population parent-selection order is unchanged."
        ),
        "ranked_candidates": scored,
    }

    out = Path(output_directory) / "p4_target_evidence"
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"target_shortlist_generation_{generation:04d}.json"
    path.write_text(json.dumps(report, indent=2, default=str) + "\n", encoding="utf-8")
    report["report_path"] = str(path)
    return report
