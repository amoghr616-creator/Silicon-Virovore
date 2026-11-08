"""Deterministic, budget-bounded refill for unavailable P4 structures.

Candidates are ordered before any structure prediction. Prediction outcomes do
not influence queue order; they only determine whether the next predeclared
candidate is attempted. Every attempt is written to an auditable JSON report.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Callable


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _candidate_refill_metadata(candidate: Any) -> dict[str, Any]:
    metadata = getattr(candidate, "metadata", None)
    if not isinstance(metadata, dict):
        raise TypeError("Every refill candidate must expose a metadata dictionary.")
    return metadata


def run_bounded_structure_refill(
    ordered_candidates: list[Any],
    *,
    required_count: int,
    max_refill_candidates: int,
    predict_one: Callable[[Any], Any],
    is_valid_structure: Callable[[Any], bool],
    output_directory: str | Path,
    generation: int,
    run_id: str,
) -> tuple[list[Any], dict[str, Any]]:
    """Try a frozen primary-plus-refill queue and retain exactly N selection slots.

    ``ordered_candidates`` must already be ranked deterministically before this
    function is called. The first ``required_count`` candidates are primary; up
    to ``max_refill_candidates`` subsequent candidates are backups. The function
    stops once N valid structures exist or the fixed queue is exhausted.

    If the queue is exhausted early, failed primary candidates reserve the
    unfilled shortlist slots. This deliberately causes downstream strict target
    evidence checks to abstain instead of silently shrinking the shortlist.
    """
    if required_count < 1:
        raise ValueError("required_count must be at least 1")
    if max_refill_candidates < 0:
        raise ValueError("max_refill_candidates must be non-negative")
    if len(ordered_candidates) < required_count:
        raise ValueError(
            "The pre-ranked refill queue must contain at least required_count candidates."
        )

    # Avoid spending repeated predictions on duplicate sequences while retaining
    # the supplied order. Duplicate sequences cannot be independent candidates.
    queue: list[Any] = []
    seen_sequences: set[str] = set()
    for candidate in ordered_candidates:
        sequence = str(getattr(candidate, "sequence", ""))
        if not sequence:
            raise ValueError("Refill candidates must have non-empty sequences.")
        if sequence in seen_sequences:
            continue
        seen_sequences.add(sequence)
        queue.append(candidate)
        if len(queue) >= required_count + max_refill_candidates:
            break

    if len(queue) < required_count:
        raise ValueError(
            "The unique pre-ranked refill queue has fewer candidates than required_count."
        )

    report_path = (
        Path(output_directory)
        / "p4_structure_refill"
        / f"generation_{int(generation):04d}.json"
    )
    report: dict[str, Any] = {
        "schema_version": 1,
        "run_id": str(run_id),
        "generation": int(generation),
        "selection_method": "precomputed_sequence_fitness_diversity_queue",
        "queue_order_frozen_before_prediction": True,
        "required_count": int(required_count),
        "queue_count": len(queue),
        "max_refill_candidates": int(max_refill_candidates),
        "max_attempts": min(
            len(queue), required_count + max_refill_candidates
        ),
        "attempted_count": 0,
        "refill_attempted_count": 0,
        "eligible_count": 0,
        "complete": False,
        "status": "running",
        "candidates": [],
        "report_path": str(report_path),
        "interpretation_limit": (
            "Structure availability handling only. Refill status is not evidence "
            "of binding affinity, structural accuracy, or efficacy."
        ),
    }

    rows: list[dict[str, Any]] = []
    for queue_rank, candidate in enumerate(queue, start=1):
        metadata = _candidate_refill_metadata(candidate)
        shortlist = metadata.setdefault("p4_structure_shortlist", {})
        if not isinstance(shortlist, dict):
            raise TypeError("p4_structure_shortlist metadata must be a dictionary.")
        shortlist.update({
            "selected": False,
            "queue_rank": queue_rank,
            "selection_role": "primary" if queue_rank <= required_count else "refill",
            "selection_status": "queued",
        })
        metadata["p4_structure_refill"] = {
            "queue_rank": queue_rank,
            "role": "primary" if queue_rank <= required_count else "refill",
            "attempted": False,
            "selected": False,
            "status": "queued",
        }
        sequence = str(candidate.sequence)
        rows.append({
            "queue_rank": queue_rank,
            "role": "primary" if queue_rank <= required_count else "refill",
            "sequence_sha256": hashlib.sha256(sequence.encode("utf-8")).hexdigest(),
            "sequence_length": len(sequence),
            "candidate_id": metadata.get("candidate_id"),
            "status": "queued",
            "attempted": False,
            "structure_available": False,
            "structure_path": None,
            "structure_confidence": None,
            "error": None,
            "selected": False,
            "selection_rank": None,
        })

    report["candidates"] = rows
    _atomic_json(report_path, report)

    eligible: list[Any] = []
    failed: list[Any] = []

    for index, candidate in enumerate(queue):
        row = rows[index]
        metadata = _candidate_refill_metadata(candidate)
        row["attempted"] = True
        row["status"] = "attempting"
        metadata["p4_structure_refill"].update({
            "attempted": True,
            "status": "attempting",
        })
        report["attempted_count"] += 1
        if index >= required_count:
            report["refill_attempted_count"] += 1
        _atomic_json(report_path, report)

        # The predictor already implements its own bounded retry/backoff policy.
        # Execute one candidate at a time so queue order and budgets stay auditable.
        predict_one(candidate)
        valid = bool(is_valid_structure(candidate))
        structure_status = str(metadata.get("structure_status", "unknown"))
        structure_error = metadata.get("structure_error")
        row.update({
            "status": "structure_available" if valid else structure_status,
            "structure_available": valid,
            "structure_path": (
                str(getattr(candidate, "structure_path", None))
                if getattr(candidate, "structure_path", None) is not None
                else None
            ),
            "structure_confidence": getattr(candidate, "structure_confidence", None),
            "error": str(structure_error) if structure_error else None,
        })
        metadata["p4_structure_refill"].update({
            "status": row["status"],
            "structure_available": valid,
            "error": row["error"],
        })

        if valid:
            eligible.append(candidate)
            if len(eligible) <= required_count:
                selection_rank = len(eligible)
                shortlist = metadata["p4_structure_shortlist"]
                shortlist.update({
                    "selected": True,
                    "rank": selection_rank,
                    "selection_role": row["role"],
                    "selection_status": (
                        "selected_primary" if row["role"] == "primary"
                        else "selected_refill"
                    ),
                })
                metadata["p4_structure_refill"].update({
                    "selected": True,
                    "selection_rank": selection_rank,
                })
                row["selected"] = True
                row["selection_rank"] = selection_rank
        else:
            failed.append(candidate)

        report["eligible_count"] = len(eligible)
        report["complete"] = len(eligible) >= required_count
        report["status"] = "complete" if report["complete"] else "running"
        _atomic_json(report_path, report)

        if report["complete"]:
            # Mark remaining pre-ranked queue entries without attempting them.
            for remaining_index in range(index + 1, len(queue)):
                remaining = queue[remaining_index]
                remaining_meta = _candidate_refill_metadata(remaining)
                remaining_meta["p4_structure_refill"].update({
                    "attempted": False,
                    "selected": False,
                    "status": "not_attempted_quota_filled",
                })
                rows[remaining_index]["status"] = "not_attempted_quota_filled"
            break

    if len(eligible) < required_count:
        # Preserve the declared shortlist width. Failed primary candidates fill
        # vacant slots as explicit ineligible placeholders; downstream evidence
        # gates then abstain rather than treating a smaller shortlist as complete.
        selected_ids = {id(candidate) for candidate in eligible}
        placeholders = [
            candidate for candidate in failed[:required_count]
            if id(candidate) not in selected_ids
        ]
        reserved_count = 0
        for candidate in placeholders:
            if len(eligible) + reserved_count >= required_count:
                break
            metadata = _candidate_refill_metadata(candidate)
            shortlist = metadata["p4_structure_shortlist"]
            selection_rank = len(eligible) + reserved_count + 1
            shortlist.update({
                "selected": True,
                "rank": selection_rank,
                "selection_status": "reserved_incomplete_slot",
            })
            metadata["p4_structure_refill"].update({
                "selected": True,
                "selection_rank": selection_rank,
                "reserved_incomplete_slot": True,
            })
            for row in rows:
                if row["sequence_sha256"] == hashlib.sha256(
                    str(candidate.sequence).encode("utf-8")
                ).hexdigest():
                    row["selected"] = True
                    row["selection_rank"] = selection_rank
                    row["status"] = "reserved_incomplete_slot"
                    break
            reserved_count += 1

    selected_scope = [
        candidate for candidate in queue
        if _candidate_refill_metadata(candidate)
        .get("p4_structure_shortlist", {}).get("selected") is True
    ]
    selected_scope.sort(
        key=lambda candidate: int(
            _candidate_refill_metadata(candidate)
            .get("p4_structure_shortlist", {}).get("rank", 10**9)
        )
    )
    if len(selected_scope) != required_count:
        raise RuntimeError(
            "Bounded refill failed to preserve the declared shortlist width: "
            f"selected_scope_count={len(selected_scope)}, required_count={required_count}."
        )
    report["eligible_count"] = len(eligible)
    report["selected_scope_count"] = len(selected_scope)
    report["complete"] = len(eligible) >= required_count
    report["status"] = (
        "complete" if report["complete"] else "incomplete_budget_exhausted"
    )
    report["unattempted_queue_count"] = sum(
        not row["attempted"] for row in rows
    )
    report["reserved_incomplete_slot_count"] = sum(
        bool(_candidate_refill_metadata(candidate)["p4_structure_refill"].get("reserved_incomplete_slot"))
        for candidate in selected_scope
    )
    _atomic_json(report_path, report)
    return selected_scope, report
