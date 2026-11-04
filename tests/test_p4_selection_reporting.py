from run_pipeline import (
    _finalize_comparable_telemetry,
    _p4_selected_winner_record,
)
from src.models import Candidate


def test_p4_selected_winner_is_distinct_from_score_winners():
    p4_winner = Candidate("A" * 21, c_score=0.2)
    score_winner = Candidate("C" * 21, c_score=0.9)

    p4_winner.overall_score = 0.2
    p4_winner.posthoc_recomputed_score = 0.6
    p4_winner.metadata.update({
        "candidate_id": "p4-selected",
        "generation": 1,
        "generation_selected_as_best": True,
        "p4_target_pareto_front": 0,
    })

    score_winner.overall_score = 0.9
    score_winner.posthoc_recomputed_score = 0.8
    score_winner.metadata.update({
        "candidate_id": "maximum-score",
        "generation": 1,
        "generation_selected_as_best": False,
        "p4_target_pareto_front": 1,
    })

    candidates = [p4_winner, score_winner]
    target_report = {
        "pareto_requested": True,
        "pareto_reordered_candidates": True,
        "reason": "Pareto selection applied.",
    }

    p4_record = _p4_selected_winner_record(
        candidates,
        target_report,
    )

    summaries = [{
        "generation": 1,
        "p4_selected_winner": p4_record,
    }]

    history = [
        {
            "candidate_id": candidate.metadata["candidate_id"],
            "generation": 1,
            "sequence": candidate.sequence,
        }
        for candidate in candidates
    ]

    _finalize_comparable_telemetry(
        summaries,
        candidates,
        history,
        score_threshold=0.7,
    )

    summary = summaries[0]
    record = summary["p4_selected_winner"]

    # The P4-selected candidate remains separately identifiable.
    assert record["sequence"] == p4_winner.sequence
    assert record["generation_score"] == 0.2
    assert record["pareto_front"] == 0
    assert record["selection_marker_verified"] is True
    assert record["posthoc_score"] == 0.6

    # Original-score maximum is a different candidate.
    assert summary["generation_local_best_sequence"] == score_winner.sequence
    assert summary["generation_local_best_score"] == 0.9

    # Post-hoc generation winner is also independently identified.
    assert summary["generation_best_sequence"] == score_winner.sequence
    assert summary["generation_best_score"] == 0.8

    # Backward-compatible aliases remain paired.
    assert summary["best_sequence"] == score_winner.sequence
    assert summary["best_overall"] == 0.8
