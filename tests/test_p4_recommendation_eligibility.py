import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from p4_integrated_gate import recommendation_eligible_candidates


class Candidate:
    def __init__(self, decision):
        self.metadata = {"p4_integrated_gate": {"decision": decision}}


def test_abstained_candidates_excluded_from_review_shortlist():
    eligible_a = Candidate("ELIGIBLE_FOR_MULTI_OBJECTIVE_REVIEW")
    abstained = Candidate("ABSTAIN")
    eligible_b = Candidate("ELIGIBLE_FOR_MULTI_OBJECTIVE_REVIEW")

    result = recommendation_eligible_candidates(
        [eligible_a, abstained, eligible_b]
    )

    assert result == [eligible_a, eligible_b]


def test_all_abstentions_produce_empty_shortlist():
    candidates = [Candidate("ABSTAIN"), Candidate("ABSTAIN")]

    assert recommendation_eligible_candidates(candidates) == []


def test_missing_gate_evidence_fails_closed():
    missing = Candidate("ABSTAIN")
    del missing.metadata["p4_integrated_gate"]
    eligible = Candidate("ELIGIBLE_FOR_MULTI_OBJECTIVE_REVIEW")

    result = recommendation_eligible_candidates([missing, eligible])

    assert result == [eligible]


def test_original_candidate_order_is_preserved():
    candidates = [
        Candidate("ELIGIBLE_FOR_MULTI_OBJECTIVE_REVIEW"),
        Candidate("ELIGIBLE_FOR_MULTI_OBJECTIVE_REVIEW"),
    ]

    result = recommendation_eligible_candidates(candidates)

    assert result[0] is candidates[0]
    assert result[1] is candidates[1]
