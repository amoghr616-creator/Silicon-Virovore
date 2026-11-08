import json
from pathlib import Path
from types import SimpleNamespace

from src.p4_structural_refill import run_bounded_structure_refill


def make_candidate(sequence):
    return SimpleNamespace(
        sequence=sequence,
        metadata={
            "p4_structure_shortlist": {
                "selected": True,
                "rank": 1,
                "selection_method": "fixture",
            }
        },
        structure_path=None,
        structure_confidence=None,
    )


def run_queue(tmp_path, candidates, outcomes, required=2, refill=2):
    attempts = []

    def predict(candidate):
        attempts.append(candidate.sequence)
        ok = outcomes[candidate.sequence]
        if ok:
            candidate.structure_path = tmp_path / f"{candidate.sequence}.pdb"
            candidate.structure_confidence = 75.0
            candidate.metadata["structure_available"] = True
            candidate.metadata["structure_status"] = "success"
            candidate.metadata.pop("structure_error", None)
        else:
            candidate.structure_path = None
            candidate.structure_confidence = None
            candidate.metadata["structure_available"] = False
            candidate.metadata["structure_status"] = "unavailable"
            candidate.metadata["structure_error"] = "test HTTP 504"

    selected, report = run_bounded_structure_refill(
        candidates,
        required_count=required,
        max_refill_candidates=refill,
        predict_one=predict,
        is_valid_structure=lambda candidate: bool(
            candidate.metadata.get("structure_available")
            and candidate.structure_path is not None
            and candidate.structure_confidence is not None
        ),
        output_directory=tmp_path,
        generation=1,
        run_id="test-run",
    )
    return selected, report, attempts


def test_no_refill_when_primary_candidates_succeed(tmp_path):
    candidates = [make_candidate(x) for x in ("A", "B", "C", "D")]
    selected, report, attempts = run_queue(
        tmp_path, candidates, {"A": True, "B": True, "C": True, "D": True}
    )

    assert [c.sequence for c in selected] == ["A", "B"]
    assert attempts == ["A", "B"]
    assert report["complete"] is True
    assert report["attempted_count"] == 2
    assert report["refill_attempted_count"] == 0
    assert report["unattempted_queue_count"] == 2
    assert report["status"] == "complete"


def test_failed_primary_is_replaced_from_frozen_queue(tmp_path):
    candidates = [make_candidate(x) for x in ("A", "B", "C", "D")]
    selected, report, attempts = run_queue(
        tmp_path,
        candidates,
        {"A": False, "B": True, "C": True, "D": True},
    )

    assert [c.sequence for c in selected] == ["B", "C"]
    assert attempts == ["A", "B", "C"]
    assert report["complete"] is True
    assert report["attempted_count"] == 3
    assert report["refill_attempted_count"] == 1
    assert candidates[0].metadata["p4_structure_shortlist"]["selected"] is False
    assert candidates[2].metadata["p4_structure_shortlist"]["selection_status"] == "selected_refill"


def test_incomplete_queue_preserves_requested_scope_with_placeholder(tmp_path):
    candidates = [make_candidate(x) for x in ("A", "B", "C", "D")]
    selected, report, attempts = run_queue(
        tmp_path,
        candidates,
        {"A": False, "B": True, "C": False, "D": False},
    )

    assert [c.sequence for c in selected] == ["B", "A"]
    assert attempts == ["A", "B", "C", "D"]
    assert report["complete"] is False
    assert report["status"] == "incomplete_budget_exhausted"
    assert report["eligible_count"] == 1
    assert report["selected_scope_count"] == 2
    assert report["reserved_incomplete_slot_count"] == 1
    assert candidates[0].metadata["p4_structure_shortlist"]["selection_status"] == "reserved_incomplete_slot"


def test_refill_cap_zero_never_attempts_backup_candidates(tmp_path):
    candidates = [make_candidate(x) for x in ("A", "B", "C")]
    selected, report, attempts = run_queue(
        tmp_path,
        candidates,
        {"A": False, "B": True, "C": True},
        required=2,
        refill=0,
    )

    assert attempts == ["A", "B"]
    assert report["max_attempts"] == 2
    assert report["complete"] is False
    assert len(selected) == 2


def test_refill_report_is_written_after_each_candidate_and_is_valid_json(tmp_path):
    candidates = [make_candidate(x) for x in ("A", "B", "C")]
    selected, report, attempts = run_queue(
        tmp_path,
        candidates,
        {"A": False, "B": True, "C": True},
        required=1,
        refill=2,
    )

    report_path = Path(report["report_path"])
    assert report_path.is_file()
    persisted = json.loads(report_path.read_text())
    assert persisted["complete"] is True
    assert persisted["attempted_count"] == 2
    assert len(selected) == 1
    assert attempts == ["A", "B"]
