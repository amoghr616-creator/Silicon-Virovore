import json
from types import SimpleNamespace

from src.p4_generation_audit import (
    write_p4_generation_assessment,
)


def make_candidate(target_status="success", calibrated=False):
    evidence = {
        "status": target_status,
        "target_geometry_score_mean": 0.62 if target_status == "success" else None,
        "hotspot_proximity_score": 0.71 if target_status == "success" else None,
    }
    reliability = {
        "status": "diagnostic_only",
        "calibrated": calibrated,
        "feedback_active": False,
    }
    return SimpleNamespace(
        sequence="EXAMPLESEQUENCE",
        rank=1,
        overall_score=0.8,
        c_score=0.7,
        metadata={
            "candidate_id": "candidate-1",
            "p4_target_evidence": evidence,
            "p4_reliability": reliability,
        },
    )


def make_target_report(tmp_path, pareto=False):
    return {
        "generation": 2,
        "candidate_count": 1,
        "target_evidence_count": 1,
        "target_evidence_coverage": 1.0,
        "pareto_requested": pareto,
        "pareto_reordered_candidates": pareto,
        "reason": "test",
        "legacy_fitness_scores_modified": False,
        "reliability_feedback_active": False,
        "reliability_calibrated": False,
        "reliability_diagnostic_candidate_count": 1,
        "manifest": {"status": "test"},
        "report_path": str(tmp_path / "target.json"),
    }


def test_assessment_writes_json_and_csv_without_changing_candidate(tmp_path):
    candidate = make_candidate()
    before = (
        candidate.rank,
        candidate.overall_score,
        candidate.c_score,
    )

    report = write_p4_generation_assessment(
        [candidate],
        make_target_report(tmp_path),
        tmp_path,
    )

    assert Path(report["report_path"]).is_file()
    assert Path(report["candidate_csv_path"]).is_file()
    assert report["target_evidence_coverage"] == 1.0
    assert report["reliability"]["calibrated"] is False
    assert report["reliability"]["feedback_active"] is False
    assert report["evsi_and_pruning"]["pruning_allowed"] is False
    assert report["candidates"][0]["candidate_decision"] == (
        "TARGET_ASSESSED_LEGACY_ORDER_RETAINED"
    )

    after = (
        candidate.rank,
        candidate.overall_score,
        candidate.c_score,
    )
    assert before == after

    saved = json.loads(Path(report["report_path"]).read_text())
    assert saved["schema_version"] == "p4-generation-assessment-v1"


def test_missing_target_evidence_is_explicit_abstention(tmp_path):
    candidate = make_candidate(target_status="missing_complex")

    report = write_p4_generation_assessment(
        [candidate],
        make_target_report(tmp_path),
        tmp_path,
    )

    assert report["target_evidence_coverage"] == 0.0
    assert report["candidates"][0]["candidate_decision"] == (
        "ABSTAIN_TARGET_EVIDENCE_INVALID"
    )


def test_epistasis_is_not_misrepresented_as_candidate_evidence(tmp_path):
    candidate = make_candidate()

    report = write_p4_generation_assessment(
        [candidate],
        make_target_report(tmp_path),
        tmp_path,
    )

    assert report["epistasis_benchmark"]["candidate_level_epistasis_used"] is False
    assert report["candidates"][0]["epistasis_status"] == (
        "no_candidate_specific_evidence"
    )


def test_calibration_flag_requires_explicit_true_on_all_candidates(tmp_path):
    first = make_candidate(calibrated=True)
    second = make_candidate(calibrated=False)

    report = write_p4_generation_assessment(
        [first, second],
        make_target_report(tmp_path),
        tmp_path,
    )

    assert report["reliability"]["calibrated_candidate_count"] == 1
    assert report["reliability"]["calibrated"] is False
    assert report["reliability"]["feedback_active"] is False


from pathlib import Path
