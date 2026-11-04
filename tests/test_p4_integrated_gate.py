import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from p4_integrated_gate import adjudicate_population


SEQ = "ACDEFGHIKLMNPQRSTVWYA"
MODEL_HASH = "a" * 64


def candidate():
    return {"candidate_id": "test_candidate", "sequence": SEQ}


def evidence():
    target = {
        "candidate_sequence": SEQ,
        "target_pdb_id": "9SYA",
        "docking_score_is_target_specific": True,
        "strongest_anchor_delta_g": 1.25,
        "hotspot_residues": [203, 205],
        "target_score": 0.75,
    }

    reliability = {
        "candidate_sequence": SEQ,
        "feature_dimension": 410,
        "uncertainty": 0.25,
        "calibration_status": "calibrated",
        "applicability_domain_status": "in_domain",
        "validation_split_unit": "family",
    }

    epistasis = {
        "candidate_sequence": SEQ,
        "predicted_delta_epsilon": -0.35,
        "epistasis_model_sha256": MODEL_HASH,
        "validation_family_count": 19,
        "validation_status": "family_held_out_validated",
    }

    return target, reliability, epistasis


def test_complete_candidate_can_enter_multiobjective_review():
    t, r, e = evidence()
    result = adjudicate_population([candidate()], [t], [r], [e])[0]
    assert result["decision"] == "ELIGIBLE_FOR_MULTI_OBJECTIVE_REVIEW"
    assert result["scalar_fitness_generated"] is False
    assert result["pruned"] is False


def test_missing_reliability_causes_abstention():
    t, _, e = evidence()
    result = adjudicate_population([candidate()], [t], [], [e])[0]
    assert result["decision"] == "ABSTAIN"
    assert not result["reliability"]["valid"]


def test_wrong_target_causes_abstention():
    t, r, e = evidence()
    t["target_pdb_id"] = "1ABC"
    result = adjudicate_population([candidate()], [t], [r], [e])[0]
    assert result["decision"] == "ABSTAIN"


def test_epistasis_requires_family_validation_provenance():
    t, r, e = evidence()
    e["validation_family_count"] = 2
    result = adjudicate_population([candidate()], [t], [r], [e])[0]
    assert result["decision"] == "ABSTAIN"


def test_sequence_mismatch_never_joins_evidence():
    t, r, e = evidence()
    t["candidate_sequence"] = "C" + SEQ[1:]
    result = adjudicate_population([candidate()], [t], [r], [e])[0]
    assert result["decision"] == "ABSTAIN"
    assert not result["target_optimization"]["valid"]
