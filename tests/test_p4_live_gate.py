from types import SimpleNamespace

from src.p4_integrated_gate import adjudicate_candidate_objects


def test_live_candidate_gate_joins_native_target_and_abstains_unvalidated_channels():
    candidate = SimpleNamespace(
        sequence="ACDEFGHIK",
        metadata={
            "candidate_id": "run:g1:c1",
            "p4_target_evidence": {
                "status": "success",
                "target_structure_id": "9SYA",
                "target_geometry_score_mean": 0.72,
                "hotspot_distances_A_mean": {
                    "A:203:ARG": 3.2,
                    "A:205:ARG": 3.6,
                },
            },
            "p4_reliability": {
                "status": "diagnostic_only",
                "calibrated": False,
                "feedback_active": False,
            },
        },
    )

    result = adjudicate_candidate_objects([candidate])[0]

    assert result["target_optimization"]["valid"] is True
    assert result["target_optimization"]["evidence_type"] == (
        "target_specific_structural_geometry"
    )
    assert result["target_optimization"]["affinity_inferred"] is False
    assert result["reliability"]["valid"] is False
    assert result["epistasis"]["valid"] is False
    assert result["decision"] == "ABSTAIN"
    assert result["scalar_fitness_generated"] is False


def test_live_candidate_gate_rejects_non_9sya_target():
    candidate = SimpleNamespace(
        sequence="ACDEFGHIK",
        metadata={
            "p4_target_evidence": {
                "status": "success",
                "target_structure_id": "OTHER",
                "target_geometry_score_mean": 0.72,
                "hotspot_distances_A_mean": {
                    "A:203:ARG": 3.2,
                    "A:205:ARG": 3.6,
                },
            },
        },
    )

    result = adjudicate_candidate_objects([candidate])[0]
    assert result["target_optimization"]["valid"] is False
    assert result["decision"] == "ABSTAIN"


def test_duplicate_population_sequences_do_not_duplicate_identical_evidence():
    metadata = {
        "p4_target_evidence": {
            "status": "success",
            "target_structure_id": "9SYA",
            "target_geometry_score_mean": 0.72,
            "hotspot_distances_A_mean": {
                "A:203:ARG": 3.2,
                "A:205:ARG": 3.6,
            },
        }
    }
    candidates = [
        SimpleNamespace(sequence="ACDEFGHIK", metadata=dict(metadata)),
        SimpleNamespace(sequence="ACDEFGHIK", metadata=dict(metadata)),
    ]

    results = adjudicate_candidate_objects(candidates)
    assert len(results) == 2
    assert all(
        result["target_optimization"]["valid"] is True
        for result in results
    )
