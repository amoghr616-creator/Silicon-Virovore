from types import SimpleNamespace

from src.p4_selection_policy import (
    p4_candidate_preference,
    p4_objective_vector,
)


def candidate(fitness, geometry, hotspot, overall=0.5):
    return SimpleNamespace(
        fitness=fitness,
        target_geometry_score_mean=geometry,
        hotspot_proximity_score=hotspot,
        overall_score=overall,
        selection_eligible=True,
    )


def test_pareto_dominating_candidate_wins():
    incumbent = candidate(0.70, 0.60, 0.60, overall=0.95)
    improved = candidate(0.75, 0.70, 0.70, overall=0.50)

    assert p4_candidate_preference(improved, incumbent) is True


def test_pareto_dominated_candidate_loses():
    incumbent = candidate(0.80, 0.80, 0.80)
    inferior = candidate(0.70, 0.70, 0.70)

    assert p4_candidate_preference(inferior, incumbent) is False


def test_incomparable_candidates_use_geometry_tiebreak():
    incumbent = candidate(0.90, 0.60, 0.90)
    candidate_with_better_geometry = candidate(0.70, 0.80, 0.70)

    assert p4_candidate_preference(
        candidate_with_better_geometry,
        incumbent,
    ) is True


def test_ineligible_candidate_has_no_p4_objectives():
    item = candidate(0.80, 0.80, 0.80)
    item.selection_eligible = False

    assert p4_objective_vector(item) is None


def test_incomplete_target_evidence_is_rejected():
    item = SimpleNamespace(
        fitness=0.80,
        overall_score=0.90,
    )

    assert p4_objective_vector(item) is None


def test_legacy_fallback_when_both_candidates_lack_p4_data():
    a = SimpleNamespace(overall_score=0.80)
    b = SimpleNamespace(overall_score=0.70)

    assert p4_candidate_preference(a, b) is None