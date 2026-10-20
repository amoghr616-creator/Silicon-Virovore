
import math

from src.mutation_credit import (
    MutationCreditARISE,
    generate_position_weighted_population,
    mutation_event_budget,
)


def test_single_mutation_records_parent_relative_reward():
    model = MutationCreditARISE(reward_source="sequence_proxy")

    accepted = model.observe_mutation(
        "AAAA", "CAAA", 0.20, 0.45, generation=1
    )

    assert accepted
    assert model.position_observations[0] == 1
    assert math.isclose(model.position_scores[0], 0.25)
    assert model.records[-1]["from_aa"] == "A"
    assert model.records[-1]["to_aa"] == "C"


def test_multi_mutation_is_not_attributed_to_one_position():
    model = MutationCreditARISE(reward_source="sequence_proxy")

    accepted = model.observe_mutation(
        "AAAA", "CCAA", 0.20, 0.60
    )

    assert not accepted
    assert not model.stats
    assert model.skipped_observations == 1


def test_probabilities_sum_to_one_and_keep_exploration_floor():
    model = MutationCreditARISE(
        reward_source="sequence_proxy",
        exploration_floor=0.25,
    )
    model.observe_mutation("AAAA", "CAAA", 0.20, 0.60)
    probabilities = model.position_probabilities(4)

    assert len(probabilities) == 4
    assert all(p >= 0.25 / 4 for p in probabilities)
    assert math.isclose(sum(probabilities), 1.0, abs_tol=1e-12)


def test_uniform_baseline_and_adaptive_use_same_event_budget():
    assert mutation_event_budget(21, 0.05) == 1

    parent = "MKLAVFALLVFFAGSSDLIRR"
    baseline = generate_position_weighted_population(
        parent, pop_size=20, mutation_rate=0.05, seed=17
    )
    adaptive = generate_position_weighted_population(
        parent,
        pop_size=20,
        mutation_rate=0.05,
        probabilities=[1.0] * len(parent),
        seed=17,
    )

    for population in (baseline, adaptive):
        assert len(population) == 20
        for child in population:
            changed = sum(a != b for a, b in zip(parent, child))
            assert changed == 1


def test_generator_is_deterministic_for_a_fixed_seed():
    parent = "MKLAVFALLVFFAGSSDLIRR"
    kwargs = dict(
        pop_size=12,
        mutation_rate=0.05,
        probabilities=[1.0] * len(parent),
        seed=991,
    )

    assert generate_position_weighted_population(
        parent, **kwargs
    ) == generate_position_weighted_population(parent, **kwargs)
def test_positive_mutation_credit_changes_position_probabilities():
    model = MutationCreditARISE(
        reward_source="sequence_proxy",
        memory_decay=0.90,
        exploration_floor=0.25,
        exploration_c=1.0,
    )

    before = model.position_probabilities(4)

    # Synthetic data tests the algorithm only; it is not protein-binding data.
    for _ in range(4):
        assert model.observe_mutation(
            "AAAA",
            "CAAA",
            parent_score=0.0,
            child_score=1.0,
            generation=1,
        )

    after = model.position_probabilities(4)

    assert after[0] > before[0]
    assert min(after) >= 0.25 / 4
    assert abs(sum(after) - 1.0) < 1e-12


def test_reset_clears_learning_and_reward_scale():
    model = MutationCreditARISE(
        reward_source="sequence_proxy",
    )

    model.observe_mutation("AAAA", "CAAA", 0.0, 1.0)
    model.reset()

    assert model.stats == {}
    assert model.records == []
    assert model.position_scores == {}
    assert model.reward_scale == model.prior_reward_scale
    assert model.position_probabilities(4) == [0.25] * 4