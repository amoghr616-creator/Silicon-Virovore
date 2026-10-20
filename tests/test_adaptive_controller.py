import pytest

from src.adaptive_controller import (
    allocate_event_budget_counts,
    resolve_adaptive_control,
    reward_weighted_position_probabilities,
    record_mutation_credit,
)


def test_budget_always_equals_population():
    for size in (0, 1, 8, 16, 101):
        counts = allocate_event_budget_counts(
            size, {1: 0.35, 2: 0.40, 3: 0.25}
        )
        assert sum(counts.values()) == size
        assert all(value >= 0 for value in counts.values())


def test_low_diversity_triggers_exploration():
    control = resolve_adaptive_control(
        base_mutation_rate=0.05,
        sequence_length=21,
        diversity=0.03,
        stagnation_generations=3,
    )
    assert control.action == "strong_exploration"
    assert control.effective_mutation_rate > 0.05
    assert control.hotspot_probability == 0.0


def test_missing_diversity_does_not_trigger_false_alarm():
    control = resolve_adaptive_control(
        base_mutation_rate=0.05,
        sequence_length=21,
        diversity=None,
        stagnation_generations=0,
    )
    assert control.action == "steady"
    assert control.effective_mutation_rate == pytest.approx(0.05)


def test_position_probabilities_are_valid():
    probabilities = reward_weighted_position_probabilities(
        [1.0, 0.0, -1.0],
        [10, 10, 10],
        temperature=0.025,
        exploration_fraction=0.25,
        prior_count=3.0,
    )
    assert sum(probabilities) == pytest.approx(1.0)
    assert all(p >= 0.25 / 3 for p in probabilities)
    assert probabilities[0] > probabilities[1] > probabilities[2]


def test_invalid_reward_counts_are_rejected():
    with pytest.raises(ValueError):
        reward_weighted_position_probabilities([1, 2], [1])

    with pytest.raises(ValueError):
        reward_weighted_position_probabilities([1, 2], [1, -1])


def test_mutation_credit_uses_parent_relative_improvement():
    reward_sums = [0.0, 0.0, 0.0]
    reward_counts = [0, 0, 0]

    record_mutation_credit(
        parent_fitness=0.60,
        child_fitness=0.66,
        mutated_positions=[0, 2],
        reward_sums=reward_sums,
        reward_counts=reward_counts,
    )

    assert reward_sums == pytest.approx([0.03, 0.0, 0.03])
    assert reward_counts == [1, 0, 1]

def test_mutation_credit_records_negative_improvement():
    reward_sums = [0.0, 0.0]
    reward_counts = [0, 0]

    record_mutation_credit(
        parent_fitness=0.70,
        child_fitness=0.68,
        mutated_positions=[1],
        reward_sums=reward_sums,
        reward_counts=reward_counts,
    )

    assert reward_sums == pytest.approx([0.0, -0.02])
    assert reward_counts == [0, 1]


def test_mutation_credit_rejects_invalid_inputs():
    sums = [0.0, 0.0]
    counts = [0, 0]

    with pytest.raises(ValueError):
        record_mutation_credit(
            parent_fitness=float("nan"),
            child_fitness=0.5,
            mutated_positions=[0],
            reward_sums=sums,
            reward_counts=counts,
        )

    with pytest.raises(ValueError):
        record_mutation_credit(
            parent_fitness=0.5,
            child_fitness=0.6,
            mutated_positions=[2],
            reward_sums=sums,
            reward_counts=counts,
        )

    with pytest.raises(ValueError):
        record_mutation_credit(
            parent_fitness=0.5,
            child_fitness=0.6,
            mutated_positions=[0],
            reward_sums=sums,
            reward_counts=[-1, 0],
        )