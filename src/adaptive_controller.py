"""Conservative diversity/stagnation controller for adaptive search.

This changes how the search explores; it does not change the fitness
function and does not establish target binding.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
import math
from typing import Mapping


@dataclass(frozen=True)
class AdaptiveControl:
    action: str
    effective_mutation_rate: float
    hotspot_probability: float
    event_budget_mix: dict[int, float]
    diversity: float | None
    stagnation_generations: int

    def as_dict(self) -> dict:
        return asdict(self)


def resolve_adaptive_control(
    *,
    base_mutation_rate: float,
    sequence_length: int,
    diversity: float | None,
    stagnation_generations: int,
    hotspot_model=None,
) -> AdaptiveControl:
    """Adjust exploration using diversity and stagnation.

    Low diversity increases the probability of multi-event mutations.
    Hotspot guidance remains disabled pending independent validation.
    """
    if (
        not math.isfinite(float(base_mutation_rate))
        or base_mutation_rate < 0
    ):
        raise ValueError("base_mutation_rate must be finite and nonnegative")

    if not isinstance(sequence_length, int) or sequence_length < 1:
        raise ValueError("sequence_length must be a positive integer")

    if (
        not isinstance(stagnation_generations, int)
        or stagnation_generations < 0
    ):
        raise ValueError(
            "stagnation_generations must be a nonnegative integer"
        )

    if diversity is not None:
        diversity = float(diversity)
        if not math.isfinite(diversity) or not 0 <= diversity <= 1:
            raise ValueError("diversity must be None or in [0, 1]")

    action = "steady"
    multiplier = 1.0
    event_mix = {1: 0.80, 2: 0.17, 3: 0.03}

    if diversity is not None and diversity < 0.05:
        action = "strong_exploration"
        multiplier = 1.50 if stagnation_generations >= 2 else 1.25
        event_mix = {1: 0.35, 2: 0.40, 3: 0.25}

    elif (
        (diversity is not None and diversity < 0.09)
        or stagnation_generations >= 3
    ):
        action = "explore"
        multiplier = 1.20
        event_mix = {1: 0.55, 2: 0.32, 3: 0.13}

    effective_rate = min(
        0.25, max(0.0, float(base_mutation_rate) * multiplier)
    )

    return AdaptiveControl(
        action=action,
        effective_mutation_rate=effective_rate,
        hotspot_probability=0.0,
        event_budget_mix=event_mix,
        diversity=diversity,
        stagnation_generations=stagnation_generations,
    )


def allocate_event_budget_counts(
    population_size: int,
    mix: Mapping[int, float],
) -> dict[int, int]:
    """Allocate integer mutation-event budgets; counts sum exactly."""
    if not isinstance(population_size, int) or population_size < 0:
        raise ValueError("population_size must be a nonnegative integer")
    if not mix:
        raise ValueError("mix must not be empty")

    weights = {}
    for events, weight in mix.items():
        if not isinstance(events, int) or events < 0:
            raise ValueError("event keys must be nonnegative integers")
        weight = float(weight)
        if not math.isfinite(weight) or weight < 0:
            raise ValueError("weights must be finite and nonnegative")
        weights[events] = weight

    total = sum(weights.values())
    if total <= 0:
        raise ValueError("weights must sum to a positive value")

    exact = {
        events: population_size * weight / total
        for events, weight in weights.items()
    }
    counts = {
        events: math.floor(value)
        for events, value in exact.items()
    }

    remainder = population_size - sum(counts.values())
    order = sorted(
        weights,
        key=lambda event: (
            exact[event] - counts[event],
            -event,
        ),
        reverse=True,
    )

    for event in order[:remainder]:
        counts[event] += 1

    return {
        event: count
        for event, count in counts.items()
        if count > 0
    }


def reward_weighted_position_probabilities(
    reward_sums,
    reward_counts,
    *,
    temperature: float = 0.025,
    exploration_fraction: float = 0.25,
    prior_count: float = 3.0,
) -> list[float]:
    """Learn mutation-position probabilities without eliminating positions.

    Rewards are shrunk toward zero. A uniform exploration component prevents
    the algorithm from prematurely committing to a few positions.
    """
    if len(reward_sums) != len(reward_counts) or not reward_sums:
        raise ValueError(
            "reward_sums and reward_counts must have equal nonzero length"
        )
    if not math.isfinite(temperature) or temperature <= 0:
        raise ValueError("temperature must be finite and positive")
    if (
        not math.isfinite(exploration_fraction)
        or not 0 <= exploration_fraction <= 1
    ):
        raise ValueError("exploration_fraction must be in [0, 1]")
    if not math.isfinite(prior_count) or prior_count < 0:
        raise ValueError("prior_count must be finite and nonnegative")

    means = []
    for total, count in zip(reward_sums, reward_counts):
        total, count = float(total), float(count)
        if (
            not math.isfinite(total)
            or not math.isfinite(count)
            or count < 0
        ):
            raise ValueError("rewards must be finite; counts nonnegative")

        denominator = count + prior_count
        means.append(total / denominator if denominator else 0.0)

    center = max(means)
    logits = [
        max(-60.0, min(0.0, (value - center) / temperature))
        for value in means
    ]

    exponentials = [math.exp(value) for value in logits]
    denominator = sum(exponentials)
    learned = [value / denominator for value in exponentials]

    uniform = 1.0 / len(learned)
    probabilities = [
        (1.0 - exploration_fraction) * probability
        + exploration_fraction * uniform
        for probability in learned
    ]

    total_probability = sum(probabilities)
    return [p / total_probability for p in probabilities]
def record_mutation_credit(
    *,
    parent_fitness: float,
    child_fitness: float,
    mutated_positions: list[int],
    reward_sums: list[float],
    reward_counts: list[int],
) -> None:
    """Credit each changed position with an equal share of parent-relative gain.

    All fitness values must use the same objective. Positions are zero-based.
    This is heuristic credit assignment: interactions between mutations mean
    the contribution of an individual position is not causally identified.
    """
    if not math.isfinite(float(parent_fitness)):
        raise ValueError("parent_fitness must be finite")
    if not math.isfinite(float(child_fitness)):
        raise ValueError("child_fitness must be finite")

    if len(reward_sums) != len(reward_counts):
        raise ValueError("reward arrays must have equal length")

    # Validate the complete arrays before modifying any entries.
    for total in reward_sums:
        if not math.isfinite(float(total)):
            raise ValueError("reward_sums must contain only finite values")

    for count in reward_counts:
        if (
            not isinstance(count, int)
            or isinstance(count, bool)
            or count < 0
        ):
            raise ValueError(
                "reward_counts must contain nonnegative integers"
            )

    positions = list(mutated_positions)

    if any(
        not isinstance(pos, int) or isinstance(pos, bool)
        for pos in positions
    ):
        raise ValueError("mutated_positions must contain integer indices")

    positions = sorted(set(positions))

    if any(pos < 0 or pos >= len(reward_sums) for pos in positions):
        raise ValueError("mutation position outside reward arrays")

    if not positions:
        return

    delta = float(child_fitness) - float(parent_fitness)
    credit = delta / len(positions)

    for pos in positions:
        reward_sums[pos] = float(reward_sums[pos]) + credit
        reward_counts[pos] += 1