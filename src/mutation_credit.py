
"""Interpretable mutation-credit learning for Silicon Virovore.

Scores must be higher-is-better and the parent and child scores must come
from the same declared objective. Target-specific runs should use only
independently validated structural scores.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import random
from typing import Any, Sequence


AA_POOL = "ACDEFGHIKLMNPQRSTVWY"


@dataclass
class PositionStats:
    observations: int = 0
    mean_delta: float = 0.0
    variance: float = 0.0


class MutationCreditARISE:
    """Learn which positions have produced beneficial mutations.

    The learned value is expected mutation reward, not biological residue
    importance. Multi-position mutations are excluded to avoid falsely
    attributing a whole-sequence score change to one position.
    """

    def __init__(
        self,
        *,
        reward_source: str = "validated_target_structural",
        memory_decay: float = 0.90,
        exploration_floor: float = 0.25,
        exploration_c: float = 1.0,
        prior_reward_scale: float = 0.01,
    ) -> None:
        if not reward_source:
            raise ValueError("reward_source must be explicit")
        if not 0.0 < memory_decay < 1.0:
            raise ValueError("memory_decay must be between 0 and 1")
        if not 0.0 <= exploration_floor <= 1.0:
            raise ValueError("exploration_floor must be in [0, 1]")
        if exploration_c < 0 or prior_reward_scale <= 0:
            raise ValueError("invalid exploration parameters")

        self.reward_source = reward_source
        self.memory_decay = memory_decay
        self.exploration_floor = exploration_floor
        self.exploration_c = exploration_c
        self.prior_reward_scale = float(prior_reward_scale)
        self.reward_scale = self.prior_reward_scale

        self.sequence_length = 0
        self.stats: dict[int, PositionStats] = {}
        self.position_scores: dict[int, float] = {}
        self.position_observations: dict[int, int] = {}
        self.records: list[dict[str, Any]] = []
        self.skipped_observations = 0

    def reset(self) -> None:
        """Clear learning state between independent runs."""
        self.sequence_length = 0
        self.stats.clear()
        self.position_scores.clear()
        self.position_observations.clear()
        self.records.clear()
        self.skipped_observations = 0
        self.reward_scale = 0.01
        self.reward_scale = self.prior_reward_scale

    @staticmethod
    def changed_positions(parent: str, child: str) -> list[int]:
        if not parent or len(parent) != len(child):
            raise ValueError("parent and child must have equal nonzero length")
        return [
            i for i, (a, b) in enumerate(zip(parent, child))
            if a != b
        ]

    def observe_mutation(
        self,
        parent_sequence: str,
        child_sequence: str,
        parent_score: float,
        child_score: float,
        *,
        generation: int | None = None,
    ) -> bool:
        """Record one eligible single-substitution outcome.

        Returns False for multi-site changes because their individual
        contributions cannot be identified from this observation alone.
        """
        if (
            not parent_sequence
            or len(parent_sequence) != len(child_sequence)
            or any(a not in AA_POOL for a in parent_sequence + child_sequence)
        ):
            raise ValueError("sequences must be equal-length canonical peptides")

        parent_score = float(parent_score)
        child_score = float(child_score)
        if not math.isfinite(parent_score) or not math.isfinite(child_score):
            raise ValueError("scores must be finite")

        changed = self.changed_positions(parent_sequence, child_sequence)
        if len(changed) != 1:
            self.skipped_observations += 1
            return False

        position = changed[0]
        delta = child_score - parent_score  # Positive means improvement.
        stat = self.stats.setdefault(position, PositionStats())

        if stat.observations == 0:
            stat.mean_delta = delta
            stat.variance = 0.0
        else:
            alpha = 1.0 - self.memory_decay
            error = delta - stat.mean_delta
            stat.mean_delta += alpha * error
            stat.variance = (
                self.memory_decay * stat.variance
                + alpha * error * (delta - stat.mean_delta)
            )
            stat.variance = max(0.0, stat.variance)

        stat.observations += 1
        self.sequence_length = len(parent_sequence)
        self.position_scores[position] = stat.mean_delta
        self.position_observations[position] = stat.observations

        alpha = 1.0 - self.memory_decay
        self.reward_scale = max(
            1e-6,
            self.memory_decay * self.reward_scale
            + alpha * abs(delta),
        )

        self.records.append({
            "generation": generation,
            "parent_sequence": parent_sequence,
            "child_sequence": child_sequence,
            "position_0_based": position,
            "position_1_based": position + 1,
            "from_aa": parent_sequence[position],
            "to_aa": child_sequence[position],
            "parent_score": parent_score,
            "child_score": child_score,
            "delta": delta,
            "position_observations": stat.observations,
            "reward_source": self.reward_source,
        })
        return True

    def observe_generation(self, candidates: list[Any]) -> None:
        """Learn only from validated, comparable parent-child rewards.

        Required metadata must come from the structural evaluation layer.
        A generic ranking-valid flag is not sufficient.
        """
        for candidate in candidates:
            metadata = getattr(candidate, "metadata", None)
            if not isinstance(metadata, dict):
                continue

            if metadata.get("ranking_observation_valid") is not True:
                continue

            if metadata.get("target_pose_passed_geometry_gate") is not True:
                continue

            if metadata.get("target_structural_score_validated") is not True:
                continue

            if metadata.get("mutation_reward_source") != self.reward_source:
                continue

            parent = metadata.get("parent_sequence")
            parent_score = metadata.get("mutation_parent_score")
            child_score = metadata.get("mutation_child_score")

            if (
                not isinstance(parent, str)
                or parent_score is None
                or child_score is None
            ):
                self.skipped_observations += 1
                continue

            try:
                accepted = self.observe_mutation(
                    parent,
                    candidate.sequence,
                    parent_score,
                    child_score,
                    generation=metadata.get("generation"),
                )
            except (TypeError, ValueError):
                self.skipped_observations += 1
                continue

            # Multi-position changes are intentionally excluded by
            # observe_mutation(), which returns False for those cases.
            if not accepted:
                continue

    def position_probabilities(self, length: int) -> list[float]:
        """Return normalized mutation-position probabilities.

        Uses a simple upper-confidence rule plus a guaranteed uniform
        exploration component. All positions remain explorable.
        """
        if length <= 0:
            raise ValueError("length must be positive")

        if not self.stats:
            return [1.0 / length] * length

        scale = max(self.reward_scale, 1e-6)
        upper_confidence: list[float] = []

        for position in range(length):
            stat = self.stats.get(position)
            if stat is None:
                count, mean = 0, 0.0
            else:
                count, mean = stat.observations, stat.mean_delta

            uncertainty_bonus = (
                self.exploration_c * scale / math.sqrt(count + 1)
            )
            upper_confidence.append(mean + uncertainty_bonus)

        minimum = min(upper_confidence)
        shifted = [
            max(0.0, value - minimum) + 1e-12
            for value in upper_confidence
        ]
        total = sum(shifted)
        exploitation = [value / total for value in shifted]

        floor = self.exploration_floor / length
        probabilities = [
            (1.0 - self.exploration_floor) * p + floor
            for p in exploitation
        ]

        # Remove floating-point summation error.
        norm = sum(probabilities)
        return [p / norm for p in probabilities]

    # Compatibility methods used by the existing population runner.

    def get_importance_array(self, length: int) -> list[float]:
        return self.position_probabilities(length)

    def importance_map(self) -> dict[int, float]:
        if not self.sequence_length:
            return {}
        probs = self.position_probabilities(self.sequence_length)
        return {i: round(p, 5) for i, p in enumerate(probs)}

    def importance_observations(self) -> dict[int, int]:
        return {
            i: self.position_observations.get(i, 0)
            for i in range(self.sequence_length)
        }

    def importance_statistics(self) -> dict[str, float | int]:
        if self.sequence_length == 0:
            return {
                "mean": 0.0,
                "nonuniformity": 0.0,
                "observed_positions": 0,
            }

        probs = self.position_probabilities(self.sequence_length)
        return {
            "mean": round(sum(probs) / len(probs), 5),
            "nonuniformity": round(max(probs) - min(probs), 5),
            "observed_positions": len(self.stats),
        }


def mutation_event_budget(length: int, mutation_rate: float) -> int:
    """Use one explicit, reproducible definition of mutation budget."""
    if length <= 0:
        raise ValueError("length must be positive")
    if not math.isfinite(mutation_rate) or not 0.0 <= mutation_rate <= 1.0:
        raise ValueError("mutation_rate must be finite and within [0, 1]")

    # Round half up instead of depending on Python's banker's rounding.
    return min(length, int(math.floor(length * mutation_rate + 0.5)))


def generate_position_weighted_population(
    parent: str,
    *,
    pop_size: int,
    mutation_rate: float,
    probabilities: Sequence[float] | None = None,
    seed: int | None = None,
) -> list[str]:
    """Generate fixed-budget mutants using weighted sampling without replacement.

    Passing probabilities=None gives the uniform-mutation baseline. Passing
    ARISE probabilities changes only which positions are selected; population
    size, mutation budget, amino-acid alphabet, and replacement rule stay fixed.
    """
    if not parent or any(aa not in AA_POOL for aa in parent):
        raise ValueError("parent must be a non-empty canonical peptide")
    if pop_size < 0:
        raise ValueError("pop_size cannot be negative")

    length = len(parent)
    budget = mutation_event_budget(length, mutation_rate)

    if probabilities is None:
        weights = [1.0] * length
    else:
        weights = [float(x) for x in probabilities]
        if len(weights) != length:
            raise ValueError("probabilities must match peptide length")
        if any(not math.isfinite(x) or x < 0 for x in weights):
            raise ValueError("probabilities must be finite and nonnegative")
        if sum(weights) <= 0:
            raise ValueError("at least one probability must be positive")

    rng = random.Random(seed) if seed is not None else random
    result: list[str] = []

    for _ in range(pop_size):
        child = list(parent)
        available = list(range(length))
        available_weights = list(weights)

        for _event in range(budget):
            total = sum(available_weights)
            if total <= 0:
                # Defensive fallback for externally supplied weights.
                available_weights = [1.0] * len(available)
                total = float(len(available))

            draw = rng.random() * total
            cumulative = 0.0
            selected_index = len(available) - 1

            for i, weight in enumerate(available_weights):
                cumulative += weight
                if draw < cumulative:
                    selected_index = i
                    break

            position = available.pop(selected_index)
            available_weights.pop(selected_index)

            choices = [aa for aa in AA_POOL if aa != child[position]]
            child[position] = rng.choice(choices)

        result.append("".join(child))

    return result
