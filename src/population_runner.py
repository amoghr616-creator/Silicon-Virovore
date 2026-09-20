"""
population_runner.py

Generates candidate peptides using the native C backend and
implements the Adaptive Recursive Intelligent Selection Engine (ARISE).

ARISE discovers local sequence windows associated with higher
computational performance. ALE can then guide future mutation positions
using a frozen hotspot model, without treating the association as causal.
"""

from __future__ import annotations

import logging

from c.bridge import (
    generate_c_population,
    generate_adaptive_population,
    generate_policy_population,
    process_candidate_peptide,
)

from src.ale_policy import resolve_mutation_policy, mutation_positions
from src.models import Candidate

from src.config import (
    ARISE_ENABLED,
    ARISE_MIN_OBSERVATIONS,
    DIVERSITY_REFILL_ATTEMPTS,
    REMOVE_DUPLICATE_CANDIDATES,
    POPULATION_SIZE,
    MUTATION_RATE,
)

logger = logging.getLogger(__name__)


# ==========================================================
# ARISE
# ==========================================================

class ARISEEngine:
    """
    Adaptive Recursive Intelligent Selection Engine.

    Estimates position-specific importance from the relationship
    between residue identity and candidate performance.
    """

    def __init__(self):
        # Running raw importance accumulated across generations.
        self.position_scores: dict[int, float] = {}

        # Current normalized importance values.
        self._importance: dict[int, float] = {}
        self.sequence_length = 0
        self.position_observations: dict[int, int] = {}

    def reset(self) -> None:
        """Reset learning state at the beginning of an experiment."""

        self.position_scores.clear()
        self._importance.clear()
        self.position_observations.clear()
        self.sequence_length = 0

    def observe_generation(
        self,
        ranked_candidates: list[Candidate],
    ) -> None:
        """
        Estimate positional importance from score differences
        associated with residue identity.

        A position receives stronger importance when different
        residues at that position are associated with meaningfully
        different candidate scores.
        """

        valid_candidates = [
            candidate
            for candidate in ranked_candidates
            if candidate.metadata.get("ranking_observation_valid") is True
        ]

        if not valid_candidates:
            return

        for candidate in valid_candidates:
            candidate.metadata["arise_observation_score"] = float(
                candidate.overall_score
            )
            candidate.metadata["arise_observation_source"] = "overall_score"

        logger.info(
            "ARISE objective: overall_score | candidates=%d",
            len(valid_candidates),
        )

        sequence_length = len(valid_candidates[0].sequence)
        self.sequence_length = sequence_length

        for position in range(sequence_length):

            residue_scores: dict[str, list[float]] = {}

            for candidate in valid_candidates:

                if position >= len(candidate.sequence):
                    continue

                residue = candidate.sequence[position]

                residue_scores.setdefault(
                    residue,
                    [],
                )

                residue_scores[residue].append(
                    float(candidate.overall_score)
                )

            total_observations = sum(
                len(scores)
                for scores in residue_scores.values()
            )
            self.position_observations[position] = total_observations

            # Need at least two distinct residue groups
            # to estimate a positional effect.
            if len(residue_scores) < 2:
                continue

            if total_observations < ARISE_MIN_OBSERVATIONS:
                continue

            residue_means = []

            for scores in residue_scores.values():

                if not scores:
                    continue

                residue_means.append(
                    sum(scores) / len(scores)
                )

            if len(residue_means) < 2:
                continue

            positional_signal = (
                max(residue_means)
                - min(residue_means)
            )

            self.position_scores[position] = (
                0.7 * self.position_scores.get(
                    position,
                    0.0,
                )
                + 0.3 * positional_signal
            )

    def update_importance(self) -> None:
        """
        Normalize accumulated positional signals to [0, 1].
        """

        if not self.position_scores:
            return

        values = list(
            self.position_scores.values()
        )

        minimum = min(values)
        maximum = max(values)

        if abs(maximum - minimum) < 1e-12:

            self._importance = {
                position: 0.5
                for position in self.position_scores
            }

            return

        self._importance = {
            position: (
                (value - minimum)
                / (maximum - minimum)
            )
            for position, value
            in self.position_scores.items()
        }

    def importance_map(self) -> dict[int, float]:
        """
        Return normalized importance for every residue position.
        """

        return {
            position: round(self._importance.get(position, 0.5), 3)
            for position in range(self.sequence_length)
        }

    def importance_observations(self) -> dict[int, int]:
        """Return the number of valid observations used per position."""

        return {
            position: self.position_observations.get(position, 0)
            for position in range(self.sequence_length)
        }

    def importance_statistics(self) -> dict[str, float | int]:
        values = self.get_importance_array(self.sequence_length)
        if not values:
            return {
                "mean": 0.0,
                "nonuniformity": 0.0,
                "observed_positions": 0,
            }

        mean = sum(values) / len(values)
        return {
            "mean": round(mean, 4),
            "nonuniformity": round(
                max(values) - min(values),
                4,
            ),
            "observed_positions": sum(
                1
                for position in range(self.sequence_length)
                if self.position_observations.get(position, 0) > 0
            ),
        }

    def get_importance_array(
        self,
        length: int,
    ) -> list[float]:
        """
        Return importance values in sequence-position order.
        """

        return [
            self._importance.get(
                position,
                0.5,
            )
            for position in range(length)
        ]


ARISE_ENGINE = ARISEEngine()

# ==========================================================
# Population generation
# ==========================================================


def generate_candidates(
    seed_sequence: str,
    population_size: int = POPULATION_SIZE,
    mutation_rate: float = MUTATION_RATE,
    arise_enabled: bool = ARISE_ENABLED,
    cache_enabled: bool = True,
    hotspot_model: dict | None = None,
    hotspot_enabled: bool = False,
    ale_enabled: bool = False,
    mutation_guidance_mode: str = "legacy",
    generation: int = 1,
    run_id: str = "",
    preserve_sequence: str | None = None,
) -> list[Candidate]:
    """
    Generate candidate peptides.

    ``legacy`` preserves the historical native generator. Controlled
    baseline/guided experiments use the policy generator so both conditions
    receive the same fixed mutation-event budget.

    Controlled baseline/guided conditions use the policy generator with a
    fixed mutation-event budget. The legacy positional-importance generator
    remains available for compatibility when ``mutation_guidance_mode`` is
    ``legacy``.

    ALE therefore changes only the mutation-position policy and does not
    perform a second Python-side mutation pass.
    """

    logger.info("Generating peptide candidates...")

    importance_map = ARISE_ENGINE.get_importance_array(len(seed_sequence))
    has_learned_importance = bool(ARISE_ENGINE.position_scores)
    logger.info(
        "ARISE legacy mode %s | importance=%s | observations=%s",
        "enabled" if arise_enabled else "disabled",
        [round(x, 3) for x in importance_map],
        ARISE_ENGINE.importance_observations(),
    )

    policy = resolve_mutation_policy(
        mutation_guidance_mode
        if mutation_guidance_mode != "legacy"
        else "uniform",
        hotspot_model,
        enabled=hotspot_enabled and ale_enabled,
    )
    if mutation_guidance_mode != "legacy":
        logger.info(
            "ALE policy | requested=%s | active=%s | guided=%s | hotspot=%s:%s",
            policy["requested_mode"],
            policy["mode"],
            policy["guided"],
            policy["hotspot_start"],
            policy["hotspot_end"],
        )

    def generate_batch(size: int) -> list[str]:
        if mutation_guidance_mode != "legacy":
            return generate_policy_population(
                seed_sequence,
                mutation_rate=mutation_rate,
                guidance_mode=policy["mode"],
                hotspot_start=policy["hotspot_start"],
                hotspot_end=policy["hotspot_end"],
                pop_size=size,
            )
        if arise_enabled and has_learned_importance:
            return generate_adaptive_population(
                seed_sequence,
                importance_map,
                pop_size=size,
            )
        return generate_c_population(
            seed_sequence,
            pop_size=size,
            mutation_rate=mutation_rate,
        )

    sequences: list[str] = []
    seen: set[str] = set()
    attempts = 0
    while len(sequences) < population_size:
        batch = generate_batch(population_size - len(sequences))
        if not REMOVE_DUPLICATE_CANDIDATES:
            sequences.extend(batch)
            break

        before = len(sequences)
        for sequence in batch:
            if sequence not in seen:
                seen.add(sequence)
                sequences.append(sequence)
                if len(sequences) >= population_size:
                    break

        attempts += 1
        if len(sequences) >= population_size:
            break
        if attempts >= DIVERSITY_REFILL_ATTEMPTS:
            logger.warning(
                "Could only generate %d/%d unique candidates after %d attempts.",
                len(sequences),
                population_size,
                attempts,
            )
            if before == len(sequences):
                break

    preserved_elite = False
    if preserve_sequence and preserve_sequence not in sequences and sequences:
        sequences[-1] = preserve_sequence
        preserved_elite = True

    candidates: list[Candidate] = []

    for sequence in sequences:

        candidate = process_candidate_peptide(
            sequence,
            cache_enabled=cache_enabled,
        )

        if candidate is None:
            continue

        changed_positions = mutation_positions(
            seed_sequence,
            sequence,
        )
        is_preserved = preserved_elite and sequence == preserve_sequence
        candidate.metadata["changed_positions_count"] = len(changed_positions)
        event_budget = 0 if is_preserved else round(
            mutation_rate * len(seed_sequence)
        )
        candidate.metadata["mutation_event_budget"] = event_budget
        candidate.metadata["mutation_count"] = (
            event_budget if mutation_guidance_mode != "legacy" else len(changed_positions)
        )
        candidate.metadata["mutation_positions"] = changed_positions
        candidate.metadata["parent_sequence"] = seed_sequence
        candidate.metadata["run_id"] = run_id
        candidate.metadata["generation"] = generation
        candidate.metadata["hotspot_guided"] = policy["guided"]
        candidate.metadata["mutation_guidance_mode"] = policy["mode"]
        candidate.metadata["requested_mutation_guidance_mode"] = (
            policy["requested_mode"]
        )
        candidate.metadata["hotspot_discovery_generation"] = policy[
            "hotspot_discovery_generation"
        ]
        candidate.metadata["hotspot_start"] = policy["hotspot_start"]
        candidate.metadata["hotspot_end"] = policy["hotspot_end"]
        candidate.metadata["is_duplicate"] = False
        candidate.metadata["preserved_elite"] = (
            is_preserved
        )
        candidates.append(candidate)

    logger.info(
        "Generated %d candidates (%d unique).",
        len(candidates),
        len({candidate.sequence for candidate in candidates}),
    )

    return candidates
# ==========================================================
# Debug
# ==========================================================


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)

    seed = "ACDEFGHIKLMNPQRSTVWYA"

    population = generate_candidates(seed)

    print()

    for candidate in population:

        print("--------------------------------------")
        print(candidate.sequence)

    print()
    print("ARISE importance map:")
    print(ARISE_ENGINE.importance_map())
