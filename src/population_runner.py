"""
population_runner.py

Generates candidate peptides using the native C backend and
implements the Adaptive Recursive Intelligent Selection Engine (ARISE).

ARISE learns which residue positions consistently contribute to
high-performing peptides and biases future generations toward
those positions while preserving diversity.
"""

from __future__ import annotations

import logging

from c.bridge import (
    generate_c_population,
    generate_adaptive_population,
    process_candidate_peptide,
)

from src.models import Candidate

from src.config import (
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

        if not ranked_candidates:
            return

        sequence_length = len(ranked_candidates[0].sequence)

        for position in range(sequence_length):

            residue_scores: dict[str, list[float]] = {}

            for candidate in ranked_candidates:

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

            # Need at least two distinct residue groups
            # to estimate a positional effect.
            if len(residue_scores) < 2:
                continue

            total_observations = sum(
                len(scores)
                for scores in residue_scores.values()
            )

            if total_observations < 5:
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
        Return normalized position importance.
        """

        return {
            position: round(
                value,
                3,
            )
            for position, value
            in sorted(
                self._importance.items()
            )
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
) -> list[Candidate]:
    """
    Generate candidate peptides.

    Generation 1 uses the standard native C mutation engine.

    Later generations use ARISE-derived positional importance
    through the native adaptive population generator.

    ARISE therefore changes the mutation distribution without
    performing a second Python-side mutation pass.
    """

    logger.info("Generating peptide candidates...")

    importance_map = ARISE_ENGINE.get_importance_array(
        len(seed_sequence)
    )

    has_learned_importance = bool(
        ARISE_ENGINE._importance
    )
    logger.info(
    "ARISE mode enabled | importance=%s",
    [
        round(x, 3)
        for x in importance_map
    ],
)
    if has_learned_importance:

        sequences = generate_adaptive_population(
            seed_sequence,
            importance_map,
            pop_size=population_size,
        )

    else:

        sequences = generate_c_population(
            seed_sequence,
            pop_size=population_size,
            mutation_rate=mutation_rate,
        )

    candidates: list[Candidate] = []

    for sequence in sequences:

        candidate = process_candidate_peptide(
            sequence
        )

        if candidate is None:
            continue

        candidates.append(candidate)

    logger.info(
        "Generated %d candidates.",
        len(candidates),
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