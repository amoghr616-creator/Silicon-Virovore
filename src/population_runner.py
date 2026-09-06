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
import random

from c.bridge import (
    generate_c_population,
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

    Learns positional importance from successful peptides.
    """

    def __init__(self):

        self.position_scores: dict[int, float] = {}

    def observe_generation(
        self,
        ranked_candidates: list[Candidate],
    ) -> None:
        """
        Learn from the best candidates in the generation.
        """

        if not ranked_candidates:
            return

        top = ranked_candidates[:10]

        for candidate in top:

            score = candidate.overall_score

            for i, residue in enumerate(candidate.sequence):

                self.position_scores.setdefault(i, 0.0)

                self.position_scores[i] += score

    def update_importance(self) -> None:
        """
        Normalize importance values.
        """

        if not self.position_scores:
            return

        maximum = max(self.position_scores.values())

        if maximum <= 0:
            return

        for position in self.position_scores:

            self.position_scores[position] /= maximum

    def importance_map(self) -> dict[int, float]:

        return {
            k: round(v, 3)
            for k, v in sorted(self.position_scores.items())
        }

    def bias_sequence(
        self,
        sequence: str,
    ) -> str:
        """
        Preserve residues at positions ARISE believes
        are highly important.
        """

        if not self.position_scores:
            return sequence

        seq = list(sequence)

        amino_acids = "ACDEFGHIKLMNPQRSTVWY"

        for i in range(len(seq)):

            importance = self.position_scores.get(i, 0)

            # Highly important positions mutate less often.
            mutation_probability = max(
                0.05,
                1.0 - importance,
            )

            if random.random() < mutation_probability:

                seq[i] = random.choice(amino_acids)

        return "".join(seq)


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
    Generate peptide candidates.

    The native C backend performs the primary mutation.

    Afterwards ARISE lightly biases candidates according to
    learned positional importance.
    """

    logger.info("Generating peptide candidates...")

    sequences = generate_c_population(
        seed_sequence,
        pop_size=population_size,
        mutation_rate=mutation_rate,
    )

    candidates: list[Candidate] = []

    for sequence in sequences:

        # Apply ARISE bias after native mutation
        sequence = ARISE_ENGINE.bias_sequence(sequence)

        candidate = process_candidate_peptide(sequence)

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