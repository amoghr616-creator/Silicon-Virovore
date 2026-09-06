"""
population_runner.py

Generates candidate peptides using the native C backend.
Returns Candidate objects.
"""

import logging

from c.bridge import (
    generate_c_population,
    process_candidate_peptide,
)

from src.config import (
    POPULATION_SIZE,
    MUTATION_RATE,
)


logger = logging.getLogger(__name__)


def generate_candidates(
    seed_sequence: str,
    population_size: int = POPULATION_SIZE,
    mutation_rate: float = MUTATION_RATE,
):
    """
    Generate a population of Candidate objects using the native
    Silicon Virovore mutation engine.
    """

    logger.info("Generating peptide candidates...")

    logger.info("Using standard random mutation engine...")

    sequences = generate_c_population(
        seed_sequence,
        pop_size=population_size,
        mutation_rate=mutation_rate,
    )

    candidates = []

    for sequence in sequences:

        candidate = process_candidate_peptide(sequence)

        if candidate is None:
            continue

        candidates.append(candidate)

    logger.info(
        "Successfully generated %d candidates.",
        len(candidates),
    )

    return candidates


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)

    seed = "ACDEFGHIKLMNPQRSTVWYA"

    population = generate_candidates(seed)

    print()

    for candidate in population:

        print("--------------------------------------")
        print(f"Sequence      : {candidate.sequence}")
        print(f"C Score       : {candidate.c_score:.4f}")

        if candidate.hydrophobic_moment is not None:
            print(f"Hydro Moment  : {candidate.hydrophobic_moment:.3f}")

        if candidate.helix_propensity is not None:
            print(f"Helix Score   : {candidate.helix_propensity:.3f}")

        if candidate.solvation_energy is not None:
            print(f"Solvation ΔG  : {candidate.solvation_energy:.3f}")