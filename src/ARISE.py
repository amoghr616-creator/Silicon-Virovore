"""
arise.py

Adaptive Residue Importance Scoring Engine (ARISE)

Phase I

ARISE observes each completed generation of peptide candidates,
learns which residue positions are consistently associated with
high-performing peptides, and produces a Residue Importance Score
(RIS) for every position.

Future phases will add:

    • Adaptive mutation memory
    • Epistasis detection
    • Evolutionary confidence
    • Realism estimation
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

from src.models import Candidate


# ============================================================
# Position Statistics
# ============================================================

@dataclass(slots=True)
class PositionStatistics:

    observations: int = 0
    cumulative_score: float = 0.0

    importance: float = 0.5

    mutation_count: int = 0
    successful_mutations: int = 0

    conservation: float = 0.0


# ============================================================
# ARISE Engine
# ============================================================

class ARISE:

    def __init__(self):

        self.reference_sequence: str | None = None

        self.sequence_length = 0

        self.generation = 0

        self.positions: dict[int, PositionStatistics] = defaultdict(
            PositionStatistics
        )

    # ========================================================
    # Observe Generation
    # ========================================================

    def observe_generation(
        self,
        candidates: list[Candidate],
    ) -> None:

        if not candidates:
            return

        self.generation += 1

        if self.reference_sequence is None:

            self.reference_sequence = candidates[0].sequence

            self.sequence_length = len(
                self.reference_sequence
            )

        best = max(
            c.overall_score
            for c in candidates
        )

        worst = min(
            c.overall_score
            for c in candidates
        )

        score_range = max(
            best - worst,
            1e-6,
        )

        for candidate in candidates:

            normalized = (

                candidate.overall_score - worst

            ) / score_range

            for position, residue in enumerate(
                candidate.sequence
            ):

                stats = self.positions[position]

                if residue == self.reference_sequence[position]:

                    stats.observations += 1
                    stats.cumulative_score += normalized

    # ========================================================
    # Compute Residue Importance
    # ========================================================

    def update_importance(self) -> None:

        values = []

        for stats in self.positions.values():

            if stats.observations == 0:

                stats.importance = 0.0

            else:

                stats.importance = (

                    stats.cumulative_score
                    / stats.observations

                )

            values.append(
                stats.importance
            )

        if not values:
            return

        minimum = min(values)
        maximum = max(values)

        scale = max(
            maximum - minimum,
            1e-6,
        )

        for stats in self.positions.values():

            stats.importance = (

                stats.importance - minimum

            ) / scale

    # ========================================================
    # Mutation Memory
    # ========================================================

    def record_mutation(
        self,
        position: int,
        improved: bool,
    ) -> None:

        stats = self.positions[position]

        stats.mutation_count += 1

        if improved:

            stats.successful_mutations += 1

    # ========================================================
    # Mutation Statistics
    # ========================================================

    def mutation_success_rate(
        self,
        position: int,
    ) -> float:

        stats = self.positions[position]

        if stats.mutation_count == 0:

            return 0.5

        return (

            stats.successful_mutations

            / stats.mutation_count

        )

    # ========================================================
    # Adaptive Mutation Policy
    # ========================================================

    def mutation_probability(
        self,
        position: int,
        base_rate: float = 0.05,
    ) -> float:

        importance = self.positions[position].importance

        probability = base_rate * (

            2.0 - importance

        )

        probability = max(

            0.01,

            min(
                0.50,
                probability,
            ),
        )

        return probability

    def mutation_rate_map(
        self,
        base_rate: float = 0.05,
    ) -> list[float]:
        """
        Returns the mutation probability for every
        residue position.

        This will be passed to the native C engine.
        """

        if self.reference_sequence is None:

            return []

        return [
        self.mutation_probability(i, base_rate)
        for i in range(self.sequence_length)
    ]

    # ========================================================
    # Queries
    # ========================================================

    def importance_map(self) -> list[float]:

        if self.reference_sequence is None:

            return []

        return [

            round(
                self.positions[i].importance,
                4,
            )

            for i in range(
                self.sequence_length
            )

        ]

    def mutation_success_rates(self):

        return {

            position: round(

                self.mutation_success_rate(position),

                3,

            )

            for position in range(
                self.sequence_length
            )

        }

    # ========================================================
    # Summary
    # ========================================================

    def summary(self):

        return {

            "generation": self.generation,

            "sequence_length": self.sequence_length,

            "importance_map":
                self.importance_map(),

            "mutation_success":
                self.mutation_success_rates(),

        }