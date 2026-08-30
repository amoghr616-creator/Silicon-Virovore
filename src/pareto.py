"""
pareto.py

Multi-objective optimization for Silicon Virovore.

Instead of selecting candidates solely by a weighted
overall score, Pareto optimization identifies peptides
that cannot be improved in one objective without
becoming worse in another.

Future versions:
    - NSGA-II ranking
    - Crowding distance
    - Hypervolume indicator
    - Pareto evolution tracking
"""

from __future__ import annotations

from dataclasses import dataclass


# ------------------------------------------------------------

@dataclass(slots=True)
class ParetoCandidate:

    candidate: object

    objectives: tuple

    dominates: int = 0

    dominated_by: int = 0

    front: int = 1


# ------------------------------------------------------------


class ParetoAnalyzer:

    """
    Computes Pareto-optimal candidate sets.
    """

    # --------------------------------------------------------

    @staticmethod
    def objective_vector(candidate):

        """
        Objective ordering:

        maximize:
            fitness
            helix
            confidence

        minimize:
            docking ΔG
            solvation penalty
        """

        return (

            candidate.c_score,

            candidate.helix_propensity,

            getattr(candidate, "confidence", 1.0),

            -candidate.strongest_anchor_delta_g,

            -candidate.solvation_energy,
        )

    # --------------------------------------------------------

    @staticmethod
    def dominates(a, b):

        """
        Returns True if candidate A dominates B.
        """

        better_or_equal = all(
            x >= y
            for x, y in zip(a, b)
        )

        strictly_better = any(
            x > y
            for x, y in zip(a, b)
        )

        return better_or_equal and strictly_better

    # --------------------------------------------------------

    def compute_front(self, candidates):

        pareto = [

            ParetoCandidate(

                candidate=c,

                objectives=self.objective_vector(c),

            )

            for c in candidates

        ]

        for i, A in enumerate(pareto):

            for j, B in enumerate(pareto):

                if i == j:
                    continue

                if self.dominates(

                    A.objectives,

                    B.objectives,

                ):

                    A.dominates += 1

                elif self.dominates(

                    B.objectives,

                    A.objectives,

                ):

                    A.dominated_by += 1

        front = [

            p

            for p in pareto

            if p.dominated_by == 0

        ]

        return front

    # --------------------------------------------------------

    def rank(self, candidates):

        front = self.compute_front(candidates)

        ranked = sorted(

            front,

            key=lambda x: (

                x.candidate.c_score,

                -x.candidate.strongest_anchor_delta_g,

                x.candidate.helix_propensity,

            ),

            reverse=True,

        )

        return ranked

    # --------------------------------------------------------

    def summary(self, candidates):

        front = self.rank(candidates)

        return {

            "num_candidates":

                len(candidates),

            "pareto_front_size":

                len(front),

            "pareto_sequences":

                [

                    p.candidate.sequence

                    for p in front

                ],

            "front":

                front,
        }