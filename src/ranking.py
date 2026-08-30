"""
ranking.py

Multi-objective ranking engine for Silicon Virovore.

Combines every computational stage into one overall score.
"""

from __future__ import annotations

from src.models import Candidate


class CandidateRanker:

    def __init__(
        self,
        fitness_weight=0.30,
        docking_weight=0.30,
        structure_weight=0.10,
        helix_weight=0.08,
        hydro_weight=0.08,
        solvation_weight=0.06,
        md_weight=0.05,
        safety_weight=0.03,
    ):

        self.weights = {
            "fitness": fitness_weight,
            "docking": docking_weight,
            "structure": structure_weight,
            "helix": helix_weight,
            "hydrophobicity": hydro_weight,
            "solvation": solvation_weight,
            "md": md_weight,
            "safety": safety_weight,
        }

    # --------------------------------------------------------

    def normalize(self, values, reverse=False):

        if not values:
            return []

        minimum = min(values)
        maximum = max(values)

        if minimum == maximum:
            return [1.0] * len(values)

        normalized = [
            (v - minimum) / (maximum - minimum)
            for v in values
        ]

        if reverse:
            normalized = [1.0 - x for x in normalized]

        return normalized

    # --------------------------------------------------------

    def rank(
        self,
        candidates: list[Candidate],
    ) -> list[Candidate]:

        fitness = [c.c_score for c in candidates]

        docking = [
            c.strongest_anchor_delta_g
            if c.strongest_anchor_delta_g is not None
            else 0
            for c in candidates
        ]

        structure = [
            c.structure_confidence or 0
            for c in candidates
        ]

        helix = [
            c.helix_propensity
            for c in candidates
        ]

        hydro = [
            c.hydrophobic_moment
            for c in candidates
        ]

        solvation = [
            c.solvation_energy
            for c in candidates
        ]

        md = [
            c.md_stability_score or 0
            for c in candidates
        ]

        safety = [
            c.toxicity_score or 0
            for c in candidates
        ]

        fitness = self.normalize(fitness)

        docking = self.normalize(
            docking,
            reverse=True,
        )

        structure = self.normalize(structure)

        helix = self.normalize(helix)

        hydro = self.normalize(hydro)

        solvation = self.normalize(
            solvation,
            reverse=True,
        )

        md = self.normalize(md)

        safety = self.normalize(
            safety,
            reverse=True,
        )

        for i, candidate in enumerate(candidates):

            breakdown = {}

            breakdown["fitness"] = (
                fitness[i] * self.weights["fitness"]
            )

            breakdown["docking"] = (
                docking[i] * self.weights["docking"]
            )

            breakdown["structure"] = (
                structure[i] * self.weights["structure"]
            )

            breakdown["helix"] = (
                helix[i] * self.weights["helix"]
            )

            breakdown["hydrophobicity"] = (
                hydro[i] * self.weights["hydrophobicity"]
            )

            breakdown["solvation"] = (
                solvation[i] * self.weights["solvation"]
            )

            breakdown["md"] = (
                md[i] * self.weights["md"]
            )

            breakdown["safety"] = (
                safety[i] * self.weights["safety"]
            )

            candidate.ranking_breakdown = breakdown

            candidate.overall_score = sum(
                breakdown.values()
            )

        candidates.sort(
            key=lambda c: c.overall_score,
            reverse=True,
        )

        for rank, candidate in enumerate(
            candidates,
            start=1,
        ):
            candidate.rank = rank

        return candidates


if __name__ == "__main__":

    candidates = []

    for i in range(5):

        c = Candidate(
            sequence=f"SEQ{i}",
            c_score=i + 5,
        )

        c.strongest_anchor_delta_g = -6 - i

        c.helix_propensity = i * 0.3

        c.hydrophobic_moment = i * 0.4

        c.solvation_energy = -10 + i

        candidates.append(c)

    ranked = CandidateRanker().rank(candidates)

    for c in ranked:

        print(
            c.rank,
            c.sequence,
            round(c.overall_score, 3),
        )