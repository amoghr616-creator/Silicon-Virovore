"""
ranking.py

Multi-objective ranking engine for Silicon Virovore.

Combines every computational stage into one overall score.
"""

from __future__ import annotations
from asyncio.log import logger

from src.models import Candidate
import logging
logger = logging.getLogger(__name__)

class CandidateRanker:

    def __init__(
        self,
        fitness_weight=0.30,
        docking_weight=0.30,
        structure_weight=0.10,
        helix_weight=0.08,
        hydro_weight=0.0,
        solvation_weight=0.06,
        md_weight=0.05,
        safety_weight=0.0,
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

    def normalize(
        self,
        values,
        reverse=False,
    ):
        """
        Min-max normalize available values to [0, 1].

        Missing values remain missing and are represented by None.
        A metric with no variation among available observations
        receives 0.5 rather than falsely becoming perfect.
        """

        if not values:
            return []

        cleaned = [
            None if value is None else float(value)
            for value in values
        ]

        available = [
            value for value in cleaned
            if value is not None
        ]

        if not available:
            return [None] * len(cleaned)

        minimum = min(available)
        maximum = max(available)

        if abs(maximum - minimum) < 1e-9:
            normalized = [
                0.5 if value is not None else None
                for value in cleaned
            ]
        else:
            normalized = [
                None
                if value is None
                else (value - minimum) / (maximum - minimum)
                for value in cleaned
            ]

        if reverse:
            normalized = [
                None
                if value is None
                else 1.0 - value
                for value in normalized
            ]

        return normalized
    # --------------------------------------------------------

    @staticmethod
    def _safe(
        value,
        default=0.0,
    ):

        if value is None:
            return default

        return value

    # --------------------------------------------------------

    @staticmethod
    def _metadata(
        candidate: Candidate,
        key: str,
        default=0.0,
    ):

        return candidate.metadata.get(
            key,
            default,
        )

    # --------------------------------------------------------

    @staticmethod
    def _confidence(
        candidate: Candidate,
    ) -> float:

        confidence = 0.0

        if candidate.structure_confidence is not None:

            confidence += (
                candidate.structure_confidence / 100.0
            ) * 0.40

        if candidate.vina_delta_g is not None:

            confidence += 0.20

        if candidate.md_stability_score is not None:

            confidence += 0.20

        if candidate.toxicity_score is not None:

            confidence += 0.20

        return round(
            min(confidence, 1.0),
            3,
        )

    def rank(
        self,
        candidates: list[Candidate],
    ) -> list[Candidate]:

        if not candidates:
            return []

        # ----------------------------------------------------
        # Extract raw metrics
        # ----------------------------------------------------

        fitness = [
            self._safe(c.c_score)
            for c in candidates
        ]

        docking = [
            self._metadata(
                c,
                "consensus_docking",
                self._safe(
                    c.strongest_anchor_delta_g
                ),
            )
            for c in candidates
        ]

        docking_consistency = [
            self._metadata(
                c,
                "docking_std",
                0.0,
            )
            for c in candidates
        ]

        structure = [
        
            c.structure_confidence
            for c in candidates
        ]

        helix = [
            self._safe(
                c.helix_propensity
            )
            for c in candidates
        ]

        hydro = [
            self._safe(
                c.hydrophobic_moment
            )
            for c in candidates
        ]

        solvation = [
            self._safe(
                c.solvation_energy
            )
            for c in candidates
        ]

        md = [
        
            c.md_stability_score
            for c in candidates
        ]

        safety = [
          
            c.toxicity_score
            for c in candidates
        ]

        # ----------------------------------------------------
        # Normalize metrics
        # ----------------------------------------------------

        fitness_norm = self.normalize(
            fitness
        )

        docking_norm = self.normalize(
            docking
        )

        structure_norm = self.normalize(
            structure
        )

        helix_norm = self.normalize(
            helix
        )

        hydro_norm = self.normalize(
            hydro
        )

        solvation_norm = self.normalize(
            solvation,
            reverse=True,
        )

        md_norm = self.normalize(
            md
        )

        safety_norm = self.normalize(
            safety,
            reverse=True,
        )

        # Lower docking standard deviation is better.
        consistency_norm = self.normalize(
            docking_consistency,
            reverse=True,
        )

        # ----------------------------------------------------
        # Candidate scoring
        # ----------------------------------------------------

        for i, candidate in enumerate(candidates):

            breakdown = {}

            def weighted_component(
                normalized_value,
                weight,
            ):
                if normalized_value is None:
                    return 0.0
                return normalized_value * weight

            breakdown["fitness"] = weighted_component(
                fitness_norm[i],
                self.weights["fitness"],
            )

            breakdown["docking"] = weighted_component(
                docking_norm[i],
                self.weights["docking"],
            )

            breakdown["structure"] = weighted_component(
                structure_norm[i],
                self.weights["structure"],
            )

            breakdown["helix"] = weighted_component(
                helix_norm[i],
                self.weights["helix"],
            )

            breakdown["hydrophobicity"] = weighted_component(
                hydro_norm[i],
                self.weights["hydrophobicity"],
            )

            breakdown["solvation"] = weighted_component(
                solvation_norm[i],
                self.weights["solvation"],
            )

            breakdown["md"] = weighted_component(
                md_norm[i],
                self.weights["md"],
            )

            breakdown["safety"] = weighted_component(
                safety_norm[i],
                self.weights["safety"],
            )

            candidate.ranking_breakdown = breakdown

            candidate.metadata["ranking_total"] = round(
                sum(breakdown.values()),
                6,
            )

            candidate.metadata["ranking_components"] = {
                key: round(value, 6)
                for key, value in breakdown.items()
            }

            available_weight = 0.0
            weighted_total = 0.0

            for key, value in breakdown.items():

                if value is None:
                    continue

                weight = self.weights[key]

                if weight <= 0.0:
                    continue

                weighted_total += value
                available_weight += weight

            if available_weight > 0.0:
                candidate.overall_score = round(
                    weighted_total / available_weight,
                    6,
                )
                logger.debug(
    "RANK %s | score=%.4f | breakdown=%s",
    candidate.sequence,
    candidate.overall_score,
    {
        k: round(v, 4)
        for k, v in breakdown.items()
    },
)
            else:
                candidate.overall_score = 0.0

            candidate.confidence = self._confidence(candidate)

        # --------------------------------------------------------
        # Sort
        # --------------------------------------------------------

        candidates.sort(
            key=lambda c: (
                c.overall_score,
                c.confidence,
            ),
            reverse=True,
        )

        # --------------------------------------------------------
        # Assign ranks
        # --------------------------------------------------------

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