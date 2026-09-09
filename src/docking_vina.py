"""
docking_vina.py

Docking stage for Silicon Virovore.

Evaluates peptide candidates using a lightweight surrogate model
and optionally escalates promising candidates to AutoDock Vina.

This module is intentionally backend-agnostic so future docking
engines (DiffDock, Boltz-2, Rosetta, GNINA, etc.) can be swapped
without changing the remainder of the pipeline.
"""

from __future__ import annotations

import logging
import statistics

from src.models import Candidate

logger = logging.getLogger(__name__)

TIER2_THRESHOLD = -7.0


# ============================================================
# Docking Backends
# ============================================================

class DockingBackend:
    """
    Abstract docking backend.
    """

    def score_fragment(self, fragment: str) -> float:
        raise NotImplementedError


class MLSurrogateBackend(DockingBackend):
    """
    Lightweight heuristic docking surrogate.

    Used during evolutionary optimization because it is
    several orders of magnitude faster than physical docking.
    """

    HYDROPHOBIC = set("AILMFWYV")
    CHARGED = set("RHKDE")

    def score_fragment(self, fragment: str) -> float:

        score = -5.0

        for aa in fragment:

            if aa in self.HYDROPHOBIC:
                score -= 0.30

            elif aa in self.CHARGED:
                score += 0.10

        return round(score, 2)


class AutoDockBackend(DockingBackend):
    """
    Placeholder AutoDock Vina backend.

    Replace score_fragment() with a real Vina invocation later.
    """

    def score_fragment(self, fragment: str) -> float:

        logger.info(
            "Running AutoDock Vina placeholder for %s",
            fragment,
        )

        #
        # TODO:
        # Launch AutoDock Vina
        # Parse docking score
        #

        return -8.40


# ============================================================
# Docking Engine
# ============================================================

class PeptideDockingScorer:
    """
    Performs fragment-based docking evaluation.

    Computes:

    • Best fragment ΔG
    • Mean ΔG
    • Docking consistency
    • Consensus docking score
    • Optional Tier-2 AutoDock validation
    """

    def __init__(
        self,
        backend: DockingBackend | None = None,
    ):

        self.backend = backend or MLSurrogateBackend()

        self.vina_backend = AutoDockBackend()

    # --------------------------------------------------------

    @staticmethod
    def _fragments(sequence: str) -> list[str]:
        """
        Generate overlapping 9-mer fragments.
        """

        if len(sequence) <= 9:
            return [sequence]

        return [
            sequence[i:i + 9]
            for i in range(len(sequence) - 8)
        ]

    # --------------------------------------------------------

    @staticmethod
    def _consensus_score(
        best_delta_g: float,
        mean_delta_g: float,
    ) -> float:
        """
        Consensus docking score.

        Combines strongest interaction with overall
        fragment performance.
        """

        return round(
            (
                0.60 * abs(best_delta_g)
                + 0.40 * abs(mean_delta_g)
            ),
            3,
        )

    # --------------------------------------------------------

    def evaluate_candidate(
        self,
        candidate: Candidate,
    ) -> Candidate:

        logger.info(
            "Docking %s",
            candidate.sequence,
        )

        fragments = self._fragments(
            candidate.sequence
        )

        scores: list[float] = []

        candidate.fragment_scores.clear()

        # ----------------------------------------------------
        # Score every fragment
        # ----------------------------------------------------

        for fragment in fragments:

            delta_g = self.backend.score_fragment(
                fragment
            )

            scores.append(delta_g)

            candidate.fragment_scores[
                fragment
            ] = delta_g

        candidate.fragments = fragments
        candidate.docking_scores = scores

        # ----------------------------------------------------
        # Summary statistics
        # ----------------------------------------------------

        candidate.mean_delta_g = round(
            sum(scores) / len(scores),
            2,
        )

        candidate.strongest_anchor_delta_g = min(scores)

        best_index = scores.index(
            candidate.strongest_anchor_delta_g
        )

        candidate.best_fragment = fragments[
            best_index
        ]

        docking_std = (
            statistics.stdev(scores)
            if len(scores) > 1
            else 0.0
        )

        candidate.metadata[
            "docking_std"
        ] = docking_std

        candidate.metadata[
            "consensus_docking"
        ] = self._consensus_score(
            candidate.strongest_anchor_delta_g,
            candidate.mean_delta_g,
        )

        # ----------------------------------------------------
        # Tier-2 Validation
        # ----------------------------------------------------

        if (
            candidate.strongest_anchor_delta_g
            <= TIER2_THRESHOLD
        ):

            logger.info(
                "Tier-2 docking triggered."
            )

            candidate.passed_tier_2 = True

            candidate.vina_delta_g = (
                self.vina_backend.score_fragment(
                    candidate.best_fragment
                )
            )

            candidate.add_note(
                "Validated using AutoDock backend."
            )

        else:

            candidate.passed_tier_2 = False
            candidate.vina_delta_g = None

        logger.info(
            (
                "Docking complete | "
                "Best %.2f | "
                "Mean %.2f | "
                "Consensus %.2f"
            ),
            candidate.strongest_anchor_delta_g,
            candidate.mean_delta_g,
            candidate.metadata[
                "consensus_docking"
            ],
        )

        return candidate

    # --------------------------------------------------------

    def evaluate(
        self,
        candidate: Candidate,
    ) -> Candidate:
        """
        Backwards-compatible alias.
        """

        return self.evaluate_candidate(
            candidate
        )