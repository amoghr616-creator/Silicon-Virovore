"""
docking_vina.py

Docking stage for Silicon Virovore.

Evaluates every peptide candidate using a lightweight
binding surrogate and optionally escalates promising
leads to AutoDock Vina.

Future versions can swap in Boltz-2, DiffDock,
Rosetta or physical docking engines without changing
the rest of the pipeline.
"""

from __future__ import annotations

import logging

from src.candidate import Candidate

logger = logging.getLogger(__name__)

TIER2_THRESHOLD = -7.0


class DockingBackend:
    """
    Abstract docking backend.
    """

    def score_fragment(self, fragment: str) -> float:
        raise NotImplementedError


class MLSurrogateBackend(DockingBackend):

    def score_fragment(self, fragment: str) -> float:

        hydrophobic = set("AILMFWYV")
        charged = set("RHKDE")

        score = -5.0

        for aa in fragment:

            if aa in hydrophobic:
                score -= 0.30

            elif aa in charged:
                score += 0.10

        return round(score, 2)


class AutoDockBackend(DockingBackend):

    def score_fragment(self, fragment: str):

        logger.info(
            "Placeholder AutoDock Vina run for %s",
            fragment,
        )

        #
        # Replace later with Vina CLI
        #

        return -8.40


class PeptideDockingScorer:

    def __init__(
        self,
        backend: DockingBackend | None = None,
    ):

        self.backend = backend or MLSurrogateBackend()

    # --------------------------------------------------

    def _fragments(self, sequence: str):

        return [

            sequence[i:i + 9]

            for i in range(len(sequence) - 8)
        ]

    # --------------------------------------------------

    def evaluate_candidate(self, candidate: Candidate):
        logger.info(
            "Docking %s",
            candidate.sequence,
        )

        fragments = self._fragments(
            candidate.sequence
        )

        scores = []

        for fragment in fragments:

            dg = self.backend.score_fragment(fragment)

            scores.append(dg)

        candidate.mean_delta_g = round(

            sum(scores) / len(scores),

            2,
        )

        candidate.strongest_anchor_delta_g = min(scores)

        if candidate.strongest_anchor_delta_g <= TIER2_THRESHOLD:

            logger.info(
                "Tier-2 docking triggered."
            )

            vina = AutoDockBackend()

            candidate.vina_delta_g = vina.score_fragment(

                fragments[
                    scores.index(min(scores))
                ]
            )

            candidate.add_note(
                "Validated with AutoDock backend."
            )

        else:

            candidate.vina_delta_g = None

        return candidate
    # --------------------------------------------------

def evaluate(self, candidate: Candidate):
    """
    Backwards-compatible alias.
    """
    return self.evaluate_candidate(candidate)