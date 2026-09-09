"""
docking_vina.py

Docking stage for Silicon Virovore.

Evaluates peptide candidates using a lightweight surrogate model
and optionally escalates promising candidates to AutoDock Vina.

This module is intentionally backend-agnostic so future docking
engines (DiffDock, Boltz-2, Rosetta, GNINA, etc.) can be swapped
without changing the remainder of the pipeline.
"""

import logging
import os
import re
import shutil
import statistics
import subprocess
import tempfile
from pathlib import Path

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

    This backend must remain fast and must not invoke
    AutoDock Vina. It is used for Tier-1 screening.
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
    Optional structure-validation backend.

    This class is intentionally disabled unless a valid executable,
    receptor, and candidate structure are available.
    """

    def __init__(self):

        self.vina = Path(
            "/Users/centurion616/Desktop/Silicon-Virovore/"
            "vina_1.2.7_mac_aarch64"
        )

        self.obabel = (
            shutil.which("obabel")
            or shutil.which("babel")
        )

        self.receptor = Path(
            "data/receptor/receptor.pdbqt"
        )

        if not self.vina.exists():
            raise RuntimeError(
                f"Vina executable not found: {self.vina}"
            )

        if not os.access(self.vina, os.X_OK):
            raise RuntimeError(
                f"Vina executable is not executable: {self.vina}"
            )

        if self.obabel is None:
            raise RuntimeError(
                "Open Babel executable not found."
            )

        if not self.receptor.exists():
            raise RuntimeError(
                f"Missing receptor: {self.receptor}"
            )

    def score_fragment(
        self,
        fragment: str,
        structure_path: str | Path | None = None,
    ) -> float:
        """
        Placeholder for structure-based validation.

        Do not silently substitute the Tier-1 surrogate here.
        """

        if structure_path is None:
            raise RuntimeError(
                "Structure path required for validation."
            )

        raise NotImplementedError(
            "Structure-based docking execution is not "
            "implemented in this backend yet."
        )
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

        try:

            self.vina_backend = AutoDockBackend()

        except Exception as exc:

            logger.warning(
                "AutoDock disabled: %s",
                exc,
            )

            self.vina_backend = None

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
        Composite surrogate docking score.

        Lower (more negative) remains better, matching the
        sign convention of the underlying docking scores.
        """

        return round(
            (
                0.60 * best_delta_g
                + 0.40 * mean_delta_g
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

        consensus_score = self._consensus_score(
            candidate.strongest_anchor_delta_g,
            candidate.mean_delta_g,
        )

        candidate.metadata[
            "consensus_docking"
        ] = consensus_score

        candidate.metadata[
            "docking_backend"
        ] = type(self.backend).__name__

        candidate.metadata[
            "docking_is_surrogate"
        ] = isinstance(
            self.backend,
            MLSurrogateBackend,
        )

        candidate.metadata[
            "docking_validated"
        ] = False

        # --------------------------------------------------------
        # Tier-2 Validation
        # --------------------------------------------------------

        if (
            self.vina_backend is not None
            and candidate.strongest_anchor_delta_g
            <= TIER2_THRESHOLD
        ):

            logger.info(
                "Tier-2 validation triggered."
            )

            if candidate.structure_path is None:

                logger.warning(
                    "No structure available for %s; "
                    "Tier-2 validation skipped.",
                    candidate.sequence,
                )

                candidate.passed_tier_2 = False
                candidate.vina_delta_g = None

                candidate.add_note(
                    "Tier-2 eligibility reached, but "
                    "validation was skipped because no "
                    "structure was available."
                )

            else:

                try:

                    candidate.vina_delta_g = (
                        self.vina_backend.score_fragment(
                            candidate.best_fragment,
                            candidate.structure_path,
                        )
                    )

                    candidate.passed_tier_2 = True

                    candidate.add_note(
                        "Tier-2 structure-based validation completed."
                    )

                except Exception as exc:

                    logger.warning(
                        "Tier-2 validation failed for %s: %s",
                        candidate.sequence,
                        exc,
                    )

                    candidate.passed_tier_2 = False
                    candidate.vina_delta_g = None

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
    consensus_score,
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