"""
analysis.py

High-level analysis stage for Silicon Virovore.

Aggregates all scientific analyses into one PipelineReport.
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path

from src.models import (
    Candidate,
    PipelineReport,
    RankedCandidate,
)

from src.diversity import DiversityAnalyzer
from src.correlation import CorrelationAnalyzer
from src.pareto import ParetoAnalyzer
from src.bootstrap import BootstrapAnalyzer
from src.validation import Validator

logger = logging.getLogger(__name__)


class AnalysisEngine:

    def analyze(
        self,
        candidates: list[Candidate],
        runtime_seconds: float,
        structure_directory: str | Path = "structures",
    ) -> PipelineReport:

        if not candidates:
            raise ValueError("No candidates supplied to AnalysisEngine.")

        logger.info("Preparing ranked candidates...")

        ranked = sorted(
            [
                RankedCandidate(
                    candidate=c,
                    overall_score=c.overall_score,
                    confidence=c.confidence,
                    rank=c.rank,
                )
                for c in candidates
            ],
            key=lambda rc: rc.rank,
        )

        # --------------------------------------------------
        # Diversity
        # --------------------------------------------------

        logger.info("Running diversity analysis...")

        sequences = [c.sequence for c in candidates]

        reference = ranked[0].candidate.sequence

        diversity = DiversityAnalyzer().summary(
            sequences=sequences,
            reference=reference,
        )

        # --------------------------------------------------
        # Correlation
        # --------------------------------------------------

        logger.info("Running correlation analysis...")

        correlations = CorrelationAnalyzer().summary(
            ranked
        )

        # --------------------------------------------------
        # Pareto
        # --------------------------------------------------

        logger.info("Computing Pareto front...")

        pareto = ParetoAnalyzer().summary(
            candidates
        )

        # --------------------------------------------------
        # Bootstrap
        # --------------------------------------------------

        logger.info("Running bootstrap analysis...")

        bootstrap = BootstrapAnalyzer().summary(
            candidates
        )

        # --------------------------------------------------
        # Validation
        # --------------------------------------------------

        logger.info("Running validation...")

        validation = Validator().validate(
            candidates,
            structure_directory,
        )

        # --------------------------------------------------
        # Report
        # --------------------------------------------------

        report = PipelineReport(
            top_candidates=ranked,
            diversity=diversity,
            correlations=correlations,
            pareto_front=pareto.get("front", []),
            bootstrap=bootstrap.get("results", []),
            validation=validation,
            runtime_seconds=runtime_seconds,
            timestamp=datetime.now().isoformat(),
        )

        logger.info("Analysis complete.")

        return report