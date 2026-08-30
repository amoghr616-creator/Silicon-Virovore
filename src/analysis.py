"""
analysis.py

High-level analysis stage for Silicon Virovore.

Aggregates all scientific analyses into one PipelineReport.
"""

from __future__ import annotations

import logging

from src.models import (
    Candidate,
    PipelineReport,
    RankedCandidate,
)

from src.diversity import analyze_diversity
from src.correlation import analyze_correlations
from src.pareto import compute_pareto_front
from src.bootstrap import bootstrap_rankings
from src.validation import validate_candidates

logger = logging.getLogger(__name__)


class AnalysisEngine:

    def analyze(
        self,
        candidates: list[Candidate],
        runtime_seconds: float,
    ) -> PipelineReport:

        logger.info("Running diversity analysis...")
        diversity = analyze_diversity(candidates)

        logger.info("Running correlation analysis...")
        correlations = analyze_correlations(candidates)

        logger.info("Computing Pareto front...")
        pareto = compute_pareto_front(candidates)

        logger.info("Running bootstrap stability...")
        bootstrap = bootstrap_rankings(candidates)

        logger.info("Running validation...")
        validation = validate_candidates(candidates)

        ranked = [
            RankedCandidate(
                candidate=c,
                overall_score=c.overall_score,
                confidence=c.confidence,
                rank=c.rank,
            )
            for c in candidates
        ]

        report = PipelineReport(
            top_candidates=ranked,
            diversity=diversity,
            correlations=correlations,
            pareto_front=pareto,
            bootstrap=bootstrap,
            validation=validation,
            runtime_seconds=runtime_seconds,
        )

        logger.info("Analysis complete.")

        return report