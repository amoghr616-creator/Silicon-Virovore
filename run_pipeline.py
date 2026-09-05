"""
run_pipeline.py

Master execution pipeline for Silicon Virovore.
"""

from __future__ import annotations

import logging
import time

from src.predict_structure import predict_population_structures
from src.docking_vina import PeptideDockingScorer
from src.ranking import CandidateRanker
from src.analysis import AnalysisEngine
from src.report import ReportGenerator
from src.plots import PlotGenerator
from src.population_runner import (
    generate_candidates,
    ARISE_ENGINE,
)
from src.config import (
    RESULTS_DIR,
    DEFAULT_SEED_SEQUENCE,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)

logger = logging.getLogger(__name__)


def run_pipeline():

    start = time.time()

    logger.info("===================================")
    logger.info("Silicon Virovore")
    logger.info("Pipeline Started")
    logger.info("===================================")

    # ------------------------------------------------------
    # Population Generation
    # ------------------------------------------------------

    logger.info("Generating candidates...")

    candidates = generate_candidates(
        DEFAULT_SEED_SEQUENCE,
    )



    # ------------------------------------------------------
    # Structure Prediction
    # ------------------------------------------------------

    logger.info("Predicting structures...")

    candidates = predict_population_structures(
        candidates
    )

    # ------------------------------------------------------
    # Docking
    # ------------------------------------------------------

    logger.info("Running docking...")

    scorer = PeptideDockingScorer()

    docked = []

    for candidate in candidates:
        docked.append(
            scorer.evaluate_candidate(candidate)
        )

    # ------------------------------------------------------
    # Ranking
    # ------------------------------------------------------

    logger.info("Ranking candidates...")

    ranked = CandidateRanker().rank(
        docked
    )

    # ------------------------------------
    # ARISE learns from this generation
    # ------------------------------------

    ARISE_ENGINE.observe_generation(ranked)
    ARISE_ENGINE.update_importance()

    logger.info("ARISE Importance Map")
    logger.info(ARISE_ENGINE.importance_map())

    logger.info("ARISE Mutation Rates")
    logger.info(ARISE_ENGINE.mutation_rate_map())
    # ------------------------------------------------------
    # Scientific Analysis
    # ------------------------------------------------------

    runtime = time.time() - start

    report = AnalysisEngine().analyze(
        ranked,
        runtime_seconds=runtime,
    )

    # ------------------------------------------------------
    # Figures
    # ------------------------------------------------------

    logger.info("Generating figures...")

    PlotGenerator(
        RESULTS_DIR,
    ).generate_all(
        ranked,
        [c.overall_score for c in ranked],
    )

    # ------------------------------------------------------
    # Report
    # ------------------------------------------------------

    logger.info("Generating report...")

    ReportGenerator(
        RESULTS_DIR,
    ).export(report)

    logger.info("===================================")
    logger.info("Pipeline Complete")
    logger.info("Runtime %.2f sec", runtime)
    logger.info("===================================")


if __name__ == "__main__":

    run_pipeline()