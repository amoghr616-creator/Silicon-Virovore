"""
run_pipeline.py

Master execution pipeline for Silicon Virovore.
"""

from __future__ import annotations

import logging
import time
from src.population_runner import (
    ARISE_ENGINE,
    generate_candidates,
)
from src.predict_structure import predict_population_structures
from src.docking_vina import PeptideDockingScorer
from src.ranking import CandidateRanker
from src.analysis import AnalysisEngine
from src.report import ReportGenerator
from src.plots import PlotGenerator
from src.config import (
    RESULTS_DIR,
    DEFAULT_SEED_SEQUENCE,
    GENERATIONS,
    TOP_K,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)

logger = logging.getLogger(__name__)


def run_pipeline():
    start = time.time()

    logger.info("=" * 50)
    logger.info("Silicon Virovore")
    logger.info("Pipeline Started")
    logger.info("=" * 50)

    seed_sequence = DEFAULT_SEED_SEQUENCE
    best_scores = []
    final_ranked = []

    scorer = PeptideDockingScorer()
    ranker = CandidateRanker()

    for generation in range(1, GENERATIONS + 1):
        logger.info("-" * 50)
        logger.info(
            "Generation %d / %d",
            generation,
            GENERATIONS,
        )
        logger.info("Seed: %s", seed_sequence)
        logger.info("-" * 50)

        ####################################################
        # Candidate Generation
        ####################################################
        candidates = generate_candidates(seed_sequence)

        if not candidates:
            logger.error("No candidates generated.")
            break

        ####################################################
        # Structure Prediction
        ####################################################
        candidates = predict_population_structures(candidates)

        ####################################################
        # Docking
        ####################################################
        docked = [
            scorer.evaluate_candidate(candidate)
            for candidate in candidates
        ]

        ####################################################
        # Ranking
        ####################################################
        ranked = ranker.rank(docked)

        if not ranked:
            logger.error("Ranking failed.")
            break

        ####################################################
        # Record Best Candidate
        ####################################################
        best = ranked[0]
        best_scores.append(best.overall_score)

        ####################################################
        # ARISE Learning
        ####################################################
        ARISE_ENGINE.observe_generation(ranked)
        ARISE_ENGINE.update_importance()

        logger.info(
            "ARISE importance: %s",
            ARISE_ENGINE.importance_map(),
        )
        logger.info("Best Candidate: %s", best.sequence)
        logger.info("Score: %.4f", best.overall_score)

        ####################################################
        # Seed Next Generation
        ####################################################
        seed_sequence = best.sequence
        final_ranked = ranked

    ########################################################
    # Final Results
    ########################################################
    if not final_ranked:
        logger.error("Pipeline produced no candidates.")
        return

    runtime = time.time() - start
    top_candidates = final_ranked[:TOP_K]

    ########################################################
    # Analysis
    ########################################################
    report = AnalysisEngine().analyze(
        top_candidates,
        runtime_seconds=runtime,
    )

    ########################################################
    # Figures
    ########################################################
    PlotGenerator(RESULTS_DIR).generate_all(
        top_candidates,
        best_scores,
    )

    ########################################################
    # Report
    ########################################################
    ReportGenerator(RESULTS_DIR).export(report)

    ########################################################
    # Finish
    ########################################################
    logger.info("=" * 50)
    logger.info("Pipeline Complete")
    logger.info("Generations: %d", len(best_scores))
    logger.info(
        "Best Overall Score: %.4f",
        final_ranked[0].overall_score,
    )
    logger.info("Runtime: %.2f sec", runtime)
    logger.info("=" * 50)


if __name__ == "__main__":
    run_pipeline()