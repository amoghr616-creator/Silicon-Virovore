"""
run_pipeline.py

Master execution pipeline for Silicon Virovore.
Prototype 1:
Computational peptide discovery + ARISE learning.
"""

from __future__ import annotations

import logging
import time

from src.population_runner import (
    generate_candidates,
    ARISE_ENGINE,
)
from src.predict_structure import (
    predict_population_structures,
)
from src.docking_vina import (
    PeptideDockingScorer,
)
from src.ranking import (
    CandidateRanker,
)
from src.analysis import (
    AnalysisEngine,
)
from src.report import (
    ReportGenerator,
)
from src.plots import (
    PlotGenerator,
)
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

    logger.info("===================================")
    logger.info("Silicon Virovore")
    logger.info("Prototype 1 Pipeline Started")
    logger.info("===================================")

    # ------------------------------------------------------
    # Evolutionary state
    # ------------------------------------------------------

    seed_sequence = DEFAULT_SEED_SEQUENCE

    best_scores = []
    final_ranked = []

    # ------------------------------------------------------
    # Multi-generation optimization
    # ------------------------------------------------------

    for generation in range(1, GENERATIONS + 1):

        logger.info("-----------------------------------")
        logger.info(
            "Generation %d / %d",
            generation,
            GENERATIONS,
        )
        logger.info(
            "Seed: %s",
            seed_sequence,
        )
        logger.info("-----------------------------------")

        # --------------------------------------------------
        # Population Generation
        # --------------------------------------------------

        logger.info("Generating candidates...")

        candidates = generate_candidates(
            seed_sequence,
        )

        if not candidates:

            logger.error(
                "No candidates generated for generation %d.",
                generation,
            )
            break

        # --------------------------------------------------
        # Structure Prediction
        # --------------------------------------------------

        logger.info(
            "Predicting structures..."
        )

        candidates = predict_population_structures(
            candidates
        )

        # --------------------------------------------------
        # Docking
        # --------------------------------------------------

        logger.info(
            "Running docking..."
        )

        scorer = PeptideDockingScorer()

        docked = []

        for candidate in candidates:

            docked.append(
                scorer.evaluate_candidate(candidate)
            )

        # --------------------------------------------------
        # Ranking
        # --------------------------------------------------

        logger.info(
            "Ranking candidates..."
        )

        ranked = CandidateRanker().rank(
            docked
        )

        if not ranked:

            logger.error(
                "No ranked candidates for generation %d.",
                generation,
            )
            break

        # --------------------------------------------------
        # ARISE learns from completed generation
        # --------------------------------------------------

        ARISE_ENGINE.observe_generation(
            ranked
        )

        ARISE_ENGINE.update_importance()

        best = ranked[0]

        best_scores.append(
            best.overall_score
        )

        logger.info(
            "Generation %d best score: %.4f",
            generation,
            best.overall_score,
        )

        logger.info(
            "ARISE importance map: %s",
            ARISE_ENGINE.importance_map(),
        )

        # --------------------------------------------------
        # Select best candidate as next-generation seed
        # --------------------------------------------------

        seed_sequence = best.sequence

        final_ranked = ranked

    # ------------------------------------------------------
    # Final Results
    # ------------------------------------------------------

    if not final_ranked:

        logger.error(
            "Pipeline produced no final candidates."
        )
        return

    logger.info("Preparing final results...")

    top_candidates = final_ranked[:TOP_K]

    # ------------------------------------------------------
    # Final Scientific Analysis
    # ------------------------------------------------------

    runtime = time.time() - start

    report = AnalysisEngine().analyze(
        top_candidates,
        runtime_seconds=runtime,
    )

    # ------------------------------------------------------
    # Figures
    # ------------------------------------------------------

    logger.info(
        "Generating figures..."
    )

    PlotGenerator(
        RESULTS_DIR,
    ).generate_all(
        top_candidates,
        best_scores,
    )

    # ------------------------------------------------------
    # Report
    # ------------------------------------------------------

    logger.info(
        "Generating report..."
    )

    ReportGenerator(
        RESULTS_DIR,
    ).export(report)

    # ------------------------------------------------------
    # Completion
    # ------------------------------------------------------

    logger.info("===================================")
    logger.info("Prototype 1 Pipeline Complete")
    logger.info(
        "Generations completed: %d",
        len(best_scores),
    )
    logger.info(
        "Final best score: %.4f",
        final_ranked[0].overall_score,
    )
    logger.info(
        "Runtime %.2f sec",
        runtime,
    )
    logger.info("===================================")


if __name__ == "__main__":
    run_pipeline()