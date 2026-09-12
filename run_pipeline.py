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
)

logging.basicConfig(
    level=logging.INFO,
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
    all_candidates = []

    run_status = {
        "structure_prediction": "not_started",
        "docking_backend": None,
        "validated_docking_count": 0,
        "surrogate_docking_count": 0,
    }
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

        for candidate in candidates:
            candidate.metadata["generation"] = generation
            candidate.metadata["generation_seed"] = seed_sequence

        structure_count = sum(
            1
            for candidate in candidates
            if candidate.structure_path is not None
        )

        if structure_count == len(candidates):
            run_status["structure_prediction"] = "complete"
        elif structure_count > 0:
            run_status["structure_prediction"] = "partial"
        else:
            run_status["structure_prediction"] = "unavailable"

        logger.info(
            "Structure prediction status: %s | %d/%d candidates have structures",
            run_status["structure_prediction"],
            structure_count,
            len(candidates),
        )

        ####################################################
        # Docking
        ####################################################
        docked = [
            scorer.evaluate_candidate(candidate)
            for candidate in candidates
        ]

        all_candidates.extend(docked)

        generation_surrogate_count = sum(
            1
            for candidate in docked
            if candidate.metadata.get(
                "docking_is_surrogate",
                False,
            )
        )

        generation_validated_count = sum(
            1
            for candidate in docked
            if candidate.metadata.get(
                "docking_validated",
                False,
            )
        )

        run_status["surrogate_docking_count"] += generation_surrogate_count
        run_status["validated_docking_count"] += generation_validated_count

        logger.info(
            "Docking status | surrogate=%d | validated=%d",
            generation_surrogate_count,
            generation_validated_count,
        )

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
    if not final_ranked or not all_candidates:
        logger.error("Pipeline produced no candidates.")
        return

    runtime = time.time() - start

    # Keep per-generation ranking for ARISE, but rank the complete run for
    # the final report.
    report_ranked = ranker.rank(all_candidates)

    ########################################################
    # Analysis
    ########################################################
    report = AnalysisEngine().analyze(
        report_ranked,
        runtime_seconds=runtime,
    )

    ########################################################
    # Figures
    ########################################################
    PlotGenerator(RESULTS_DIR).generate_all(
        all_candidates,
        best_scores,
    )

    ########################################################
    # Report
    ########################################################
    ReportGenerator(RESULTS_DIR).export(report)

    best_overall = report.top_candidates[0]
    best_candidate = best_overall.candidate

    logger.info(
        (
            "Final best overall | sequence=%s | overall=%.4f | "
            "confidence=%.3f | pLDDT=%s | surrogate_best=%s | "
            "surrogate_mean=%s | surrogate_consensus=%s | "
            "Vina=%s | tier2=%s"
        ),
        best_candidate.sequence,
        best_overall.overall_score,
        best_overall.confidence,
        best_candidate.structure_confidence,
        best_candidate.strongest_anchor_delta_g,
        best_candidate.mean_delta_g,
        best_candidate.metadata.get("consensus_docking"),
        best_candidate.vina_delta_g,
        best_candidate.passed_tier_2,
    )

    ########################################################
    # Finish
    ########################################################
    logger.info(
        "Evidence summary | structure=%s | surrogate docking=%d | validated docking=%d",
        run_status["structure_prediction"],
        run_status["surrogate_docking_count"],
        run_status["validated_docking_count"],
    )
    logger.info("=" * 50)
    logger.info("Pipeline Complete")
    logger.info("Generations: %d", len(best_scores))
    logger.info(
        "Best Overall Score: %.4f",
        report_ranked[0].overall_score,
    )
    logger.info("Runtime: %.2f sec", runtime)
    logger.info("=" * 50)


if __name__ == "__main__":
    run_pipeline()
