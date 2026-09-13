"""Master execution pipeline for Silicon Virovore."""

from __future__ import annotations

import json
import logging
import random
import statistics
import time

from c.bridge import seed_native_random
from src.analysis import AnalysisEngine
from src.config import (
    ARISE_ENABLED,
    BOOTSTRAP_ITERATIONS,
    DEFAULT_SEED_SEQUENCE,
    ELITE_COUNT,
    GENERATIONS,
    MUTATION_RATE,
    PEPTIDE_FRAGMENT_SIZE,
    POPULATION_SIZE,
    RANDOM_SEED,
    RESULTS_DIR,
    STRUCTURE_DIR,
    TIER2_DOCKING_THRESHOLD,
    validate_configuration,
)
from src.diversity import DiversityAnalyzer
from src.docking_vina import PeptideDockingScorer
from src.models import (
    Candidate,
    candidate_has_valid_structure,
    mark_structure_unavailable,
)
from src.plots import PlotGenerator
from src.population_runner import (
    ARISE_ENGINE,
    generate_candidates,
)
from src.predict_structure import predict_population_structures
from src.ranking import CandidateRanker
from src.report import ReportGenerator

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def seed_experiment(seed: int | None) -> None:
    """Seed Python, optional NumPy, and native C stochastic components."""

    if seed is None:
        return

    random.seed(seed)
    try:
        import numpy as np
    except ImportError:
        pass
    else:
        np.random.seed(seed)

    seed_native_random(seed)


def _synchronize_structure_state(candidates: list[Candidate]) -> None:
    """Downgrade stale structure paths before ranking or evidence counts."""

    for candidate in candidates:
        if candidate.metadata.get("structure_available") is True:
            if not candidate_has_valid_structure(candidate):
                mark_structure_unavailable(
                    candidate,
                    status="unavailable",
                    error="Structure file or pLDDT evidence is missing.",
                )


def _mean(values: list[float]) -> float | None:
    return round(statistics.mean(values), 4) if values else None


def evidence_summary(candidates: list[Candidate]) -> dict[str, int]:
    """Count evidence states directly from candidate records."""

    return {
        "total_candidates": len(candidates),
        "structure_available": sum(
            candidate_has_valid_structure(candidate)
            for candidate in candidates
        ),
        "surrogate_available": sum(
            candidate.metadata.get("docking_status") == "complete"
            for candidate in candidates
        ),
        "surrogate_incomplete": sum(
            candidate.metadata.get("docking_status") == "incomplete"
            for candidate in candidates
        ),
        "tier2_eligible": sum(
            candidate.metadata.get("tier2_eligible", False)
            for candidate in candidates
        ),
        "tier2_attempted": sum(
            candidate.metadata.get("tier2_attempted", False)
            for candidate in candidates
        ),
        "tier2_validated": sum(
            candidate.metadata.get("tier2_validated", False)
            for candidate in candidates
        ),
        "vina_success": sum(
            candidate.metadata.get("vina_status") == "success"
            for candidate in candidates
        ),
        "vina_missing": sum(
            candidate.vina_delta_g is None
            for candidate in candidates
        ),
    }


def _candidate_provenance(candidate: Candidate) -> dict:
    """Serialize final-candidate provenance from the Candidate itself."""

    return {
        "sequence": candidate.sequence,
        "generation_first_observed": candidate.metadata.get("generation"),
        "generation_rank": candidate.metadata.get("generation_rank"),
        "generation_selected_as_best": candidate.metadata.get(
            "generation_selected_as_best", False
        ),
        "generation_seed": candidate.metadata.get("generation_seed"),
        "structure_status": candidate.metadata.get("structure_status"),
        "structure_available": candidate_has_valid_structure(candidate),
        "pLDDT": candidate.structure_confidence,
        "surrogate_best": candidate.strongest_anchor_delta_g,
        "surrogate_mean": candidate.mean_delta_g,
        "surrogate_consensus": candidate.metadata.get("consensus_docking"),
        "tier2_eligible": candidate.metadata.get("tier2_eligible", False),
        "tier2_attempted": candidate.metadata.get("tier2_attempted", False),
        "tier2_validated": candidate.metadata.get("tier2_validated", False),
        "tier2_status": candidate.metadata.get("tier2_status"),
        "vina_score": candidate.vina_delta_g,
        "vina_status": candidate.metadata.get("vina_status"),
        "overall_score": candidate.overall_score,
        "confidence": candidate.confidence,
        "generation_overall_score": candidate.metadata.get(
            "generation_overall_score"
        ),
        "generation_confidence": candidate.metadata.get(
            "generation_confidence"
        ),
        "generation_ranking_components": candidate.metadata.get(
            "generation_ranking_components", {}
        ),
        "ranking_components": candidate.metadata.get(
            "ranking_components", {}
        ),
    }


def _generation_summary(
    generation: int,
    candidates: list[Candidate],
    arise_stats: dict,
    arise_importance: dict[int, float] | None = None,
    arise_observations: dict[int, int] | None = None,
    generation_seed: str | None = None,
) -> dict:
    docking = [
        c.strongest_anchor_delta_g
        for c in candidates
        if c.strongest_anchor_delta_g is not None
    ]
    consensus = [
        c.metadata["consensus_docking"]
        for c in candidates
        if c.metadata.get("consensus_docking") is not None
    ]
    plddt = [
        c.structure_confidence
        for c in candidates
        if candidate_has_valid_structure(c)
    ]
    sequences = [c.sequence for c in candidates]
    diversity = DiversityAnalyzer().summary(
        sequences=sequences,
        reference=sequences[0],
    ) if sequences else {}
    sequence_length = (
        len(generation_seed)
        if generation_seed is not None
        else (len(sequences[0]) if sequences else 0)
    )
    mutation_frequency = {
        str(position): round(
            sum(
                candidate.sequence[position] != generation_seed[position]
                for candidate in candidates
                if generation_seed is not None
                and position < len(candidate.sequence)
                and position < len(generation_seed)
            ) / len(candidates),
            4,
        ) if candidates and generation_seed is not None else None
        for position in range(sequence_length)
    }

    surrogate_available = sum(
        c.metadata.get("docking_status") == "complete"
        for c in candidates
    )
    surrogate_incomplete = sum(
        c.metadata.get("docking_status") == "incomplete"
        for c in candidates
    )
    tier2_eligible = sum(
        c.metadata.get("tier2_eligible", False)
        for c in candidates
    )
    tier2_attempted = sum(
        c.metadata.get("tier2_attempted", False)
        for c in candidates
    )
    tier2_validated = sum(
        c.metadata.get("tier2_validated", False)
        for c in candidates
    )
    tier2_statuses = {
        c.metadata.get("tier2_status", "unknown")
        for c in candidates
    }
    if tier2_validated:
        tier2_evidence_status = "validated"
    elif "unavailable" in tier2_statuses or "eligible_no_structure" in tier2_statuses:
        tier2_evidence_status = "unavailable"
    else:
        tier2_evidence_status = "not_eligible"

    return {
        "generation": generation,
        "best_sequence": candidates[0].sequence if candidates else None,
        "best_overall": max(
            (c.overall_score for c in candidates),
            default=None,
        ),
        "mean_overall": _mean([c.overall_score for c in candidates]),
        "best_docking": min(docking) if docking else None,
        "mean_docking": _mean(docking),
        "best_consensus_docking": min(consensus) if consensus else None,
        "mean_consensus_docking": _mean(consensus),
        "mean_pLDDT": _mean(plddt),
        "best_pLDDT": max(plddt) if plddt else None,
        "structure_available": sum(
            candidate_has_valid_structure(c)
            for c in candidates
        ),
        "surrogate_available": surrogate_available,
        "surrogate_incomplete": surrogate_incomplete,
        "tier2_eligible": tier2_eligible,
        "tier2_attempted": tier2_attempted,
        "tier2_validated": tier2_validated,
        "vina_success": sum(
            c.metadata.get("vina_status") == "success"
            for c in candidates
        ),
        "tier2_evidence_status": tier2_evidence_status,
        "unique_sequences": len(set(sequences)),
        "duplicate_sequences": len(sequences) - len(set(sequences)),
        "mean_mutation_count": _mean([
            float(c.metadata.get("mutation_count", 0))
            for c in candidates
        ]),
        "mean_pairwise_identity": diversity.get(
            "mean_pairwise_identity"
        ),
        "mutation_frequency_by_position": mutation_frequency,
        "total_candidates": len(candidates),
        "diversity": diversity,
        "ARISE_mean_importance": arise_stats.get("mean", 0.0),
        "ARISE_nonuniformity": arise_stats.get("nonuniformity", 0.0),
        "ARISE_observed_positions": arise_stats.get("observed_positions", 0),
        "ARISE_importance": arise_importance or {},
        "ARISE_observations": arise_observations or {},
        "ARISE_objective": "overall_score",
    }


def _write_audit(audit: dict) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    (RESULTS_DIR / "experiment_audit.json").write_text(
        json.dumps(audit, indent=2, default=str) + "\n"
    )


def run_pipeline():
    validate_configuration()
    seed_experiment(RANDOM_SEED)
    ARISE_ENGINE.reset()

    start = time.time()
    logger.info("=" * 50)
    logger.info("Silicon Virovore")
    logger.info("Pipeline Started")
    logger.info("=" * 50)

    seed_sequence = DEFAULT_SEED_SEQUENCE
    best_scores: list[float] = []
    final_ranked: list[Candidate] = []
    all_candidates: list[Candidate] = []
    global_best: Candidate | None = None
    global_best_score: float | None = None
    global_best_sequence: str | None = None
    generation_summaries: list[dict] = []

    run_status = {
        "structure_success_count": 0,
        "structure_available_count": 0,
        "structure_cached_count": 0,
        "structure_unavailable_count": 0,
        "surrogate_docking_count": 0,
        "surrogate_incomplete_count": 0,
        "tier2_eligible_count": 0,
        "tier2_attempted_count": 0,
        "tier2_validated_count": 0,
        "total_candidates": 0,
    }
    scorer = PeptideDockingScorer()
    ranker = CandidateRanker()

    for generation in range(1, GENERATIONS + 1):
        logger.info("-" * 50)
        logger.info("Generation %d / %d", generation, GENERATIONS)
        logger.info("Seed: %s", seed_sequence)
        logger.info("-" * 50)

        candidates = generate_candidates(
            seed_sequence,
            population_size=POPULATION_SIZE,
            mutation_rate=MUTATION_RATE,
            arise_enabled=ARISE_ENABLED,
        )
        if not candidates:
            logger.error("No candidates generated.")
            break

        candidates = predict_population_structures(candidates)
        for candidate in candidates:
            candidate.metadata["generation"] = generation
            candidate.metadata["generation_seed"] = seed_sequence
        _synchronize_structure_state(candidates)

        docked = [scorer.evaluate_candidate(candidate) for candidate in candidates]
        all_candidates.extend(docked)

        structure_count = sum(
            candidate_has_valid_structure(candidate)
            for candidate in docked
        )
        cached_count = sum(
            candidate.metadata.get("structure_status") == "cached"
            for candidate in docked
        )
        surrogate_count = sum(
            c.metadata.get("docking_status") == "complete"
            for c in docked
        )
        surrogate_incomplete_count = sum(
            c.metadata.get("docking_status") == "incomplete"
            for c in docked
        )
        tier2_eligible_count = sum(
            c.metadata.get("tier2_eligible", False)
            for c in docked
        )
        tier2_attempted_count = sum(
            c.metadata.get("tier2_attempted", False)
            for c in docked
        )
        tier2_validated_count = sum(
            c.metadata.get("tier2_validated", False)
            for c in docked
        )
        run_status["structure_success_count"] += structure_count - cached_count
        run_status["structure_available_count"] += structure_count
        run_status["structure_cached_count"] += cached_count
        run_status["structure_unavailable_count"] += (
            len(docked) - structure_count
        )
        run_status["surrogate_docking_count"] += surrogate_count
        run_status["surrogate_incomplete_count"] += surrogate_incomplete_count
        run_status["tier2_eligible_count"] += tier2_eligible_count
        run_status["tier2_attempted_count"] += tier2_attempted_count
        run_status["tier2_validated_count"] += tier2_validated_count
        run_status["total_candidates"] += len(docked)

        ranked = ranker.rank(docked)
        if not ranked:
            logger.error("Ranking failed.")
            break

        current_best = ranked[0]
        for generation_rank, candidate in enumerate(ranked, start=1):
            candidate.metadata["generation_rank"] = generation_rank
            candidate.metadata["generation_selected_as_best"] = (
                generation_rank == 1
            )
            candidate.metadata["generation_overall_score"] = candidate.overall_score
            candidate.metadata["generation_confidence"] = candidate.confidence
            candidate.metadata["generation_ranking_components"] = dict(
                candidate.metadata.get("ranking_components", {})
            )
        best_scores.append(current_best.overall_score)
        if global_best is None or current_best.overall_score > global_best.overall_score:
            global_best = current_best
            global_best_score = current_best.overall_score
            global_best_sequence = current_best.sequence

        ARISE_ENGINE.observe_generation(ranked)
        ARISE_ENGINE.update_importance()
        arise_stats = ARISE_ENGINE.importance_statistics()
        summary = _generation_summary(
            generation,
            ranked,
            arise_stats,
            arise_importance=ARISE_ENGINE.importance_map(),
            arise_observations=ARISE_ENGINE.importance_observations(),
            generation_seed=seed_sequence,
        )
        generation_summaries.append(summary)
        logger.info("Generation %d summary: %s", generation, summary)
        logger.info("ARISE importance: %s", ARISE_ENGINE.importance_map())
        logger.info(
            "ARISE observations: %s",
            ARISE_ENGINE.importance_observations(),
        )
        logger.info("Best Candidate: %s", current_best.sequence)
        logger.info("Score: %.4f", current_best.overall_score)
        logger.info(
            "Selected objective contributions: %s",
            current_best.metadata.get("ranking_components", {}),
        )

        seed_sequence = current_best.sequence
        final_ranked = ranked

    if not final_ranked or not all_candidates:
        logger.error("Pipeline produced no candidates.")
        return

    runtime = time.time() - start
    report_ranked = ranker.rank(all_candidates)
    report = AnalysisEngine().analyze(
        report_ranked,
        runtime_seconds=runtime,
        structure_directory=STRUCTURE_DIR,
    )
    report_evidence = evidence_summary(all_candidates)
    final_candidate_provenance = _candidate_provenance(report_ranked[0])

    audit = {
        "configuration": {
            "population_size": POPULATION_SIZE,
            "generations": GENERATIONS,
            "mutation_rate": MUTATION_RATE,
            "elite_count": ELITE_COUNT,
            "random_seed": RANDOM_SEED,
            "tier2_threshold": TIER2_DOCKING_THRESHOLD,
            "fragment_size": PEPTIDE_FRAGMENT_SIZE,
            "arise_enabled": ARISE_ENABLED,
            "bootstrap_iterations": BOOTSTRAP_ITERATIONS,
        },
        "execution": run_status | {
            "runtime_seconds": round(runtime, 3),
        },
        "optimization": {
            "best_current_generation_score": final_ranked[0].overall_score,
            "best_global_score": report_ranked[0].overall_score,
            "best_global_sequence": report_ranked[0].sequence,
            "best_global_observed_score": (
                global_best_score
            ),
            "best_global_observed_sequence": (
                global_best_sequence
            ),
            "generation_summaries": generation_summaries,
            "ARISE_importance": ARISE_ENGINE.importance_map(),
            "ARISE_observations": ARISE_ENGINE.importance_observations(),
            "ARISE_statistics": ARISE_ENGINE.importance_statistics(),
            "final_candidate_provenance": final_candidate_provenance,
        },
        "evidence": {
            "structure_status": "available" if run_status[
                "structure_available_count"
            ] else "unavailable",
            "structure_cached_count": run_status[
                "structure_cached_count"
            ],
            "surrogate_status": "available" if run_status[
                "surrogate_docking_count"
            ] else "unavailable",
            "tier2_status": (
                "validated"
                if run_status["tier2_validated_count"]
                else (
                    "unavailable_or_failed"
                    if run_status["tier2_eligible_count"]
                    else "not_eligible"
                )
            ),
            "missing_structure_count": run_status[
                "structure_unavailable_count"
            ],
            "candidate_record_counts": report_evidence,
        },
    }
    report.audit = audit
    _write_audit(audit)

    PlotGenerator(RESULTS_DIR).generate_all(all_candidates, best_scores)
    ReportGenerator(RESULTS_DIR).export(report)

    best_overall = report.top_candidates[0]
    best_candidate = best_overall.candidate
    logger.info(
        "Final best overall | sequence=%s | overall=%.4f | confidence=%.3f | "
        "pLDDT=%s | surrogate_best=%s | surrogate_mean=%s | "
        "surrogate_consensus=%s | Vina=%s | tier2=%s",
        best_candidate.sequence,
        best_overall.overall_score,
        best_overall.confidence,
        best_candidate.structure_confidence,
        best_candidate.strongest_anchor_delta_g,
        best_candidate.mean_delta_g,
        best_candidate.metadata.get("consensus_docking"),
        best_candidate.vina_delta_g,
        best_candidate.metadata.get("tier2_validated", False),
    )
    logger.info(
        "Final Vina state | score=%s | status=%s | attempted=%s | validated=%s",
        best_candidate.vina_delta_g,
        best_candidate.metadata.get("vina_status"),
        best_candidate.metadata.get("tier2_attempted", False),
        best_candidate.metadata.get("tier2_validated", False),
    )
    logger.info(
        "Final objective contributions: %s",
        best_candidate.metadata.get("ranking_components", {}),
    )
    logger.info(
        "Evidence summary | structures=%d/%d | surrogate_available=%d/%d | "
        "tier2_eligible=%d/%d | tier2_attempted=%d/%d | "
        "tier2_validated=%d/%d",
        run_status["structure_available_count"],
        run_status["total_candidates"],
        run_status["surrogate_docking_count"],
        run_status["total_candidates"],
        run_status["tier2_eligible_count"],
        run_status["total_candidates"],
        run_status["tier2_attempted_count"],
        run_status["total_candidates"],
        run_status["tier2_validated_count"],
        run_status["total_candidates"],
    )
    logger.info("=" * 50)
    logger.info("Pipeline Complete")
    logger.info("Generations: %d", len(best_scores))
    logger.info("Best Overall Score: %.4f", report_ranked[0].overall_score)
    logger.info("Runtime: %.2f sec", runtime)
    logger.info("=" * 50)


if __name__ == "__main__":
    run_pipeline()
