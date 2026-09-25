"""Master execution pipeline for Silicon Virovore."""

from __future__ import annotations

import argparse
import copy
from html import parser
import json
import logging
import random
import secrets
import statistics
import time
from dataclasses import replace
from pathlib import Path
from src import config
from src.arise_memory import (
    load_arise_memory,
    save_arise_memory,
    get_memory_history,
    upsert_arise_run,
    aggregate_memory_positional_signal,
)
from c.bridge import (
    clear_evaluation_cache,
    evaluation_cache_stats,
    seed_native_random,
)
from src.analysis import AnalysisEngine
from src.config import (
    BOOTSTRAP_ITERATIONS,
    PEPTIDE_FRAGMENT_SIZE,
    PipelineSettings,
    PILOT_SEED_SEQUENCE,
    validate_configuration,
)
from src.diversity import DiversityAnalyzer
from src.docking_vina import PeptideDockingScorer
from src.hotspot import (
    discover_hotspot_model,
    derive_positional_signal,
)
from src.models import (
    Candidate,
    candidate_has_valid_structure,
    mark_structure_unavailable,
)
from src.plots import PlotGenerator
from src.population_runner import (
    generate_candidates,
)
from src.predict_structure import predict_population_structures
from src.ranking import CandidateRanker
from src.report import ReportGenerator

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def resolve_experiment_seed(seed: int | None) -> int:
    """Return an explicit seed, generating one when the user did not supply it."""

    if seed is not None:
        return int(seed)

    # Fresh OS-generated seed for a new stochastic run.
    # The generated value is logged so the run can be reproduced later.
    return secrets.randbelow(2**32)

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


def _historical_occurrences(
    sequence: str,
    candidate_history: list[dict] | None,
) -> list[dict]:
    """Return every observed local score for a sequence with provenance."""

    if not candidate_history:
        return []
    occurrences = []
    for record in candidate_history:
        if record.get("sequence") != sequence:
            continue
        occurrences.append({
            "candidate_id": record.get("candidate_id"),
            "generation": record.get("generation"),
            "generation_local_score": record.get(
                "generation_local_score",
                record.get("score"),
            ),
            "generation_local_score_source": record.get(
                "generation_local_score_source",
                "generation_rank.current_candidate_set",
            ),
        })
    return sorted(
        occurrences,
        key=lambda item: (
            item.get("generation") or 0,
            item.get("candidate_id") or "",
        ),
    )


def _candidate_provenance(
    candidate: Candidate,
    candidate_history: list[dict] | None = None,
) -> dict:
    """Serialize final-candidate provenance from the Candidate itself."""

    return {
        "candidate_id": candidate.metadata.get("candidate_id"),
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
        "generation_local_score": candidate.metadata.get(
            "generation_local_score",
            candidate.overall_score,
        ),
        "generation_local_score_generation": candidate.metadata.get(
            "generation_local_score_generation",
            candidate.metadata.get("generation"),
        ),
        "generation_local_score_source": candidate.metadata.get(
            "generation_local_score_source",
            "generation_rank.current_candidate_set",
        ),
        "posthoc_recomputed_score": candidate.posthoc_recomputed_score,
        "posthoc_run_normalized_score": candidate.posthoc_recomputed_score,
        "posthoc_run_normalized_score_source": "full_run_candidate_set",
        "score_scope": candidate.metadata.get("score_scope", "generation_local"),
        "historical_occurrences": _historical_occurrences(
            candidate.sequence,
            candidate_history or [],
        ),
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
        "parent_sequence": candidate.metadata.get("parent_sequence"),
        "hotspot_guided": candidate.metadata.get("hotspot_guided", False),
        "mutation_guidance_mode": candidate.metadata.get(
            "mutation_guidance_mode"
        ),
        "mutation_positions": candidate.metadata.get("mutation_positions", []),
        "raw_objectives": candidate.metadata.get("raw_objectives", {}),
        "normalized_objectives": candidate.metadata.get(
            "normalized_objectives",
            {},
        ),
        "posthoc_normalized_objectives": candidate.metadata.get(
            "posthoc_normalized_objectives",
            {},
        ),
    }


def _candidate_history_record(candidate: Candidate) -> dict:
    """Serialize the complete per-candidate history without model objects."""

    return {
        "candidate_id": candidate.metadata.get("candidate_id"),
        "sequence": candidate.sequence,
        "generation": candidate.metadata.get("generation"),
        "run_id": candidate.metadata.get("run_id"),
        "score": candidate.overall_score,
        "generation_local_score": candidate.overall_score,
        "generation_local_score_generation": candidate.metadata.get(
            "generation"
        ),
        "generation_local_score_source": (
            "generation_rank.current_candidate_set"
        ),
        "score_scope": candidate.metadata.get(
            "score_scope",
            "generation_local",
        ),
        "arise_comparable_score": candidate.metadata.get(
            "arise_comparable_score"
        ),
        "rank": candidate.metadata.get("generation_rank"),
        "native_fitness": candidate.c_score,
        "consensus_docking": candidate.metadata.get("consensus_docking"),
        "vina_score": candidate.vina_delta_g,
        "structure_confidence": candidate.structure_confidence,
        "parent_sequence": candidate.metadata.get("parent_sequence"),
        "mutation_count": candidate.metadata.get("mutation_count", 0),
        "changed_positions_count": candidate.metadata.get(
            "changed_positions_count",
            0,
        ),
        "mutation_event_budget": candidate.metadata.get(
            "mutation_event_budget",
            0,
        ),
        "mutation_positions": candidate.metadata.get("mutation_positions", []),
        "hotspot_guided": candidate.metadata.get("hotspot_guided", False),
        "mutation_guidance_mode": candidate.metadata.get(
            "mutation_guidance_mode"
        ),
        "raw_objectives": candidate.metadata.get("raw_objectives", {}),
        "generation_normalized_objectives": candidate.metadata.get(
            "normalized_objectives",
            {},
        ),
        "ranking_missing_evidence": candidate.metadata.get(
            "ranking_missing_evidence",
            [],
        ),
        "available_weight": candidate.metadata.get("available_weight"),
        "hotspot_discovery_generation": candidate.metadata.get(
            "hotspot_discovery_generation"
        ),
        "hotspot_start": candidate.metadata.get("hotspot_start"),
        "hotspot_end": candidate.metadata.get("hotspot_end"),
    }


def _generation_summary(
    generation: int,
    candidates: list[Candidate],
    generation_seed: str | None = None,
    hotspot_state: dict | None = None,
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
        "mutation_count": sum(
            int(c.metadata.get("mutation_count", 0))
            for c in candidates
        ),
        "changed_positions_count": sum(
            int(c.metadata.get("changed_positions_count", 0))
            for c in candidates
        ),
        "mean_pairwise_identity": diversity.get(
            "mean_pairwise_identity"
        ),
        "mutation_frequency_by_position": mutation_frequency,
        "total_candidates": len(candidates),
        "diversity": diversity,
        "hotspot_state": hotspot_state,
        "ARISE_objective": "pooled_history_overall_score",
        "arise_score_source": (hotspot_state or {}).get(
            "arise_score_source",
            "not_used",
        ),
        "arise_score_comparable_across_generations": (hotspot_state or {}).get(
            "arise_score_comparable_across_generations",
            False,
        ),
    }


def _write_audit(audit: dict, output_directory: Path) -> None:
    output_directory.mkdir(parents=True, exist_ok=True)
    (output_directory / "experiment_audit.json").write_text(
        json.dumps(audit, indent=2, default=str) + "\n"
    )


def _write_json(payload: object, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2, default=str) + "\n")


def _posthoc_rank_without_overwriting_generation_scores(
    ranker: CandidateRanker,
    candidates: list[Candidate],
) -> tuple[list[Candidate], dict[str, Candidate]]:
    """Recompute full-run normalized scores on copies only."""

    copies = copy.deepcopy(candidates)
    ranked = ranker.rank(
        copies,
        score_field="posthoc_recomputed_score",
        metadata_prefix="posthoc",
        normalization_scope="full_run",
    )
    originals = {
        candidate.metadata.get("candidate_id"): candidate
        for candidate in candidates
    }
    for candidate in ranked:
        candidate.metadata["score_scope"] = "posthoc_full_run_normalized"
        candidate.metadata["generation_local_score"] = candidate.overall_score
        candidate.metadata["generation_local_score_generation"] = (
            candidate.metadata.get("generation")
        )
        candidate.metadata["generation_local_score_source"] = (
            "generation_rank.current_candidate_set"
        )
        candidate.metadata["posthoc_run_normalized_score"] = (
            candidate.posthoc_recomputed_score
        )
        candidate.metadata["posthoc_run_normalized_score_source"] = (
            "full_run_candidate_set"
        )
        candidate.overall_score = candidate.posthoc_recomputed_score or 0.0
        candidate.rank = candidate.metadata["posthoc_rank"]
        original = originals.get(candidate.metadata.get("candidate_id"))
        if original is not None:
            original.posthoc_recomputed_score = candidate.posthoc_recomputed_score
            for key, value in candidate.metadata.items():
                if key.startswith("posthoc_"):
                    original.metadata[key] = value
            original.metadata["posthoc_recomputed_score"] = (
                candidate.posthoc_recomputed_score
            )
            original.metadata["generation_local_score_generation"] = (
                candidate.metadata.get("generation")
            )
            original.metadata["generation_local_score_source"] = (
                "generation_rank.current_candidate_set"
            )
            original.metadata["posthoc_run_normalized_score"] = (
                candidate.posthoc_recomputed_score
            )
            original.metadata["posthoc_run_normalized_score_source"] = (
                "full_run_candidate_set"
            )
    return ranked, originals


def _finalize_comparable_telemetry(
    generation_summaries: list[dict],
    candidates: list[Candidate],
    candidate_history: list[dict],
    score_threshold: float | None,
) -> dict:
    """Add run-normalized trajectories without changing the search behavior."""

    history_by_id = {
        record.get("candidate_id"): record
        for record in candidate_history
    }
    for candidate in candidates:
        record = history_by_id.get(candidate.metadata.get("candidate_id"))
        if record is None:
            continue
        record.update({
            "generation_local_score": candidate.overall_score,
            "generation_local_score_generation": candidate.metadata.get(
                "generation"
            ),
            "generation_local_score_source": (
                "generation_rank.current_candidate_set"
            ),
            "posthoc_recomputed_score": candidate.posthoc_recomputed_score,
            "posthoc_run_normalized_score": candidate.posthoc_recomputed_score,
            "posthoc_run_normalized_score_source": (
                "full_run_candidate_set"
            ),
            "run_normalized_score": candidate.posthoc_recomputed_score,
            "raw_objectives": candidate.metadata.get("raw_objectives", {}),
            "generation_normalized_objectives": candidate.metadata.get(
                "normalized_objectives",
                {},
            ),
            "posthoc_normalized_objectives": candidate.metadata.get(
                "posthoc_normalized_objectives",
                {},
            ),
        })

    best_so_far_score = None
    best_so_far_sequence = None
    trajectory = []
    for summary in generation_summaries:
        generation_candidates = [
            candidate
            for candidate in candidates
            if candidate.metadata.get("generation") == summary["generation"]
        ]
        if not generation_candidates:
            continue
        local_best = max(
            generation_candidates,
            key=lambda candidate: candidate.overall_score,
        )
        comparable_best = max(
            generation_candidates,
            key=lambda candidate: candidate.posthoc_recomputed_score or 0.0,
        )
        comparable_score = comparable_best.posthoc_recomputed_score or 0.0
        if best_so_far_score is None or comparable_score > best_so_far_score:
            best_so_far_score = comparable_score
            best_so_far_sequence = comparable_best.sequence
        comparable_scores = [
            candidate.posthoc_recomputed_score or 0.0
            for candidate in generation_candidates
        ]
        summary.update({
            "generation_local_best_score": local_best.overall_score,
            "generation_local_best_sequence": local_best.sequence,
            "generation_best_score": comparable_score,
            "generation_best_sequence": comparable_best.sequence,
            "best_so_far_score": best_so_far_score,
            "best_so_far_sequence": best_so_far_sequence,
            "population_mean_score": round(statistics.mean(comparable_scores), 6),
            "population_median_score": round(statistics.median(comparable_scores), 6),
            "score_normalization_scope": "full_run",
            # Backward-compatible alias; new consumers should use
            # generation_best_score explicitly.
            "best_overall": comparable_score,
        })
        trajectory.append({
            "generation": summary["generation"],
            "generation_best_score": comparable_score,
            "generation_best_sequence": comparable_best.sequence,
            "best_so_far_score": best_so_far_score,
            "best_so_far_sequence": best_so_far_sequence,
            "mean_score": round(statistics.mean(comparable_scores), 6),
            "median_score": round(statistics.median(comparable_scores), 6),
            "diversity": summary.get("diversity"),
            "unique_sequences": summary.get("unique_sequences"),
            "candidate_count": summary.get("total_candidates"),
            "mutation_count": summary.get("mutation_count"),
            "runtime_seconds": summary.get("runtime_seconds"),
        })

    reached = [
        candidate
        for candidate in sorted(
            candidates,
            key=lambda item: item.metadata.get("evaluation_index", 0),
        )
        if score_threshold is not None
        and (candidate.posthoc_recomputed_score or 0.0) >= score_threshold
    ]
    first = reached[0] if reached else None
    return {
        "trajectory": trajectory,
        "candidates_to_threshold": (
            first.metadata.get("evaluation_index") if first else None
        ),
        "generation_to_threshold": (
            first.metadata.get("generation") if first else None
        ),
        "threshold": score_threshold,
    }


def _run_directory(settings: PipelineSettings) -> Path:
    output_directory = Path(settings.output_directory)
    if settings.run_id:
        output_directory = output_directory / settings.run_id
    return output_directory


def run_pipeline(settings: PipelineSettings | None = None):
    settings = settings or PipelineSettings()

    # Explicit seed = deterministic/reproducible run.
    # No seed = fresh stochastic run.
    effective_seed = resolve_experiment_seed(settings.random_seed)

    effective_run_id = (
        settings.run_id
        or f"{settings.condition}-seed{effective_seed}"
    )

    # Do not mutate the caller's PipelineSettings object.
    settings = replace(
        settings,
        random_seed=effective_seed,
        run_id=effective_run_id,
    )

    validate_configuration(settings)
    output_directory = _run_directory(settings)
    structure_directory = output_directory / "structures"
    docking_directory = output_directory / "docking"
    seed_experiment(settings.random_seed)
    clear_evaluation_cache()

    start = time.time()
    start_perf = time.perf_counter()
    logger.info("=" * 50)
    logger.info(
        "Silicon Virovore | run_id=%s | condition=%s | seed=%s | population=%d | generations=%d",
        settings.run_id or "default",
        settings.condition,
        settings.random_seed,
        settings.population_size,
        settings.generations,
    )
    logger.info("Pipeline Started")
    logger.info("=" * 50)

    seed_sequence = settings.seed_sequence
    best_scores: list[float] = []
    final_ranked: list[Candidate] = []
    all_candidates: list[Candidate] = []
    global_best: Candidate | None = None
    global_best_score: float | None = None
    global_best_sequence: str | None = None
    generation_summaries: list[dict] = []
    candidate_history: list[dict] = []
    hotspot_models: list[dict] = []
    hotspot_model: dict | None = None
    previous_best_sequence: str | None = None

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
        "generation_runtime_seconds": [],
        "structure_runtime_seconds": [],
        "docking_runtime_seconds": [],
        "ranking_runtime_seconds": [],
    }
    scorer = PeptideDockingScorer(
        docking_directory=docking_directory,
        workers=settings.docking_workers,
        cache_enabled=settings.cache_enabled,
        tier2_threshold=settings.docking_threshold,
    )
    ranker = CandidateRanker(
    fitness_weight=config.FITNESS_WEIGHT,
    docking_weight=config.DOCKING_WEIGHT,
    helix_weight=config.HELIX_WEIGHT,
    solvation_weight=config.SOLVATION_WEIGHT,
)

    arise_memory = load_arise_memory(
        settings.arise_memory_path,
        max_records=settings.arise_memory_max_records,
        max_runs=settings.arise_memory_max_runs,
    )
    prior_memory_history = get_memory_history(arise_memory)
    use_prior_memory = (
        settings.arise_memory_enabled
        and settings.arise_memory_use_history
    )

    logger.info(
        "ARISE memory | enabled=%s | use_history=%s | prior_runs=%d | prior_records=%d",
        settings.arise_memory_enabled,
        use_prior_memory,
        len(arise_memory.get("runs", [])),
        len(prior_memory_history),
    )

    for generation in range(1, settings.generations + 1):
        generation_start = time.perf_counter()
        logger.info("-" * 50)
        logger.info("Generation %d / %d", generation, settings.generations)
        logger.info("Seed: %s", seed_sequence)
        logger.info("-" * 50)

        candidates = generate_candidates(
            seed_sequence,
            population_size=settings.population_size,
            mutation_rate=settings.mutation_rate,
            arise_enabled=(
                settings.adaptive_enabled
                and settings.arise_enabled
                and settings.ale_enabled
                and settings.mutation_guidance_mode == "legacy"
            ),
            cache_enabled=settings.cache_enabled,
            hotspot_model=hotspot_model,
            hotspot_enabled=settings.hotspot_enabled,
            ale_enabled=(
                settings.adaptive_enabled
                and settings.ale_enabled
            ),
            mutation_guidance_mode=settings.mutation_guidance_mode,
            generation=generation,
            run_id=settings.run_id,
            preserve_sequence=(
                previous_best_sequence
                if settings.preserve_best
                else None
            ),
        )
        if not candidates:
            logger.error("No candidates generated.")
            break

        structure_start = time.perf_counter()
        candidates = predict_population_structures(
            candidates,
            structure_directory=structure_directory,
            workers=settings.structure_workers,
        )
        structure_runtime = time.perf_counter() - structure_start
        for candidate in candidates:
            candidate.metadata["generation"] = generation
            candidate.metadata["generation_seed"] = seed_sequence
        _synchronize_structure_state(candidates)

        docking_start = time.perf_counter()
        docked = scorer.evaluate_candidates(candidates)
        docking_runtime = time.perf_counter() - docking_start
        evaluation_start = len(all_candidates) + 1
        for offset, candidate in enumerate(docked):
            candidate.metadata["evaluation_index"] = evaluation_start + offset
            candidate.metadata["candidate_id"] = (
                f"{settings.run_id or 'run'}:g{generation}:c{offset + 1}"
            )
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

        ranking_start = time.perf_counter()
        ranked = ranker.rank(docked)
        ranking_runtime = time.perf_counter() - ranking_start
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

        candidate_history.extend(
            _candidate_history_record(candidate)
            for candidate in ranked
        )

        if (
            settings.hotspot_enabled
            and settings.arise_enabled
            and generation % settings.hotspot_update_interval == 0
        ):
            discovery_history = (
                prior_memory_history + candidate_history
                if use_prior_memory
                else candidate_history
            )

            hotspot_model = discover_hotspot_model(
                discovery_history,
                discovery_generation=generation,
                run_id=settings.run_id,
                top_quantile=settings.hotspot_top_quantile,
                window_min=settings.hotspot_discovery_window_min,
                window_max=settings.hotspot_discovery_window_max,
                bootstrap_iterations=settings.hotspot_bootstrap_iterations,
                random_seed=(settings.random_seed or 0) + generation,
                hotspot_start=settings.hotspot_start,
                hotspot_end=settings.hotspot_end,
            )
            hotspot_model["derived_positional_signal"] = derive_positional_signal(
                hotspot_model,
                len(seed_sequence),
            )
            hotspot_model["derived_signal_source"] = (
                "arise-hotspot-v2.selection_score"
            )
            hotspot_model["derived_signal_used_for_mutation"] = False
            hotspot_model["memory_used"] = use_prior_memory
            hotspot_model["memory_prior_run_count"] = len(
                arise_memory.get("runs", [])
            )
            hotspot_model["memory_prior_record_count"] = len(
                prior_memory_history
            )
            logger.info(
                "ARISE positional signal | values=%s | used_for_mutation=%s",
                hotspot_model["derived_positional_signal"],
                hotspot_model["derived_signal_used_for_mutation"],
            )
            hotspot_models.append(hotspot_model)
            logger.info(
                "ARISE score source | source=%s | comparable_across_generations=%s",
                hotspot_model.get("arise_score_source"),
                hotspot_model.get(
                    "arise_score_comparable_across_generations",
                    False,
                ),
            )
            selected = hotspot_model.get("selected_hotspot")
            if selected:
                logger.info(
                    "ARISE hotspot | generation=%d | positions=%d:%d | "
                    "length=%d | train=%.4f | heldout=%s | stability=%.3f",
                    generation,
                    selected["start_position"],
                    selected["end_position"],
                    selected["window_length"],
                    selected["hotspot_score"],
                    selected.get("held_out_association"),
                    selected.get("bootstrap_stability", 0.0),
                )
            else:
                logger.info(
                    "ARISE hotspot unavailable after generation %d: %s",
                    generation,
                    hotspot_model.get("status"),
                )

        summary = _generation_summary(
            generation,
            ranked,
            generation_seed=seed_sequence,
            hotspot_state=hotspot_model,
        )
        generation_runtime = time.perf_counter() - generation_start
        summary.update({
            "median_overall": statistics.median(
                candidate.overall_score for candidate in ranked
            ),
            "runtime_seconds": round(generation_runtime, 4),
            "cumulative_runtime_seconds": round(
                time.perf_counter() - start_perf,
                4,
            ),
            "stage_runtime_seconds": {
                "structure": round(structure_runtime, 4),
                "docking": round(docking_runtime, 4),
                "ranking": round(ranking_runtime, 4),
            },
            "cache_stats": scorer.cache_stats() | {
                "c_fitness": evaluation_cache_stats(),
            },
        })
        generation_summaries.append(summary)
        logger.info("Generation %d summary: %s", generation, summary)
        logger.info("Generation-local best candidate: %s", current_best.sequence)
        logger.info("Generation-local score: %.4f", current_best.overall_score)
        logger.info(
            "Selected objective contributions: %s",
            current_best.metadata.get("ranking_components", {}),
        )
        logger.info(
            "Experiment telemetry | generation=%d | best=%.4f | median=%.4f | "
            "unique=%d/%d | evaluated=%d | runtime=%.3fs | cumulative=%.3fs",
            generation,
            current_best.overall_score,
            summary["median_overall"],
            summary["unique_sequences"],
            summary["total_candidates"],
            summary["total_candidates"],
            summary["runtime_seconds"],
            summary["cumulative_runtime_seconds"],
        )

        run_status["generation_runtime_seconds"].append(
            round(generation_runtime, 4)
        )
        run_status["structure_runtime_seconds"].append(
            round(structure_runtime, 4)
        )
        run_status["docking_runtime_seconds"].append(
            round(docking_runtime, 4)
        )
        run_status["ranking_runtime_seconds"].append(
            round(ranking_runtime, 4)
        )

        seed_sequence = current_best.sequence
        previous_best_sequence = current_best.sequence
        final_ranked = ranked

    if not final_ranked or not all_candidates:
        logger.error("Pipeline produced no candidates.")
        return

    runtime = time.time() - start
    report_ranked, _original_candidates = (
        _posthoc_rank_without_overwriting_generation_scores(
            ranker,
            all_candidates,
        )
    )
    comparable_metrics = _finalize_comparable_telemetry(
        generation_summaries,
        all_candidates,
        candidate_history,
        settings.score_threshold,
    )
    if settings.arise_memory_enabled:
        final_signal = (
        derive_positional_signal(
            hotspot_model,
            len(settings.seed_sequence),
        )
        if hotspot_model is not None
        else [0.0] * len(settings.seed_sequence)
    )

    upsert_arise_run(
        arise_memory,
        run_id=settings.run_id or f"run-{settings.random_seed}",
        condition=settings.condition,
        random_seed=settings.random_seed,
        candidate_history=candidate_history,
        hotspot_models=hotspot_models,
        derived_positional_signal=final_signal,
        max_records=settings.arise_memory_max_records,
        max_runs=settings.arise_memory_max_runs,
    )

    save_arise_memory(
        settings.arise_memory_path,
        arise_memory,
    )

    aggregate_signal = (
        aggregate_memory_positional_signal(
            arise_memory,
            len(settings.seed_sequence),
        )
    )

    logger.info(
        "ARISE memory saved | runs=%d | records=%d | "
        "aggregate_signal=%s",
        aggregate_signal["run_count"],
        len(arise_memory.get("candidate_history", [])),
        aggregate_signal["mean_signal"],
    )
    best_scores = [
        summary["generation_best_score"]
        for summary in generation_summaries
    ]
    final_generation_number = max(
        candidate.metadata.get("generation", 0)
        for candidate in all_candidates
    )
    final_population_candidates = [
        candidate
        for candidate in all_candidates
        if candidate.metadata.get("generation") == final_generation_number
    ]
    final_population_best = max(
        final_population_candidates,
        key=lambda candidate: candidate.posthoc_recomputed_score or 0.0,
    )
    final_generation_local_best = max(
        final_population_candidates,
        key=lambda candidate: candidate.overall_score,
    )
    final_generation_summary = generation_summaries[-1]
    report = AnalysisEngine().analyze(
        report_ranked,
        runtime_seconds=runtime,
        structure_directory=structure_directory,
    )
    report_evidence = evidence_summary(all_candidates)
    final_candidate_provenance = _candidate_provenance(
        report_ranked[0],
        candidate_history,
    )

    audit = {
        "configuration": {
            "population_size": settings.population_size,
            "generations": settings.generations,
            "mutation_rate": settings.mutation_rate,
            "elite_count": settings.elite_count,
            "tournament_size": settings.tournament_size,
            "random_seed": settings.random_seed,
            "seed_sequence": settings.seed_sequence,
            "adaptive_enabled": settings.adaptive_enabled,
            "arise_enabled": settings.arise_enabled,
            "ale_enabled": settings.ale_enabled,
            "condition": settings.condition,
            "run_id": settings.run_id,
            "output_directory": str(output_directory),
            "structure_workers": settings.structure_workers,
            "docking_workers": settings.docking_workers,
            "cache_enabled": settings.cache_enabled,
            "generate_plots": settings.generate_plots,
            "hotspot_enabled": settings.hotspot_enabled,
            "hotspot_start": settings.hotspot_start,
            "hotspot_end": settings.hotspot_end,
            "hotspot_discovery_window_min": settings.hotspot_discovery_window_min,
            "hotspot_discovery_window_max": settings.hotspot_discovery_window_max,
            "hotspot_update_interval": settings.hotspot_update_interval,
            "hotspot_top_quantile": settings.hotspot_top_quantile,
            "hotspot_bootstrap_iterations": settings.hotspot_bootstrap_iterations,
            "mutation_guidance_mode": settings.mutation_guidance_mode,
            "preserve_best": settings.preserve_best,
            "score_threshold": settings.score_threshold,
            "legacy_residue_importance_status": (
                "active"
                if settings.mutation_guidance_mode == "legacy"
                and settings.arise_enabled
                else "inactive"
            ),
            "arise_hotspot_discovery_status": (
                "active"
                if settings.hotspot_enabled and settings.arise_enabled
                else "inactive"
            ),
            "arise_score_source": (
                hotspot_model.get("arise_score_source", "not_used")
                if hotspot_model is not None
                else "not_used"
            ),
            "arise_score_comparable_across_generations": (
                hotspot_model.get(
                    "arise_score_comparable_across_generations",
                    False,
                )
                if hotspot_model is not None
                else False
            ),
            "ale_hotspot_guidance_status": (
                "active"
                if settings.hotspot_enabled and settings.ale_enabled
                else "inactive"
            ),
            "tier2_threshold": settings.docking_threshold,
            "fragment_size": PEPTIDE_FRAGMENT_SIZE,
            "bootstrap_iterations": BOOTSTRAP_ITERATIONS,
            "selection_implementation": (
                "python_best_candidate_carry_forward; "
                "elite/tournament apply to c/ga_loop only"
            ),
        },
        "execution": run_status | {
            "runtime_seconds": round(runtime, 3),
            "candidate_history_file": "candidate_history.json",
            "cache_stats": scorer.cache_stats() | {
                "c_fitness": evaluation_cache_stats(),
            },
        },
        "optimization": {
            "generation_best_score": final_generation_summary.get(
                "generation_best_score"
            ),
            "generation_best_sequence": final_generation_summary.get(
                "generation_best_sequence"
            ),
            "generation_best_generation": final_generation_number,
            "generation_best_score_source": "full_run_candidate_set",
            "generation_local_best_score": final_generation_summary.get(
                "generation_local_best_score"
            ),
            "generation_local_best_sequence": final_generation_summary.get(
                "generation_local_best_sequence"
            ),
            "generation_local_best_generation": final_generation_number,
            "generation_local_best_score_source": (
                "generation_rank.current_candidate_set"
            ),
            "best_so_far_score": comparable_metrics["trajectory"][-1].get(
                "best_so_far_score"
            ),
            "best_so_far_sequence": comparable_metrics["trajectory"][-1].get(
                "best_so_far_sequence"
            ),
            "best_so_far_score_source": "full_run_candidate_set",
            "final_population_best_score": final_population_best.posthoc_recomputed_score,
            "final_population_best_sequence": final_population_best.sequence,
            "final_population_generation": final_generation_number,
            "final_population_best_score_source": "full_run_candidate_set",
            "posthoc_recomputed_best_score": report_ranked[0].overall_score,
            "posthoc_recomputed_best_sequence": report_ranked[0].sequence,
            "posthoc_run_normalized_best_score": report_ranked[0].overall_score,
            "posthoc_run_normalized_best_sequence": report_ranked[0].sequence,
            "posthoc_run_normalized_best_generation": report_ranked[0].metadata.get(
                "generation"
            ),
            "posthoc_run_normalized_best_score_source": (
                "full_run_candidate_set"
            ),
            "score_normalization_scope": "full_run_for_comparable_telemetry",
            "trajectory": comparable_metrics["trajectory"],
            "threshold_metrics": {
                "threshold": comparable_metrics["threshold"],
                "candidates_to_threshold": comparable_metrics[
                    "candidates_to_threshold"
                ],
                "generation_to_threshold": comparable_metrics[
                    "generation_to_threshold"
                ],
            },
            "generation_summaries": generation_summaries,
            "hotspot_models": hotspot_models,
            "final_hotspot_model": hotspot_model,
            "candidate_history_count": len(candidate_history),
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
    _write_audit(audit, output_directory)
    _write_json(candidate_history, output_directory / "candidate_history.json")

    if settings.generate_plots:
        try:
            PlotGenerator(output_directory).generate_all(
                all_candidates,
                best_scores,
            )
        except Exception as exc:
            logger.warning("Plot generation failed; reports remain available: %s", exc)
    ReportGenerator(output_directory).export(report)

    best_overall = report.top_candidates[0]
    best_candidate = best_overall.candidate
    logger.info(
        "Posthoc run-normalized best | sequence=%s | score=%.4f | confidence=%.3f | "
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
        "Final-generation score comparison | sequence=%s | "
        "generation_local_score=%s | generation_local_score_generation=%s | "
        "generation_local_score_source=%s | "
        "posthoc_run_normalized_score=%s | "
        "posthoc_run_normalized_score_source=%s",
        final_generation_local_best.sequence,
        final_generation_local_best.overall_score,
        final_generation_local_best.metadata.get("generation"),
        "generation_rank.current_candidate_set",
        final_generation_local_best.posthoc_recomputed_score,
        "full_run_candidate_set",
    )
    logger.info(
        "Posthoc-global candidate score provenance | sequence=%s | "
        "generation_local_score=%s | generation_local_score_generation=%s | "
        "generation_local_score_source=%s | "
        "posthoc_run_normalized_score=%s | "
        "posthoc_run_normalized_score_source=%s | historical_occurrences=%d",
        best_candidate.sequence,
        best_candidate.metadata.get("generation_local_score"),
        best_candidate.metadata.get("generation_local_score_generation"),
        best_candidate.metadata.get(
            "generation_local_score_source",
            "generation_rank.current_candidate_set",
        ),
        best_candidate.posthoc_recomputed_score,
        "full_run_candidate_set",
        len(_historical_occurrences(best_candidate.sequence, candidate_history)),
    )
    logger.info(
        "Best definitions | final_population=%s (%s) | posthoc_global=%s (%s) | "
        "threshold_candidates=%s | threshold_generation=%s",
        final_population_best.posthoc_recomputed_score,
        final_population_best.sequence,
        report_ranked[0].overall_score,
        report_ranked[0].sequence,
        comparable_metrics["candidates_to_threshold"],
        comparable_metrics["generation_to_threshold"],
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


def parse_settings(argv: list[str] | None = None) -> PipelineSettings:
    """Build an experiment configuration without changing source defaults."""

    parser = argparse.ArgumentParser(
        description="Run a reproducible Silicon Virovore experiment."
    )
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--population-size", type=int)
    parser.add_argument("--generations", type=int)
    parser.add_argument("--mutation-rate", type=float)
    parser.add_argument("--elite-count", type=int)
    parser.add_argument("--tournament-size", type=int)
    parser.add_argument(
        "--seed",
        type=int,
        help="Explicit RNG seed. Omit for a fresh randomly generated seed.",
    )
    parser.add_argument("--seed-sequence")
    parser.add_argument(
        "--condition",
        choices=("baseline", "adaptive", "pilot"),
    )
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--run-id")
    parser.add_argument(
        "--adaptive-enabled",
        action=argparse.BooleanOptionalAction,
        default=None,
    )
    parser.add_argument(
        "--arise-enabled",
        action=argparse.BooleanOptionalAction,
        default=None,
    )
    parser.add_argument(
        "--ale-enabled",
        action=argparse.BooleanOptionalAction,
        default=None,
    )
    parser.add_argument(
        "--arise-memory-enabled",
        action=argparse.BooleanOptionalAction,
        default=None,
    )
    parser.add_argument(
        "--arise-memory-use-history",
        action=argparse.BooleanOptionalAction,
        default=None,
    )
    parser.add_argument(
        "--arise-memory-max-records",
        type=int,
    )
    parser.add_argument(
        "--arise-memory-max-runs",
        type=int,
    )
    parser.add_argument("--structure-workers", type=int)
    parser.add_argument("--docking-workers", type=int)
    parser.add_argument(
        "--hotspot-enabled",
        action=argparse.BooleanOptionalAction,
        default=None,
    )
    parser.add_argument("--hotspot-start", type=int)
    parser.add_argument("--hotspot-end", type=int)
    parser.add_argument("--hotspot-discovery-window-min", type=int)
    parser.add_argument("--hotspot-discovery-window-max", type=int)
    parser.add_argument("--hotspot-update-interval", type=int)
    parser.add_argument("--hotspot-top-quantile", type=float)
    parser.add_argument("--hotspot-bootstrap-iterations", type=int)
    parser.add_argument(
        "--preserve-best",
        action=argparse.BooleanOptionalAction,
        default=None,
    )
    parser.add_argument("--score-threshold", type=float)
    parser.add_argument(
        "--mutation-guidance-mode",
        choices=("uniform", "hotspot_only", "hotspot_biased"),
    )
    parser.add_argument(
        "--cache-enabled",
        action=argparse.BooleanOptionalAction,
        default=None,
    )
    parser.add_argument(
        "--generate-plots",
        action=argparse.BooleanOptionalAction,
        default=None,
    )
    args = parser.parse_args(argv)

    values: dict[str, object] = {}
    for field_name, argument_name in (
        ("population_size", "population_size"),
        ("generations", "generations"),
        ("mutation_rate", "mutation_rate"),
        ("elite_count", "elite_count"),
        ("tournament_size", "tournament_size"),
        ("random_seed", "seed"),
        ("seed_sequence", "seed_sequence"),
        ("output_directory", "output_dir"),
        ("run_id", "run_id"),
        ("structure_workers", "structure_workers"),
        ("docking_workers", "docking_workers"),
        ("hotspot_start", "hotspot_start"),
        ("hotspot_end", "hotspot_end"),
        ("hotspot_discovery_window_min", "hotspot_discovery_window_min"),
        ("hotspot_discovery_window_max", "hotspot_discovery_window_max"),
        ("hotspot_update_interval", "hotspot_update_interval"),
        ("hotspot_top_quantile", "hotspot_top_quantile"),
        ("hotspot_bootstrap_iterations", "hotspot_bootstrap_iterations"),
        ("mutation_guidance_mode", "mutation_guidance_mode"),
        ("score_threshold", "score_threshold"),
        ("arise_memory_max_records", "arise_memory_max_records"),
        ("arise_memory_max_runs", "arise_memory_max_runs"),
    ):
        value = getattr(args, argument_name)
        if value is not None:
            values[field_name] = value

    if args.condition is not None:
        values["condition"] = args.condition
    if args.adaptive_enabled is not None:
        values["adaptive_enabled"] = args.adaptive_enabled
    if args.arise_enabled is not None:
        values["arise_enabled"] = args.arise_enabled
    if args.ale_enabled is not None:
        values["ale_enabled"] = args.ale_enabled
    if args.arise_memory_enabled is not None:
        values["arise_memory_enabled"] = args.arise_memory_enabled
    if args.arise_memory_use_history is not None:
        values["arise_memory_use_history"] = args.arise_memory_use_history
    if args.cache_enabled is not None:
        values["cache_enabled"] = args.cache_enabled
    if args.generate_plots is not None:
        values["generate_plots"] = args.generate_plots
    if args.hotspot_enabled is not None:
        values["hotspot_enabled"] = args.hotspot_enabled
    if args.preserve_best is not None:
        values["preserve_best"] = args.preserve_best

    if args.pilot:
        values.update({
            "population_size": 16,
            "generations": 3,
            "random_seed": 1001,
            "condition": "pilot",
            "seed_sequence": values.get(
                "seed_sequence",
                PILOT_SEED_SEQUENCE,
            ),
        })
        values.setdefault("run_id", "pilot-1001")

    condition = values.get("condition")
    if condition == "baseline":
        values.setdefault("adaptive_enabled", False)
        values.setdefault("arise_enabled", False)
        values.setdefault("ale_enabled", False)
        values.setdefault("hotspot_enabled", False)
        values.setdefault("mutation_guidance_mode", "uniform")
    elif condition == "adaptive":
        values.setdefault("adaptive_enabled", True)
        values.setdefault("arise_enabled", True)
        values.setdefault("ale_enabled", True)
        values.setdefault("hotspot_enabled", True)
        values.setdefault("mutation_guidance_mode", "hotspot_biased")
    elif condition == "pilot":
        values.setdefault("hotspot_enabled", True)
        values.setdefault("mutation_guidance_mode", "hotspot_biased")

    return PipelineSettings(**values)


if __name__ == "__main__":
    run_pipeline(parse_settings())
