"""
population_runner.py

Generates candidate peptides using the native C backend or controlled
Python mutation policies.

ARISE's mutation-credit mode learns from eligible, explicitly recorded
parent-child score changes. A missing or invalid structural reward is never
replaced silently with sequence-level fitness.

Modes:
    legacy                  Existing native C mutation generator.
    uniform_matched         Fixed-budget uniform mutation baseline.
    mutation_credit         Fixed-budget mutation-credit adaptive search.
    uniform                 Existing ALE uniform policy.
    hotspot_only            Existing ALE hotspot-only policy.
    hotspot_biased          Existing ALE hotspot-biased policy.
    hotspot_biased_annealed Existing ALE annealed hotspot policy.
"""

from __future__ import annotations

import logging
import math
from collections import defaultdict
from typing import Any

from c.bridge import (
    generate_c_population,
    generate_policy_population,
    process_candidate_peptide,
)

from src.ale_policy import resolve_mutation_policy, mutation_positions
from src.config import (
    ARISE_ENABLED,
    DIVERSITY_REFILL_ATTEMPTS,
    GENERATIONS,
    MUTATION_RATE,
    POPULATION_SIZE,
    REMOVE_DUPLICATE_CANDIDATES,
)
from src.models import Candidate
from src.mutation_credit import (
    MutationCreditARISE,
    generate_position_weighted_population,
    mutation_event_budget,
)

logger = logging.getLogger(__name__)


# ==========================================================
# Legacy diagnostic: sequence-position association
# ==========================================================

def estimate_position_effects(
    sequences: list[str],
    scores: list[float],
    *,
    min_group_size: int = 3,
    shrinkage: float = 5.0,
) -> dict[int, float]:
    """
    Estimate the strength of the association between residue identity
    at each position and scores in the supplied candidate set.

    This is an observational association, not a causal mutation effect
    and not evidence of target binding. It is retained for compatibility
    and analysis; mutation_credit mode does not use it as its reward.
    """
    if len(sequences) != len(scores):
        raise ValueError("sequences and scores must have equal length")

    if not sequences:
        return {}

    if min_group_size < 2 or shrinkage < 0:
        raise ValueError("invalid estimator parameters")

    try:
        numeric_scores = [float(score) for score in scores]
    except (TypeError, ValueError) as exc:
        raise ValueError("scores must be numeric and finite") from exc

    if any(not math.isfinite(score) for score in numeric_scores):
        raise ValueError("scores must be finite")

    if any(not isinstance(seq, str) or not seq for seq in sequences):
        raise ValueError("sequences must be nonempty strings")

    normalized_sequences = [seq.upper() for seq in sequences]
    length = len(normalized_sequences[0])

    if any(len(seq) != length for seq in normalized_sequences):
        raise ValueError("all sequences must have equal length")

    ordered = sorted(numeric_scores)

    def quantile(q: float) -> float:
        index = q * (len(ordered) - 1)
        lo = int(index)
        hi = min(lo + 1, len(ordered) - 1)
        fraction = index - lo
        return (
            ordered[lo] * (1.0 - fraction)
            + ordered[hi] * fraction
        )

    median = quantile(0.5)
    scale = max(quantile(0.75) - quantile(0.25), 1e-8)

    standardized = [
        max(-3.0, min(3.0, (score - median) / scale))
        for score in numeric_scores
    ]

    overall = sum(standardized) / len(standardized)
    total_ss = sum((score - overall) ** 2 for score in standardized)

    if total_ss <= 1e-12:
        return {position: 0.0 for position in range(length)}

    effects: dict[int, float] = {}

    for position in range(length):
        groups: dict[str, list[float]] = defaultdict(list)

        for sequence, score in zip(normalized_sequences, standardized):
            groups[sequence[position]].append(score)

        eligible = [
            values
            for values in groups.values()
            if len(values) >= min_group_size
        ]

        if len(eligible) < 2:
            effects[position] = 0.0
            continue

        between_ss = 0.0

        for values in eligible:
            count = len(values)
            group_mean = sum(values) / count
            reliability = count / (count + shrinkage)
            shrunk_mean = overall + reliability * (group_mean - overall)
            between_ss += count * (shrunk_mean - overall) ** 2

        effects[position] = math.sqrt(
            max(0.0, min(1.0, between_ss / total_ss))
        )

    return effects


# ==========================================================
# ARISE mutation-credit engine
# ==========================================================

# The public ARISE_ENGINE name points to the mutation-credit learner.
# Its observe_generation() method accepts only candidates carrying
# explicitly validated, comparable reward metadata.
ARISE_ENGINE = MutationCreditARISE(
    reward_source="validated_target_structural",
    memory_decay=0.90,
    exploration_floor=0.25,
    exploration_c=1.0,
)
# Backward-compatible import name used by existing tests and modules.
ARISEEngine = MutationCreditARISE


# ==========================================================
# Population generation
# ==========================================================

CONTROLLED_MODES = {
    "uniform_matched",
    "mutation_credit",
}

VALID_MODES = {
    "legacy",
    "uniform_matched",
    "mutation_credit",
    "uniform",
    "hotspot_only",
    "hotspot_biased",
    "hotspot_biased_annealed",
    "hotspot_biased_annealed_exploratory",
}


def generate_candidates(
    seed_sequence: str,
    population_size: int = POPULATION_SIZE,
    mutation_rate: float = MUTATION_RATE,
    arise_enabled: bool = ARISE_ENABLED,
    cache_enabled: bool = True,
    hotspot_model: dict | None = None,
    hotspot_enabled: bool = False,
    ale_enabled: bool = False,
    mutation_guidance_mode: str = "legacy",
    generation: int = 1,
    run_id: str = "",
    preserve_sequence: str | None = None,
    total_generations: int | None = None,
    random_seed: int | None = None,
) -> list[Candidate]:
    """
    Generate and evaluate a population of candidate peptides.

    Controlled comparisons:
        uniform_matched:
            Uses uniform position probabilities with a fixed event budget.

        mutation_credit:
            Uses ARISE's learned position probabilities with the same
            fixed event budget, alphabet, and replacement rule.

    Both controlled modes differ only in the distribution used to select
    mutation positions. The mutation-credit mode learns only when valid,
    comparable parent-child structural reward metadata is supplied by the
    downstream scoring/evaluation layer.

    random_seed:
        Optional deterministic seed. Use the same seed for the baseline and
        adaptive arms of each paired experiment.
    """
    if not isinstance(seed_sequence, str) or not seed_sequence:
        raise ValueError("seed_sequence must be a nonempty string")

    seed_sequence = seed_sequence.upper()

    if population_size <= 0:
        raise ValueError("population_size must be positive")

    if not math.isfinite(float(mutation_rate)):
        raise ValueError("mutation_rate must be finite")

    if not 0.0 <= mutation_rate <= 1.0:
        raise ValueError("mutation_rate must be between 0 and 1")

    if mutation_guidance_mode not in {
        "uniform",
        "hotspot_only",
        "hotspot_biased",
        "hotspot_biased_annealed",
        "hotspot_biased_annealed_exploratory",
        "uniform_matched",
        "mutation_credit",
    }:
        raise ValueError(
            "mutation_guidance_mode must be one of: "
            "uniform, hotspot_only, hotspot_biased, "
            "hotspot_biased_annealed, "
            "hotspot_biased_annealed_exploratory, "
            "uniform_matched, mutation_credit"
        )

    if (
        mutation_guidance_mode == "mutation_credit"
        and not arise_enabled
    ):
        raise ValueError(
            "mutation_credit mode requires arise_enabled=True"
        )

    total_generations = total_generations or GENERATIONS

    if total_generations <= 0:
        raise ValueError("total_generations must be positive")

    if generation <= 0:
        raise ValueError("generation must be positive")

    mode = mutation_guidance_mode

    logger.info(
        "Generating candidates | mode=%s | generation=%d/%d | "
        "population_size=%d | mutation_rate=%.4f | random_seed=%s",
        mode,
        generation,
        total_generations,
        population_size,
        mutation_rate,
        random_seed,
    )

    # These two modes bypass ALE's hotspot resolver: they are the controlled
    # baseline/adaptive pair for the P2 experiment.
    if mode in CONTROLLED_MODES:
        policy: dict[str, Any] = {
            "requested_mode": mode,
            "mode": mode,
            "guided": mode == "mutation_credit",
            "hotspot_start": None,
            "hotspot_end": None,
            "hotspot_discovery_generation": None,
            "hotspot_probability": 0.0,
        }

        if mode == "uniform_matched":
            probabilities = None
        else:
            probabilities = ARISE_ENGINE.position_probabilities(
                len(seed_sequence)
            )

        logger.info(
            "Controlled mutation policy | mode=%s | budget=%d | "
            "observed_positions=%d | probability_range=[%.5f, %.5f]",
            mode,
            mutation_event_budget(len(seed_sequence), mutation_rate),
            len(ARISE_ENGINE.stats),
            min(probabilities) if probabilities is not None else (
                1.0 / len(seed_sequence)
            ),
            max(probabilities) if probabilities is not None else (
                1.0 / len(seed_sequence)
            ),
        )
    else:
        probabilities = None

        # Keep the legacy mode on its original native C population generator.
        # Other established ALE modes use the existing policy resolver.
        policy = resolve_mutation_policy(
            "uniform" if mode == "legacy" else mode,
            hotspot_model,
            enabled=hotspot_enabled and ale_enabled,
            generation=generation,
            total_generations=total_generations,
        )

        if mode != "legacy":
            logger.info(
                "ALE policy | requested=%s | active=%s | guided=%s | "
                "hotspot=%s:%s | hotspot_probability=%.3f",
                policy["requested_mode"],
                policy["mode"],
                policy["guided"],
                policy["hotspot_start"],
                policy["hotspot_end"],
                policy["hotspot_probability"],
            )

    # Candidate generation can be retried when duplicate removal is active.
    attempts = 0

    def generate_batch(size: int) -> list[str]:
        nonlocal attempts

        if mode in CONTROLLED_MODES:
            batch_seed = (
                None
                if random_seed is None
                else (
                    int(random_seed)
                    + generation * 1_000_003
                    + attempts
                )
            )

            return generate_position_weighted_population(
                seed_sequence,
                pop_size=size,
                mutation_rate=mutation_rate,
                probabilities=probabilities,
                seed=batch_seed,
            )

        if mode == "legacy":
            return generate_c_population(
                seed_sequence,
                pop_size=size,
                mutation_rate=mutation_rate,
            )

        return generate_policy_population(
            seed_sequence,
            mutation_rate=mutation_rate,
            guidance_mode=policy["mode"],
            hotspot_start=policy["hotspot_start"],
            hotspot_end=policy["hotspot_end"],
            pop_size=size,
            hotspot_probability=policy["hotspot_probability"],
        )

    # ------------------------------------------------------
    # Generate population and remove duplicates if required
    # ------------------------------------------------------

    sequences: list[str] = []
    seen: set[str] = set()

    while len(sequences) < population_size:
        requested_size = population_size - len(sequences)
        batch = generate_batch(requested_size)

        if not REMOVE_DUPLICATE_CANDIDATES:
            sequences.extend(batch)
            break

        before = len(sequences)

        for sequence in batch:
            if sequence not in seen:
                seen.add(sequence)
                sequences.append(sequence)

                if len(sequences) >= population_size:
                    break

        attempts += 1

        if len(sequences) >= population_size:
            break

        if attempts >= DIVERSITY_REFILL_ATTEMPTS:
            logger.warning(
                "Could only generate %d/%d unique candidates after %d "
                "attempts.",
                len(sequences),
                population_size,
                attempts,
            )
            break

        # Avoid spinning indefinitely if a generator returns no new sequence.
        if before == len(sequences):
            logger.debug(
                "Duplicate refill made no progress on attempt %d",
                attempts,
            )

    # ------------------------------------------------------
    # Preserve the requested elite, if applicable
    # ------------------------------------------------------

    preserved_elite = False

    if (
        preserve_sequence
        and preserve_sequence not in sequences
        and sequences
    ):
        sequences[-1] = preserve_sequence.upper()
        preserved_elite = True

    # ------------------------------------------------------
    # Evaluate candidates and attach mutation metadata
    # ------------------------------------------------------

    candidates: list[Candidate] = []

    controlled_budget = mutation_event_budget(
        len(seed_sequence),
        mutation_rate,
    )

    for sequence in sequences:
        sequence = sequence.upper()

        candidate = process_candidate_peptide(
            sequence,
            cache_enabled=cache_enabled,
        )

        if candidate is None:
            logger.warning(
                "Candidate processing returned None; candidate skipped."
            )
            continue

        if not isinstance(candidate.metadata, dict):
            candidate.metadata = {}

        changed_positions = mutation_positions(seed_sequence, sequence)

        is_preserved = (
            preserved_elite
            and sequence == (preserve_sequence or "").upper()
        )

        if is_preserved:
            event_budget = 0
        elif mode in CONTROLLED_MODES:
            event_budget = controlled_budget
        else:
            # Preserve the recorded convention used by the existing modes.
            event_budget = round(mutation_rate * len(seed_sequence))

        candidate.metadata["changed_positions_count"] = len(
            changed_positions
        )
        candidate.metadata["mutation_event_budget"] = event_budget

        candidate.metadata["mutation_event_count"] = (
            event_budget
            if mode in CONTROLLED_MODES or mode != "legacy"
            else len(changed_positions)
        )

        candidate.metadata["mutation_count"] = len(changed_positions)
        candidate.metadata["mutation_positions"] = changed_positions

        # The seed sequence is the direct parent for this population.
        # Parent and child structural scores are intentionally not fabricated
        # here; the downstream evaluator must add them when valid.
        candidate.metadata["parent_sequence"] = seed_sequence
        candidate.metadata["run_id"] = run_id
        candidate.metadata["generation"] = generation

        candidate.metadata["mutation_guidance_mode"] = policy["mode"]
        candidate.metadata["requested_mutation_guidance_mode"] = (
            policy["requested_mode"]
        )
        candidate.metadata["hotspot_guided"] = policy["guided"]
        candidate.metadata["hotspot_probability"] = policy[
            "hotspot_probability"
        ]
        candidate.metadata["hotspot_discovery_generation"] = policy[
            "hotspot_discovery_generation"
        ]
        candidate.metadata["hotspot_start"] = policy["hotspot_start"]
        candidate.metadata["hotspot_end"] = policy["hotspot_end"]

        candidate.metadata["is_duplicate"] = False
        candidate.metadata["preserved_elite"] = is_preserved

        candidate.metadata["mutation_reward_source"] = None
        candidate.metadata["mutation_parent_score"] = None
        candidate.metadata["mutation_child_score"] = None
        candidate.metadata["ranking_observation_valid"] = False

        candidates.append(candidate)

    logger.info(
        "Generated %d candidates (%d unique) | mode=%s | budget=%d",
        len(candidates),
        len({candidate.sequence for candidate in candidates}),
        mode,
        controlled_budget,
    )

    return candidates


# ==========================================================
# Debug
# ==========================================================

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    seed = "ACDEFGHIKLMNPQRSTVWYA"

    population = generate_candidates(
        seed,
        mutation_guidance_mode="uniform_matched",
        random_seed=1234,
        generation=1,
    )

    for candidate in population:
        print("-" * 40)
        print(candidate.sequence)
        print(candidate.metadata)

    print("\nARISE mutation-credit position probabilities:")
    print(ARISE_ENGINE.importance_map())