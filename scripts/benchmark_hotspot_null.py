"""Null calibration and power benchmark for the ARISE hotspot permutation test.

The structured null includes shared run/generation effects and inherited
lineages. Lineage score effects are generated independently of sequence
content, so there is no causal position-specific sequence effect; however,
lineage is a common cause of sequence pattern and score. This intentionally
challenges the detector's within-run/generation exchangeability assumption.

Synthetic results assess statistical behavior only, not biological validity.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import json
import math
from pathlib import Path
import random
import statistics
import sys
import time

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from src.hotspot import (
    _split_history,
    _window_association,
    discover_hotspot_model,
)


AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
TARGET_WINDOW = (4, 7)
TARGET_PATTERN = "KLM"
SCENARIOS = (
    "exchangeable_null",
    "structured_null",
)
DEFAULT_EFFECT_SIZES = (0.4, 0.8, 1.2)
PERMUTATION_METHODS = (
    "candidate_shuffle",
    "lineage_generation_block_shuffle",
)
LINEAGE_BLOCK_FIELDS = (
    "run_id",
    "generation",
    "synthetic_lineage_id",
)


def generate_synthetic_history(
    scenario: str,
    *,
    seed: int,
    effect_size: float = 0.0,
    run_count: int = 8,
    generations: int = 4,
    candidates_per_generation: int = 12,
    sequence_length: int = 12,
) -> tuple[list[dict], tuple[int, int] | None]:
    """Generate production-schema hotspot records under a named scenario.

    ``exchangeable_null``: sequences are IID across positions and scores are
    independently sampled inside each run/generation stratum.

    ``structured_null``: records share inherited founder sequences by lineage;
    run effects, generation trends, and lineage-level score offsets create
    dependence. The offsets are sampled independently of sequence content, so
    sequence positions have no causal effect, but lineage can confound their
    association with score.

    ``injected_hotspot``: sequences are otherwise IID, with a controlled
    motif at ``TARGET_WINDOW`` in a subset of candidates and an additive score
    effect for that motif. ``effect_size`` controls that effect.
    """
    allowed = {*SCENARIOS, "injected_hotspot"}
    if scenario not in allowed:
        raise ValueError(f"unknown scenario: {scenario}")
    if run_count < 2 or generations < 1 or candidates_per_generation < 4:
        raise ValueError("need >=2 runs, >=1 generation, and >=4 candidates per generation")
    if sequence_length < TARGET_WINDOW[1]:
        raise ValueError("sequence_length is too short for the target window")
    if not math.isfinite(effect_size) or effect_size < 0:
        raise ValueError("effect_size must be finite and nonnegative")

    rng = random.Random(seed)
    records: list[dict] = []
    lineage_count = 4
    run_offsets = {
        run_index: rng.gauss(0.0, 0.8)
        for run_index in range(run_count)
    }
    lineage_offsets = {
        (run_index, lineage): rng.gauss(0.0, 0.65)
        for run_index in range(run_count)
        for lineage in range(lineage_count)
    }
    founders = {
        (run_index, lineage): [
            rng.choice(AMINO_ACIDS) for _ in range(sequence_length)
        ]
        for run_index in range(run_count)
        for lineage in range(lineage_count)
    }

    for run_index in range(run_count):
        lineage_ancestors = {
            lineage: list(founders[(run_index, lineage)])
            for lineage in range(lineage_count)
        }
        for generation in range(1, generations + 1):
            generation_offset = 0.35 * (generation - 1)
            next_ancestors = dict(lineage_ancestors)
            for candidate_index in range(candidates_per_generation):
                lineage = candidate_index % lineage_count
                parent_sequence = "".join(lineage_ancestors[lineage])
                if scenario == "structured_null":
                    sequence = list(lineage_ancestors[lineage])
                    # Descendants inherit the previous generation's lineage
                    # ancestor. The score equation below never uses sequence.
                    mutation_probability = 0.02 * generation
                    for position in range(sequence_length):
                        if rng.random() < mutation_probability:
                            sequence[position] = rng.choice(AMINO_ACIDS)
                    next_ancestors[lineage] = sequence
                    score = (
                        run_offsets[run_index]
                        + generation_offset
                        + lineage_offsets[(run_index, lineage)]
                        + rng.gauss(0.0, 0.35)
                    )
                else:
                    sequence = [rng.choice(AMINO_ACIDS) for _ in range(sequence_length)]
                    if scenario == "exchangeable_null":
                        # Stratum offsets are shared; residual scores are
                        # exchangeable among records in the same run/gen.
                        score = (
                            run_offsets[run_index]
                            + generation_offset
                            + rng.gauss(0.0, 1.0)
                        )
                    else:
                        has_motif = rng.random() < 0.5
                        if has_motif:
                            sequence[TARGET_WINDOW[0]:TARGET_WINDOW[1]] = TARGET_PATTERN
                        else:
                            # Use the same motif alphabet but not the target
                            # triplet, yielding a controlled local signal.
                            sequence[TARGET_WINDOW[0]:TARGET_WINDOW[1]] = list("AAA")
                        score = (
                            run_offsets[run_index]
                            + generation_offset
                            + effect_size * has_motif
                            + rng.gauss(0.0, 1.0)
                        )

                records.append({
                    "candidate_id": (
                        f"run-{run_index:03d}-gen-{generation:03d}-"
                        f"candidate-{candidate_index:03d}"
                    ),
                    "sequence": "".join(sequence),
                    "parent_sequence": (
                        parent_sequence if scenario == "structured_null" else None
                    ),
                    "score": float(score),
                    "generation": generation,
                    "run_id": f"synthetic-run-{run_index:03d}",
                    "arise_score_comparable_across_generations": True,
                    "synthetic_lineage_id": f"run-{run_index:03d}-lineage-{lineage}",
                    "synthetic_candidate_id": (
                        f"run-{run_index:03d}-gen-{generation:03d}-"
                        f"candidate-{candidate_index:03d}"
                    ),
                })
            lineage_ancestors = next_ancestors

    target_window = TARGET_WINDOW if scenario == "injected_hotspot" else None
    return records, target_window


def wilson_interval(successes: int, trials: int, z: float = 1.959963984540054) -> tuple[float, float]:
    """Return the two-sided Wilson score interval for a binomial rate."""
    if trials < 1 or not 0 <= successes <= trials:
        raise ValueError("require trials >= 1 and 0 <= successes <= trials")
    proportion = successes / trials
    denominator = 1.0 + z * z / trials
    center = (proportion + z * z / (2.0 * trials)) / denominator
    half_width = (
        z
        * math.sqrt(
            proportion * (1.0 - proportion) / trials
            + z * z / (4.0 * trials * trials)
        )
        / denominator
    )
    return max(0.0, center - half_width), min(1.0, center + half_width)


def interval_overlap(
    selected: tuple[int, int] | None,
    target: tuple[int, int] | None,
) -> float | None:
    """Calculate interval intersection-over-union for half-open windows."""
    if selected is None or target is None:
        return None
    intersection = max(0, min(selected[1], target[1]) - max(selected[0], target[0]))
    union = max(selected[1], target[1]) - min(selected[0], target[0])
    return intersection / union if union else 0.0


def _selected_window(model: dict) -> tuple[int, int] | None:
    selected = model.get("selected_hotspot")
    if not selected:
        return None
    return int(selected["start_position"]), int(selected["end_position"])


def _permute_scores_within_blocks(
    records: list[dict],
    *,
    seed: int,
    block_fields: tuple[str, ...],
) -> list[dict]:
    """Shuffle scores only within declared exchangeability blocks."""
    blocks: dict[tuple[str, ...], list[int]] = defaultdict(list)
    for index, record in enumerate(records):
        key = tuple(str(record.get(field)) for field in block_fields)
        blocks.setdefault(key, []).append(index)

    rng = random.Random(seed)
    permuted = [dict(record) for record in records]
    for indices in blocks.values():
        scores = [records[index]["score"] for index in indices]
        rng.shuffle(scores)
        for index, score in zip(indices, scores):
            permuted[index]["score"] = score
    return permuted


def _blocked_max_statistic_p_value(
    records: list[dict],
    *,
    iterations: int,
    seed: int,
    block_fields: tuple[str, ...],
    window_min: int = 2,
    window_max: int = 6,
    top_quantile: float = 0.25,
) -> float | None:
    """Maximum-window permutation conditional on fixed block membership.

    For lineage blocks this tests for within-lineage/generation sequence-score
    association, conditional on each block's score distribution. It is valid
    only if candidate residual scores are exchangeable within every declared
    block; it does not test marginal between-lineage association.
    """
    if iterations < 1:
        raise ValueError("iterations must be >= 1")
    if not records or any(
        any(record.get(field) is None for field in block_fields)
        for record in records
    ):
        return None

    sequence_length = len(records[0]["sequence"])
    lengths = range(max(2, window_min), min(sequence_length, window_max) + 1)
    starts = range(0, max(0, sequence_length - max(2, window_min) + 1))

    def maximum_score(sample: list[dict]) -> float | None:
        values = [
            association["hotspot_score"]
            for length in lengths
            for start in starts
            if (
                association := _window_association(
                    sample,
                    start,
                    length,
                    top_quantile,
                )
            ) is not None
        ]
        return max(values) if values else None

    observed = maximum_score(records)
    if observed is None:
        return None

    exceedances = 0
    valid_permutations = 0
    for iteration in range(iterations):
        permuted = _permute_scores_within_blocks(
            records,
            seed=seed + iteration,
            block_fields=block_fields,
        )
        permuted_maximum = maximum_score(permuted)
        if permuted_maximum is not None:
            valid_permutations += 1
            exceedances += permuted_maximum >= observed
    if valid_permutations == 0:
        return None
    return (1 + exceedances) / (1 + valid_permutations)


def _run_trial(
    scenario: str,
    *,
    trial_seed: int,
    permutation_iterations: int,
    effect_size: float = 0.0,
    run_count: int = 8,
    generations: int = 4,
    candidates_per_generation: int = 12,
    sequence_length: int = 12,
) -> dict:
    records, target = generate_synthetic_history(
        scenario,
        seed=trial_seed,
        effect_size=effect_size,
        run_count=run_count,
        generations=generations,
        candidates_per_generation=candidates_per_generation,
        sequence_length=sequence_length,
    )
    model = discover_hotspot_model(
        records,
        discovery_generation=generations,
        run_id="synthetic-run-000",
        bootstrap_iterations=0,
        random_seed=trial_seed,
        permutation_iterations=permutation_iterations,
        permutation_alpha=0.05,
        require_significant_hotspot=False,
    )
    discovery_records, _holdout_records = _split_history(records, trial_seed)
    lineage_p_value = _blocked_max_statistic_p_value(
        discovery_records,
        iterations=permutation_iterations,
        seed=trial_seed + 7919,
        block_fields=LINEAGE_BLOCK_FIELDS,
    )
    candidate_p_value = model.get("hotspot_permutation_p_value")
    selected = _selected_window(model)
    overlap = interval_overlap(selected, target)
    return {
        "scenario": scenario,
        "trial_seed": trial_seed,
        "effect_size": effect_size,
        "candidate_shuffle_p_value": candidate_p_value,
        "candidate_shuffle_significant": (
            candidate_p_value is not None and candidate_p_value < 0.05
        ),
        "lineage_generation_block_shuffle_p_value": lineage_p_value,
        "lineage_generation_block_shuffle_significant": (
            lineage_p_value is not None and lineage_p_value < 0.05
        ),
        "selected_window": selected,
        "target_window": target,
        "overlap_iou": overlap,
        "localized": overlap is not None and overlap >= 0.5,
        "status": model.get("status"),
        "discovery_size": model.get("discovery_sample_size_train"),
        "holdout_size": model.get("holdout_sample_size"),
    }


def _summarize_trials(
    trials: list[dict],
    *,
    target_present: bool,
    permutation_method: str,
) -> dict:
    total = len(trials)
    significance_field = f"{permutation_method}_significant"
    significant_count = sum(trial[significance_field] for trial in trials)
    interval = wilson_interval(significant_count, total)
    result = {
        "trials": total,
        "significant_count": significant_count,
        "rate": significant_count / total,
        "rate_ci95": interval,
    }
    if target_present:
        localized_count = sum(
            trial[significance_field] and trial["localized"]
            for trial in trials
        )
        power_ci = wilson_interval(localized_count, total)
        overlaps = [trial["overlap_iou"] for trial in trials]
        result.update({
            "detection_power": significant_count / total,
            "detection_power_ci95": interval,
            "localized_detection_count": localized_count,
            "localized_detection_power": localized_count / total,
            "localized_detection_power_ci95": power_ci,
            "mean_selected_window_iou": statistics.mean(overlaps),
            "mean_iou_among_significant": (
                statistics.mean(
                    trial["overlap_iou"]
                    for trial in trials
                    if trial[significance_field]
                )
                if significant_count
                else None
            ),
            "localization_rate_given_detection": (
                localized_count / significant_count
                if significant_count
                else None
            ),
        })
    return result


def run_benchmark(
    *,
    mode: str = "quick",
    trials: int | None = None,
    permutations: int | None = None,
    seed: int | None = None,
) -> dict:
    """Run reproducible null calibration and injected-effect power tests."""
    if mode not in {"quick", "full"}:
        raise ValueError("mode must be 'quick' or 'full'")
    seed = seed if seed is not None else (
        20261015 if mode == "quick" else 20261016
    )
    trial_count = trials if trials is not None else (8 if mode == "quick" else 40)
    permutation_count = permutations if permutations is not None else (24 if mode == "quick" else 199)
    if trial_count < 1 or permutation_count < 1:
        raise ValueError("trials and permutations must be >= 1")

    dimensions = (
        {
            "run_count": 6,
            "generations": 3,
            "candidates_per_generation": 8,
            "sequence_length": 10,
        }
        if mode == "quick"
        else {
            "run_count": 8,
            "generations": 4,
            "candidates_per_generation": 12,
            "sequence_length": 12,
        }
    )

    started = time.perf_counter()
    result_groups: list[dict] = []
    trial_rows: list[dict] = []
    scenarios: list[tuple[str, float]] = [
        ("exchangeable_null", 0.0),
        ("structured_null", 0.0),
        *(('injected_hotspot', effect) for effect in DEFAULT_EFFECT_SIZES),
    ]
    for scenario_index, (scenario, effect_size) in enumerate(scenarios):
        scenario_trials = []
        for trial_index in range(trial_count):
            trial_seed = seed + scenario_index * 1_000_000 + trial_index * 10_007
            trial = _run_trial(
                scenario,
                trial_seed=trial_seed,
                permutation_iterations=permutation_count,
                effect_size=effect_size,
                **dimensions,
            )
            scenario_trials.append(trial)
            trial_rows.append(trial)
        summary = {
            "scenario": scenario,
            "effect_size": effect_size,
        }
        result_groups.extend(
            {
                **summary,
                "permutation_method": permutation_method,
                **_summarize_trials(
                    scenario_trials,
                    target_present=scenario == "injected_hotspot",
                    permutation_method=permutation_method,
                ),
            }
            for permutation_method in PERMUTATION_METHODS
        )

    runtime_seconds = time.perf_counter() - started
    return {
        "benchmark": "hotspot-null-calibration-and-power",
        "benchmark_version": 1,
        "mode": mode,
        "base_seed": seed,
        "trial_count_per_scenario": trial_count,
        "permutation_iterations_per_trial": permutation_count,
        "permutation_alpha": 0.05,
        "permutation_methods": list(PERMUTATION_METHODS),
        "lineage_block_fields": list(LINEAGE_BLOCK_FIELDS),
        "lineage_block_assumption": (
            "Scores are exchangeable among candidate residuals within each run/generation/lineage block under the conditional null; this tests within-lineage association and preserves between-lineage score/sequence structure."
        ),
        "production_lineage_schema_status": (
            "Current history stores candidate_id and parent_sequence but no unambiguous lineage/root identifier; the synthetic lineage-block method cannot be applied to production histories until lineage provenance is recorded."
        ),
        "confidence_interval": "95% Wilson score interval",
        "target_window_half_open": TARGET_WINDOW,
        "target_pattern": TARGET_PATTERN,
        **dimensions,
        "trial_seed_schedule": (
            "seed + scenario_index * 1,000,000 + trial_index * 10,007"
        ),
        "runtime_seconds": runtime_seconds,
        "scenario_generation_notes": {
            "exchangeable_null": (
                "IID sequences; independent score residuals within each run/generation stratum, with shared stratum offsets."
            ),
            "structured_null": (
                "Inherited lineage ancestor sequences, correlated descendants, run offsets, generation trend, and lineage score effects sampled independently of sequence content. No sequence position causes score; lineage is a common-cause confounder."
            ),
            "injected_hotspot": (
                "IID background sequences with the KLM motif injected at positions 4:7 in half of records; motif presence adds the configured score effect."
            ),
            "interpretation": (
                "Passing the exchangeable null alone does not establish calibration under dependence. These are synthetic statistical diagnostics, not biological validation."
            ),
        },
        "summaries": result_groups,
        "trials": trial_rows,
    }


def write_outputs(payload: dict, output_dir: Path) -> None:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"refusing to overwrite non-empty directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "raw_results.json").write_text(json.dumps(payload, indent=2) + "\n")
    with (output_dir / "summary.csv").open("w", newline="") as handle:
        rows = payload["summaries"]
        fieldnames = list(dict.fromkeys(
            field
            for row in rows
            for field in row
        ))
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "# Hotspot permutation null calibration and power",
        "",
        f"Mode: {payload['mode']}; datasets per scenario: {payload['trial_count_per_scenario']}; "
        f"permutations per dataset: {payload['permutation_iterations_per_trial']}; "
        f"base seed: {payload['base_seed']}; runtime: {payload['runtime_seconds']:.2f}s.",
        f"Synthetic dimensions: {payload['run_count']} runs, {payload['generations']} generations, "
        f"{payload['candidates_per_generation']} candidates/generation, sequence length {payload['sequence_length']}.",
        "False-positive rate and detection power intervals are 95% Wilson score intervals.",
        "",
        "| Scenario | Permutation | Effect | Significant | Rate / power | 95% CI | Localized power | Mean window IoU |",
        "|---|---|---:|---:|---:|---|---:|---:|",
    ]
    for row in payload["summaries"]:
        rate_name = "detection_power" if "detection_power" in row else "rate"
        interval_name = "detection_power_ci95" if "detection_power_ci95" in row else "rate_ci95"
        ci = row[interval_name]
        localized = row.get("localized_detection_power")
        mean_iou = row.get("mean_selected_window_iou")
        lines.append(
            f"| {row['scenario']} | {row['permutation_method']} "
            f"| {row['effect_size']:.2f} "
            f"| {row['significant_count']}/{row['trials']} "
            f"| {row[rate_name]:.3f} "
            f"| [{ci[0]:.3f}, {ci[1]:.3f}] "
            f"| {'n/a' if localized is None else f'{localized:.3f}'} "
            f"| {'n/a' if mean_iou is None else f'{mean_iou:.3f}'} |"
        )
    lines.extend([
        "",
        "Structured-null scores have run, generation, and lineage components; lineage sequence patterns are inherited while lineage score effects are generated independently of sequence content. Thus there is no causal sequence-position effect, but observational lineage confounding violates candidate-level exchangeability.",
        "The candidate shuffle destroys sequence-lineage/score-lineage alignment. The restricted shuffle fixes run/generation/lineage membership and permutes scores only among siblings; this is conditional calibration, requiring exchangeable within-block residuals, and does not test between-lineage association.",
        "The production history schema currently lacks an unambiguous lineage/root identifier, so this synthetic block test is not currently deployable on production histories.",
        "A high structured-null rate is evidence against calibration for dependent histories. A low rate in the exchangeable null alone is not sufficient evidence of calibration.",
        "No result here establishes biological validity or peptide-binding performance.",
    ])
    (output_dir / "report.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("quick", "full"), default="quick")
    parser.add_argument("--trials", type=int)
    parser.add_argument("--permutations", type=int)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()

    payload = run_benchmark(
        mode=args.mode,
        trials=args.trials,
        permutations=args.permutations,
        seed=args.seed,
    )
    if args.output_dir:
        write_outputs(payload, args.output_dir)
        print(f"Wrote {args.output_dir / 'raw_results.json'}")
        print(f"Wrote {args.output_dir / 'summary.csv'}")
        print(f"Wrote {args.output_dir / 'report.md'}")
        print((args.output_dir / "report.md").read_text(), end="")
    else:
        print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
