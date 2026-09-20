"""Small deterministic before/after benchmark for pipeline optimizations.

This benchmark keeps the generated candidate sequences identical between the
two cases. The baseline disables deterministic caches and workers; the
optimized case enables caches and bounded candidate workers. It intentionally
does not fabricate ESMFold or Vina timings: those external stages are reported
by the normal experiment audit when enabled.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import statistics
import sys
import time

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from c.bridge import (
    clear_evaluation_cache,
    evaluation_cache_stats,
    generate_c_population,
    process_candidate_peptide,
    seed_native_random,
)
from src.docking_vina import MLSurrogateBackend, PeptideDockingScorer


def make_inputs(
    seed_sequence: str,
    population_size: int,
    generations: int,
    mutation_rate: float,
    random_seed: int,
) -> list[list[str]]:
    seed_native_random(random_seed)
    inputs = []
    current_seed = seed_sequence
    for _generation in range(generations):
        population = generate_c_population(
            current_seed,
            pop_size=population_size,
            mutation_rate=mutation_rate,
        )
        inputs.append(population)
        # Keep input construction identical for both benchmark cases.
        current_seed = population[0]
    return inputs


def run_case(
    inputs: list[list[str]],
    cache_enabled: bool,
    workers: int,
) -> dict:
    clear_evaluation_cache()
    backend = MLSurrogateBackend(cache_enabled=cache_enabled)
    scorer = PeptideDockingScorer(
        backend=backend,
        workers=workers,
        cache_enabled=cache_enabled,
    )
    scorer.vina_backend = None
    generation_times = []
    candidate_count = 0
    start = time.perf_counter()
    cpu_start = time.process_time()

    for sequences in inputs:
        generation_start = time.perf_counter()
        candidates = [
            process_candidate_peptide(
                sequence,
                cache_enabled=cache_enabled,
            )
            for sequence in sequences
        ]
        scorer.evaluate_candidates(candidates)
        generation_times.append(time.perf_counter() - generation_start)
        candidate_count += len(candidates)

    wall = time.perf_counter() - start
    cpu = time.process_time() - cpu_start
    return {
        "cache_enabled": cache_enabled,
        "workers": workers,
        "candidate_count": candidate_count,
        "total_runtime_seconds": round(wall, 6),
        "runtime_per_candidate_seconds": round(wall / candidate_count, 8),
        "runtime_per_generation_seconds": [
            round(value, 6) for value in generation_times
        ],
        "median_generation_seconds": round(
            statistics.median(generation_times),
            6,
        ),
        "cpu_seconds_since_process_start": round(cpu, 6),
        "cpu_utilization_estimate": round(cpu / wall, 3) if wall else None,
        "c_fitness_cache": evaluation_cache_stats(),
        "docking_cache": scorer.cache_stats(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed-sequence", required=True)
    parser.add_argument("--population-size", type=int, default=16)
    parser.add_argument("--generations", type=int, default=3)
    parser.add_argument("--mutation-rate", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=1001)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    inputs = make_inputs(
        args.seed_sequence,
        args.population_size,
        args.generations,
        args.mutation_rate,
        args.seed,
    )
    result = {
        "configuration": {
            "population_size": args.population_size,
            "generations": args.generations,
            "mutation_rate": args.mutation_rate,
            "random_seed": args.seed,
            "same_inputs": True,
        },
        "baseline": run_case(inputs, cache_enabled=False, workers=1),
        "optimized": run_case(
            inputs,
            cache_enabled=True,
            workers=max(1, args.workers),
        ),
    }
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
