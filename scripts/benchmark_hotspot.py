"""Synthetic ARISE hotspot benchmark with no project-specific biology."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from src.hotspot import discover_hotspot_model


def synthetic_history(
    population_size: int,
    generations: int,
    hotspot_start: int,
    hotspot_pattern: str,
) -> list[dict]:
    history = []
    for generation in range(1, generations + 1):
        for index in range(population_size):
            sequence = list("A" * 21)
            if index < population_size // 2:
                sequence[hotspot_start:hotspot_start + len(hotspot_pattern)] = list(
                    hotspot_pattern
                )
                score = 1.0 + (index / max(1, population_size))
            else:
                sequence[hotspot_start:hotspot_start + len(hotspot_pattern)] = list(
                    "C" * len(hotspot_pattern)
                )
                score = 0.1 + (index / max(1, population_size))
            history.append({
                "sequence": "".join(sequence),
                "score": score,
                "generation": generation,
                "run_id": "synthetic-hotspot-benchmark",
                "arise_score_comparable_across_generations": True,
            })
    return history


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--population-size", type=int, default=32)
    parser.add_argument("--generations", type=int, default=4)
    parser.add_argument("--hotspot-start", type=int, default=8)
    parser.add_argument("--pattern", default="KLMN")
    parser.add_argument("--bootstrap-iterations", type=int, default=100)
    parser.add_argument("--seed", type=int, default=616)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    history = synthetic_history(
        args.population_size,
        args.generations,
        args.hotspot_start,
        args.pattern,
    )
    model = discover_hotspot_model(
        history,
        discovery_generation=args.generations,
        run_id="synthetic-hotspot-benchmark",
        bootstrap_iterations=args.bootstrap_iterations,
        random_seed=args.seed,
    )
    result = {
        "configuration": vars(args),
        "selected_hotspot": model.get("selected_hotspot"),
        "correlation_length": model.get("correlation_length"),
        "model": model,
    }
    text = json.dumps(result, indent=2, default=str) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
