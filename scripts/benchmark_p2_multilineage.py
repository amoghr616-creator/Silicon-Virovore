
from __future__ import annotations

import csv
import json
import random
import statistics
from pathlib import Path

import numpy as np


DIMENSIONS = 16
BUDGET = 400
INSTANCES = 20
SEEDS = 5
POOL_SIZE = 8


def make_landscape(kind: str, instance_seed: int):
    """Create a reproducible synthetic maximization problem."""
    rng = np.random.default_rng(instance_seed)
    weights = rng.uniform(0.2, 1.0, DIMENSIONS)
    pairs = [
        (i, j, float(rng.uniform(0.2, 0.8)))
        for i in range(DIMENSIONS)
        for j in range(i + 1, DIMENSIONS)
        if rng.random() < 0.12
    ]

    def objective(x):
        x = np.asarray(x, dtype=np.int8)

        if kind == "smooth":
            return float(np.dot(weights, x))

        if kind == "deceptive":
            # Four deceptive blocks of four bits.
            total = 0.0
            for start in range(0, DIMENSIONS, 4):
                ones = int(x[start:start + 4].sum())
                total += 4.0 if ones == 4 else 3.0 - ones
            return total

        if kind == "interactions":
            base = float(np.dot(weights, x))
            synergy = sum(
                strength * int(x[i]) * int(x[j])
                for i, j, strength in pairs
            )
            return base + synergy

        if kind == "noisy":
            # The true landscape is smooth; observations are noisy.
            return float(np.dot(weights, x))

        raise ValueError(f"Unknown landscape: {kind}")

    def evaluate(x, observation_rng):
        true_score = objective(x)
        observed = true_score
        if kind == "noisy":
            observed += float(observation_rng.normal(0, 1.0))
        return observed, true_score

    return evaluate


def random_state(rng):
    return tuple(rng.randrange(2) for _ in range(DIMENSIONS))


def mutate(state, rng):
    """One reversible bit flip."""
    child = list(state)
    position = rng.randrange(DIMENSIONS)
    child[position] = 1 - child[position]
    return tuple(child)


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def run_single_parent(evaluate, seed):
    """Baseline: keep and mutate one current state."""
    rng = random.Random(seed)
    obs_rng = np.random.default_rng(seed + 100_000)

    current = random_state(rng)
    current_score, _ = evaluate(current, obs_rng)
    best_state = current
    best_observed = current_score

    evaluations = 1

    while evaluations < BUDGET:
        candidate = mutate(current, rng)
        score, _ = evaluate(candidate, obs_rng)
        evaluations += 1

        if score >= current_score:
            current = candidate
            current_score = score

        if score > best_observed:
            best_state = candidate
            best_observed = score

    # In a simulation, inspect the known objective without charging
    # another search evaluation. For noisy landscapes, this is the
    # true score of the state selected using noisy observations.
    _, best_true = evaluate(
        best_state, np.random.default_rng(seed + 200_000)
    )

    return {
        "best_observed": best_observed,
        "best_true": best_true,
        "evaluations": evaluations,
    }


def run_multilineage(evaluate, seed):
    """Maintain a small pool of distinct candidate states."""
    rng = random.Random(seed)
    obs_rng = np.random.default_rng(seed + 100_000)

    pool = []
    seen = set()
    evaluations = 0

    # Initial pool evaluations count against the same total budget.
    while len(pool) < POOL_SIZE and evaluations < BUDGET:
        state = random_state(rng)
        if state in seen:
            continue
        seen.add(state)

        score, _ = evaluate(state, obs_rng)
        pool.append((state, score))
        evaluations += 1

    best_state, best_observed = max(pool, key=lambda item: item[1])

    while evaluations < BUDGET:
        # Prefer stronger candidates, but select among several
        # distinct parents rather than only the current winner.
        ranked = sorted(pool, key=lambda item: item[1], reverse=True)
        parent_choices = ranked[:max(2, len(ranked) // 2)]
        parent, _ = rng.choice(parent_choices)

        child = mutate(parent, rng)

        score, _ = evaluate(child, obs_rng)
        evaluations += 1

        if score > best_observed:
            best_state = child
            best_observed = score

        if child not in {state for state, _ in pool}:
            distances = [
                hamming(child, state) for state, _ in pool
            ]
            sufficiently_distinct = min(distances) >= 2
            worst_index = min(
                range(len(pool)), key=lambda i: pool[i][1]
            )
            pool_scores = [item[1] for item in pool]
            at_least_median = score >= statistics.median(pool_scores)

            if score > pool[worst_index][1]:
                pool[worst_index] = (child, score)
            elif sufficiently_distinct and at_least_median:
                pool[worst_index] = (child, score)

    _, best_true = evaluate(
        best_state, np.random.default_rng(seed + 200_000)
    )

    return {
        "best_observed": best_observed,
        "best_true": best_true,
        "evaluations": evaluations,
    }


def mean_ci(values):
    """Approximate 95% CI across independent landscape instances."""
    values = list(values)
    mean = statistics.mean(values)
    if len(values) < 2:
        return mean, None, None

    # Normal approximation; report the small-sample limitation.
    se = statistics.stdev(values) / np.sqrt(len(values))
    return mean, mean - 1.96 * se, mean + 1.96 * se


def main():
    output = Path("results/p2_multilineage_v1")
    output.mkdir(parents=True, exist_ok=True)

    families = ["smooth", "deceptive", "interactions", "noisy"]
    rows = []

    for family_index, family in enumerate(families):
        for instance in range(INSTANCES):
            instance_seed = 10_000 + family_index * 1_000 + instance
            evaluate = make_landscape(family, instance_seed)

            for seed_index in range(SEEDS):
                seed = 50_000 + instance * 100 + seed_index

                for method, runner in [
                    ("single_parent", run_single_parent),
                    ("multi_lineage", run_multilineage),
                ]:
                    result = runner(evaluate, seed)
                    rows.append({
                        "family": family,
                        "instance": instance,
                        "seed": seed,
                        "method": method,
                        **result,
                    })

    csv_path = output / "raw_results.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    summaries = []

    for family in families:
        family_rows = [r for r in rows if r["family"] == family]

        # Average seeds within each instance first. The independent
        # landscape instance—not each seed—is the replication unit.
        instance_differences = []

        for instance in range(INSTANCES):
            baseline = [
                r["best_true"] for r in family_rows
                if r["instance"] == instance
                and r["method"] == "single_parent"
            ]
            multi = [
                r["best_true"] for r in family_rows
                if r["instance"] == instance
                and r["method"] == "multi_lineage"
            ]
            instance_differences.append(
                statistics.mean(multi) - statistics.mean(baseline)
            )

        mean, low, high = mean_ci(instance_differences)

        summaries.append({
            "family": family,
            "instances": INSTANCES,
            "seeds_per_instance": SEEDS,
            "budget_per_run": BUDGET,
            "mean_paired_true_score_difference": mean,
            "approx_ci95_low": low,
            "approx_ci95_high": high,
            "instances_multi_lineage_won": sum(
                d > 0 for d in instance_differences
            ),
            "instances_tied": sum(
                d == 0 for d in instance_differences
            ),
        })

    with (output / "summary.json").open("w") as f:
        json.dump(summaries, f, indent=2)

    print(json.dumps(summaries, indent=2))
    print(f"\nRaw results: {csv_path}")
    print(f"Summary: {output / 'summary.json'}")


if __name__ == "__main__":
    main()
