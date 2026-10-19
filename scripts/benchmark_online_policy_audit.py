"""Synthetic integration benchmark for online Idea 3 matched audits.

Methods:
- baseline
- fixed_adaptive
- online_audit

The signal is selected from development landscapes only. This isolates
online policy auditing; independent signal validation (Idea 1) is not
part of this comparison yet.
Every test run includes its allocated share of discovery evaluations,
main-search evaluations, and online-audit evaluations in one fixed budget.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import random
import statistics
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.benchmark_combined_adaptive_search import (
    FAMILY_OFFSETS,
    SyntheticLandscape,
    make_landscape,
)
from src.online_policy_audit import run_online_audit


MODES = ("baseline", "fixed_adaptive", "online_audit")


def discover_signal(
    landscapes: list[SyntheticLandscape],
    *,
    seed: int,
) -> tuple[int, int, dict[str, float]]:
    """Choose the best average one-bit signal using development landscapes."""
    if not landscapes:
        raise ValueError("at least one development landscape is required")

    feature_count = landscapes[0].feature_count
    if any(landscape.feature_count != feature_count for landscape in landscapes):
        raise ValueError("development landscapes must share a feature count")

    rng = random.Random(seed)
    totals = [0.0] * feature_count
    for landscape in landscapes:
        for position in range(feature_count):
            state = [0] * feature_count
            state[position] = 1
            totals[position] += landscape.evaluate(tuple(state), rng)

    scores = {
        f"position:{position}": total / len(landscapes)
        for position, total in enumerate(totals)
    }
    selected = max(range(feature_count), key=lambda position: totals[position])
    return selected, len(landscapes) * feature_count, scores


def _proposal_position(
    state: tuple[int, ...],
    rng: random.Random,
    preferred: int | None,
    adaptive: bool,
) -> int:
    guidance_draw = rng.random()
    fallback_draw = rng.random()
    if (
        adaptive
        and preferred is not None
        and state[preferred] == 0
        and guidance_draw >= 0.25
    ):
        return preferred

    alternatives = [
        position for position in range(len(state))
        if not (adaptive and preferred is not None and position == preferred)
    ]
    if not alternatives:
        return preferred if preferred is not None else 0
    index = min(int(fallback_draw * len(alternatives)), len(alternatives) - 1)
    return alternatives[index]


def _run_case(
    landscape: SyntheticLandscape,
    *,
    mode: str,
    seed: int,
    total_budget: int,
    discovery_allocation: int,
    preferred_position: int,
    branch_budget: int,
    audit_replications: int,
    audit_interval: int,
) -> dict:
    audit_cost = 2 * audit_replications * branch_budget
    available_after_discovery = total_budget - discovery_allocation
    planned_audits = (
        available_after_discovery // (audit_interval + audit_cost)
        if mode == "online_audit"
        else 0
    )
    online_audit_evaluations = planned_audits * audit_cost
    search_budget = available_after_discovery - online_audit_evaluations
    if search_budget < 1:
        raise ValueError("total budget leaves no evaluations for main search")

    search_rng = random.Random(seed)
    objective_rng = random.Random(seed + 100_000)
    state = (0,) * landscape.feature_count
    current_score: float | None = None
    best_observed = float("-inf")
    best_expected = float("-inf")
    audit_count = 0
    adaptive_enabled = mode == "fixed_adaptive"
    audits: list[dict] = []

    for search_evaluation in range(1, search_budget + 1):
        if search_evaluation == 1:
            candidate = state
            mutation_position = None
        else:
            mutation_position = _proposal_position(
                state,
                search_rng,
                preferred_position,
                adaptive_enabled and mode != "baseline",
            )
            proposed = list(state)
            proposed[mutation_position] = 1 - proposed[mutation_position]
            candidate = tuple(proposed)

        outcome = float(landscape.evaluate(candidate, objective_rng))
        if not math.isfinite(outcome):
            raise ValueError("landscape returned a non-finite objective")
        best_observed = max(best_observed, outcome)
        best_expected = max(
            best_expected,
            landscape.noiseless_value(candidate),
        )

        if current_score is None or outcome > current_score:
            state = candidate
            current_score = outcome

        if (
            mode == "online_audit"
            and audit_count < planned_audits
            and search_evaluation % audit_interval == 0
        ):
            audit = run_online_audit(
                start_state=state,
                objective=landscape.evaluate,
                preferred_position=preferred_position,
                branch_budget=branch_budget,
                replications=audit_replications,
                seed=seed + 1_000_003 + audit_count * 10_007,
            )
            audit_count += 1
            adaptive_enabled = audit.selected_policy == "adaptive"
            audits.append(audit.to_dict())
            for trial in audit.trials:
                for branch in (trial.baseline, trial.adaptive):
                    best_expected = max(
                        best_expected,
                        *(landscape.noiseless_value(visited)
                          for visited in branch.evaluated_states),
                    )
                    best_observed = max(
                        best_observed,
                        branch.best_score,
                    )

    if audit_count != planned_audits:
        raise AssertionError(
            f"scheduled {planned_audits} audits but ran {audit_count}"
        )

    optimum = landscape.optimum()
    total_evaluations = (
        discovery_allocation
        + search_budget
        + audit_count * audit_cost
    )
    return {
        "mode": mode,
        "seed": seed,
        "search_evaluations": search_budget,
        "online_audit_evaluations": audit_count * audit_cost,
        "allocated_discovery_evaluations": discovery_allocation,
        "audit_count": audit_count,
        "total_evaluations": total_evaluations,
        "best_objective": best_observed,
        "best_expected_objective": best_expected,
        "known_optimum": optimum,
        "regret": optimum - best_expected,
        "audits": audits,
    }


def run_benchmark(
    *,
    development_instances: int = 20,
    final_instances: int = 20,
    seeds: list[int] | None = None,
    total_budget: int = 400,
    branch_budget: int = 8,
    audit_replications: int = 4,
    audit_interval: int = 20,
) -> dict:
    seeds = list(seeds) if seeds is not None else [11, 23, 37, 41, 53]
    if development_instances < 1:
        raise ValueError("development_instances must be >= 1")
    if final_instances < 1:
        raise ValueError("final_instances must be >= 1")
    if not seeds:
        raise ValueError("at least one seed is required")
    if total_budget < 1 or branch_budget < 2:
        raise ValueError("total_budget must be positive and branch_budget >= 2")
    if audit_replications < 1 or audit_interval < 1:
        raise ValueError("audit_replications and audit_interval must be >= 1")

    families = ("smooth_additive", "deceptive", "non_additive", "noisy_misleading")
    runs: list[dict] = []
    discoveries: dict[str, dict] = {}
    summaries: list[dict] = []

    for family in families:
        development = [
            make_landscape(family, index)
            for index in range(development_instances)
        ]
        preferred, discovery_cost, signal_scores = discover_signal(
            development,
            seed=17 + FAMILY_OFFSETS[family],
        )
        test_cases = [
            (instance_id, seed, make_landscape(family, instance_id))
            for instance_id in range(
                development_instances,
                development_instances + final_instances,
            )
            for seed in seeds
        ]
        quotient, remainder = divmod(discovery_cost, len(test_cases))
        discoveries[family] = {
            "selected_signal": f"position:{preferred}",
            "selected_position": preferred,
            "evaluations": discovery_cost,
            "scores": signal_scores,
        }

        for mode in MODES:
            mode_runs = []
            for run_index, (instance_id, seed, landscape) in enumerate(test_cases):
                allocation = quotient + (1 if run_index < remainder else 0)
                result = _run_case(
                    landscape,
                    mode=mode,
                    seed=seed,
                    total_budget=total_budget,
                    discovery_allocation=allocation,
                    preferred_position=preferred,
                    branch_budget=branch_budget,
                    audit_replications=audit_replications,
                    audit_interval=audit_interval,
                )
                run = {
                    "family": family,
                    "instance_id": instance_id,
                    **result,
                }
                runs.append(run)
                mode_runs.append(run)

            if sum(run["allocated_discovery_evaluations"] for run in mode_runs) != discovery_cost:
                raise AssertionError(f"discovery allocation mismatch for {family}/{mode}")

    for family in families:
        for mode in MODES:
            selected = [
                run for run in runs
                if run["family"] == family and run["mode"] == mode
            ]
            summaries.append({
                "family": family,
                "mode": mode,
                "n_runs": len(selected),
                "best_expected_objective_mean": statistics.mean(
                    run["best_expected_objective"] for run in selected
                ),
                "regret_mean": statistics.mean(run["regret"] for run in selected),
                "total_evaluations_equal": all(
                    run["total_evaluations"] == total_budget for run in selected
                ),
            })

    return {
        "benchmark_version": "p2-online-policy-audit-v1",
        "methods": list(MODES),
        "landscape_families": list(families),
        "development_instances_per_family": development_instances,
        "final_instances_per_family": final_instances,
        "seeds": seeds,
        "total_evaluation_budget": total_budget,
        "branch_budget": branch_budget,
        "audit_replications": audit_replications,
        "audit_interval": audit_interval,
        "discoveries": discoveries,
        "runs": runs,
        "summary": summaries,
        "scientific_status": "synthetic benchmark only; no biological performance claim",
    }


def write_outputs(payload: dict, output_dir: Path) -> None:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"refusing to overwrite non-empty directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    (output_dir / "raw_results.json").write_text(json.dumps(payload, indent=2) + "\n")
    config = {key: value for key, value in payload.items() if key not in {"runs", "summary", "discoveries"}}
    (output_dir / "config.json").write_text(json.dumps(config, indent=2) + "\n")
    with (output_dir / "summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=payload["summary"][0].keys())
        writer.writeheader()
        writer.writerows(payload["summary"])

    lines = [
        "# P2 online policy audit benchmark",
        "",
        "Synthetic, non-biological optimization only. Larger objective values are better.",
        "Discovery and online-audit costs are included in every held-out run's fixed budget.",
        "Audit branch states are included in the best-visited noiseless objective.",
        "",
        "| Family | Method | Runs | Best noiseless mean | Regret mean | Equal total budget |",
        "|---|---|---:|---:|---:|:---:|",
    ]
    for row in payload["summary"]:
        lines.append(
            f"| {row['family']} | {row['mode']} | {row['n_runs']} "
            f"| {row['best_expected_objective_mean']:.3f} "
            f"| {row['regret_mean']:.3f} "
            f"| {row['total_evaluations_equal']} |"
        )
    lines.extend([
        "",
        "Online policy audits compare baseline and adaptive branches from a saved state.",
        "Greedy acceptance uses one noisy observation; repeated evaluation of candidate states remains a prespecified future comparison.",
        "These results do not establish peptide binding, biological activity, or algorithmic novelty.",
    ])
    (output_dir / "report.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--development-instances", type=int, default=20)
    parser.add_argument("--final-instances", type=int, default=20)
    parser.add_argument("--total-budget", type=int, default=400)
    parser.add_argument("--branch-budget", type=int, default=8)
    parser.add_argument("--audit-replications", type=int, default=4)
    parser.add_argument("--audit-interval", type=int, default=20)
    parser.add_argument("--seeds", nargs="+", type=int, default=[11, 23, 37, 41, 53])
    args = parser.parse_args()

    payload = run_benchmark(
        development_instances=args.development_instances,
        final_instances=args.final_instances,
        seeds=args.seeds,
        total_budget=args.total_budget,
        branch_budget=args.branch_budget,
        audit_replications=args.audit_replications,
        audit_interval=args.audit_interval,
    )
    write_outputs(payload, args.output_dir)
    print(f"Wrote {args.output_dir / 'raw_results.json'}")
    print(f"Wrote {args.output_dir / 'summary.csv'}")
    print(f"Wrote {args.output_dir / 'report.md'}")


if __name__ == "__main__":
    main()
