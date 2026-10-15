"""Reproducible synthetic benchmark for P2's combined adaptive-search gate.

The benchmark is deliberately independent of Silicon Virovore's biological
pipeline. It spends a fixed audit budget on discovery, grouped validation, and
matched policy trials for every method, then gives every method the same
remaining search budget.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
import itertools
import json
import math
from pathlib import Path
import random
import statistics
import sys
import time

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from src.p2_adaptive_search import (
    CombinedAdaptiveController,
    IndependentSignalValidator,
    MatchedPolicyTrialRunner,
    SignalObservation,
    asdict_records,
)


MODES = (
    "baseline",
    "adaptive_unvalidated",
    "signal_validation_only",
    "matched_policy_only",
    "combined",
)
SIGNALS = tuple(f"position:{index}" for index in range(8))
FAMILY_OFFSETS = {
    "smooth_additive": 101,
    "deceptive": 211,
    "non_additive": 307,
    "noisy_misleading": 401,
}


@dataclass(frozen=True)
class SyntheticLandscape:
    family: str
    instance_id: int
    feature_count: int
    weights: tuple[float, ...]
    pairs: tuple[tuple[int, int, float], ...] = ()
    noise_std: float = 0.0
    deceptive: bool = False

    def noiseless_value(self, state: tuple[int, ...]) -> float:
        if self.deceptive:
            value = 0.0
            for start in range(0, self.feature_count, 4):
                count = sum(state[start:start + 4])
                value += 5.0 if count == 4 else 4.0 - count
            return value
        value = sum(weight for bit, weight in zip(state, self.weights) if bit)
        value += sum(effect for left, right, effect in self.pairs if state[left] and state[right])
        return value

    def evaluate(self, state: tuple[int, ...], rng: random.Random) -> float:
        value = self.noiseless_value(state)
        return value + (rng.gauss(0.0, self.noise_std) if self.noise_std else 0.0)

    def optimum(self) -> float:
        return max(self.noiseless_value(state) for state in itertools.product((0, 1), repeat=self.feature_count))

    def true_best_signals(self) -> set[str]:
        zero = (0,) * self.feature_count
        values = {}
        for position in range(self.feature_count):
            state = list(zero)
            state[position] = 1
            values[f"position:{position}"] = self.noiseless_value(tuple(state))
        maximum = max(values.values())
        return {signal for signal, value in values.items() if value == maximum}


def make_landscape(family: str, instance_id: int) -> SyntheticLandscape:
    rng = random.Random(20_000 + instance_id * 31 + FAMILY_OFFSETS[family])
    weights = tuple(round(rng.uniform(0.4, 1.4), 4) for _ in range(8))
    if family == "smooth_additive":
        # A stable, independently replicated signal is intentionally present
        # here: position 0 is consistently best, with small instance-level
        # perturbations. The noisy family tests the opposite case.
        base = (1.40, 1.20, 1.00, 0.80, 0.60, 0.45, 0.30, 0.20)
        perturbed = tuple(round(value + rng.uniform(-0.02, 0.02), 4) for value in base)
        return SyntheticLandscape(family, instance_id, 8, perturbed)
    if family == "deceptive":
        return SyntheticLandscape(family, instance_id, 8, weights, deceptive=True)
    if family == "non_additive":
        pairs = ((0, 1, 2.5), (2, 3, 2.0), (4, 5, 1.8), (6, 7, -0.7))
        return SyntheticLandscape(family, instance_id, 8, weights, pairs=pairs)
    if family == "noisy_misleading":
        pairs = ((0, 1, -1.4), (2, 3, 1.8), (5, 6, 1.2))
        return SyntheticLandscape(family, instance_id, 8, weights, pairs=pairs, noise_std=0.8)
    raise ValueError(f"unknown landscape family: {family}")


def _flip(state: tuple[int, ...], position: int) -> tuple[int, ...]:
    result = list(state)
    result[position] = 1 - result[position]
    return tuple(result)


def _set_position(state: tuple[int, ...], position: int) -> tuple[int, ...]:
    result = list(state)
    result[position] = 1
    return tuple(result)


def _signal_discovery_records(
    landscapes: list[SyntheticLandscape],
    *,
    split: str,
    seed: int,
) -> tuple[list[SignalObservation], int]:
    records: list[SignalObservation] = []
    evaluations = 0
    for landscape in landscapes:
        rng = random.Random(seed + landscape.instance_id * 101)
        zero = (0,) * landscape.feature_count
        for position in range(landscape.feature_count):
            outcome = landscape.evaluate(_flip(zero, position), rng)
            evaluations += 1
            records.append(SignalObservation(
                observation_id=f"{split}-{landscape.family}-{landscape.instance_id}-{position}",
                independent_group_id=f"{landscape.family}-instance-{landscape.instance_id}",
                run_id=f"{split}-run-{landscape.instance_id}",
                lineage_id=f"{landscape.family}-lineage-{landscape.instance_id}",
                context_id=landscape.family,
                signal=f"position:{position}",
                outcome=outcome,
                split=split,
                observation_order=evaluations,
            ))
    return records, evaluations


def _preferred_position(signal: str | None) -> int | None:
    if not signal or not signal.startswith("position:"):
        return None
    return int(signal.split(":", 1)[1])


def _bootstrap_ci(values: list[float], groups: list[str], seed: int = 71) -> tuple[float, float] | None:
    grouped: dict[str, list[float]] = {}
    for group, value in zip(groups, values):
        grouped.setdefault(group, []).append(value)
    means = [statistics.mean(group_values) for group_values in grouped.values()]
    if len(means) < 2:
        return None
    rng = random.Random(seed)
    boot = [statistics.mean(rng.choices(means, k=len(means))) for _ in range(2000)]
    quantiles = statistics.quantiles(boot, n=100)
    return quantiles[1], quantiles[-2]


def _run_search(
    landscape: SyntheticLandscape,
    *,
    mode: str,
    seed: int,
    budget: int,
    controller: CombinedAdaptiveController,
    checkpoint_interval: int,
) -> dict:
    started = time.perf_counter()
    rng = random.Random(seed)
    objective_rng = random.Random(seed + 100_000)
    state = (0,) * landscape.feature_count
    best = float("-inf")
    best_expected = float("-inf")
    events: list[dict] = []
    decisions: list[dict] = []
    failures = 0
    preferred = _preferred_position(controller.signal_result.selected_signal)
    for evaluation in range(1, budget + 1):
        if evaluation == 1 or (evaluation - 1) % checkpoint_interval == 0:
            decision = controller.decide(checkpoint=f"evaluation:{evaluation}", mode=mode)
            decisions.append(asdict(decision))
        try:
            outcome = landscape.evaluate(state, objective_rng)
            if not math.isfinite(outcome):
                raise ValueError("non-finite objective")
        except Exception as exc:
            failures += 1
            events.append({"evaluation": evaluation, "state": list(state), "failure": type(exc).__name__})
            continue
        best = max(best, outcome)
        best_expected = max(best_expected, landscape.noiseless_value(state))
        selected = decisions[-1]["adaptive_selected"] if decisions else False
        available = [index for index, value in enumerate(state) if not value]
        if available and selected and preferred in available and rng.random() >= 0.25:
            position = preferred
        elif available:
            position = rng.choice(available)
        else:
            position = rng.randrange(landscape.feature_count)
        events.append({
            "evaluation": evaluation,
            "state": list(state),
            "outcome": outcome,
            "adaptive_selected": selected,
            "mutation_position": position,
        })
        state = _set_position(state, position)
    return {
        "best_objective": best,
        "best_expected_objective": best_expected,
        "regret": landscape.optimum() - best_expected,
        "search_evaluations": budget,
        "failures": failures,
        "failure_rate": failures / budget,
        "runtime_seconds": time.perf_counter() - started,
        "adaptive_selected_count": sum(decision["adaptive_selected"] for decision in decisions),
        "fallback_count": sum(not decision["adaptive_selected"] for decision in decisions),
        "events": events,
        "decisions": decisions,
    }


def _prepare_controller(
    family: str,
    discovery_instances: list[SyntheticLandscape],
    validation_instances: list[SyntheticLandscape],
    *,
    seed: int,
    branch_budget: int,
    thresholds: dict,
) -> tuple[CombinedAdaptiveController, dict, int]:
    discovery, discovery_cost = _signal_discovery_records(discovery_instances, split="discovery", seed=seed)
    validation, validation_cost = _signal_discovery_records(validation_instances, split="validation", seed=seed + 1_000)
    validator = IndependentSignalValidator(
        alpha=thresholds["signal_alpha"],
        minimum_groups=thresholds["signal_minimum_groups"],
        minimum_effect=thresholds["signal_minimum_effect"],
        permutations=thresholds["signal_permutations"],
        seed=seed,
    )
    signal_result = validator.fit_and_validate(discovery, validation)
    runner = MatchedPolicyTrialRunner(branch_budget=branch_budget)
    preferred = _preferred_position(signal_result.selected_signal)
    matched_trials = []
    for landscape in validation_instances:
        matched_trials.append(runner.run_trial(
            trial_id=f"{family}-matched-{landscape.instance_id}",
            independent_group_id=f"{family}-instance-{landscape.instance_id}",
            seed=seed + landscape.instance_id * 31,
            start_state=(0,) * landscape.feature_count,
            objective=landscape.evaluate,
            preferred_position=preferred,
        ))
    policy_result = runner.compare(
        matched_trials,
        minimum_trials=thresholds["policy_minimum_trials"],
        minimum_improvement=thresholds["policy_minimum_improvement"],
        required_win_rate=thresholds["policy_required_win_rate"],
    )
    controller = CombinedAdaptiveController(signal_result=signal_result, policy_result=policy_result)
    audit = {
        "family": family,
        "discovery_observations": asdict_records(discovery),
        "validation_observations": asdict_records(validation),
        "signal_result": asdict(signal_result),
        "matched_trials": asdict_records(matched_trials),
        "policy_result": asdict(policy_result),
        "selected_signal_true_in_discovery": signal_result.selected_signal,
    }
    return controller, audit, discovery_cost + validation_cost + 2 * branch_budget * len(validation_instances)


def run_benchmark(
    *,
    development_instances: int = 20,
    final_instances: int = 20,
    seeds: list[int] | None = None,
    total_budget: int = 400,
    branch_budget: int = 8,
    checkpoint_interval: int = 20,
) -> dict:
    seeds = seeds or [11, 23, 37, 41, 53]
    thresholds = {
        "signal_alpha": 0.05,
        "signal_minimum_groups": 4,
        "signal_minimum_effect": 0.05,
        "signal_permutations": 999,
        "policy_minimum_trials": 4,
        "policy_minimum_improvement": 0.05,
        "policy_required_win_rate": 0.60,
    }
    families = ("smooth_additive", "deceptive", "non_additive", "noisy_misleading")
    runs: list[dict] = []
    audits: dict[str, dict] = {}
    for family in families:
        dev = [make_landscape(family, index) for index in range(development_instances)]
        split = max(1, development_instances // 2)
        discovery_instances = dev[:split]
        validation_instances = dev[split:]
        if len(discovery_instances) < thresholds["signal_minimum_groups"] or len(validation_instances) < thresholds["signal_minimum_groups"]:
            raise ValueError("development_instances must provide at least four discovery and four validation instances")
        controller, audit, audit_cost = _prepare_controller(
            family, discovery_instances, validation_instances, seed=17, branch_budget=branch_budget, thresholds=thresholds
        )
        search_budget = total_budget - audit_cost
        if search_budget <= 0:
            raise ValueError(f"total_budget {total_budget} is too small for audit cost {audit_cost}")
        audits[family] = {"audit_cost": audit_cost, **audit}
        for instance_id in range(development_instances, development_instances + final_instances):
            landscape = make_landscape(family, instance_id)
            true_signals = landscape.true_best_signals()
            for seed in seeds:
                for mode in MODES:
                    result = _run_search(
                        landscape, mode=mode, seed=seed, budget=search_budget,
                        controller=controller, checkpoint_interval=checkpoint_interval,
                    )
                    runs.append({
                        "family": family,
                        "instance_id": instance_id,
                        "seed": seed,
                        "mode": mode,
                        "known_optimum": landscape.optimum(),
                        "audit_evaluations": audit_cost,
                        "total_evaluations": audit_cost + result["search_evaluations"],
                        "signal_selected": controller.signal_result.selected_signal,
                        "signal_validation_passed": controller.signal_result.passes,
                        "signal_selected_is_true_best_single": controller.signal_result.selected_signal in true_signals,
                        "policy_comparison_passed": controller.policy_result.passes,
                        "audit_cost_included": True,
                        "controller_thresholds_frozen_before_final_test": True,
                        **{key: value for key, value in result.items() if key not in {"events", "decisions"}},
                        "event_log": result["events"],
                        "policy_decisions": result["decisions"],
                    })

    paired = []
    keyed = {(run["family"], run["instance_id"], run["seed"], run["mode"]): run for run in runs}
    for family in families:
        for instance_id in range(development_instances, development_instances + final_instances):
            for seed in seeds:
                baseline = keyed[(family, instance_id, seed, "baseline")]
                for mode in MODES[1:]:
                    candidate = keyed[(family, instance_id, seed, mode)]
                    paired.append({
                        "family": family,
                        "instance_id": instance_id,
                        "seed": seed,
                        "mode": mode,
                        "best_objective_difference": candidate["best_objective"] - baseline["best_objective"],
                        "regret_difference": candidate["regret"] - baseline["regret"],
                        "adaptive_won": candidate["best_objective"] > baseline["best_objective"],
                    })
    summary = _summarize(runs, paired, families)
    return {
        "benchmark_version": "p2-combined-adaptive-search-v1",
        "methods": list(MODES),
        "landscape_families": list(families),
        "development_instances_per_family": development_instances,
        "final_instances_per_family": final_instances,
        "seeds": seeds,
        "total_evaluation_budget": total_budget,
        "matched_branch_budget": branch_budget,
        "checkpoint_interval": checkpoint_interval,
        "thresholds_frozen_before_final_test": thresholds,
        "audits": audits,
        "runs": runs,
        "paired_results": paired,
        "summary": summary,
        "scientific_status": "synthetic benchmark only; no biological performance or novelty claim",
    }


def _summarize(runs: list[dict], paired: list[dict], families: tuple[str, ...]) -> list[dict]:
    rows = []
    for family in families:
        family_runs = [run for run in runs if run["family"] == family]
        for mode in MODES:
            selected = [run for run in family_runs if run["mode"] == mode]
            differences = [row["best_objective_difference"] for row in paired if row["family"] == family and row["mode"] == mode]
            groups = [f"{row['instance_id']}" for row in paired if row["family"] == family and row["mode"] == mode]
            best_values = [run["best_objective"] for run in selected]
            paired_win = sum(value > 0 for value in differences) / len(differences) if differences else None
            rows.append({
                "family": family,
                "mode": mode,
                "n_instance_seed_runs": len(selected),
                "n_independent_instances": len({run["instance_id"] for run in selected}),
                "best_mean": statistics.mean(best_values),
                "best_median": statistics.median(best_values),
                "best_worst_decile": sorted(best_values)[max(0, math.ceil(len(best_values) * 0.1) - 1)],
                "regret_mean": statistics.mean(run["regret"] for run in selected),
                "paired_difference_mean": statistics.mean(differences) if differences else None,
                "paired_difference_median": statistics.median(differences) if differences else None,
                "paired_difference_ci95": _bootstrap_ci(differences, groups),
                "paired_win_rate": paired_win,
                "failure_rate_mean": statistics.mean(run["failure_rate"] for run in selected),
                "signal_validation_rate": statistics.mean(run["signal_validation_passed"] for run in selected),
                "false_discovery_rate": statistics.mean(not run["signal_selected_is_true_best_single"] for run in selected),
                "fallback_frequency": statistics.mean(run["fallback_count"] for run in selected),
                "adaptive_selection_frequency": statistics.mean(run["adaptive_selected_count"] for run in selected),
                "adaptive_selected_and_improved": statistics.mean(
                    run["adaptive_selected_count"] > 0 and run["best_objective"] > next(
                        baseline["best_objective"] for baseline in family_runs
                        if baseline["instance_id"] == run["instance_id"] and baseline["seed"] == run["seed"] and baseline["mode"] == "baseline"
                    ) for run in selected
                ),
                "total_evaluations_equal": len({run["total_evaluations"] for run in selected}) == 1,
                "runtime_seconds_mean": statistics.mean(run["runtime_seconds"] for run in selected),
            })
    return rows


def write_outputs(payload: dict, output_dir: Path) -> None:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"refusing to overwrite non-empty directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "config.json").write_text(json.dumps({key: payload[key] for key in payload if key not in {"runs", "paired_results", "summary", "audits"}}, indent=2) + "\n")
    (output_dir / "raw_results.json").write_text(json.dumps(payload, indent=2) + "\n")
    with (output_dir / "summary.csv").open("w", newline="") as handle:
        rows = payload["summary"]
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    lines = [
        "# P2 combined adaptive-search benchmark",
        "",
        "Synthetic, non-biological optimization only. Larger objective values are better.",
        "Audit discovery, validation, and matched-trial costs are included in every method's fixed total evaluation budget.",
        "",
        "| Family | Method | Instances | Best mean | Regret mean | Paired difference | 95% CI | Win rate | Signal validation | Fallback frequency | Equal evaluations |",
        "|---|---|---:|---:|---:|---:|---|---:|---:|---:|:---:|",
    ]
    for row in payload["summary"]:
        ci = row["paired_difference_ci95"]
        ci_text = "n/a" if ci is None else f"[{ci[0]:.3f}, {ci[1]:.3f}]"
        lines.append(
            f"| {row['family']} | {row['mode']} | {row['n_independent_instances']} | {row['best_mean']:.3f} | "
            f"{row['regret_mean']:.3f} | {row['paired_difference_mean'] if row['paired_difference_mean'] is not None else float('nan'):.3f} | "
            f"{ci_text} | {row['paired_win_rate'] if row['paired_win_rate'] is not None else float('nan'):.3f} | "
            f"{row['signal_validation_rate']:.3f} | {row['fallback_frequency']:.2f} | {row['total_evaluations_equal']} |"
        )
    lines += [
        "",
        "The combined gate requires both independent signal validation and a matched-policy pass.",
        "Unvalidated adaptive, signal-only, and matched-policy-only modes are comparators, not state-of-the-art claims.",
        "",
        "Final-test results do not change thresholds or policy decisions.",
    ]
    (output_dir / "report.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--development-instances", type=int, default=20)
    parser.add_argument("--final-instances", type=int, default=20)
    parser.add_argument("--total-budget", type=int, default=400)
    parser.add_argument("--branch-budget", type=int, default=8)
    parser.add_argument("--seeds", nargs="+", type=int, default=[11, 23, 37, 41, 53])
    args = parser.parse_args()
    payload = run_benchmark(
        development_instances=args.development_instances,
        final_instances=args.final_instances,
        seeds=args.seeds,
        total_budget=args.total_budget,
        branch_budget=args.branch_budget,
    )
    write_outputs(payload, args.output_dir)
    print(f"Wrote {args.output_dir / 'raw_results.json'}")
    print(f"Wrote {args.output_dir / 'summary.csv'}")
    print(f"Wrote {args.output_dir / 'report.md'}")


if __name__ == "__main__":
    main()
