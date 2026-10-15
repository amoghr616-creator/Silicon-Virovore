"""Audit-budget efficiency study for the isolated P2 adaptive mechanism.

This study does not alter the completed P2 benchmark. It uses new development
and confirmatory instance IDs, keeps each method at the same total objective
budget, and treats the 0/5/10/20/80 percent allocations as predeclared audit
plans. Thresholds are frozen before the confirmatory split is evaluated.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import json
import math
from pathlib import Path
import random
import statistics
import sys
import time

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.benchmark_combined_adaptive_search import (
    MODES,
    _run_search,
    _signal_discovery_records,
    make_landscape,
)
from src.p2_adaptive_search import (
    CombinedAdaptiveController,
    IndependentSignalValidator,
    MatchedPolicyTrialRunner,
    PolicyComparisonResult,
    SignalValidationResult,
    asdict_records,
)


FAMILIES = ("smooth_additive", "deceptive", "non_additive", "noisy_misleading")
PLAN_DEFINITIONS = {
    0: {"discovery_groups": 0, "validation_groups": 0, "matched_trials": 0, "branch_budget": 0},
    5: {"discovery_groups": 1, "validation_groups": 1, "matched_trials": 0, "branch_budget": 2},
    10: {"discovery_groups": 2, "validation_groups": 2, "matched_trials": 2, "branch_budget": 2},
    20: {"discovery_groups": 4, "validation_groups": 4, "matched_trials": 4, "branch_budget": 2},
    80: {"discovery_groups": 10, "validation_groups": 10, "matched_trials": 10, "branch_budget": 8},
}

# Only smooth_additive has a predeclared, stable one-bit signal. The other
# families are deliberately signal-negative for this one-bit discovery design.
TRUSTED_SINGLE_SIGNAL = {"smooth_additive": "position:0"}


def _empty_signal() -> SignalValidationResult:
    return SignalValidationResult(
        selected_signal=None, discovery_groups=(), validation_groups=(), candidate_signals=(),
        discovery_effects={}, adjusted_p_values={}, selected_discovery_effect=None,
        validation_effect=None, validation_ci95=None, validation_win_rate=None,
        passes=False, reason="zero_audit_budget", permutations=0, minimum_groups=4,
        alpha=0.05, minimum_effect=0.05, signal_count=0,
    )


def _empty_policy() -> PolicyComparisonResult:
    return PolicyComparisonResult(
        trials=(), independent_trials=0, mean_difference=None, median_difference=None,
        ci95=None, win_rate=None, passes=False, reason="zero_audit_budget",
        minimum_trials=4, minimum_improvement=0.05, required_win_rate=0.60,
    )


def _sign_permutation_p(values: list[float], *, seed: int, permutations: int = 4095) -> float | None:
    if not values:
        return None
    observed = statistics.mean(values)
    rng = random.Random(seed)
    exceedances = 0
    for _ in range(permutations):
        null = statistics.mean(value if rng.random() < 0.5 else -value for value in values)
        exceedances += null >= observed
    return (1.0 + exceedances) / (permutations + 1.0)


def _sequential_audit(
    *,
    discovery_instances,
    validation_instances,
    matched_trials,
    family: str,
    seed: int,
    branch_budget: int,
    thresholds: dict,
) -> dict:
    """Analyze fixed looks with Bonferroni alpha allocation.

    Looks are fixed in advance at 4, 6, 8, and 10 independent groups. The
    signal max-statistic test uses alpha/K at every look. The policy comparison
    uses a one-sided paired sign permutation test with the same alpha/K; this
    avoids repeatedly inspecting ordinary confidence intervals.
    """
    look_counts = (4, 6, 8, 10)
    look_alpha = thresholds["signal_alpha"] / len(look_counts)
    discovery, _ = _signal_discovery_records(discovery_instances, split="discovery", seed=seed)
    validation, _ = _signal_discovery_records(validation_instances, split="validation", seed=seed + 1000)
    observations = []
    stop = None
    for look in look_counts:
        if len(discovery_instances) < look or len(validation_instances) < look or len(matched_trials) < look:
            observations.append({"look_groups": look, "available": False, "passed": False, "reason": "insufficient_planned_groups"})
            continue
        signal = IndependentSignalValidator(
            alpha=look_alpha,
            minimum_groups=thresholds["signal_minimum_groups"],
            minimum_effect=thresholds["signal_minimum_effect"],
            permutations=thresholds["signal_permutations"],
            seed=seed,
        ).fit_and_validate(
            [record for record in discovery if record.independent_group_id in {f"{family}-instance-{index}" for index in range(100, 100 + look)}],
            [record for record in validation if record.independent_group_id in {f"{family}-instance-{index}" for index in range(110, 110 + look)}],
        )
        differences = [trial.paired_difference for trial in matched_trials[:look] if math.isfinite(trial.paired_difference)]
        policy_p = _sign_permutation_p(differences, seed=seed + look)
        policy_pass = (
            len(differences) >= thresholds["policy_minimum_trials"]
            and policy_p is not None and policy_p <= look_alpha
            and statistics.mean(differences) >= thresholds["policy_minimum_improvement"]
            and statistics.mean(value > 0 for value in differences) >= thresholds["policy_required_win_rate"]
        )
        passed = signal.passes and policy_pass
        observations.append({
            "look_groups": look,
            "available": True,
            "signal_passed": signal.passes,
            "signal_reason": signal.reason,
            "policy_p_value": policy_p,
            "policy_passed": policy_pass,
            "passed": passed,
            "alpha_per_look": look_alpha,
            "audit_cost_if_stopped": look * (16 + 2 * branch_budget),
        })
        if passed and stop is None:
            stop = observations[-1]
    return {
        "look_schedule": look_counts,
        "alpha_spending": "Bonferroni alpha/K at each of K=4 fixed looks",
        "observations": observations,
        "stopped_at_groups": stop["look_groups"] if stop else None,
        "sequential_gate_passed": stop is not None,
        "sequential_audit_cost": stop["audit_cost_if_stopped"] if stop else None,
    }


def _prepare_plan(family: str, allocation: int, *, seed: int, thresholds: dict) -> tuple[CombinedAdaptiveController, dict, int]:
    plan = PLAN_DEFINITIONS[allocation]
    all_development = [make_landscape(family, 100 + index) for index in range(20)]
    discovery_instances = all_development[:plan["discovery_groups"]]
    validation_instances = all_development[10:10 + plan["validation_groups"]]
    if allocation == 0:
        signal_result = _empty_signal()
        policy_result = _empty_policy()
        return CombinedAdaptiveController(signal_result=signal_result, policy_result=policy_result), {
            "allocation": allocation, "audit_cost": 0, "signal_result": asdict(signal_result),
            "policy_result": asdict(policy_result), "sequential": {"sequential_gate_passed": False, "reason": "zero_audit_budget"},
        }, 0

    discovery, discovery_cost = _signal_discovery_records(discovery_instances, split="discovery", seed=seed)
    validation, validation_cost = _signal_discovery_records(validation_instances, split="validation", seed=seed + 1000)
    validator = IndependentSignalValidator(
        alpha=thresholds["signal_alpha"], minimum_groups=thresholds["signal_minimum_groups"],
        minimum_effect=thresholds["signal_minimum_effect"], permutations=thresholds["signal_permutations"], seed=seed,
    )
    signal_result = validator.fit_and_validate(discovery, validation)
    runner = MatchedPolicyTrialRunner(branch_budget=plan["branch_budget"])
    preferred = None
    if signal_result.selected_signal and signal_result.selected_signal.startswith("position:"):
        preferred = int(signal_result.selected_signal.split(":", 1)[1])
    matched_trials = [
        runner.run_trial(
            trial_id=f"{family}-allocation-{allocation}-matched-{index}",
            independent_group_id=f"{family}-instance-{100 + 10 + index}",
            seed=seed + index * 31,
            start_state=(0,) * 8,
            objective=make_landscape(family, 110 + index).evaluate,
            preferred_position=preferred,
        )
        for index in range(plan["matched_trials"])
    ]
    policy_result = runner.compare(
        matched_trials,
        minimum_trials=thresholds["policy_minimum_trials"],
        minimum_improvement=thresholds["policy_minimum_improvement"],
        required_win_rate=thresholds["policy_required_win_rate"],
    )
    controller = CombinedAdaptiveController(signal_result=signal_result, policy_result=policy_result)
    audit_cost = discovery_cost + validation_cost + 2 * plan["branch_budget"] * len(matched_trials)
    sequential = _sequential_audit(
        discovery_instances=discovery_instances,
        validation_instances=validation_instances,
        matched_trials=matched_trials,
        family=family,
        seed=seed,
        branch_budget=plan["branch_budget"],
        thresholds=thresholds,
    ) if len(discovery_instances) >= 4 and len(validation_instances) >= 4 else {"sequential_gate_passed": False, "reason": "insufficient_planned_groups"}
    audit = {
        "allocation": allocation,
        "plan": plan,
        "audit_cost": audit_cost,
        "discovery_observations": asdict_records(discovery),
        "validation_observations": asdict_records(validation),
        "signal_result": asdict(signal_result),
        "matched_trials": asdict_records(matched_trials),
        "policy_result": asdict(policy_result),
        "trusted_single_signal": TRUSTED_SINGLE_SIGNAL.get(family),
        "signal_false_acceptance": bool(signal_result.passes and signal_result.selected_signal != TRUSTED_SINGLE_SIGNAL.get(family)),
        "signal_false_rejection": bool(TRUSTED_SINGLE_SIGNAL.get(family) and not signal_result.passes),
        "sequential": sequential,
    }
    return controller, audit, audit_cost


def _paired_ci(values: list[float], seed: int = 700) -> tuple[float, float] | None:
    if len(values) < 2:
        return None
    rng = random.Random(seed)
    means = [statistics.mean(rng.choices(values, k=len(values))) for _ in range(2000)]
    q = statistics.quantiles(means, n=100)
    return q[1], q[-2]


def _summary(runs: list[dict], paired: list[dict]) -> list[dict]:
    rows = []
    keys = sorted({(run["evaluation_split"], run["family"], run["allocation_percent"], run["mode"]) for run in runs})
    for evaluation_split, family, allocation, mode in keys:
        selected = [run for run in runs if (run["evaluation_split"], run["family"], run["allocation_percent"], run["mode"]) == (evaluation_split, family, allocation, mode)]
        diffs = [row["best_objective_difference"] for row in paired if (row["evaluation_split"], row["family"], row["allocation_percent"], row["mode"]) == (evaluation_split, family, allocation, mode)]
        rows.append({
            "evaluation_split": evaluation_split,
            "family": family,
            "allocation_percent": allocation,
            "mode": mode,
            "independent_instances": len({run["instance_id"] for run in selected}),
            "instance_seed_runs": len(selected),
            "best_mean": statistics.mean(run["best_objective"] for run in selected),
            "best_median": statistics.median(run["best_objective"] for run in selected),
            "best_worst_decile": sorted(run["best_objective"] for run in selected)[max(0, math.ceil(len(selected) * 0.1) - 1)],
            "regret_mean": statistics.mean(run["regret"] for run in selected),
            "paired_difference_mean": statistics.mean(diffs) if diffs else None,
            "paired_difference_median": statistics.median(diffs) if diffs else None,
            "paired_difference_ci95": _paired_ci(diffs),
            "paired_win_rate": sum(value > 0 for value in diffs) / len(diffs) if diffs else None,
            "audit_cost": selected[0]["audit_evaluations"],
            "total_evaluations": selected[0]["total_evaluations"],
            "runtime_mean": statistics.mean(run["runtime_seconds"] for run in selected),
            "failure_rate": statistics.mean(run["failure_rate"] for run in selected),
            "signal_acceptance_rate": statistics.mean(run["signal_validation_passed"] for run in selected),
            "signal_false_acceptance_rate": statistics.mean(run["signal_false_acceptance"] for run in selected),
            "signal_false_rejection_rate": statistics.mean(run["signal_false_rejection"] for run in selected),
            "false_policy_selection_rate": statistics.mean(run["false_policy_selection"] for run in selected),
            "missed_positive_policy_rate": statistics.mean(run["missed_positive_policy"] for run in selected),
            "fallback_frequency": statistics.mean(run["fallback_count"] for run in selected),
            "total_evaluations_equal": len({run["total_evaluations"] for run in selected}) == 1,
        })
    return rows


def run_efficiency_study(*, total_budget: int = 400, development_instances: int = 20, final_instances: int = 20, seeds: list[int] | None = None) -> dict:
    seeds = seeds or [11, 23, 37, 41, 53]
    thresholds = {
        "signal_alpha": 0.05, "signal_minimum_groups": 4, "signal_minimum_effect": 0.05,
        "signal_permutations": 999, "policy_minimum_trials": 4,
        "policy_minimum_improvement": 0.05, "policy_required_win_rate": 0.60,
    }
    runs = []
    paired = []
    audits = {}
    for family in FAMILIES:
        for allocation in PLAN_DEFINITIONS:
            controller, audit, audit_cost = _prepare_plan(family, allocation, seed=17, thresholds=thresholds)
            search_budget = total_budget - audit_cost
            if search_budget <= 0:
                raise ValueError(f"allocation {allocation} exceeds total budget")
            audits[f"{family}:{allocation}"] = audit
            for evaluation_split, start, count in (("development", 200, development_instances), ("confirmatory", 1000, final_instances)):
                for instance_id in range(start, start + count):
                    landscape = make_landscape(family, instance_id)
                    for seed in seeds:
                        baseline_result = None
                        mode_results = {}
                        for mode in MODES:
                            run_controller = CombinedAdaptiveController(
                                signal_result=controller.signal_result,
                                policy_result=controller.policy_result,
                            )
                            result = _run_search(
                                landscape, mode=mode, seed=seed, budget=search_budget,
                                controller=run_controller, checkpoint_interval=20,
                            )
                            mode_results[mode] = result
                            if mode == "baseline":
                                baseline_result = result
                            runs.append({
                                "evaluation_split": evaluation_split, "family": family,
                                "allocation_percent": allocation, "instance_id": instance_id,
                                "seed": seed, "mode": mode,
                                "audit_evaluations": audit_cost,
                                "search_evaluations": search_budget,
                                "total_evaluations": audit_cost + search_budget,
                                "known_optimum": landscape.optimum(),
                                "signal_validation_passed": controller.signal_result.passes,
                                "policy_comparison_passed": controller.policy_result.passes,
                                "signal_false_acceptance": audit.get("signal_false_acceptance", False),
                                "signal_false_rejection": audit.get("signal_false_rejection", False),
                                "audit_plan_valid": controller.signal_result.passes or allocation == 0 or allocation >= 20,
                                **{key: value for key, value in result.items() if key not in {"events", "decisions"}},
                                "event_log": result["events"], "policy_decisions": result["decisions"],
                            })
                        raw = mode_results["adaptive_unvalidated"]
                        for mode in MODES[1:]:
                            candidate = mode_results[mode]
                            paired.append({
                                "evaluation_split": evaluation_split, "family": family,
                                "allocation_percent": allocation, "instance_id": instance_id,
                                "seed": seed, "mode": mode,
                                "best_objective_difference": candidate["best_objective"] - baseline_result["best_objective"],
                                "regret_difference": candidate["regret"] - baseline_result["regret"],
                                "adaptive_won": candidate["best_objective"] > baseline_result["best_objective"],
                                "false_policy_selection": candidate["adaptive_selected_count"] > 0 and candidate["best_objective"] <= baseline_result["best_objective"],
                                "missed_positive_policy": candidate["adaptive_selected_count"] == 0 and raw["best_objective"] > baseline_result["best_objective"],
                            })
    paired_by_key = {
        (
            row["evaluation_split"], row["family"], row["allocation_percent"],
            row["instance_id"], row["seed"], row["mode"],
        ): row
        for row in paired
    }
    for run in runs:
        if run["mode"] == "baseline":
            run["false_policy_selection"] = False
            run["missed_positive_policy"] = False
        else:
            pair = paired_by_key[(run["evaluation_split"], run["family"], run["allocation_percent"], run["instance_id"], run["seed"], run["mode"])]
            run["false_policy_selection"] = pair["false_policy_selection"]
            run["missed_positive_policy"] = pair["missed_positive_policy"]
    return {
        "benchmark_version": "p2-audit-efficiency-v1",
        "methods": list(MODES), "allocations_percent": list(PLAN_DEFINITIONS),
        "landscape_families": list(FAMILIES), "total_evaluation_budget": total_budget,
        "development_instances_per_family": development_instances,
        "confirmatory_instances_per_family": final_instances, "seeds": seeds,
        "thresholds_frozen_before_confirmatory_test": thresholds,
        "audit_plans": PLAN_DEFINITIONS, "audits": audits,
        "runs": runs, "paired_results": paired, "summary": _summary(runs, paired),
        "prior_final_results_reused": False,
        "scientific_status": "synthetic audit-efficiency study; no biological performance or novelty claim",
    }


def write_outputs(payload: dict, output_dir: Path) -> None:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"refusing to overwrite non-empty directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "config.json").write_text(json.dumps({key: payload[key] for key in payload if key not in {"runs", "paired_results", "summary", "audits"}}, indent=2) + "\n")
    (output_dir / "raw_results.json").write_text(json.dumps(payload, indent=2) + "\n")
    with (output_dir / "summary.csv").open("w", newline="") as handle:
        rows = payload["summary"]
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader(); writer.writerows(rows)
    lines = [
        "# P2 audit-budget efficiency study", "",
        "Synthetic, non-biological study. Total objective evaluations are fixed within every allocation.",
        "The confirmatory split uses new instance IDs and does not update thresholds or policy decisions.", "",
        "See `summary.csv` for allocation, decision-quality, and optimizer-performance metrics.", "",
        "Sequential looks use four predeclared group counts and Bonferroni alpha/K spending; ordinary CIs are not used for opportunistic stopping.",
    ]
    (output_dir / "report.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--total-budget", type=int, default=400)
    parser.add_argument("--development-instances", type=int, default=20)
    parser.add_argument("--confirmatory-instances", type=int, default=20)
    parser.add_argument("--seeds", nargs="+", type=int, default=[11, 23, 37, 41, 53])
    args = parser.parse_args()
    payload = run_efficiency_study(
        total_budget=args.total_budget, development_instances=args.development_instances,
        final_instances=args.confirmatory_instances, seeds=args.seeds,
    )
    write_outputs(payload, args.output_dir)
    print(f"Wrote {args.output_dir / 'raw_results.json'}")
    print(f"Wrote {args.output_dir / 'summary.csv'}")
    print(f"Wrote {args.output_dir / 'report.md'}")


if __name__ == "__main__":
    main()
