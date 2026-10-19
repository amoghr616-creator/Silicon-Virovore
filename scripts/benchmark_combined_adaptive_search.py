"""Reproducible synthetic benchmark for P2 adaptive search.

This benchmark is independent of the biological pipeline.

It compares:
- baseline search;
- adaptive search without evidence gates;
- signal-validation-only search;
- matched-policy-only search;
- combined evidence-gated search.

Audit evaluations are charged to the corresponding method's total budget.
Signal discovery, signal validation, and matched policy trials use disjoint
development groups. Final test instances are not used to set thresholds.

Synthetic results do not establish biological performance.
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
    deceptive_blocks: tuple[tuple[int, ...], ...] = ()
    deceptive_target: tuple[int, ...] = ()
    deceptive_baselines: tuple[float, ...] = ()
    deceptive_rewards: tuple[float, ...] = ()

    def noiseless_value(self, state: tuple[int, ...]) -> float:
        if len(state) != self.feature_count:
            raise ValueError("state length does not match feature_count")
        if any(bit not in (0, 1) for bit in state):
            raise ValueError("state must contain only binary values")

        if self.deceptive:
            if (
                len(self.deceptive_target) != self.feature_count
                or not self.deceptive_blocks
                or len(self.deceptive_baselines) != len(self.deceptive_blocks)
                or len(self.deceptive_rewards) != len(self.deceptive_blocks)
            ):
                raise ValueError("invalid deceptive landscape configuration")

            covered = [
                position
                for block in self.deceptive_blocks
                for position in block
            ]
            if sorted(covered) != list(range(self.feature_count)):
                raise ValueError(
                    "deceptive blocks must partition all positions exactly once"
                )

            value = 0.0
            for block_index, block in enumerate(self.deceptive_blocks):
                matches = sum(
                    state[position] == self.deceptive_target[position]
                    for position in block
                )
                baseline = self.deceptive_baselines[block_index]
                reward = self.deceptive_rewards[block_index]

                if reward <= baseline:
                    raise ValueError(
                        "each deceptive trap reward must exceed its baseline"
                    )

                # The target gets a reward. Partial progress toward the
                # target decreases the score, creating a local trap at
                # the blockwise complement of the target.
                if matches == len(block):
                    value += reward
                else:
                    value += baseline - matches

            return value

        value = sum(
            weight
            for bit, weight in zip(state, self.weights)
            if bit
        )
        value += sum(
            effect
            for left, right, effect in self.pairs
            if state[left] and state[right]
        )
        return value

    def evaluate(
        self,
        state: tuple[int, ...],
        rng: random.Random,
    ) -> float:
        value = self.noiseless_value(state)
        if self.noise_std:
            value += rng.gauss(0.0, self.noise_std)
        return value

    def optimum(self) -> float:
        """Calculate the exact optimum by enumerating all binary states."""
        return max(
            self.noiseless_value(state)
            for state in itertools.product(
                (0, 1), repeat=self.feature_count
            )
        )

    def true_best_signals(self) -> set[str]:
        """Return best one-bit signals relative to the all-zero state.

        This is a narrow reference label, not a claim that a signal is
        globally useful once interactions are present.
        """
        zero = (0,) * self.feature_count
        values: dict[str, float] = {}

        for position in range(self.feature_count):
            state = list(zero)
            state[position] = 1
            values[f"position:{position}"] = self.noiseless_value(
                tuple(state)
            )

        maximum = max(values.values())
        return {
            signal
            for signal, value in values.items()
            if math.isclose(value, maximum, rel_tol=0.0, abs_tol=1e-12)
        }


def make_landscape(
    family: str,
    instance_id: int,
) -> SyntheticLandscape:
    if family not in FAMILY_OFFSETS:
        raise ValueError(f"unknown landscape family: {family}")

    rng = random.Random(
        20_000 + instance_id * 31 + FAMILY_OFFSETS[family]
    )

    weights = tuple(
        round(rng.uniform(0.4, 1.4), 4)
        for _ in range(8)
    )

    if family == "smooth_additive":
        # A prespecified stable signal: position 0 has the highest
        # additive value across instances.
        base = (1.40, 1.20, 1.00, 0.80, 0.60, 0.45, 0.30, 0.20)
        perturbed = tuple(
            round(value + rng.uniform(-0.02, 0.02), 4)
            for value in base
        )
        return SyntheticLandscape(
            family=family,
            instance_id=instance_id,
            feature_count=8,
            weights=perturbed,
        )

    if family == "deceptive":
        positions = list(range(8))
        rng.shuffle(positions)

        blocks = (
            tuple(positions[:4]),
            tuple(positions[4:]),
        )
        target = tuple(rng.randrange(2) for _ in range(8))
        baselines = tuple(
            round(rng.uniform(3.2, 4.2), 4)
            for _ in blocks
        )
        rewards = tuple(
            round(baseline + rng.uniform(1.0, 2.0), 4)
            for baseline in baselines
        )

        return SyntheticLandscape(
            family=family,
            instance_id=instance_id,
            feature_count=8,
            weights=weights,
            deceptive=True,
            deceptive_blocks=blocks,
            deceptive_target=target,
            deceptive_baselines=baselines,
            deceptive_rewards=rewards,
        )

    if family == "non_additive":
        # Interaction coefficients vary by instance. The optimum is
        # calculated exactly rather than assumed to be all ones.
        pairs = (
            (0, 1, round(rng.uniform(1.5, 3.0), 4)),
            (2, 3, round(rng.uniform(-2.5, 2.5), 4)),
            (4, 5, round(rng.uniform(1.0, 2.5), 4)),
            (6, 7, round(rng.uniform(-2.0, 1.0), 4)),
        )
        return SyntheticLandscape(
            family=family,
            instance_id=instance_id,
            feature_count=8,
            weights=weights,
            pairs=pairs,
        )

    if family == "noisy_misleading":
        # Mixed interactions and observation noise test whether a
        # signal chosen from a few noisy observations generalizes.
        pairs = (
            (0, 1, round(rng.uniform(-2.5, 1.5), 4)),
            (2, 3, round(rng.uniform(1.0, 2.5), 4)),
            (5, 6, round(rng.uniform(-1.8, 2.0), 4)),
        )
        return SyntheticLandscape(
            family=family,
            instance_id=instance_id,
            feature_count=8,
            weights=weights,
            pairs=pairs,
            noise_std=round(rng.uniform(0.6, 1.0), 4),
        )

    raise ValueError(f"unknown landscape family: {family}")


def _flip(
    state: tuple[int, ...],
    position: int,
) -> tuple[int, ...]:
    result = list(state)
    result[position] = 1 - result[position]
    return tuple(result)


def _set_position(
    state: tuple[int, ...],
    position: int,
) -> tuple[int, ...]:
    """Compatibility helper for older tests/scripts."""
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
            outcome = landscape.evaluate(
                _flip(zero, position),
                rng,
            )
            evaluations += 1

            records.append(
                SignalObservation(
                    observation_id=(
                        f"{split}-{landscape.family}-"
                        f"{landscape.instance_id}-{position}"
                    ),
                    independent_group_id=(
                        f"{landscape.family}-instance-{landscape.instance_id}"
                    ),
                    run_id=f"{split}-run-{landscape.instance_id}",
                    lineage_id=(
                        f"{landscape.family}-lineage-{landscape.instance_id}"
                    ),
                    context_id=landscape.family,
                    signal=f"position:{position}",
                    outcome=outcome,
                    split=split,
                    observation_order=evaluations,
                )
            )

    return records, evaluations


def _preferred_position(signal: str | None) -> int | None:
    if not signal or not signal.startswith("position:"):
        return None
    try:
        return int(signal.split(":", 1)[1])
    except ValueError:
        return None


def _bootstrap_ci(
    values: list[float],
    groups: list[str],
    seed: int = 71,
) -> tuple[float, float] | None:
    """Bootstrap independent-group means rather than treating seeds as groups."""
    if len(values) != len(groups):
        raise ValueError("values and groups must have equal length")

    grouped: dict[str, list[float]] = {}
    for group, value in zip(groups, values):
        grouped.setdefault(group, []).append(value)

    group_means = [
        statistics.mean(group_values)
        for group_values in grouped.values()
    ]
    if len(group_means) < 2:
        return None

    rng = random.Random(seed)
    bootstrap_means = [
        statistics.mean(
            rng.choices(group_means, k=len(group_means))
        )
        for _ in range(2000)
    ]
    quantiles = statistics.quantiles(bootstrap_means, n=100)
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
    if budget < 1:
        raise ValueError("budget must be >= 1")
    if checkpoint_interval < 1:
        raise ValueError("checkpoint_interval must be >= 1")

    started = time.perf_counter()
    policy_rng = random.Random(seed)
    objective_rng = random.Random(seed + 100_000)

    state = (0,) * landscape.feature_count
    current_score: float | None = None
    best_observed = float("-inf")
    best_noiseless = float("-inf")

    events: list[dict] = []
    decisions: list[dict] = []
    failures = 0

    preferred = _preferred_position(
        controller.signal_result.selected_signal
    )

    for evaluation in range(1, budget + 1):
        if (
            evaluation == 1
            or (evaluation - 1) % checkpoint_interval == 0
        ):
            decision = controller.decide(
                checkpoint=f"evaluation:{evaluation}",
                mode=mode,
            )
            decisions.append(asdict(decision))

        selected = decisions[-1]["adaptive_selected"]

        # Evaluate the initial state before making any moves.
        if current_score is None:
            try:
                current_score = float(
                    landscape.evaluate(state, objective_rng)
                )
                if not math.isfinite(current_score):
                    raise ValueError("non-finite objective")
            except Exception as exc:
                failures += 1
                events.append({
                    "evaluation": evaluation,
                    "evaluation_type": "initial",
                    "state": list(state),
                    "failure": type(exc).__name__,
                })
                continue

            true_value = landscape.noiseless_value(state)
            best_observed = max(best_observed, current_score)
            best_noiseless = max(best_noiseless, true_value)

            events.append({
                "evaluation": evaluation,
                "evaluation_type": "initial",
                "state": list(state),
                "candidate_state": list(state),
                "state_after": list(state),
                "outcome": current_score,
                "noiseless_value": true_value,
                "adaptive_selected": selected,
                "mutation_position": None,
            })
            continue

        # Propose a reversible single-bit mutation.
        try:
            guidance_draw = policy_rng.random()
            fallback_draw = policy_rng.random()

            position = min(
                int(fallback_draw * landscape.feature_count),
                landscape.feature_count - 1,
            )

            if (
                selected
                and preferred is not None
                and 0 <= preferred < landscape.feature_count
            ):
                alternatives = [
                    index
                    for index in range(landscape.feature_count)
                    if index != preferred
                ]

                if (
                    state[preferred] == 0
                    and guidance_draw >= 0.25
                ):
                    # Bias toward activating the selected signal.
                    position = preferred
                elif alternatives:
                    # Do not intentionally deactivate the selected bit.
                    alternative_index = min(
                        int(fallback_draw * len(alternatives)),
                        len(alternatives) - 1,
                    )
                    position = alternatives[alternative_index]

            candidate_state = _flip(state, position)
            proposal_score = float(
                landscape.evaluate(candidate_state, objective_rng)
            )

            if not math.isfinite(proposal_score):
                raise ValueError("non-finite objective")

        except Exception as exc:
            failures += 1
            events.append({
                "evaluation": evaluation,
                "evaluation_type": "proposal",
                "state": list(state),
                "candidate_state": (
                    list(candidate_state)
                    if "candidate_state" in locals()
                    else None
                ),
                "mutation_position": position if "position" in locals() else None,
                "failure": type(exc).__name__,
            })
            continue

        # Record both the proposed state and the decision separately.
        current_state_before = state
        current_score_before = current_score
        accepted = proposal_score > current_score_before

        true_value = landscape.noiseless_value(candidate_state)
        best_observed = max(best_observed, proposal_score)
        best_noiseless = max(best_noiseless, true_value)

        if accepted:
            state = candidate_state
            current_score = proposal_score

        events.append({
            "evaluation": evaluation,
            "evaluation_type": "proposal",
            "current_state_before": list(current_state_before),
            "candidate_state": list(candidate_state),
            "state_after": list(state),
            # Preserve the old field for existing diagnostic scripts.
            "state": list(state),
            "outcome": proposal_score,
            "current_score_before": current_score_before,
            "noiseless_value": true_value,
            "mutation_position": position,
            "accepted": accepted,
            "adaptive_selected": selected,
        })

    optimum = landscape.optimum()

    return {
        "best_objective": best_observed,
        "best_expected_objective": best_noiseless,
        "regret": optimum - best_noiseless,
        "search_evaluations": budget,
        "failures": failures,
        "failure_rate": failures / budget,
        "runtime_seconds": time.perf_counter() - started,
        "adaptive_selected_count": sum(
            decision["adaptive_selected"]
            for decision in decisions
        ),
        "fallback_count": sum(
            not decision["adaptive_selected"]
            for decision in decisions
        ),
        "events": events,
        "decisions": decisions,
    }

def _prepare_controller(
    family: str,
    discovery_instances: list[SyntheticLandscape],
    validation_instances: list[SyntheticLandscape],
    policy_trial_instances: list[SyntheticLandscape],
    *,
    seed: int,
    branch_budget: int,
    thresholds: dict,
) -> tuple[CombinedAdaptiveController, dict, dict[str, int]]:
    """Prepare signal and policy evidence using disjoint groups."""

    if len(discovery_instances) < thresholds["signal_minimum_groups"]:
        raise ValueError("not enough independent discovery groups")
    if len(validation_instances) < thresholds["signal_minimum_groups"]:
        raise ValueError("not enough independent signal-validation groups")
    if len(policy_trial_instances) < thresholds["policy_minimum_trials"]:
        raise ValueError("not enough independent matched-policy groups")

    group_sets = [
        {landscape.instance_id for landscape in discovery_instances},
        {landscape.instance_id for landscape in validation_instances},
        {landscape.instance_id for landscape in policy_trial_instances},
    ]
    if (
        group_sets[0] & group_sets[1]
        or group_sets[0] & group_sets[2]
        or group_sets[1] & group_sets[2]
    ):
        raise ValueError(
            "discovery, validation, and policy trial groups must be disjoint"
        )

    discovery, discovery_cost = _signal_discovery_records(
        discovery_instances,
        split="discovery",
        seed=seed,
    )
    validation, validation_cost = _signal_discovery_records(
        validation_instances,
        split="validation",
        seed=seed + 1_000,
    )

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

    for landscape in policy_trial_instances:
        matched_trials.append(
            runner.run_trial(
                trial_id=f"{family}-matched-{landscape.instance_id}",
                independent_group_id=(
                    f"{family}-instance-{landscape.instance_id}"
                ),
                seed=seed + landscape.instance_id * 31,
                start_state=(0,) * landscape.feature_count,
                objective=landscape.evaluate,
                preferred_position=preferred,
            )
        )

    policy_result = runner.compare(
        matched_trials,
        minimum_trials=thresholds["policy_minimum_trials"],
        minimum_improvement=thresholds["policy_minimum_improvement"],
        required_win_rate=thresholds["policy_required_win_rate"],
    )

    controller = CombinedAdaptiveController(
        signal_result=signal_result,
        policy_result=policy_result,
    )

    # Each discovery/validation evaluation and each branch evaluation counts.
    # The policy-only mode needs discovery to know which signal to test.
    matched_policy_cost = 2 * branch_budget * len(matched_trials)

    audit_costs = {
        "baseline": 0,
        "adaptive_unvalidated": discovery_cost,
        "signal_validation_only": discovery_cost + validation_cost,
        "matched_policy_only": discovery_cost + matched_policy_cost,
        "combined": (
            discovery_cost + validation_cost + matched_policy_cost
        ),
    }

    audit = {
        "family": family,
        "discovery_observations": asdict_records(discovery),
        "validation_observations": asdict_records(validation),
        "signal_result": asdict(signal_result),
        "matched_trials": asdict_records(matched_trials),
        "discovery_evaluations": discovery_cost,
        "signal_validation_evaluations": validation_cost,
        "matched_policy_evaluations": matched_policy_cost,
        "audit_cost_by_mode": audit_costs,
        "development_group_ids": {
            "discovery": [
                landscape.instance_id for landscape in discovery_instances
            ],
            "signal_validation": [
                landscape.instance_id for landscape in validation_instances
            ],
            "policy_trials": [
                landscape.instance_id for landscape in policy_trial_instances
            ],
        },
    }

    return controller, audit, audit_costs


def _summarize(
    runs: list[dict],
    paired: list[dict],
    families: tuple[str, ...],
) -> list[dict]:
    rows = []

    for family in families:
        family_runs = [
            run for run in runs if run["family"] == family
        ]

        for mode in MODES:
            selected = [
                run for run in family_runs if run["mode"] == mode
            ]
            paired_rows = [
                row
                for row in paired
                if row["family"] == family and row["mode"] == mode
            ]

            # Paired differences are grouped by independent test instance.
            per_instance: dict[int, list[dict]] = {}
            for row in paired_rows:
                per_instance.setdefault(row["instance_id"], []).append(row)

            instance_differences = [
                statistics.mean(
                    row["best_expected_objective_difference"]
                    for row in instance_rows
                )
                for instance_rows in per_instance.values()
            ]
            instance_noisy_differences = [
                statistics.mean(
                    row["best_objective_difference"]
                    for row in instance_rows
                )
                for instance_rows in per_instance.values()
            ]
            instance_regret_differences = [
                statistics.mean(
                    row["regret_difference"]
                    for row in instance_rows
                )
                for instance_rows in per_instance.values()
            ]

            groups = [str(instance_id) for instance_id in per_instance]
            ci_noiseless = _bootstrap_ci(
                instance_differences,
                groups,
                seed=71,
            ) if instance_differences else None

            ci_noisy = _bootstrap_ci(
                instance_noisy_differences,
                groups,
                seed=72,
            ) if instance_noisy_differences else None

            best_values = [
                run["best_objective"] for run in selected
            ]
            noiseless_values = [
                run["best_expected_objective"] for run in selected
            ]

            signal_correctness = [
                run["signal_selected_is_true_best_single"]
                for run in selected
            ]

            rows.append({
                "family": family,
                "mode": mode,
                "n_instance_seed_runs": len(selected),
                "n_independent_instances": len({
                    run["instance_id"] for run in selected
                }),
                "best_mean": statistics.mean(best_values),
                "best_median": statistics.median(best_values),
                "best_noiseless_mean": statistics.mean(noiseless_values),
                "best_worst_decile": sorted(best_values)[
                    max(0, math.ceil(len(best_values) * 0.1) - 1)
                ],
                "regret_mean": statistics.mean(
                    run["regret"] for run in selected
                ),
                "paired_noiseless_difference_mean": (
                    statistics.mean(instance_differences)
                    if instance_differences else None
                ),
                "paired_noiseless_difference_ci95": ci_noiseless,
                "paired_difference_mean": (
                    statistics.mean(instance_noisy_differences)
                    if instance_noisy_differences else None
                ),
                "paired_difference_median": (
                    statistics.median(instance_noisy_differences)
                    if instance_noisy_differences else None
                ),
                "paired_difference_ci95": ci_noisy,
                "paired_regret_difference_mean": (
                    statistics.mean(instance_regret_differences)
                    if instance_regret_differences else None
                ),
                "paired_win_rate": (
                    sum(value > 0 for value in instance_differences)
                    / len(instance_differences)
                    if instance_differences else None
                ),
                "failure_rate_mean": statistics.mean(
                    run["failure_rate"] for run in selected
                ),
                "signal_validation_rate": statistics.mean(
                    run["signal_validation_passed"] for run in selected
                ),
                "false_discovery_rate": statistics.mean(
                    not run["signal_selected_is_true_best_single"]
                    for run in selected
                ),
                "fallback_frequency": statistics.mean(
                    run["fallback_count"] for run in selected
                ),
                "adaptive_selection_frequency": statistics.mean(
                    run["adaptive_selected_count"] for run in selected
                ),
                "adaptive_selected_and_improved": statistics.mean(
                    run["adaptive_selected_count"] > 0
                    and run["best_expected_objective"] >
                    next(
                        baseline["best_expected_objective"]
                        for baseline in family_runs
                        if (
                            baseline["instance_id"] == run["instance_id"]
                            and baseline["seed"] == run["seed"]
                            and baseline["mode"] == "baseline"
                        )
                    )
                    for run in selected
                ),
                "total_evaluations_equal": (
                    len({run["total_evaluations"] for run in selected}) == 1
                ),
                "runtime_seconds_mean": statistics.mean(
                    run["runtime_seconds"] for run in selected
                ),
            })

    return rows


def run_benchmark(
    *,
    development_instances: int = 20,
    final_instances: int = 20,
    seeds: list[int] | None = None,
    total_budget: int = 400,
    branch_budget: int = 8,
    checkpoint_interval: int = 20,
) -> dict:
    seeds = list(seeds) if seeds is not None else [11, 23, 37, 41, 53]

    if development_instances < 12:
        raise ValueError(
            "development_instances must be at least 12: "
            "each of the three independent development stages needs "
            "at least four instances"
        )
    if final_instances < 1:
        raise ValueError("final_instances must be >= 1")
    if not seeds:
        raise ValueError("at least one seed is required")
    if total_budget < 1 or branch_budget < 2:
        raise ValueError("total_budget must be positive and branch_budget >= 2")
    if checkpoint_interval < 1:
        raise ValueError("checkpoint_interval must be >= 1")

    thresholds = {
        "signal_alpha": 0.05,
        "signal_minimum_groups": 4,
        "signal_minimum_effect": 0.05,
        "signal_permutations": 999,
        "policy_minimum_trials": 4,
        "policy_minimum_improvement": 0.05,
        "policy_required_win_rate": 0.60,
    }
    families = (
        "smooth_additive",
        "deceptive",
        "non_additive",
        "noisy_misleading",
    )

    runs: list[dict] = []
    audits: dict[str, dict] = {}

    for family in families:
        development = [
            make_landscape(family, index)
            for index in range(development_instances)
        ]

        # Three non-overlapping sets. Do not reuse a group across stages.
        discovery_count = development_instances // 3
        validation_count = development_instances // 3

        discovery_instances = development[:discovery_count]
        validation_instances = development[
            discovery_count:discovery_count + validation_count
        ]
        policy_trial_instances = development[
            discovery_count + validation_count:
        ]

        controller, audit, audit_costs = _prepare_controller(
            family,
            discovery_instances,
            validation_instances,
            policy_trial_instances,
            seed=17,
            branch_budget=branch_budget,
            thresholds=thresholds,
        )

        audits[family] = {
            "audit_cost_by_mode": audit_costs,
            **audit,
        }

        # The controller's development audit is performed once per family.
        # Allocate its cost across held-out runs without double-counting it.
        test_cases = []

        for instance_id in range(
            development_instances,
            development_instances + final_instances,
        ):
            landscape = make_landscape(family, instance_id)
            true_signals = landscape.true_best_signals()

            for seed in seeds:
                test_cases.append(
                    (instance_id, seed, landscape, true_signals)
                )

        number_of_test_runs = len(test_cases)
        if number_of_test_runs == 0:
            raise ValueError("no held-out test runs were created")

        for mode in MODES:
            one_time_audit_cost = audit_costs[mode]
            quotient, remainder = divmod(
                one_time_audit_cost,
                number_of_test_runs,
            )

            for run_index, (
                instance_id,
                seed,
                landscape,
                true_signals,
            ) in enumerate(test_cases):

                allocated_audit_cost = quotient + (
                    1 if run_index < remainder else 0
                )
                search_budget = total_budget - allocated_audit_cost

                if search_budget < 1:
                    raise ValueError(
                        f"total budget too small for {mode}"
                    )

                result = _run_search(
                    landscape,
                    mode=mode,
                    seed=seed,
                    budget=search_budget,
                    controller=controller,
                    checkpoint_interval=checkpoint_interval,
                )

                runs.append({
                    "family": family,
                    "instance_id": instance_id,
                    "seed": seed,
                    "mode": mode,
                    "known_optimum": landscape.optimum(),
                    "audit_evaluations": allocated_audit_cost,
                    "one_time_audit_evaluations": one_time_audit_cost,
                    "search_evaluations": search_budget,
                    "total_evaluations": (
                        allocated_audit_cost + search_budget
                    ),
                    "audit_cost_included": True,
                    "signal_selected": (
                        controller.signal_result.selected_signal
                    ),
                    "signal_validation_passed": (
                        controller.signal_result.passes
                    ),
                    "signal_selected_is_true_best_single": (
                        controller.signal_result.selected_signal
                        in true_signals
                    ),
                    "policy_comparison_passed": (
                        controller.policy_result.passes
                    ),
                    "controller_thresholds_frozen_before_final_test": True,
                    **{
                        key: value
                        for key, value in result.items()
                        if key not in {"events", "decisions"}
                    },
                    "event_log": result["events"],
                    "policy_decisions": result["decisions"],
                })

        # Check per-run totals and verify that audit costs sum correctly.
        for mode in MODES:
            family_mode_runs = [
                run for run in runs
                if (
                    run["family"] == family
                    and run["mode"] == mode
                    and run["instance_id"] >= development_instances
                )
            ]

            if any(
                run["total_evaluations"] != total_budget
                for run in family_mode_runs
            ):
                raise AssertionError(
                    f"per-run budget mismatch: {family}/{mode}"
                )

            allocated_total = sum(
                run["audit_evaluations"]
                for run in family_mode_runs
            )
            if allocated_total != audit_costs[mode]:
                raise AssertionError(
                    f"audit allocation mismatch: {family}/{mode}: "
                    f"allocated={allocated_total}, "
                    f"expected={audit_costs[mode]}"
                )

    keyed = {
        (
            run["family"],
            run["instance_id"],
            run["seed"],
            run["mode"],
        ): run
        for run in runs
    }

    paired: list[dict] = []

    for family in families:
        for instance_id in range(
            development_instances,
            development_instances + final_instances,
        ):
            for seed in seeds:
                baseline = keyed[(family, instance_id, seed, "baseline")]
                for mode in MODES[1:]:
                    candidate = keyed[(family, instance_id, seed, mode)]
                    paired.append({
                        "family": family,
                        "instance_id": instance_id,
                        "seed": seed,
                        "mode": mode,
                        "best_objective_difference": (
                            candidate["best_objective"]
                            - baseline["best_objective"]
                        ),
                        "best_expected_objective_difference": (
                            candidate["best_expected_objective"]
                            - baseline["best_expected_objective"]
                        ),
                        "regret_difference": (
                            candidate["regret"] - baseline["regret"]
                        ),
                        "adaptive_won": (
                            candidate["best_expected_objective"]
                            > baseline["best_expected_objective"]
                        ),
                    })

    summary = _summarize(runs, paired, families)

    return {
        "benchmark_version": "p2-combined-adaptive-search-v2",
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
        "scientific_status": (
            "synthetic benchmark only; no biological performance "
            "or novelty claim"
        ),
    }


def write_outputs(payload: dict, output_dir: Path) -> None:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(
            f"refusing to overwrite non-empty directory: {output_dir}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    config = {
        key: value
        for key, value in payload.items()
        if key not in {"runs", "paired_results", "summary", "audits"}
    }
    (output_dir / "config.json").write_text(
        json.dumps(config, indent=2) + "\n"
    )
    (output_dir / "raw_results.json").write_text(
        json.dumps(payload, indent=2) + "\n"
    )

    with (output_dir / "summary.csv").open("w", newline="") as handle:
        rows = payload["summary"]
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "# P2 combined adaptive-search benchmark",
        "",
        "Synthetic, non-biological optimization only. Larger objective values are better.",
        "One-time development audit costs are amortized across held-out runs; each held-out run has the same total evaluation budget.",
        "Primary quality comparison: best noiseless objective visited and regret.",
        "",
        "| Family | Method | Independent instances | Best noiseless mean | Regret mean | Paired noiseless difference | 95% CI | Win rate | Signal validation | Equal total budget |",
        "|---|---|---:|---:|---:|---:|---|---:|:---:|:---:|",
    ]

    for row in payload["summary"]:
        ci = row["paired_noiseless_difference_ci95"]
        ci_text = (
            "n/a"
            if ci is None
            else f"[{ci[0]:.3f}, {ci[1]:.3f}]"
        )
        paired_mean = row["paired_noiseless_difference_mean"]
        paired_text = (
            "n/a" if paired_mean is None else f"{paired_mean:.3f}"
        )

        lines.append(
            f"| {row['family']} | {row['mode']} "
            f"| {row['n_independent_instances']} "
            f"| {row['best_noiseless_mean']:.3f} "
            f"| {row['regret_mean']:.3f} "
            f"| {paired_text} "
            f"| {ci_text} "
            f"| {row['paired_win_rate'] if row['paired_win_rate'] is not None else float('nan'):.3f} "
            f"| {row['signal_validation_rate']:.3f} "
            f"| {row['total_evaluations_equal']} |"
        )

    lines += [
        "",
        "The combined gate requires both independent signal validation and a matched-policy pass.",
        "The noisy best-observed score is reported separately from the noiseless objective.",
        "Confidence intervals resample independent landscape instances, not individual seeds.",
        "These results do not establish peptide binding, biological activity, or algorithmic novelty.",
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
    parser.add_argument(
        "--seeds",
        nargs="+",
        type=int,
        default=[11, 23, 37, 41, 53],
    )
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