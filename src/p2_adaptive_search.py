"""Isolated P2 adaptive-search validation components.

This module is intentionally general-purpose and is not imported by the
biological runner. It contains:

* grouped independent-signal validation with a maximum-statistic permutation
  correction; and
* matched, equal-budget policy trials plus a frozen combined decision gate.

The public classes accept only synthetic-search observations and callables.
They do not know about peptides, docking, receptors, or ARISE.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import itertools
import math
import random
import statistics
from typing import Callable, Iterable, Sequence


@dataclass(frozen=True)
class SignalObservation:
    observation_id: str
    independent_group_id: str
    run_id: str
    lineage_id: str
    context_id: str
    signal: str
    outcome: float
    split: str
    observation_order: int


@dataclass(frozen=True)
class SignalValidationResult:
    selected_signal: str | None
    discovery_groups: tuple[str, ...]
    validation_groups: tuple[str, ...]
    candidate_signals: tuple[str, ...]
    discovery_effects: dict[str, float]
    adjusted_p_values: dict[str, float]
    selected_discovery_effect: float | None
    validation_effect: float | None
    validation_ci95: tuple[float, float] | None
    validation_win_rate: float | None
    passes: bool
    reason: str
    permutations: int
    minimum_groups: int
    alpha: float
    minimum_effect: float
    signal_count: int


class IndependentSignalValidator:
    """Validate a selected signal using independent grouped observations."""

    def __init__(
        self,
        *,
        alpha: float = 0.05,
        minimum_groups: int = 4,
        minimum_effect: float = 0.05,
        permutations: int = 999,
        seed: int = 0,
    ) -> None:
        if not 0.0 < alpha < 1.0:
            raise ValueError("alpha must be between 0 and 1")
        if minimum_groups < 2 or permutations < 1:
            raise ValueError("minimum_groups and permutations must be positive")
        if minimum_effect < 0.0 or not math.isfinite(minimum_effect):
            raise ValueError("minimum_effect must be finite and non-negative")
        self.alpha = alpha
        self.minimum_groups = minimum_groups
        self.minimum_effect = minimum_effect
        self.permutations = permutations
        self.seed = seed
        self.result: SignalValidationResult | None = None

    @staticmethod
    def _validate_records(records: Sequence[SignalObservation], expected_split: str) -> None:
        if not records:
            raise ValueError("signal validation requires observations")
        ids = [record.observation_id for record in records]
        if len(ids) != len(set(ids)):
            raise ValueError("duplicate signal observation_id")
        if any(record.split != expected_split for record in records):
            raise ValueError(f"all observations must be assigned to {expected_split}")
        if any(not math.isfinite(record.outcome) for record in records):
            raise ValueError("signal outcomes must be finite")

    @staticmethod
    def _matrix(records: Sequence[SignalObservation]) -> tuple[tuple[str, ...], tuple[str, ...], list[list[float]]]:
        groups = tuple(sorted({record.independent_group_id for record in records}))
        signals = tuple(sorted({record.signal for record in records}))
        values = {(record.independent_group_id, record.signal): record.outcome for record in records}
        if len(values) != len(records):
            raise ValueError("duplicate signal/group observations are not comparable")
        if any((group, signal) not in values for group in groups for signal in signals):
            raise ValueError("signal validation requires a complete group-by-signal matrix")
        return groups, signals, [[values[(group, signal)] for signal in signals] for group in groups]

    @staticmethod
    def _effects(matrix: Sequence[Sequence[float]]) -> list[float]:
        if not matrix:
            return []
        return [
            statistics.mean(row[index] - statistics.mean(row) for row in matrix)
            for index in range(len(matrix[0]))
        ]

    def _adjusted_p_values(self, matrix: list[list[float]], observed: list[float]) -> list[float]:
        rng = random.Random(self.seed)
        null_maxima: list[float] = []
        for _ in range(self.permutations):
            permuted = [rng.sample(row, len(row)) for row in matrix]
            effects = self._effects(permuted)
            null_maxima.append(max(abs(effect) for effect in effects))
        return [
            (1.0 + sum(null >= abs(effect) for null in null_maxima)) / (self.permutations + 1.0)
            for effect in observed
        ]

    @staticmethod
    def _bootstrap_ci(values: Sequence[float], seed: int, iterations: int = 2000) -> tuple[float, float] | None:
        if len(values) < 2:
            return None
        rng = random.Random(seed)
        means = [statistics.mean(rng.choices(list(values), k=len(values))) for _ in range(iterations)]
        return (statistics.quantiles(means, n=100)[1], statistics.quantiles(means, n=100)[-2])

    def fit_and_validate(
        self,
        discovery: Sequence[SignalObservation],
        validation: Sequence[SignalObservation],
    ) -> SignalValidationResult:
        """Select from discovery, then test only that selection on new groups."""
        self._validate_records(discovery, "discovery")
        self._validate_records(validation, "validation")
        discovery_groups, signals, matrix = self._matrix(discovery)
        validation_groups, validation_signals, validation_matrix = self._matrix(validation)
        if set(discovery_groups) & set(validation_groups):
            raise ValueError("independent groups cannot cross discovery and validation")
        if signals != validation_signals:
            raise ValueError("discovery and validation signal sets must match")
        effects = self._effects(matrix)
        adjusted = self._adjusted_p_values(matrix, effects)
        selected_index = max(range(len(signals)), key=lambda index: (effects[index], signals[index]))
        selected = signals[selected_index]
        validation_index = signals.index(selected)
        validation_effects = self._effects(validation_matrix)
        per_group_effects = [
            row[validation_index] - statistics.mean(row)
            for row in validation_matrix
        ]
        validation_effect = statistics.mean(per_group_effects) if per_group_effects else None
        ci = self._bootstrap_ci(per_group_effects, self.seed + 1)
        win_rate = (
            sum(effect > 0.0 for effect in per_group_effects) / len(per_group_effects)
            if per_group_effects else None
        )
        reasons: list[str] = []
        if len(discovery_groups) < self.minimum_groups:
            reasons.append("insufficient_discovery_groups")
        if len(validation_groups) < self.minimum_groups:
            reasons.append("insufficient_validation_groups")
        if adjusted[selected_index] > self.alpha:
            reasons.append("maximum_statistic_not_significant")
        if validation_effect is None or validation_effect < self.minimum_effect:
            reasons.append("validation_effect_below_threshold")
        passes = not reasons
        self.result = SignalValidationResult(
            # Preserve the discovery winner even when validation fails. The
            # unvalidated comparator may use it; gated modes must not.
            selected,
            discovery_groups,
            validation_groups,
            signals,
            {signal: effect for signal, effect in zip(signals, effects)},
            {signal: p_value for signal, p_value in zip(signals, adjusted)},
            effects[selected_index], validation_effect, ci, win_rate, passes,
            ";".join(reasons) if reasons else "passed",
            self.permutations, self.minimum_groups, self.alpha, self.minimum_effect, len(signals),
        )
        return self.result


@dataclass(frozen=True)
class MatchedTrialResult:
    trial_id: str
    independent_group_id: str
    seed: int
    start_state: tuple[int, ...]
    baseline_best: float
    adaptive_best: float
    baseline_evaluations: int
    adaptive_evaluations: int
    paired_difference: float
    baseline_failure: str | None
    adaptive_failure: str | None


@dataclass(frozen=True)
class PolicyComparisonResult:
    trials: tuple[MatchedTrialResult, ...]
    independent_trials: int
    mean_difference: float | None
    median_difference: float | None
    ci95: tuple[float, float] | None
    win_rate: float | None
    passes: bool
    reason: str
    minimum_trials: int
    minimum_improvement: float
    required_win_rate: float


class MatchedPolicyTrialRunner:
    """Run baseline and adaptive policies from identical states and budgets."""

    def __init__(self, *, branch_budget: int = 20) -> None:
        if branch_budget < 2:
            raise ValueError("branch_budget must be at least 2")
        self.branch_budget = branch_budget

    @staticmethod
    def _choose(rng: random.Random, state: tuple[int, ...], preferred: int | None) -> int:
        available = [index for index, value in enumerate(state) if not value]
        if not available:
            return rng.randrange(len(state))
        if preferred is not None and preferred in available and rng.random() >= 0.25:
            return preferred
        return rng.choice(available)

    def run_trial(
        self,
        *,
        trial_id: str,
        independent_group_id: str,
        seed: int,
        start_state: tuple[int, ...],
        objective: Callable[[tuple[int, ...], random.Random], float],
        preferred_position: int | None,
    ) -> MatchedTrialResult:
        def run(preferred: int | None, branch_seed: int) -> tuple[float, int, str | None]:
            rng = random.Random(branch_seed)
            state = tuple(start_state)
            best = float("-inf")
            failure = None
            evaluations = 0
            for _ in range(self.branch_budget):
                try:
                    outcome = float(objective(state, rng))
                    if not math.isfinite(outcome):
                        raise ValueError("non-finite objective")
                except Exception as exc:  # record branch failure, do not hide it
                    failure = type(exc).__name__
                    break
                evaluations += 1
                best = max(best, outcome)
                position = self._choose(rng, state, preferred)
                next_state = list(state)
                next_state[position] = 1
                state = tuple(next_state)
            return best, evaluations, failure

        baseline, baseline_evals, baseline_failure = run(None, seed)
        adaptive, adaptive_evals, adaptive_failure = run(preferred_position, seed)
        if baseline_failure or adaptive_failure:
            difference = float("nan")
        else:
            difference = adaptive - baseline
        return MatchedTrialResult(
            trial_id, independent_group_id, seed, tuple(start_state), baseline, adaptive,
            baseline_evals, adaptive_evals, difference, baseline_failure, adaptive_failure,
        )

    def compare(
        self,
        trials: Sequence[MatchedTrialResult],
        *,
        minimum_trials: int = 4,
        minimum_improvement: float = 0.05,
        required_win_rate: float = 0.60,
    ) -> PolicyComparisonResult:
        valid = [trial for trial in trials if math.isfinite(trial.paired_difference)]
        differences = [trial.paired_difference for trial in valid]
        ci = IndependentSignalValidator._bootstrap_ci(differences, 17) if differences else None
        mean_difference = statistics.mean(differences) if differences else None
        median_difference = statistics.median(differences) if differences else None
        win_rate = sum(value > 0.0 for value in differences) / len(differences) if differences else None
        reasons: list[str] = []
        if len(valid) < minimum_trials:
            reasons.append("insufficient_independent_trials")
        if mean_difference is None or mean_difference < minimum_improvement:
            reasons.append("mean_improvement_below_threshold")
        if win_rate is None or win_rate < required_win_rate:
            reasons.append("win_rate_below_threshold")
        return PolicyComparisonResult(
            tuple(trials), len(valid), mean_difference, median_difference, ci, win_rate,
            not reasons, ";".join(reasons) if reasons else "passed",
            minimum_trials, minimum_improvement, required_win_rate,
        )


@dataclass(frozen=True)
class PolicyDecision:
    checkpoint: str
    mode: str
    adaptive_selected: bool
    signal_passed: bool
    policy_passed: bool
    rationale: str


class CombinedAdaptiveController:
    """Freeze development evidence and make auditable policy decisions."""

    def __init__(
        self,
        *,
        signal_result: SignalValidationResult,
        policy_result: PolicyComparisonResult,
    ) -> None:
        self.signal_result = signal_result
        self.policy_result = policy_result
        self.decisions: list[PolicyDecision] = []

    def decide(self, *, checkpoint: str, mode: str) -> PolicyDecision:
        signal_passed = self.signal_result.passes
        policy_passed = self.policy_result.passes
        if mode == "baseline":
            selected = False
            rationale = "baseline_control"
        elif mode == "adaptive_unvalidated":
            selected = self.signal_result.selected_signal is not None
            rationale = "raw_discovery_signal_without_gate" if selected else "no_discovery_signal"
        elif mode == "signal_validation_only":
            selected = signal_passed
            rationale = "signal_validation_passed" if selected else f"signal_validation_failed:{self.signal_result.reason}"
        elif mode == "matched_policy_only":
            selected = policy_passed and self.signal_result.selected_signal is not None
            rationale = "matched_policy_passed" if selected else f"matched_policy_failed:{self.policy_result.reason}"
        elif mode == "combined":
            selected = signal_passed and policy_passed
            rationale = "signal_and_matched_policy_passed" if selected else f"combined_gate_failed:signal={signal_passed};policy={policy_passed}"
        else:
            raise ValueError(f"unknown policy mode: {mode}")
        decision = PolicyDecision(checkpoint, mode, selected, signal_passed, policy_passed, rationale)
        self.decisions.append(decision)
        return decision


def asdict_records(records: Iterable[object]) -> list[dict]:
    """Serialize dataclass records for machine-readable event logs."""
    return [asdict(record) for record in records]
