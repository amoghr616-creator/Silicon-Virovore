"""Online matched-policy audits for synthetic P2 search.

At a saved state, run baseline and adaptive branches using:
- identical starting states;
- paired random seeds;
- equal objective-evaluation budgets.

Recommend adaptive guidance only when its paired gains meet
predeclared effect-size and win-rate thresholds.

This is an experimental controller, not a validated statistical test.
Repeated audits need a prespecified schedule and multiplicity control.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
import random
import statistics
from typing import Callable


State = tuple[int, ...]
Objective = Callable[[State, random.Random], float]


@dataclass(frozen=True)
class BranchOutcome:
    policy: str
    start_state: State
    start_score: float
    best_state: State
    best_score: float
    final_state: State
    final_score: float
    evaluations: int
    evaluated_states: tuple[State, ...]
    failed: bool = False
    failure_reason: str | None = None

    @property
    def best_gain(self) -> float:
        return self.best_score - self.start_score


@dataclass(frozen=True)
class PairedTrial:
    replicate: int
    seed: int
    baseline: BranchOutcome
    adaptive: BranchOutcome
    paired_advantage: float


@dataclass(frozen=True)
class OnlineAuditResult:
    start_state: State
    preferred_position: int | None
    trials: tuple[PairedTrial, ...]
    mean_baseline_gain: float
    mean_adaptive_gain: float
    mean_paired_advantage: float
    adaptive_win_rate: float
    selected_policy: str
    total_objective_evaluations: int
    branch_budget: int
    replications: int
    minimum_mean_advantage: float
    required_win_rate: float
    rationale: str

    def to_dict(self) -> dict:
        """Return a JSON-serializable result."""
        return asdict(self)


def _validate_state(state: State) -> None:
    if not state:
        raise ValueError("state must not be empty")
    if any(bit not in (0, 1) for bit in state):
        raise ValueError("state must contain only binary values")


def _choose_position(
    state: State,
    preferred_position: int | None,
    rng: random.Random,
    *,
    guidance_probability: float,
) -> int:
    """Choose a reversible bit flip for one branch step.

    Baseline samples a position uniformly.
    Adaptive guidance favors activating the preferred bit when it is off.
    Once active, the preferred bit is not intentionally switched off.
    """

    guidance_draw = rng.random()
    fallback_draw = rng.random()

    fallback_position = min(
        int(fallback_draw * len(state)),
        len(state) - 1,
    )

    if (
        preferred_position is None
        or not 0 <= preferred_position < len(state)
    ):
        return fallback_position

    if (
        state[preferred_position] == 0
        and guidance_draw >= 1.0 - guidance_probability
    ):
        return preferred_position

    alternatives = [
        index
        for index in range(len(state))
        if index != preferred_position
    ]

    if not alternatives:
        return preferred_position

    alternative_index = min(
        int(fallback_draw * len(alternatives)),
        len(alternatives) - 1,
    )
    return alternatives[alternative_index]


def _run_branch(
    *,
    policy: str,
    start_state: State,
    objective: Objective,
    budget: int,
    seed: int,
    preferred_position: int | None,
    guidance_probability: float,
) -> BranchOutcome:
    """Run one branch; the starting-state evaluation counts toward budget."""

    policy_rng = random.Random(seed)
    objective_rng = random.Random(seed + 1_000_003)

    state = tuple(start_state)
    start_score: float | None = None
    current_score: float | None = None
    best_score = float("-inf")
    best_state = state
    evaluations = 0
    visited_states: list[State] = []

    for step in range(budget):
        if step == 0:
            candidate = state
        else:
            position = _choose_position(
                state,
                preferred_position if policy == "adaptive" else None,
                policy_rng,
                guidance_probability=guidance_probability,
            )
            mutable = list(state)
            mutable[position] = 1 - mutable[position]
            candidate = tuple(mutable)

        score = float(objective(candidate, objective_rng))
        if not math.isfinite(score):
            raise ValueError("objective returned a non-finite score")

        evaluations += 1
        visited_states.append(candidate)

        if start_score is None:
            start_score = score
            current_score = score

        if score > best_score:
            best_score = score
            best_state = candidate

        # Fitness-guided acceptance is identical in both branches.
        if current_score is None or score > current_score:
            state = candidate
            current_score = score

    assert start_score is not None
    assert current_score is not None

    return BranchOutcome(
        policy=policy,
        start_state=tuple(start_state),
        start_score=start_score,
        best_state=best_state,
        best_score=best_score,
        final_state=state,
        final_score=current_score,
        evaluations=evaluations,
        evaluated_states=tuple(visited_states),
    )


def run_online_audit(
    *,
    start_state: State,
    objective: Objective,
    preferred_position: int | None,
    branch_budget: int = 8,
    replications: int = 4,
    seed: int = 0,
    guidance_probability: float = 0.75,
    minimum_mean_advantage: float = 0.05,
    required_win_rate: float = 0.60,
) -> OnlineAuditResult:
    """Compare policies from the same saved state before choosing a policy.

    Every replication runs one baseline branch and one adaptive branch
    with matching seeds and equal objective-evaluation budgets.

    The comparison is conditional on this starting state. Replications
    use different random seeds; they are not independent landscapes.
    """

    start_state = tuple(start_state)
    _validate_state(start_state)

    if branch_budget < 2:
        raise ValueError("branch_budget must be at least 2")
    if replications < 1:
        raise ValueError("replications must be at least 1")
    if not 0.0 < guidance_probability <= 1.0:
        raise ValueError("guidance_probability must be in (0, 1]")
    if not math.isfinite(minimum_mean_advantage):
        raise ValueError("minimum_mean_advantage must be finite")
    if not 0.0 <= required_win_rate <= 1.0:
        raise ValueError("required_win_rate must be between 0 and 1")
    if (
        preferred_position is not None
        and not 0 <= preferred_position < len(start_state)
    ):
        raise ValueError("preferred_position is outside the state")

    trials: list[PairedTrial] = []

    for replicate in range(replications):
        trial_seed = seed + replicate * 10_007

        # Both branches receive the exact same state and matched seeds.
        baseline = _run_branch(
            policy="baseline",
            start_state=start_state,
            objective=objective,
            budget=branch_budget,
            seed=trial_seed,
            preferred_position=None,
            guidance_probability=guidance_probability,
        )

        adaptive = _run_branch(
            policy="adaptive",
            start_state=start_state,
            objective=objective,
            budget=branch_budget,
            seed=trial_seed,
            preferred_position=preferred_position,
            guidance_probability=guidance_probability,
        )

        trials.append(
            PairedTrial(
                replicate=replicate,
                seed=trial_seed,
                baseline=baseline,
                adaptive=adaptive,
                paired_advantage=(
                    adaptive.best_gain - baseline.best_gain
                ),
            )
        )

    baseline_gains = [trial.baseline.best_gain for trial in trials]
    adaptive_gains = [trial.adaptive.best_gain for trial in trials]
    differences = [trial.paired_advantage for trial in trials]

    mean_baseline_gain = statistics.mean(baseline_gains)
    mean_adaptive_gain = statistics.mean(adaptive_gains)
    mean_advantage = statistics.mean(differences)
    win_rate = sum(value > 0 for value in differences) / len(differences)

    # Conservative tie behavior: adaptive must pass both thresholds.
    use_adaptive = (
        mean_advantage >= minimum_mean_advantage
        and win_rate >= required_win_rate
    )

    if use_adaptive:
        selected_policy = "adaptive"
        rationale = "paired_gain_and_win_rate_thresholds_passed"
    else:
        selected_policy = "baseline"
        rationale = "insufficient_paired_evidence_for_adaptive"

    return OnlineAuditResult(
        start_state=start_state,
        preferred_position=preferred_position,
        trials=tuple(trials),
        mean_baseline_gain=mean_baseline_gain,
        mean_adaptive_gain=mean_adaptive_gain,
        mean_paired_advantage=mean_advantage,
        adaptive_win_rate=win_rate,
        selected_policy=selected_policy,
        total_objective_evaluations=2 * replications * branch_budget,
        branch_budget=branch_budget,
        replications=replications,
        minimum_mean_advantage=minimum_mean_advantage,
        required_win_rate=required_win_rate,
        rationale=rationale,
    )