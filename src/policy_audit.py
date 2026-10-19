"""Matched, equal-budget policy trials for adaptive search.

The objective is maximized. Each objective call counts as one evaluation.
Both policies receive identical starting states and budgets.

This module measures paired policy advantage; it does not by itself
establish that adaptive guidance should be trusted in production.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence
import random
import statistics


State = tuple[int, ...]
Objective = Callable[[State], float]
Policy = Callable[
    [State, Objective, random.Random],
    State,
]


@dataclass(frozen=True)
class BranchResult:
    start_score: float
    baseline_score: float
    adaptive_score: float
    paired_difference: float
    evaluations_per_branch: int


@dataclass(frozen=True)
class AuditResult:
    trials: int
    mean_advantage: float
    ci95_low: float
    ci95_high: float
    adaptive_wins: int
    baseline_wins: int
    ties: int
    adaptive_earned_trust: bool


def run_branch(
    start: State,
    objective: Objective,
    policy: Policy,
    *,
    evaluations: int,
    seed: int,
) -> float:
    """Run one policy from a copied starting state.

    The policy must evaluate candidates only through the supplied
    objective. The current state is evaluated once before searching.
    """
    if evaluations < 1:
        raise ValueError("evaluations must be >= 1")

    rng = random.Random(seed)
    state = tuple(start)
    best_score = objective(state)

    for _ in range(evaluations - 1):
        candidate = tuple(policy(state, objective, rng))
        score = objective(candidate)

        if score > best_score:
            state = candidate
            best_score = score

    return best_score


def matched_branch_trial(
    start: State,
    objective: Objective,
    baseline_policy: Policy,
    adaptive_policy: Policy,
    *,
    evaluations_per_branch: int = 40,
    seed: int = 0,
) -> BranchResult:
    """Compare policies from the same state with equal budgets."""
    start = tuple(start)
    start_score = objective(start)

    baseline_score = run_branch(
        start,
        objective,
        baseline_policy,
        evaluations=evaluations_per_branch,
        seed=seed + 1,
    )
    adaptive_score = run_branch(
        start,
        objective,
        adaptive_policy,
        evaluations=evaluations_per_branch,
        seed=seed + 2,
    )

    return BranchResult(
        start_score=start_score,
        baseline_score=baseline_score,
        adaptive_score=adaptive_score,
        paired_difference=adaptive_score - baseline_score,
        evaluations_per_branch=evaluations_per_branch,
    )


def audit_policies(
    starts: Sequence[State],
    objectives: Sequence[Objective],
    baseline_policy: Policy,
    adaptive_policy: Policy,
    *,
    evaluations_per_branch: int = 40,
    seed: int = 0,
    min_trials: int = 20,
    min_mean_advantage: float = 0.01,
) -> tuple[list[BranchResult], AuditResult]:
    """Estimate paired advantage across independent benchmark cases.

    Each starts[i] must correspond to objectives[i]. Use independent
    landscape instances as the primary unit of replication.
    """
    if len(starts) != len(objectives):
        raise ValueError("starts and objectives must have equal length")
    if len(starts) < min_trials:
        raise ValueError(f"Need at least {min_trials} matched trials")

    results = [
        matched_branch_trial(
            starts[i],
            objectives[i],
            baseline_policy,
            adaptive_policy,
            evaluations_per_branch=evaluations_per_branch,
            seed=seed + 1009 * i,
        )
        for i in range(len(starts))
    ]

    differences = [r.paired_difference for r in results]
    mean_diff = statistics.mean(differences)
    sd = statistics.stdev(differences)
    # Approximate paired t interval. For a final report, also calculate
    # a bootstrap interval resampling independent landscape instances.
    margin = 1.96 * sd / (len(differences) ** 0.5)

    wins = sum(d > 0 for d in differences)
    losses = sum(d < 0 for d in differences)
    ties = sum(d == 0 for d in differences)

    # This is an initial decision rule, not a validated production gate.
    earned = (
        len(results) >= min_trials
        and mean_diff > min_mean_advantage
        and mean_diff - margin > 0
    )

    summary = AuditResult(
        trials=len(results),
        mean_advantage=mean_diff,
        ci95_low=mean_diff - margin,
        ci95_high=mean_diff + margin,
        adaptive_wins=wins,
        baseline_wins=losses,
        ties=ties,
        adaptive_earned_trust=earned,
    )
    return results, summary