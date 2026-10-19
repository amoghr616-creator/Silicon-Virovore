"""Context-dependent selection between baseline and adaptive search.

This module is intended for controlled synthetic-search experiments first.
It does not assume that hotspot association implies biological relevance.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import random
from typing import Callable

from src.online_policy_audit import (
    OnlineAuditResult,
    run_online_audit,
)


State = tuple[int, ...]
Objective = Callable[[State, random.Random], float]


@dataclass(frozen=True)
class SearchContext:
    hotspot_p_value: float | None
    hotspot_stability: float
    stagnation_generations: int
    diversity: float
    hotspot_significance_validated: bool = False


@dataclass(frozen=True)
class PolicyDecision:
    action: str
    reason: str
    audit_cost: int
    remaining_evaluations: int


@dataclass(frozen=True)
class ContextualStepResult:
    action: str
    selected_policy: str
    reason: str
    next_state: State
    evaluations_spent: int
    remaining_evaluations: int
    audit_result: OnlineAuditResult | None


def matched_audit_cost(
    *,
    branch_budget: int,
    replications: int,
) -> int:
    """Count every baseline and adaptive branch evaluation."""
    if branch_budget < 2:
        raise ValueError("branch_budget must be at least 2")
    if replications < 1:
        raise ValueError("replications must be at least 1")
    return 2 * branch_budget * replications


def decide_contextual_policy(
    context: SearchContext,
    *,
    remaining_evaluations: int,
    branch_budget: int = 8,
    replications: int = 4,
    minimum_search_budget_after_audit: int = 16,
    hotspot_alpha: float = 0.05,
    minimum_hotspot_stability: float = 0.60,
    stagnation_threshold: int = 3,
    minimum_diversity: float = 0.35,
) -> PolicyDecision:
    """Choose a policy from evidence quality and current search conditions."""
    audit_cost = matched_audit_cost(
        branch_budget=branch_budget,
        replications=replications,
    )
    if remaining_evaluations < 0:
        raise ValueError("remaining_evaluations must be nonnegative")
    if minimum_search_budget_after_audit < 0:
        raise ValueError("minimum_search_budget_after_audit must be nonnegative")
    if not 0.0 < hotspot_alpha < 1.0:
        raise ValueError("hotspot_alpha must be between 0 and 1")
    if not 0.0 <= minimum_hotspot_stability <= 1.0:
        raise ValueError("minimum_hotspot_stability must be between 0 and 1")
    if stagnation_threshold < 1:
        raise ValueError("stagnation_threshold must be at least 1")
    if not 0.0 <= minimum_diversity <= 1.0:
        raise ValueError("minimum_diversity must be between 0 and 1")
    if context.hotspot_p_value is not None and (
        not math.isfinite(context.hotspot_p_value)
        or not 0.0 <= context.hotspot_p_value <= 1.0
    ):
        raise ValueError("hotspot_p_value must be between 0 and 1")
    if not math.isfinite(context.hotspot_stability) or not (
        0.0 <= context.hotspot_stability <= 1.0
    ):
        raise ValueError("hotspot_stability must be between 0 and 1")
    if context.stagnation_generations < 0:
        raise ValueError("stagnation_generations must be nonnegative")
    if not math.isfinite(context.diversity) or not 0.0 <= context.diversity <= 1.0:
        raise ValueError("diversity must be between 0 and 1")

    reliable_hotspot = (
        context.hotspot_significance_validated
        and
        context.hotspot_p_value is not None
        and context.hotspot_p_value < hotspot_alpha
        and context.hotspot_stability >= minimum_hotspot_stability
    )
    prolonged_stagnation = (
        context.stagnation_generations >= stagnation_threshold
    )
    low_diversity = context.diversity < minimum_diversity

    if not reliable_hotspot:
        action = "baseline"
        reason = (
            "hotspot_significance_unvalidated"
            if not context.hotspot_significance_validated
            else "hotspot_evidence_not_reliable"
        )
    elif prolonged_stagnation and low_diversity:
        if remaining_evaluations >= audit_cost + minimum_search_budget_after_audit:
            action = "audit"
            reason = "reliable_hotspot_stagnation_and_low_diversity"
        else:
            action = "baseline"
            reason = "audit_unaffordable_with_search_reserve"
    elif (
        not prolonged_stagnation
        and context.diversity >= minimum_diversity
    ):
        action = "adaptive"
        reason = "reliable_hotspot_with_diversity_and_no_prolonged_stagnation"
    else:
        action = "baseline"
        reason = "context_does_not_support_adaptive_or_audit"

    return PolicyDecision(
        action=action,
        reason=reason,
        audit_cost=audit_cost,
        remaining_evaluations=remaining_evaluations,
    )


def run_contextual_step(
    *,
    context: SearchContext,
    start_state: State,
    objective: Objective,
    preferred_position: int | None,
    remaining_evaluations: int,
    branch_budget: int = 8,
    replications: int = 4,
    seed: int = 0,
    minimum_search_budget_after_audit: int = 16,
    hotspot_alpha: float = 0.05,
    minimum_hotspot_stability: float = 0.60,
    stagnation_threshold: int = 3,
    minimum_diversity: float = 0.35,
) -> ContextualStepResult:
    """Apply a context decision, running and charging an audit if selected.

    For baseline and adaptive decisions, the caller performs its ordinary
    move and this function leaves the state unchanged.
    """
    start_state = tuple(start_state)
    if not start_state:
        raise ValueError("start_state must not be empty")
    if any(bit not in (0, 1) for bit in start_state):
        raise ValueError("start_state must contain only binary values")

    decision = decide_contextual_policy(
        context,
        remaining_evaluations=remaining_evaluations,
        branch_budget=branch_budget,
        replications=replications,
        minimum_search_budget_after_audit=minimum_search_budget_after_audit,
        hotspot_alpha=hotspot_alpha,
        minimum_hotspot_stability=minimum_hotspot_stability,
        stagnation_threshold=stagnation_threshold,
        minimum_diversity=minimum_diversity,
    )
    if decision.action != "audit":
        return ContextualStepResult(
            action=decision.action,
            selected_policy=decision.action,
            reason=decision.reason,
            next_state=start_state,
            evaluations_spent=0,
            remaining_evaluations=remaining_evaluations,
            audit_result=None,
        )

    audit = run_online_audit(
        start_state=start_state,
        objective=objective,
        preferred_position=preferred_position,
        branch_budget=branch_budget,
        replications=replications,
        seed=seed,
    )
    selected_trials = [
        trial.baseline if audit.selected_policy == "baseline" else trial.adaptive
        for trial in audit.trials
    ]
    best_branch = max(selected_trials, key=lambda branch: branch.best_score)
    evaluations_spent = audit.total_objective_evaluations
    return ContextualStepResult(
        action=decision.action,
        selected_policy=audit.selected_policy,
        reason=audit.rationale,
        next_state=best_branch.best_state,
        evaluations_spent=evaluations_spent,
        remaining_evaluations=remaining_evaluations - evaluations_spent,
        audit_result=audit,
    )