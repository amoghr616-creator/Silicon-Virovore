"""Fail-closed gate for target-context ranking.

The gate does not claim that PEPBI predicts HERV-K Env binding. It requires a
predeclared, group-held-out external transfer experiment to show that adding
complex-specific structural context improves prediction of direct experimental
measurements before those target-context scores may alter parent ordering.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SUMMARY = (
    ROOT / "prototype4/results/target_context_external_validation_fresh/summary.json"
)


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def evaluate_external_validation_summary(
    summary: dict[str, Any] | None,
    *,
    summary_path: str | Path | None = None,
    min_independent_groups: int = 10,
    alpha: float = 0.05,
) -> dict[str, Any]:
    """Return explicit authorization for target-aware candidate reordering.

    Authorization requires all of: at least 10 held-out binding groups, a
    positive P3-minus-P4 macro-MAE benefit, a group-bootstrap 95% CI wholly
    above zero, one-sided exact paired sign-flip p < alpha, and more groups
    improved than worsened. A missing/failed condition blocks ranking.
    """
    source_path = Path(summary_path).expanduser() if summary_path else None
    source_sha256 = None
    if source_path is not None and source_path.is_file():
        digest = hashlib.sha256(source_path.read_bytes()).hexdigest()
        source_sha256 = digest

    if not isinstance(summary, dict):
        return {
            "status": "missing_external_validation",
            "allow_target_selection": False,
            "reason": "External validation summary is absent or unreadable.",
            "summary_path": str(source_path) if source_path else None,
            "summary_sha256": source_sha256,
        }

    effect = summary.get("paired_primary_effect") or {}
    benefit = _finite(effect.get("macro_mae_benefit_kcal_mol"))
    ci = effect.get("group_bootstrap_95ci") or []
    ci_low = _finite(ci[0]) if len(ci) >= 2 else None
    ci_high = _finite(ci[1]) if len(ci) >= 2 else None
    p_improvement = _finite(effect.get("exact_one_sided_signflip_p"))
    groups = int(effect.get("n_independent_groups", summary.get("independent_binding_groups", 0)) or 0)
    improved = int(effect.get("groups_improved", 0) or 0)
    worsened = int(effect.get("groups_worsened", 0) or 0)

    checks = {
        "minimum_independent_groups": groups >= min_independent_groups,
        "positive_macro_mae_benefit": benefit is not None and benefit > 0,
        "bootstrap_ci_excludes_zero_in_benefit_direction": ci_low is not None and ci_low > 0,
        "one_sided_signflip_p_below_alpha": p_improvement is not None and p_improvement < alpha,
        "more_groups_improved_than_worsened": improved > worsened,
    }
    allowed = all(checks.values())
    failures = [name for name, passed in checks.items() if not passed]

    return {
        "status": "external_transfer_gate_passed" if allowed else "external_transfer_gate_failed",
        "allow_target_selection": bool(allowed),
        "reason": (
            "External group-held-out benchmark supports adding target-context evidence."
            if allowed else
            "Target-context evidence is not authorized to change selection; failed checks: "
            + ", ".join(failures)
        ),
        "summary_path": str(source_path) if source_path else None,
        "summary_sha256": source_sha256,
        "source_sha256": summary.get("source_sha256"),
        "outcome": summary.get("outcome"),
        "n_rows": int(summary.get("rows", 0) or 0),
        "n_independent_groups": groups,
        "macro_mae_benefit_kcal_mol": benefit,
        "group_bootstrap_95ci": [ci_low, ci_high],
        "one_sided_signflip_p": p_improvement,
        "groups_improved": improved,
        "groups_worsened": worsened,
        "alpha": float(alpha),
        "minimum_independent_groups": int(min_independent_groups),
        "checks": checks,
        "transfer_scope_note": (
            "PEPBI is an external peptide-protein benchmark. Passing this gate would support "
            "transfer of target-context features in that benchmark; it would not validate "
            "HERV-K Env 9SYA binding or therapeutic activity."
        ),
    }


def load_external_validation_gate(
    summary_path: str | Path | None = None,
) -> dict[str, Any]:
    """Load the frozen external validation summary and evaluate it."""
    configured = (
        summary_path
        or os.environ.get("SILICON_VIROVORE_P4_EXTERNAL_VALIDATION_SUMMARY")
        or DEFAULT_SUMMARY
    )
    path = Path(configured).expanduser()
    if not path.is_file():
        return evaluate_external_validation_summary(None, summary_path=path)
    try:
        summary = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return evaluate_external_validation_summary(None, summary_path=path)
    return evaluate_external_validation_summary(summary, summary_path=path)
