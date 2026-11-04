"""Family-matched risk--coverage analysis for nested LOFO predictions.

This analysis consumes out-of-fold predictions from prototype4/p4.py validate.
It uses ``heldout_group`` as the independent family unit, keeps at least one
record per held-out family at each requested coverage, and compares selective
retention with random retention using the same per-family retention counts.
It evaluates prediction error on an external PPI mutation-cycle benchmark; it
is not evidence of HERV-K Env binding or peptide activity.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


REQUIRED_COLUMNS = {
    "feature_uid",
    "heldout_group",
    "split_unit",
    "y_true",
    "y_pred",
    "calibrated_abs_error",
}


def _finite_column(frame: pd.DataFrame, name: str) -> np.ndarray:
    values = pd.to_numeric(frame[name], errors="coerce").to_numpy(dtype=float)
    if values.ndim != 1 or not np.isfinite(values).all():
        raise ValueError(f"Column {name!r} contains missing or non-finite values.")
    return values


def analyze_selective_abstention(
    predictions: pd.DataFrame,
    *,
    random_repetitions: int = 1000,
    family_bootstrap_repetitions: int = 10_000,
    seed: int = 20_261_103,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Compare selective and family-matched random retention on LOFO OOF data."""
    missing = REQUIRED_COLUMNS - set(predictions.columns)
    if missing:
        raise ValueError(f"OOF prediction table is missing columns: {sorted(missing)}")
    if random_repetitions < 100:
        raise ValueError("random_repetitions must be at least 100.")
    if family_bootstrap_repetitions < 1000:
        raise ValueError("family_bootstrap_repetitions must be at least 1000.")
    if predictions.empty:
        raise ValueError("OOF prediction table is empty.")
    if predictions["feature_uid"].duplicated().any():
        raise ValueError("Each feature_uid must have exactly one outer-LOFO prediction.")
    split_units = set(predictions["split_unit"].astype(str))
    if split_units != {"family"}:
        raise ValueError(
            "This analysis requires family-held-out predictions: split_unit must be 'family'."
        )

    frame = predictions.copy()
    frame["heldout_group"] = frame["heldout_group"].astype(str)
    if (frame["heldout_group"].str.len() == 0).any():
        raise ValueError("heldout_group values must be nonempty.")

    y_true = _finite_column(frame, "y_true")
    y_pred = _finite_column(frame, "y_pred")
    uncertainty = _finite_column(frame, "calibrated_abs_error")
    if np.any(uncertainty < 0):
        raise ValueError("calibrated_abs_error must be nonnegative.")
    abs_error = np.abs(y_true - y_pred)

    groups = sorted(frame["heldout_group"].unique().tolist())
    if len(groups) < 5:
        raise ValueError(
            f"Only {len(groups)} held-out families; at least five are required."
        )

    # Use stable input row indices; retain an identical per-family quota in the
    # selective and random arms so family composition cannot explain the result.
    group_indices = [
        np.flatnonzero(frame["heldout_group"].to_numpy() == group)
        for group in groups
    ]
    n = len(frame)
    n_unique_uncertainty = int(np.unique(np.round(uncertainty, 10)).size)
    rankable = n_unique_uncertainty >= max(4, int(math.ceil(0.10 * n)))

    rng = np.random.default_rng(seed)
    rows: list[dict[str, Any]] = []
    coverage_specs = [x / 10 for x in range(1, 11)]

    for target_coverage in coverage_specs:
        quotas = [
            n_g if target_coverage >= 1.0 else min(
                n_g, max(1, int(math.ceil(target_coverage * n_g)))
            )
            for n_g in map(len, group_indices)
        ]
        selected_group_mae = []
        random_group_mae = np.empty((random_repetitions, len(groups)), dtype=float)

        for group_position, (indices, quota) in enumerate(zip(group_indices, quotas)):
            # Stable tie breaking by original row order. Calibration was fitted
            # only with inner OOF data in the corresponding outer training fold.
            order = np.argsort(uncertainty[indices], kind="mergesort")
            selected = indices[order[:quota]]
            selected_group_mae.append(float(abs_error[selected].mean()))

            for repetition in range(random_repetitions):
                random_selected = rng.choice(indices, size=quota, replace=False)
                random_group_mae[repetition, group_position] = float(
                    abs_error[random_selected].mean()
                )

        selected_group_mae_array = np.asarray(selected_group_mae, dtype=float)
        selected_macro_mae = float(selected_group_mae_array.mean())
        random_macro_mae = random_group_mae.mean(axis=1)
        random_macro_mean = float(random_macro_mae.mean())
        improvement = random_macro_mean - selected_macro_mae
        realized_n = int(sum(quotas))
        realized_coverage = realized_n / n

        # Family-cluster bootstrap: resample entire held-out families, never rows.
        # For each bootstrap family draw, also sample one matched random-retention
        # replicate per family to propagate finite Monte Carlo-control variability.
        boot_selected = np.empty(family_bootstrap_repetitions, dtype=float)
        boot_delta = np.empty(family_bootstrap_repetitions, dtype=float)
        group_count = len(groups)
        for b in range(family_bootstrap_repetitions):
            sampled_group_positions = rng.integers(0, group_count, size=group_count)
            selected_values = selected_group_mae_array[sampled_group_positions]
            random_rep_for_group = rng.integers(
                0, random_repetitions, size=group_count
            )
            random_values = random_group_mae[
                random_rep_for_group, sampled_group_positions
            ]
            boot_selected[b] = float(selected_values.mean())
            boot_delta[b] = float((random_values - selected_values).mean())

        # One-sided Monte Carlo test: fraction of matched random selections with
        # a macro-family MAE no worse than the uncertainty-selected set.
        random_at_or_better = int(np.count_nonzero(random_macro_mae <= selected_macro_mae))
        randomization_p = (random_at_or_better + 1) / (random_repetitions + 1)

        rows.append({
            "target_coverage": target_coverage,
            "realized_coverage": realized_coverage,
            "retained_n": realized_n,
            "n_families": group_count,
            "selective_macro_family_mae": selected_macro_mae,
            "selective_macro_family_mae_ci95_low": float(np.quantile(boot_selected, 0.025)),
            "selective_macro_family_mae_ci95_high": float(np.quantile(boot_selected, 0.975)),
            "random_matched_macro_family_mae_mean": random_macro_mean,
            "random_matched_macro_family_mae_ci95_low": float(np.quantile(random_macro_mae, 0.025)),
            "random_matched_macro_family_mae_ci95_high": float(np.quantile(random_macro_mae, 0.975)),
            "random_minus_selective_mae": improvement,
            "random_minus_selective_family_bootstrap_ci95_low": float(np.quantile(boot_delta, 0.025)),
            "random_minus_selective_family_bootstrap_ci95_high": float(np.quantile(boot_delta, 0.975)),
            "randomization_p_random_at_or_better": randomization_p,
            "uncertainty_rankable": bool(rankable),
        })

    table = pd.DataFrame(rows)
    x = table["realized_coverage"].to_numpy(dtype=float)
    selective_y = table["selective_macro_family_mae"].to_numpy(dtype=float)
    random_y = table["random_matched_macro_family_mae_mean"].to_numpy(dtype=float)
    # Sort by realized coverage; quota rounding can cause ties between points.
    order = np.argsort(x, kind="mergesort")
    selective_aurc = float(np.trapezoid(selective_y[order], x[order]))
    random_aurc = float(np.trapezoid(random_y[order], x[order]))

    summary = {
        "status": (
            "uncertainty_not_sufficiently_variable_to_rank"
            if not rankable else
            "exploratory_fewer_than_10_heldout_families"
            if len(groups) < 10 else
            "family_matched_selective_analysis_complete"
        ),
        "validation_strength": (
            "exploratory_lt_10_families"
            if len(groups) < 10 else
            "family_count_at_least_10_still_requires_external_replication"
        ),
        "n_oof_predictions": n,
        "n_heldout_families": len(groups),
        "split_unit": "family",
        "selection_score": "nested_outer-fold calibrated expected absolute error",
        "random_control": "within-family random retention with identical per-family quotas",
        "risk_metric": "macro-family mean absolute error (kcal/mol)",
        "family_bootstrap_repetitions": family_bootstrap_repetitions,
        "random_control_repetitions": random_repetitions,
        "selective_aurc_macro_family_mae": selective_aurc if rankable else None,
        "random_matched_aurc_macro_family_mae": random_aurc,
        "aurc_improvement_vs_random": (
            random_aurc - selective_aurc if rankable else None
        ),
        "interpretation_limit": (
            "This tests selective prediction on family-held-out SKEMPI-derived "
            "protein-protein mutation-cycle data. It does not validate HERV-K Env "
            "peptide binding, target-pose accuracy, or therapeutic activity."
        ),
    }
    return table, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--predictions",
        type=Path,
        default=Path("prototype4/results/p4_validation/loco_oof_predictions_with_errors.csv"),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("prototype4/results/p4_validation/selective_abstention_family_matched"),
    )
    parser.add_argument("--random-repetitions", type=int, default=1000)
    parser.add_argument("--family-bootstrap-repetitions", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=20_261_103)
    args = parser.parse_args()

    predictions = pd.read_csv(args.predictions)
    table, summary = analyze_selective_abstention(
        predictions,
        random_repetitions=args.random_repetitions,
        family_bootstrap_repetitions=args.family_bootstrap_repetitions,
        seed=args.seed,
    )
    args.out.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out / "selective_abstention_family_matched.csv", index=False)
    (args.out / "selective_abstention_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    print(f"Wrote {args.out / 'selective_abstention_family_matched.csv'}")
    print(f"Wrote {args.out / 'selective_abstention_summary.json'}")


if __name__ == "__main__":
    main()
