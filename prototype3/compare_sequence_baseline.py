#!/usr/bin/env python3
"""Compare fixed out-of-fold predictions using a target-cluster bootstrap."""

from pathlib import Path
from datetime import datetime, timezone
import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "prototype3/results/sequence_only_mutation_baseline"
runs = sorted(
    p for p in BASE.glob("run_*")
    if (p / "out_of_fold_predictions.csv").exists()
)
if not runs:
    raise SystemExit("ERROR: no out-of-fold prediction file found.")

run_dir = runs[-1]
df = pd.read_csv(run_dir / "out_of_fold_predictions.csv")
required = {
    "target_group_id", "delta_pKd",
    "ridge_prediction", "train_mean_prediction",
}
missing = required - set(df.columns)
if missing:
    raise SystemExit(f"ERROR: missing columns: {sorted(missing)}")

def rho(y, pred):
    if len(y) < 3 or np.unique(y).size < 2 or np.unique(pred).size < 2:
        return None
    value = float(spearmanr(y, pred).statistic)
    return value if np.isfinite(value) else None

def direction_score(y, pred):
    informative = np.abs(y) >= 0.3
    if not informative.any():
        return None
    y_i, p_i = y[informative], pred[informative]
    # Neutral predictions get half credit; they do not count as wrong directions.
    correct = (np.sign(y_i) == np.sign(p_i)).astype(float)
    correct[np.abs(p_i) <= 1e-10] = 0.5
    return float(correct.mean())

def metrics(frame):
    y = frame["delta_pKd"].to_numpy(float)
    ridge = frame["ridge_prediction"].to_numpy(float)
    mean = frame["train_mean_prediction"].to_numpy(float)

    mae_ridge = float(np.mean(np.abs(y - ridge)))
    mae_mean = float(np.mean(np.abs(y - mean)))
    rmse_ridge = float(np.sqrt(np.mean((y - ridge) ** 2)))
    rmse_mean = float(np.sqrt(np.mean((y - mean) ** 2)))

    rho_ridge = rho(y, ridge)
    rho_mean = rho(y, mean)
    ds_ridge = direction_score(y, ridge)
    ds_mean = direction_score(y, mean)

    return {
        "n_pairs": int(len(frame)),
        "n_targets": int(frame["target_group_id"].nunique()),
        "ridge_mae": mae_ridge,
        "fold_mean_mae": mae_mean,
        "mae_delta_ridge_minus_fold_mean": mae_ridge - mae_mean,
        "ridge_rmse": rmse_ridge,
        "fold_mean_rmse": rmse_mean,
        "rmse_delta_ridge_minus_fold_mean": rmse_ridge - rmse_mean,
        "ridge_spearman": rho_ridge,
        "fold_mean_spearman": rho_mean,
        "spearman_delta_ridge_minus_fold_mean": (
            rho_ridge - rho_mean
            if rho_ridge is not None and rho_mean is not None else None
        ),
        "ridge_direction_score_ties_half": ds_ridge,
        "fold_mean_direction_score_ties_half": ds_mean,
        "direction_delta_ridge_minus_fold_mean": (
            ds_ridge - ds_mean
            if ds_ridge is not None and ds_mean is not None else None
        ),
    }

observed = metrics(df)
rng = np.random.default_rng(20261027)
by_target = {k: v for k, v in df.groupby("target_group_id")}
target_ids = np.array(list(by_target))
replicates = 2000

delta_keys = [
    "mae_delta_ridge_minus_fold_mean",
    "rmse_delta_ridge_minus_fold_mean",
    "spearman_delta_ridge_minus_fold_mean",
    "direction_delta_ridge_minus_fold_mean",
]
boot = {key: [] for key in delta_keys}

for _ in range(replicates):
    sampled = rng.choice(target_ids, size=len(target_ids), replace=True)
    sample_df = pd.concat([by_target[t] for t in sampled], ignore_index=True)
    result = metrics(sample_df)
    for key in delta_keys:
        value = result[key]
        if value is not None and np.isfinite(value):
            boot[key].append(float(value))

def ci(values):
    if len(values) < 100:
        return None
    return [
        float(np.quantile(values, 0.025)),
        float(np.quantile(values, 0.975)),
    ]

# Negative deltas favor Ridge for MAE/RMSE. Positive deltas favor Ridge for
# Spearman and direction score.
favorable_sign = {
    "mae_delta_ridge_minus_fold_mean": -1,
    "rmse_delta_ridge_minus_fold_mean": -1,
    "spearman_delta_ridge_minus_fold_mean": 1,
    "direction_delta_ridge_minus_fold_mean": 1,
}
comparison = {}
for key in delta_keys:
    values = boot[key]
    comparison[key] = {
        "observed_difference": observed[key],
        "target_cluster_bootstrap_95pct_ci": ci(values),
        "valid_replicates": len(values),
        "fraction_bootstrap_draws_favoring_ridge": (
            float(np.mean(np.asarray(values) * favorable_sign[key] > 0))
            if values else None
        ),
    }

result = {
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "prediction_source": str(run_dir),
    "bootstrap": {
        "method": "2,000 target-cluster resamples of fixed out-of-fold predictions",
        "seed": 20261027,
        "note": (
            "Intervals condition on the saved out-of-fold predictions; "
            "they do not include variation from refitting the model."
        ),
    },
    "observed_metrics": observed,
    "paired_differences": comparison,
    "interpretation_rules": {
        "mae_rmse": "Negative Ridge-minus-mean difference favors Ridge.",
        "spearman_direction": "Positive Ridge-minus-mean difference favors Ridge.",
        "uncertainty": "If the 95% interval for a difference includes zero, "
                       "the direction of that difference remains uncertain.",
    },
}

out_json = run_dir / "paired_model_comparison.json"
out_md = run_dir / "paired_model_comparison.md"
out_json.write_text(json.dumps(result, indent=2, allow_nan=False))

lines = [
    "# Paired sequence-model comparison",
    "",
    f"- Predictions: `{run_dir}`",
    "- Resampling unit: target group, not individual pair.",
    "- Bootstrap: 2,000 resamples; intervals condition on fixed out-of-fold predictions.",
    "",
    "## Observed performance differences",
    "",
    "| Difference (Ridge minus fold-mean baseline) | Estimate | 95% target-bootstrap interval |",
    "|---|---:|---:|",
]
for key in delta_keys:
    item = comparison[key]
    estimate = item["observed_difference"]
    interval = item["target_cluster_bootstrap_95pct_ci"]
    est_text = "NA" if estimate is None else f"{estimate:.4f}"
    ci_text = "NA" if interval is None else f"[{interval[0]:.4f}, {interval[1]:.4f}]"
    lines.append(f"| {key} | {est_text} | {ci_text} |")

lines.extend([
    "",
    "## Reading the differences",
    "",
    "- Negative MAE or RMSE differences favor Ridge.",
    "- Positive Spearman or direction-score differences favor Ridge.",
    "- If an interval contains zero, this comparison does not clearly establish a difference.",
    "- This compares a sequence-only comparator with a fold-mean baseline. It does not validate `c_score` or the structural score.",
    "",
])
out_md.write_text("\n".join(lines))
print("PREDICTION SOURCE:", run_dir)
print(out_md.read_text())
print("JSON:", out_json)
