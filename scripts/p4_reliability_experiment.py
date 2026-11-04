#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sqlite3
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import (
    mean_absolute_error,
    mean_squared_error,
    r2_score,
)
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


TARGET_ALIASES = [
    "delta_epsilon",
    "delta_eps",
    "delta_delta_g",
    "experimental_ddg",
    "measured_ddg",
    "ddg",
    "y_true",
    "target_value",
]

GROUP_ALIASES = [
    "interaction_family",
    "family_id",
    "family",
    "protein_family",
    "complex_family",
    "pdb_family",
]

EXCLUDE_TOKENS = (
    "ddg", "delta_eps", "delta_delta_g",
    "experimental", "measured", "observed",
    "affinity", "binding_energy", "outcome",
    "label", "prediction", "predicted", "residual",
    "error", "fold", "split", "family", "pdb",
    "complex_id", "sequence", "mutation",
    "residue", "chain", "seed", "replicate",
    "run_id", "source",
)


def normalized(name):
    return str(name).strip().lower().replace("-", "_").replace(" ", "_")


def resolve_column(columns, requested, aliases):
    normalized_columns = {normalized(c): c for c in columns}

    if requested:
        key = normalized(requested)
        if key not in normalized_columns:
            raise ValueError(
                f"Requested column {requested!r} not found. "
                f"Available columns: {list(columns)}"
            )
        return normalized_columns[key]

    for alias in aliases:
        if alias in normalized_columns:
            return normalized_columns[alias]

    for column in columns:
        key = normalized(column)
        if any(alias in key for alias in aliases):
            return column

    return None


def load_database(path, requested_table=None):
    conn = sqlite3.connect(path)
    tables = pd.read_sql_query(
        "SELECT name FROM sqlite_master WHERE type='table'",
        conn,
    )["name"].tolist()

    candidates = []

    for table in tables:
        if requested_table and table != requested_table:
            continue

        quoted = '"' + table.replace('"', '""') + '"'
        df = pd.read_sql_query(f"SELECT * FROM {quoted}", conn)

        if len(df) < 20:
            continue

        candidates.append((table, df))

    conn.close()

    if not candidates:
        raise RuntimeError(
            f"No suitable tables found in {path}. Tables: {tables}"
        )

    return candidates


def discover_sources(root):
    found = []

    for path in root.rglob("*"):
        if not path.is_file():
            continue

        parts = {part.lower() for part in path.parts}
        if parts & {"backup", "backups", "patch_backups", "__pycache__"}:
            continue

        suffix = path.suffix.lower()
        if suffix in {".db", ".sqlite", ".sqlite3"}:
            found.append(path)

    return sorted(found)


def find_dataset(args):
    if args.input:
        path = Path(args.input)
        if not path.exists():
            raise FileNotFoundError(path)

        if path.suffix.lower() in {".db", ".sqlite", ".sqlite3"}:
            return load_database(path, args.table)

        return [(path.name, pd.read_csv(path))]

    sources = discover_sources(Path("prototype4"))
    candidates = []

    for path in sources:
        try:
            tables = load_database(path, args.table)
        except Exception:
            continue

        for table, df in tables:
            target = resolve_column(df.columns, args.target, TARGET_ALIASES)
            group = resolve_column(df.columns, args.group, GROUP_ALIASES)

            if target is None or group is None:
                continue

            score = (
                100 if "cycle_features" in table.lower() else 0
            ) + min(len(df), 1000) / 1000

            candidates.append((score, path, table, df))

    if not candidates:
        raise RuntimeError(
            "Could not automatically locate a dataset containing an "
            "experimental target and an interaction-family column. "
            "Run the schema-inspection command below and pass --input, "
            "--target, and --group explicitly."
        )

    candidates.sort(key=lambda x: x[0], reverse=True)
    _, path, table, df = candidates[0]

    print(f"Selected database: {path}")
    print(f"Selected table: {table}")

    return [(table, df)]


def make_model(alpha):
    return make_pipeline(
        SimpleImputer(strategy="median"),
        StandardScaler(),
        Ridge(alpha=alpha),
    )


def safe_spearman(a, b):
    result = spearmanr(a, b)
    if not np.isfinite(result.statistic):
        return None
    return float(result.statistic)


def group_bootstrap_ci(frame, metric, repetitions=2000, seed=731):
    rng = np.random.default_rng(seed)
    groups = frame["group"].unique()

    if len(groups) < 3:
        return None

    group_values = {
        group: metric(frame.loc[frame["group"] == group])
        for group in groups
    }

    values = np.array(list(group_values.values()), dtype=float)
    values = values[np.isfinite(values)]

    if len(values) < 3:
        return None

    samples = rng.choice(
        values,
        size=(repetitions, len(values)),
        replace=True,
    ).mean(axis=1)

    return [
        float(np.quantile(samples, 0.025)),
        float(np.quantile(samples, 0.975)),
    ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=None)
    parser.add_argument("--table", default=None)
    parser.add_argument("--target", default=None)
    parser.add_argument("--group", default=None)
    parser.add_argument("--alpha", type=float, default=1.0)
    parser.add_argument("--bootstraps", type=int, default=60)
    parser.add_argument("--seed", type=int, default=731)
    parser.add_argument(
        "--out",
        default="prototype4/results/reliability_experiment_20261103",
    )
    args = parser.parse_args()

    if args.bootstraps < 10:
        raise ValueError("--bootstraps must be at least 10")

    tables = find_dataset(args)

    suitable = []

    for table_name, df in tables:
        target_col = resolve_column(
            df.columns, args.target, TARGET_ALIASES
        )
        group_col = resolve_column(
            df.columns, args.group, GROUP_ALIASES
        )

        if target_col is not None and group_col is not None:
            suitable.append((table_name, df, target_col, group_col))

    if not suitable:
        for table_name, df in tables:
            print(f"\nTABLE: {table_name}")
            print("COLUMNS:", list(df.columns))
        raise RuntimeError(
            "No table contained both the requested outcome and group columns."
        )

    table_name, df, target_col, group_col = max(
        suitable, key=lambda item: len(item[1])
    )

    df = df.copy()
    df["_target"] = pd.to_numeric(df[target_col], errors="coerce")
    df["_group"] = df[group_col].astype(str)

    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["_target", "_group"])

    feature_cols = []

    for column in df.columns:
        if column in {target_col, group_col, "_target", "_group"}:
            continue

        key = normalized(column)

        if any(token in key for token in EXCLUDE_TOKENS):
            continue

        numeric = pd.to_numeric(df[column], errors="coerce")

        if numeric.notna().mean() < 0.80:
            continue

        if numeric.nunique(dropna=True) < 3:
            continue

        feature_cols.append(column)
        df[column] = numeric

    if len(feature_cols) < 2:
        raise RuntimeError(
            "Fewer than two usable numeric predictors were found. "
            f"Columns: {list(df.columns)}"
        )

    df = df.dropna(subset=["_target", "_group"])
    df = df.reset_index(drop=True)

    groups = df["_group"].to_numpy()
    y = df["_target"].to_numpy(dtype=float)
    X = df[feature_cols].copy()

    unique_groups = np.unique(groups)

    if len(unique_groups) < 4:
        raise RuntimeError(
            f"Only {len(unique_groups)} groups. "
            "At least four are required for family-held-out validation."
        )

    print("\nFRESH FAMILY-HELD-OUT EXPERIMENT")
    print(f"Table: {table_name}")
    print(f"Rows: {len(df)}")
    print(f"Families/groups: {len(unique_groups)}")
    print(f"Target: {target_col}")
    print(f"Grouping column: {group_col}")
    print(f"Features: {len(feature_cols)}")
    print(f"Bootstrap models per fold: {args.bootstraps}")

    rng = np.random.default_rng(args.seed)
    records = []

    for fold, held_group in enumerate(unique_groups, start=1):
        test_mask = groups == held_group
        train_mask = ~test_mask

        X_train = X.loc[train_mask]
        X_test = X.loc[test_mask]
        y_train = y[train_mask]
        y_test = y[test_mask]
        train_groups = groups[train_mask]

        if len(y_train) < 10:
            warnings.warn(f"Skipping small training fold: {held_group}")
            continue

        model = make_model(args.alpha)
        model.fit(X_train, y_train)
        ridge_pred = model.predict(X_test)

        group_names = np.unique(train_groups)
        bootstrap_predictions = []

        for _ in range(args.bootstraps):
            sampled_groups = rng.choice(
                group_names,
                size=len(group_names),
                replace=True,
            )

            bootstrap_indices = np.concatenate([
                np.flatnonzero(train_groups == group)
                for group in sampled_groups
            ])

            boot_model = make_model(args.alpha)
            boot_model.fit(
                X_train.iloc[bootstrap_indices],
                y_train[bootstrap_indices],
            )

            bootstrap_predictions.append(
                boot_model.predict(X_test)
            )

        bootstrap_predictions = np.asarray(bootstrap_predictions)

        ensemble_mean = bootstrap_predictions.mean(axis=0)
        ensemble_std = bootstrap_predictions.std(axis=0, ddof=1)

        fold_baseline = float(np.mean(y_train))

        for i, original_index in enumerate(np.flatnonzero(test_mask)):
            records.append({
                "row_index": int(original_index),
                "group": str(held_group),
                "y_true": float(y_test[i]),
                "ridge_prediction": float(ridge_pred[i]),
                "ensemble_prediction": float(ensemble_mean[i]),
                "ensemble_disagreement": float(ensemble_std[i]),
                "fold_mean_prediction": fold_baseline,
                "absolute_error": float(
                    abs(y_test[i] - ensemble_mean[i])
                ),
            })

        print(
            f"Fold {fold}/{len(unique_groups)} "
            f"group={held_group} n_test={len(y_test)} complete"
        )

    pred = pd.DataFrame(records)

    if len(pred) < 20:
        raise RuntimeError(
            f"Only {len(pred)} out-of-family predictions were produced."
        )

    y_true = pred["y_true"].to_numpy()
    ensemble_pred = pred["ensemble_prediction"].to_numpy()
    fold_mean_pred = pred["fold_mean_prediction"].to_numpy()
    errors = pred["absolute_error"].to_numpy()
    uncertainty = pred["ensemble_disagreement"].to_numpy()

    def metrics(actual, predicted):
        return {
            "n": int(len(actual)),
            "mae": float(mean_absolute_error(actual, predicted)),
            "rmse": float(
                np.sqrt(mean_squared_error(actual, predicted))
            ),
            "r2": (
                float(r2_score(actual, predicted))
                if len(actual) >= 2 else None
            ),
            "spearman": safe_spearman(actual, predicted),
        }

    results = {
        "experiment": "fresh_family_held_out_reliability_validation",
        "source_table": table_name,
        "target_column": target_col,
        "group_column": group_col,
        "n_rows": int(len(df)),
        "n_families": int(len(unique_groups)),
        "n_features": int(len(feature_cols)),
        "feature_columns": feature_cols,
        "ridge_alpha": args.alpha,
        "bootstrap_models_per_fold": args.bootstraps,
        "seed": args.seed,
        "validation": "leave_one_interaction_family_out",
        "metrics": {
            "bootstrap_ridge": metrics(y_true, ensemble_pred),
            "ridge_single_fit": metrics(
                y_true, pred["ridge_prediction"].to_numpy()
            ),
            "fold_mean_baseline": metrics(y_true, fold_mean_pred),
            "zero_baseline": metrics(
                y_true, np.zeros_like(y_true)
            ),
            "direction_accuracy": float(np.mean(
                np.sign(y_true) == np.sign(ensemble_pred)
            )),
            "uncertainty_vs_absolute_error_spearman": safe_spearman(
                uncertainty, errors
            ),
            "mae_group_bootstrap_95pct_ci": group_bootstrap_ci(
                pred,
                lambda part: mean_absolute_error(
                    part["y_true"], part["ensemble_prediction"]
                ),
            ),
        },
        "selective_prediction": {},
        "interpretation_limits": [
            "Family-held-out prediction tests generalization to unseen families.",
            "Bootstrap disagreement is an uncertainty diagnostic, not a calibrated probability.",
            "Selective MAE measures error among retained predictions.",
            "This experiment does not establish HERV-K-specific binding affinity.",
            "This experiment does not by itself test the HADDOCK repeatability gate."
        ],
    }

    order = np.argsort(uncertainty)

    for coverage in (1.0, 0.8, 0.6):
        count = max(1, int(np.ceil(len(pred) * coverage)))
        kept = pred.iloc[order[:count]]

        results["selective_prediction"][str(coverage)] = {
            "coverage": float(len(kept) / len(pred)),
            "retained_n": int(len(kept)),
            "mae": float(mean_absolute_error(
                kept["y_true"], kept["ensemble_prediction"]
            )),
            "fold_mean_baseline_mae": float(mean_absolute_error(
                kept["y_true"], kept["fold_mean_prediction"]
            )),
            "mean_uncertainty": float(
                kept["ensemble_disagreement"].mean()
            ),
        }

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    pred.to_csv(out / "oof_predictions.csv", index=False)

    (out / "summary.json").write_text(
        json.dumps(results, indent=2) + "\n"
    )

    print("\n========== FRESH RESULTS ==========")
    print(json.dumps(results["metrics"], indent=2))
    print("\n========== SELECTIVE MAE ==========")
    print(json.dumps(results["selective_prediction"], indent=2))
    print(f"\nSaved results to: {out.resolve()}")


if __name__ == "__main__":
    main()
