#!/usr/bin/env python3
"""Exploratory target-held-out baseline for experimental single-mutant effects.

This is a sequence-only comparator. It does NOT evaluate c_score or the
target-aware structural score. It never pads/truncates sequences to fit c_score.
"""

from __future__ import annotations

import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "prototype3/results/ppikb_kd_benchmark/filtered_kd_measurements.csv"
OUTPUT_ROOT = ROOT / "prototype3/results/sequence_only_mutation_baseline"

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {aa: i for i, aa in enumerate(AA)}

# Fixed, documented properties for a small sequence-only comparison.
# These are approximate features, not physical binding energies.
HYDROPATHY = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
    "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
    "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
    "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}
CHARGE = {
    aa: (-1.0 if aa in "DE" else 1.0 if aa in "KR" else 0.1 if aa == "H" else 0.0)
    for aa in AA
}
AROMATIC = {aa: float(aa in "FWY") for aa in AA}
SULFUR = {aa: float(aa in "CM") for aa in AA}
SPECIAL = {aa: float(aa in "GP") for aa in AA}


def norm_col(value: object) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def find_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    lookup = {norm_col(col): col for col in df.columns}
    for candidate in candidates:
        if norm_col(candidate) in lookup:
            return lookup[norm_col(candidate)]
    return None


def clean_seq(value: object) -> str:
    if pd.isna(value):
        return ""
    return re.sub(r"\s+", "", str(value)).upper()


def stable_id(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()[:16]


def safe_spearman(y: np.ndarray, pred: np.ndarray) -> float | None:
    if len(y) < 3 or np.unique(y).size < 2 or np.unique(pred).size < 2:
        return None
    value = float(spearmanr(y, pred).statistic)
    return value if np.isfinite(value) else None


def metrics_for(y: np.ndarray, pred: np.ndarray, informative_cutoff: float = 0.3) -> dict:
    y = np.asarray(y, dtype=float)
    pred = np.asarray(pred, dtype=float)
    informative = np.abs(y) >= informative_cutoff

    result = {
        "n": int(len(y)),
        "mae": float(np.mean(np.abs(y - pred))) if len(y) else None,
        "rmse": float(np.sqrt(np.mean((y - pred) ** 2))) if len(y) else None,
        "spearman_rho": safe_spearman(y, pred),
        "informative_cutoff_abs_delta_pKd": float(informative_cutoff),
        "n_informative": int(informative.sum()),
    }

    if informative.any():
        observed = y[informative]
        predicted = pred[informative]
        has_direction = np.abs(predicted) > 1e-10

        result["informative_direction_coverage"] = float(has_direction.mean())
        result["informative_direction_accuracy"] = (
            float(np.mean(
                np.sign(observed[has_direction])
                == np.sign(predicted[has_direction])
            ))
            if has_direction.any() else None
        )
        result["informative_direction_accuracy_ties_half"] = float(np.mean(
            np.where(
                has_direction,
                (np.sign(observed) == np.sign(predicted)).astype(float),
                0.5,
            )
        ))
    else:
        result["informative_direction_coverage"] = None
        result["informative_direction_accuracy"] = None
        result["informative_direction_accuracy_ties_half"] = None

    return result


def encode_pair(seq_a: str, seq_b: str, max_len: int) -> np.ndarray:
    """Encode a directed single-substitution pair as fixed sequence features."""
    if len(seq_a) != len(seq_b):
        raise ValueError(f"Sequences have different lengths: {seq_a} vs {seq_b}")

    positions = [i for i, (a, b) in enumerate(zip(seq_a, seq_b)) if a != b]
    if len(positions) != 1:
        raise ValueError(
            f"Expected exactly one substitution, found {len(positions)}: "
            f"{seq_a} -> {seq_b}"
        )
    if len(seq_a) > max_len:
        raise ValueError(f"Sequence length {len(seq_a)} exceeds max_len={max_len}")

    pos = positions[0]
    aa_from, aa_to = seq_a[pos], seq_b[pos]
    if aa_from not in AA_INDEX or aa_to not in AA_INDEX:
        raise ValueError(f"Unsupported amino acid in pair: {seq_a} -> {seq_b}")

    x = np.zeros(400 + 2 * max_len + 5, dtype=float)
    x[AA_INDEX[aa_from] * 20 + AA_INDEX[aa_to]] = 1.0
    x[400 + pos] = 1.0
    x[400 + max_len + len(seq_a) - 1] = 1.0

    deltas = [
        HYDROPATHY[aa_to] - HYDROPATHY[aa_from],
        CHARGE[aa_to] - CHARGE[aa_from],
        AROMATIC[aa_to] - AROMATIC[aa_from],
        SULFUR[aa_to] - SULFUR[aa_from],
        SPECIAL[aa_to] - SPECIAL[aa_from],
    ]
    x[400 + 2 * max_len:] = deltas
    return x


def cluster_bootstrap(
    oof: pd.DataFrame,
    n_boot: int = 2000,
    seed: int = 20261026,
) -> dict:
    """Resample entire target groups, not individual correlated pairs."""
    rng = np.random.default_rng(seed)
    target_ids = oof["target_group_id"].drop_duplicates().to_numpy()
    by_target = {key: frame for key, frame in oof.groupby("target_group_id")}
    rho_values = []
    direction_values = []

    for _ in range(n_boot):
        sampled_targets = rng.choice(target_ids, size=len(target_ids), replace=True)
        sampled = pd.concat([by_target[t] for t in sampled_targets], ignore_index=True)
        y = sampled["delta_pKd"].to_numpy(float)
        pred = sampled["ridge_prediction"].to_numpy(float)

        rho = safe_spearman(y, pred)
        if rho is not None:
            rho_values.append(rho)

        informative = np.abs(y) >= 0.3
        if informative.any():
            direction_values.append(float(np.mean(
                np.sign(y[informative]) == np.sign(pred[informative])
            )))

    def ci(values):
        if len(values) < 100:
            return None
        return [
            float(np.quantile(values, 0.025)),
            float(np.quantile(values, 0.975)),
        ]

    return {
        "method": "target-cluster bootstrap; 2000 replicates",
        "seed": seed,
        "spearman_rho_95pct_ci": ci(rho_values),
        "informative_direction_accuracy_95pct_ci": ci(direction_values),
        "valid_rho_replicates": int(len(rho_values)),
        "valid_direction_replicates": int(len(direction_values)),
    }


def main() -> int:
    if not INPUT.exists():
        print(f"ERROR: Input file not found: {INPUT}", file=sys.stderr)
        return 2

    raw = pd.read_csv(INPUT, low_memory=False)
    peptide_col = find_col(raw, [
        "peptide_sequence", "peptide_seq", "peptide", "sequence",
    ])
    target_col = find_col(raw, [
        "target_sequence", "protein_sequence", "target_protein_sequence",
        "protein", "target", "protein_id", "target_id", "uniprot_id",
        "target_group_id",
    ])
    affinity_col = find_col(raw, [
        "p_affinity", "p_affinity_value", "pKd", "pKD", "p_kd", "paffinity",
    ])

    if not all([peptide_col, target_col, affinity_col]):
        print("ERROR: Could not identify the required columns.", file=sys.stderr)
        print(f"Detected peptide={peptide_col}, target={target_col}, "
              f"p-affinity={affinity_col}", file=sys.stderr)
        print("Actual columns:", list(raw.columns), file=sys.stderr)
        return 2

    # Preserve experimental context. Splitting more conservatively is preferable
    # to combining measurements from unlike assays without warning.
    possible_context_cols = [
        "assay", "assay_id", "assay_name", "assay_method",
        "assay_type", "method", "affinity_type", "affinity_unit", "source",
    ]
    context_cols = []
    for name in possible_context_cols:
        found = find_col(raw, [name])
        if found and found not in context_cols:
            context_cols.append(found)

    if not context_cols:
        print("ERROR: No assay/context field found. Refusing to combine all assays.",
              file=sys.stderr)
        print("Actual columns:", list(raw.columns), file=sys.stderr)
        return 2

    df = pd.DataFrame({
        "sequence": raw[peptide_col].map(clean_seq),
        "target_value": raw[target_col].fillna("").astype(str).str.strip(),
        "pKd": pd.to_numeric(raw[affinity_col], errors="coerce"),
    })
    for col in context_cols:
        df[f"context::{col}"] = raw[col].fillna("MISSING").astype(str).str.strip()

    context_feature_cols = [c for c in df.columns if c.startswith("context::")]
    df["context_key"] = df[context_feature_cols].astype(str).agg(" | ".join, axis=1)

    valid_seq = df["sequence"].str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]{5,50}")
    valid = (
        valid_seq
        & df["target_value"].ne("")
        & df["pKd"].notna()
        & np.isfinite(df["pKd"])
    )
    n_raw = int(len(df))
    df = df.loc[valid].copy()

    if df.empty:
        print("ERROR: No valid rows after sequence/target/affinity filtering.",
              file=sys.stderr)
        return 2

    df["target_group_id"] = df["target_value"].map(stable_id)
    df["seq_len"] = df["sequence"].str.len()

    # Repeated measurements of the same peptide in the same target/context are
    # collapsed to a median before pairs are generated.
    unique_measurements = (
        df.groupby(
            ["target_group_id", "target_value", "context_key", "seq_len", "sequence"],
            as_index=False,
        )["pKd"].median()
    )

    pair_rows = []
    grouping = unique_measurements.groupby(
        ["target_group_id", "target_value", "context_key", "seq_len"],
        sort=False,
    )

    for (target_id, target_value, context_key, seq_len), frame in grouping:
        seq_to_y = dict(zip(frame["sequence"], frame["pKd"]))
        sequences = sorted(seq_to_y)
        for seq_a in sequences:
            chars = list(seq_a)
            for pos, aa_from in enumerate(chars):
                for aa_to in AA:
                    if aa_to == aa_from:
                        continue
                    mutated = chars.copy()
                    mutated[pos] = aa_to
                    seq_b = "".join(mutated)

                    # Lexicographic orientation makes every pair occur exactly once.
                    if seq_b not in seq_to_y or seq_a >= seq_b:
                        continue

                    y_a = float(seq_to_y[seq_a])
                    y_b = float(seq_to_y[seq_b])
                    pair_rows.append({
                        "target_group_id": target_id,
                        "target_value": target_value,
                        "context_key": context_key,
                        "seq_len": int(seq_len),
                        "sequence_a": seq_a,
                        "sequence_b": seq_b,
                        "pKd_a": y_a,
                        "pKd_b": y_b,
                        "delta_pKd": y_b - y_a,
                        "abs_delta_pKd": abs(y_b - y_a),
                    })

    pairs = pd.DataFrame(pair_rows)
    if pairs.empty:
        print("ERROR: No same-target, same-context, same-length "
              "single-substitution pairs found.", file=sys.stderr)
        return 2

    # A target is the split unit. No target appears in both train and test
    # within a fold, even if it has multiple assays or many variants.
    target_count = pairs["target_group_id"].nunique()
    if target_count < 5:
        print(f"ERROR: Only {target_count} target groups; five-fold "
              "target-held-out evaluation is not supported.", file=sys.stderr)
        return 2

    max_len = int(pairs["seq_len"].max())
    X = np.vstack([
        encode_pair(a, b, max_len)
        for a, b in zip(pairs["sequence_a"], pairs["sequence_b"])
    ])
    y = pairs["delta_pKd"].to_numpy(float)
    groups = pairs["target_group_id"].to_numpy()

    # Fixed alpha chosen before evaluation; no tuning on held-out folds.
    cv = GroupKFold(n_splits=5)
    oof_rows = []
    fold_sizes = []

    for fold, (train_idx, test_idx) in enumerate(cv.split(X, y, groups), start=1):
        train_groups = set(groups[train_idx])
        test_groups = set(groups[test_idx])
        assert train_groups.isdisjoint(test_groups), "Target leakage between folds."

        model = Ridge(alpha=10.0)
        model.fit(X[train_idx], y[train_idx])
        ridge_pred = model.predict(X[test_idx])
        train_mean = float(np.mean(y[train_idx]))

        fold_sizes.append({
            "fold": fold,
            "n_train_pairs": int(len(train_idx)),
            "n_test_pairs": int(len(test_idx)),
            "n_train_targets": int(len(train_groups)),
            "n_test_targets": int(len(test_groups)),
            "target_overlap": 0,
        })

        for j, pred in zip(test_idx, ridge_pred):
            oof_rows.append({
                "pair_index": int(j),
                "fold": fold,
                "target_group_id": pairs.iloc[j]["target_group_id"],
                "context_key": pairs.iloc[j]["context_key"],
                "sequence_a": pairs.iloc[j]["sequence_a"],
                "sequence_b": pairs.iloc[j]["sequence_b"],
                "pKd_a": float(pairs.iloc[j]["pKd_a"]),
                "pKd_b": float(pairs.iloc[j]["pKd_b"]),
                "delta_pKd": float(y[j]),
                "zero_change_prediction": 0.0,
                "train_mean_prediction": train_mean,
                "ridge_prediction": float(pred),
            })

    oof = pd.DataFrame(oof_rows).sort_values("pair_index").reset_index(drop=True)
    assert len(oof) == len(pairs), "Some pairs did not receive an out-of-fold prediction."
    assert oof["pair_index"].is_unique, "An experimental pair was predicted more than once."

    y_oof = oof["delta_pKd"].to_numpy(float)
    metrics = {
        "experiment": "Exploratory sequence-only mutation-effect comparator",
        "input_csv": str(INPUT),
        "detected_columns": {
            "peptide": peptide_col,
            "target_grouping": target_col,
            "p_affinity": affinity_col,
            "context": context_cols,
        },
        "row_counts": {
            "input_rows": n_raw,
            "valid_rows": int(len(df)),
            "unique_measurements_after_median_aggregation": int(len(unique_measurements)),
            "single_substitution_pairs": int(len(pairs)),
            "target_groups": int(target_count),
            "context_groups": int(pairs["context_key"].nunique()),
            "informative_pairs_abs_delta_pKd_ge_0_3": int((pairs["abs_delta_pKd"] >= 0.3).sum()),
        },
        "split": {
            "method": "5-fold GroupKFold, grouped by target",
            "all_targets_disjoint_between_train_and_test_each_fold": True,
            "folds": fold_sizes,
        },
        "model": {
            "type": "Ridge regression",
            "features": "directed amino-acid substitution + position + length + simple property deltas",
            "alpha": 10.0,
            "alpha_tuned_on_test_set": False,
        },
        "out_of_fold_metrics": {
            "zero_change_baseline": metrics_for(y_oof, oof["zero_change_prediction"].to_numpy()),
            "training_fold_mean_baseline": metrics_for(y_oof, oof["train_mean_prediction"].to_numpy()),
            "sequence_only_ridge": metrics_for(y_oof, oof["ridge_prediction"].to_numpy()),
        },
        "cluster_bootstrap": cluster_bootstrap(oof),
        "interpretation_limits": [
            "This is a sequence-only comparator, not the project c_score.",
            "It does not validate the target-aware structural score.",
            "Affinity differences are compared only within the same target and detected assay context.",
            "Results are retrospective and exploratory, not prospective biological validation.",
        ],
    }

    run_id = datetime.now(timezone.utc).strftime("run_%Y%m%dT%H%M%SZ")
    outdir = OUTPUT_ROOT / run_id
    outdir.mkdir(parents=True, exist_ok=False)

    pairs.to_csv(outdir / "single_substitution_pairs.csv", index=False)
    oof.to_csv(outdir / "out_of_fold_predictions.csv", index=False)
    with (outdir / "metrics.json").open("w") as f:
        json.dump(metrics, f, indent=2, allow_nan=False)

    def show_metric(label: str, item: dict) -> str:
        rho = item.get("spearman_rho")
        coverage = item.get("informative_direction_coverage")
        accuracy = item.get("informative_direction_accuracy")

        rho_text = "NA" if rho is None else f"{rho:.4f}"
        coverage_text = "NA" if coverage is None else f"{coverage:.3f}"
        accuracy_text = "NA" if accuracy is None else f"{accuracy:.3f}"

        return (
            f"| {label} | {item['n']} | {item['mae']:.4f} | "
            f"{item['rmse']:.4f} | {rho_text} | "
            f"{coverage_text} | {accuracy_text} |"
        )

    sections = [
        "# Sequence-only mutation-effect benchmark",
        "",
        f"- Run: `{run_id}`",
        f"- Input: `{INPUT}`",
        f"- Peptide column: `{peptide_col}`",
        f"- Target grouping column: `{target_col}`",
        f"- Affinity column: `{affinity_col}`",
        f"- Context columns: `{', '.join(context_cols)}`",
        "- Split: five-fold GroupKFold grouped by target; no target overlap across train/test folds.",
        "- Model: fixed Ridge regression (`alpha=10.0`); no hyperparameter selection on the held-out data.",
        "",
        "## Data counts",
        "",
        f"- Input rows: {n_raw}",
        f"- Valid rows after filtering: {len(df)}",
        f"- Unique target/context/peptide measurements after median aggregation: {len(unique_measurements)}",
        f"- Matched single-substitution pairs: {len(pairs)}",
        f"- Distinct target groups: {target_count}",
        f"- Pairs with absolute experimental change >= 0.3 pKd: {(pairs['abs_delta_pKd'] >= 0.3).sum()}",
        "",
        "## Out-of-fold results",
        "",
        "| Method | Pairs | MAE (pKd) | RMSE (pKd) | Spearman rho | Direction coverage | Accuracy among nonzero directions |",
        "|---|---:|---:|---:|---:|---:|",
        show_metric("Predict no change", metrics["out_of_fold_metrics"]["zero_change_baseline"]),
        show_metric("Training-fold mean", metrics["out_of_fold_metrics"]["training_fold_mean_baseline"]),
        show_metric("Sequence-only Ridge", metrics["out_of_fold_metrics"]["sequence_only_ridge"]),
        "",
        "## Cluster-bootstrap uncertainty",
        "",
        f"- Spearman rho 95% CI: `{metrics['cluster_bootstrap']['spearman_rho_95pct_ci']}`",
        f"- Informative direction accuracy 95% CI: `{metrics['cluster_bootstrap']['informative_direction_accuracy_95pct_ci']}`",
        "",
        "## Interpretation",
        "",
        "This benchmark tests a simple sequence-only comparator on retrospective experimental mutation pairs. "
        "It does not use or validate the native `c_score`, and it does not validate the target-aware structural score. "
        "The target-held-out split reduces target leakage, but assay heterogeneity, data selection, and limited "
        "generalization remain relevant. Do not claim real-world binding or efficacy from this benchmark.",
        "",
    ]
    (outdir / "REPORT.md").write_text("\n".join(sections))

    print(f"\nWROTE: {outdir}")
    print(f"Pairs: {len(pairs)} | Target groups: {target_count}")
    print("\nOut-of-fold metrics:")
    for name, values in metrics["out_of_fold_metrics"].items():
        print(name, json.dumps(values, indent=2))
    print("\nCluster bootstrap:")
    print(json.dumps(metrics["cluster_bootstrap"], indent=2))
    print("\nReview REPORT.md and metrics.json before using any result in the logbook.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
