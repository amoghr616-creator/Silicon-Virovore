#!/usr/bin/env python3
"""Leakage-resistant, target-held-out validation of a separate sequence-only pair model.

This script does NOT evaluate the native c_score, the adaptive optimizer, or the
structural scorer. It builds experimental single-substitution pairs from the
filtered PPIKB table and evaluates a separate fixed-alpha Ridge regression on
out-of-fold affinity differences. Entire target sequences and repeated exact
peptide-pair signatures are grouped into the same CV component.

The score is a research comparator only; it is not used by Silicon Virovore's
optimizer. Model hyperparameters are fixed (Ridge alpha=10) before CV and are
not tuned on the held-out folds.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {aa: i for i, aa in enumerate(AA)}
CANONICAL_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
DEFAULT_INPUT = Path("prototype3/results/ppikb_kd_benchmark/filtered_kd_measurements.csv")
DEFAULT_OUTPUT = Path("prototype3/results/sequence_only_grouped_target_validation")
RIDGE_ALPHA = 10.0  # Fixed in advance; no hyperparameter search on held-out data.
INFORMATIVE_CUTOFF = 0.3
DIRECTION_EPSILON = 0.05


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


class UnionFind:
    def __init__(self, n: int):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


def _resolve_col(df: pd.DataFrame, aliases: Iterable[str], label: str) -> str:
    normalized = {str(c).strip().casefold(): str(c) for c in df.columns}
    for candidate in aliases:
        if candidate.casefold() in normalized:
            return normalized[candidate.casefold()]
    raise ValueError(
        f"Cannot identify {label} column. Tried aliases {list(aliases)!r}. "
        f"Available columns: {list(df.columns)!r}"
    )


def prepare_measurements(path: str | Path) -> pd.DataFrame:
    """Load and collapse canonical, exact affinity rows to target/assay/peptide."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(
            f"Filtered PPIKB input not found: {path}. Run prototype3/scripts/ppikb_kd_benchmark.py first."
        )
    df = pd.read_csv(path)
    pep_col = _resolve_col(df, ["_peptide", "Peptide_Sequence", "peptide_sequence", "peptide"], "peptide sequence")
    target_col = _resolve_col(df, ["_protein", "Protein_Sequence", "protein_sequence", "target_sequence"], "target protein sequence")
    assay_col = _resolve_col(df, ["_assay", "Affinity_Assay", "assay", "assay_name"], "assay")
    pkd_col = _resolve_col(df, ["_p_affinity", "pKd", "pkd", "p_affinity", "pKd_value"], "exact pKd")

    out = pd.DataFrame({
        "peptide": df[pep_col].fillna("").astype(str).str.replace(r"\s+", "", regex=True),
        "target_sequence": df[target_col].fillna("").astype(str).str.upper().str.replace(r"\s+", "", regex=True),
        "assay": df[assay_col].fillna("").astype(str).str.strip(),
        "pKd": pd.to_numeric(df[pkd_col], errors="coerce"),
    })
    valid = (
        out["peptide"].str.fullmatch(CANONICAL_RE).fillna(False)
        & out["target_sequence"].str.fullmatch(CANONICAL_RE).fillna(False)
        & out["assay"].ne("")
        & out["pKd"].notna()
        & np.isfinite(out["pKd"])
        & out["peptide"].str.len().between(9, 20)
    )
    out = out.loc[valid].copy()
    if out.empty:
        raise ValueError("No eligible canonical 9-20 aa peptide rows with exact numeric pKd were found.")
    # Aggregate technical/biological replicate rows before pair construction.
    collapsed = (
        out.groupby(["target_sequence", "assay", "peptide"], as_index=False, sort=True)
        .agg(pKd=("pKd", "median"), source_measurement_rows=("pKd", "size"))
    )
    return collapsed


def build_single_substitution_pairs(measurements: pd.DataFrame) -> pd.DataFrame:
    """Enumerate unordered pairs within each same-target/same-assay group."""
    rows: list[dict[str, Any]] = []
    for (target, assay), block in measurements.groupby(["target_sequence", "assay"], sort=True):
        records = block[["peptide", "pKd"]].drop_duplicates("peptide").to_dict("records")
        records.sort(key=lambda r: r["peptide"])
        for a, b in itertools.combinations(records, 2):
            seq_a, seq_b = a["peptide"], b["peptide"]
            if len(seq_a) != len(seq_b):
                continue
            diffs = [i for i, (aa, bb) in enumerate(zip(seq_a, seq_b)) if aa != bb]
            if len(diffs) != 1:
                continue
            if not (math.isfinite(float(a["pKd"])) and math.isfinite(float(b["pKd"]))):
                continue
            pos = diffs[0]
            pair_key = f"{seq_a}|{seq_b}"
            rows.append({
                "pair_id": _sha(f"{target}\t{assay}\t{pair_key}")[:16],
                "target_group_id": _sha(target)[:16],
                "target_sequence": target,
                "assay": assay,
                "peptide_a": seq_a,
                "peptide_b": seq_b,
                "pKd_a": float(a["pKd"]),
                "pKd_b": float(b["pKd"]),
                "delta_pKd": float(b["pKd"] - a["pKd"]),
                "mutation_position_zero_based": pos,
                "residue_a": seq_a[pos],
                "residue_b": seq_b[pos],
                "pair_signature": _sha(pair_key),
                "sequence_length": len(seq_a),
            })
    result = pd.DataFrame(rows)
    if result.empty:
        raise ValueError("No same-length, exactly-single-substitution peptide pairs were found.")
    if result["pair_id"].duplicated().any():
        # Duplicate target/assay/sequence pairs should be impossible after collapsing.
        raise AssertionError("Duplicate pair IDs indicate an input aggregation problem.")
    return result


def make_leakage_components(pairs: pd.DataFrame) -> np.ndarray:
    """Connect rows sharing a target sequence, peptide sequence, or exact pair.

    These connected components keep all assays/pairs for an exact target
    together and prevent any individual peptide or unordered peptide pair from
    appearing on both sides of the train/test split.
    """
    uf = UnionFind(len(pairs))
    first_seen: dict[tuple[str, str], int] = {}
    target_values = pairs["target_sequence"].astype(str).tolist()
    pair_values = pairs["pair_signature"].astype(str).tolist()
    peptide_a = pairs["peptide_a"].astype(str).tolist() if "peptide_a" in pairs else [""] * len(pairs)
    peptide_b = pairs["peptide_b"].astype(str).tolist() if "peptide_b" in pairs else [""] * len(pairs)
    for i, (target, pair_sig, pep_a, pep_b) in enumerate(zip(target_values, pair_values, peptide_a, peptide_b)):
        tokens = [("target", target), ("pair", pair_sig)]
        tokens.extend(("peptide", pep) for pep in {pep_a, pep_b} if pep)
        for token in tokens:
            if token in first_seen:
                uf.union(i, first_seen[token])
            else:
                first_seen[token] = i
    roots = [uf.find(i) for i in range(len(pairs))]
    root_to_group: dict[int, int] = {}
    return np.asarray([root_to_group.setdefault(root, len(root_to_group)) for root in roots], dtype=int)


def encode_pair(peptide_a: str, peptide_b: str) -> np.ndarray:
    """Fixed, target-free directional features for an exactly one-substitution pair.

    Feature blocks: 20x20 directed substitution identity; one-hot peptide length
    (9..20); normalized mutation position; one-hot left/right flanking residues.
    Pair orientation is lexical, and the target is pKd_b - pKd_a accordingly.
    No target ID, assay ID, affinity, native c_score, structural score, or fold
    assignment enters the predictor.
    """
    a, b = str(peptide_a), str(peptide_b)
    if len(a) != len(b) or not 9 <= len(a) <= 20:
        raise ValueError("Pair sequences must have equal lengths between 9 and 20 residues.")
    diffs = [i for i, (aa, bb) in enumerate(zip(a, b)) if aa != bb]
    if len(diffs) != 1:
        raise ValueError("Pair must differ at exactly one position.")
    pos = diffs[0]
    if any(aa not in AA_INDEX for aa in a + b):
        raise ValueError("Noncanonical amino acid in pair.")

    n_transition = 20 * 20
    n_length = 12
    n_position = 1
    n_flanks = 20 * 2
    feat = np.zeros(n_transition + n_length + n_position + n_flanks, dtype=float)
    old, new = AA_INDEX[a[pos]], AA_INDEX[b[pos]]
    feat[old * 20 + new] = 1.0

    length_start = n_transition
    feat[length_start + (len(a) - 9)] = 1.0
    position_idx = length_start + n_length
    feat[position_idx] = pos / max(1, len(a) - 1)
    flank_start = position_idx + 1
    if pos > 0:
        feat[flank_start + AA_INDEX[a[pos - 1]]] = 1.0
    if pos < len(a) - 1:
        feat[flank_start + 20 + AA_INDEX[a[pos + 1]]] = 1.0
    return feat


def metric_dict(y: np.ndarray, pred: np.ndarray) -> dict[str, Any]:
    y = np.asarray(y, dtype=float)
    pred = np.asarray(pred, dtype=float)
    finite = np.isfinite(y) & np.isfinite(pred)
    y, pred = y[finite], pred[finite]
    if not len(y):
        return {"n": 0, "mae": None, "rmse": None, "spearman_rho": None,
                "n_informative": 0, "informative_direction_coverage": None,
                "informative_direction_accuracy": None,
                "informative_direction_accuracy_ties_half": None}
    mae = float(np.mean(np.abs(y - pred)))
    rmse = float(np.sqrt(np.mean((y - pred) ** 2)))
    rho = None
    if len(y) >= 3 and np.unique(y).size > 1 and np.unique(pred).size > 1:
        rho = float(spearmanr(y, pred).statistic)
    informative = np.abs(y) >= INFORMATIVE_CUTOFF
    n_inf = int(informative.sum())
    if n_inf:
        yi, pi = y[informative], pred[informative]
        covered = np.abs(pi) >= DIRECTION_EPSILON
        coverage = float(covered.mean())
        accuracy = float(np.mean(np.sign(pi[covered]) == np.sign(yi[covered]))) if covered.any() else None
        tie_half = float(np.mean(np.where(np.abs(pi) < DIRECTION_EPSILON, 0.5, (np.sign(pi) == np.sign(yi)).astype(float))))
    else:
        coverage, accuracy, tie_half = None, None, None
    return {"n": int(len(y)), "mae": mae, "rmse": rmse, "spearman_rho": rho,
            "n_informative": n_inf, "informative_direction_coverage": coverage,
            "informative_direction_accuracy": accuracy,
            "informative_direction_accuracy_ties_half": tie_half}


def _cluster_bootstrap_differences(
    frame: pd.DataFrame, y_col: str, pred_col_a: str, pred_col_b: str,
    cluster_col: str, replicates: int, seed: int,
) -> dict[str, Any]:
    clusters = frame[cluster_col].drop_duplicates().tolist()
    if len(clusters) < 3:
        return {"cluster_count": len(clusters), "replicates_requested": replicates,
                "valid_replicates": 0, "note": "Fewer than 3 independent leakage components; CI suppressed."}
    rng = np.random.default_rng(seed)
    boot_mae_diff, boot_rho_diff, boot_direction_diff = [], [], []
    by_cluster = {g: frame.loc[frame[cluster_col].eq(g)] for g in clusters}
    for _ in range(replicates):
        sampled = rng.choice(clusters, size=len(clusters), replace=True)
        sample = pd.concat([by_cluster[g] for g in sampled], ignore_index=True)
        y = sample[y_col].to_numpy(float)
        a = sample[pred_col_a].to_numpy(float)
        b = sample[pred_col_b].to_numpy(float)
        ma, mb = metric_dict(y, a), metric_dict(y, b)
        if ma["mae"] is not None and mb["mae"] is not None:
            boot_mae_diff.append(float(ma["mae"] - mb["mae"]))
        if ma["spearman_rho"] is not None and mb["spearman_rho"] is not None:
            boot_rho_diff.append(float(ma["spearman_rho"] - mb["spearman_rho"]))
        da, db = ma["informative_direction_accuracy_ties_half"], mb["informative_direction_accuracy_ties_half"]
        if da is not None and db is not None:
            boot_direction_diff.append(float(da - db))
    def interval(values: list[float]) -> dict[str, float] | None:
        if not values:
            return None
        return {"ci95_low": float(np.quantile(values, 0.025)), "ci95_high": float(np.quantile(values, 0.975))}
    return {
        "cluster_count": len(clusters), "replicates_requested": int(replicates),
        "valid_replicates_mae": len(boot_mae_diff), "mae_difference_a_minus_b_ci95": interval(boot_mae_diff),
        "valid_replicates_spearman": len(boot_rho_diff), "spearman_difference_a_minus_b_ci95": interval(boot_rho_diff),
        "valid_replicates_direction_tie_half": len(boot_direction_diff), "direction_tie_half_difference_a_minus_b_ci95": interval(boot_direction_diff),
        "bootstrap_unit": "connected leakage component, resampled as a cluster",
    }


def run_validation(input_csv: str | Path, output_root: str | Path, n_splits: int = 5,
                   bootstrap_replicates: int = 2000, seed: int = 20261026) -> Path:
    measurements = prepare_measurements(input_csv)
    pairs = build_single_substitution_pairs(measurements)
    groups = make_leakage_components(pairs)
    pairs["leakage_component_id"] = groups
    n_components = int(np.unique(groups).size)
    if n_components < 3:
        sizes = pairs["leakage_component_id"].value_counts().to_dict()
        raise ValueError(
            f"Only {n_components} independent leakage components remain after grouping targets, individual peptides, and exact pairs. "
            f"Cannot make a credible 3+ fold held-out evaluation. Component sizes: {sizes}. "
            "Do not fall back to random row splits; obtain more independent target data or report this limitation."
        )
    folds = min(int(n_splits), n_components)
    if folds < 3:
        raise ValueError("At least 3 folds are required for this validation script.")

    X = np.vstack([encode_pair(a, b) for a, b in zip(pairs["peptide_a"], pairs["peptide_b"])])
    y = pairs["delta_pKd"].to_numpy(dtype=float)
    cv = GroupKFold(n_splits=folds)
    pred_ridge = np.full(len(pairs), np.nan)
    pred_zero = np.zeros(len(pairs), dtype=float)
    pred_fold_mean = np.full(len(pairs), np.nan)
    fold_ids = np.full(len(pairs), -1, dtype=int)
    split_rows = []
    for fold, (train_idx, test_idx) in enumerate(cv.split(X, y, groups=groups), start=1):
        train_groups = set(groups[train_idx].tolist())
        test_groups = set(groups[test_idx].tolist())
        if train_groups & test_groups:
            raise AssertionError(f"Leakage component overlap in fold {fold}.")
        train_targets = set(pairs.iloc[train_idx]["target_sequence"].astype(str))
        test_targets = set(pairs.iloc[test_idx]["target_sequence"].astype(str))
        train_pairs = set(pairs.iloc[train_idx]["pair_signature"].astype(str))
        test_pairs = set(pairs.iloc[test_idx]["pair_signature"].astype(str))
        train_peptides = set(pairs.iloc[train_idx]["peptide_a"].astype(str)) | set(pairs.iloc[train_idx]["peptide_b"].astype(str))
        test_peptides = set(pairs.iloc[test_idx]["peptide_a"].astype(str)) | set(pairs.iloc[test_idx]["peptide_b"].astype(str))
        if train_targets & test_targets:
            raise AssertionError(f"Target sequence leakage in fold {fold}: {len(train_targets & test_targets)} shared targets.")
        if train_pairs & test_pairs:
            raise AssertionError(f"Exact peptide-pair leakage in fold {fold}: {len(train_pairs & test_pairs)} shared pairs.")
        if train_peptides & test_peptides:
            raise AssertionError(f"Peptide-sequence leakage in fold {fold}: {len(train_peptides & test_peptides)} shared peptide sequences.")
        y_train = y[train_idx]
        pred_fold_mean[test_idx] = float(np.mean(y_train))
        model = Ridge(alpha=RIDGE_ALPHA)
        model.fit(X[train_idx], y_train)
        pred_ridge[test_idx] = model.predict(X[test_idx])
        fold_ids[test_idx] = fold
        split_rows.append({
            "fold": fold, "n_train_pairs": int(len(train_idx)), "n_test_pairs": int(len(test_idx)),
            "n_train_targets": int(len(train_targets)), "n_test_targets": int(len(test_targets)),
            "n_train_components": int(len(train_groups)), "n_test_components": int(len(test_groups)),
            "target_sequence_overlap": int(len(train_targets & test_targets)),
            "exact_pair_signature_overlap": int(len(train_pairs & test_pairs)),
            "peptide_sequence_overlap": int(len(train_peptides & test_peptides)),
            "train_target_sha256": sorted(_sha(s)[:16] for s in train_targets),
            "test_target_sha256": sorted(_sha(s)[:16] for s in test_targets),
        })
    if np.isnan(pred_ridge).any() or np.isnan(pred_fold_mean).any() or (fold_ids < 1).any():
        raise AssertionError("Every pair must receive exactly one held-out prediction.")

    out = pairs.copy()
    out["fold"] = fold_ids
    out["prediction_zero_change"] = pred_zero
    out["prediction_training_fold_mean"] = pred_fold_mean
    out["prediction_sequence_only_ridge"] = pred_ridge
    out["prediction_source"] = "grouped_target_and_pair_component_cross_validation"

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    root = Path(output_root)
    run_dir = root / f"run_{stamp}"
    run_dir.mkdir(parents=True, exist_ok=False)
    out.to_csv(run_dir / "out_of_fold_predictions.csv", index=False)
    pd.DataFrame(split_rows).to_json(run_dir / "fold_split_audit.json", orient="records", indent=2)
    # Store the group assignment and hashes needed to audit split identity.
    split_manifest = out[["pair_id", "fold", "leakage_component_id", "target_group_id", "pair_signature"]].copy()
    split_manifest["target_group_id"] = split_manifest["target_group_id"].astype(str)
    split_manifest.to_csv(run_dir / "split_manifest.csv", index=False)

    metrics = {
        "benchmark": "experimental_single_substitution_pKd_difference_sequence_only_grouped_validation",
        "claim_scope": "separate sequence-only comparator; not native c_score, not the adaptive optimizer, and not target-aware structural score",
        "input_csv": str(Path(input_csv).resolve()),
        "input_sha256": hashlib.sha256(Path(input_csv).read_bytes()).hexdigest(),
        "collapsed_target_assay_peptide_records": int(len(measurements)),
        "single_substitution_pairs": int(len(pairs)),
        "unique_target_sequences": int(pairs["target_sequence"].nunique()),
        "unique_exact_pair_signatures": int(pairs["pair_signature"].nunique()),
        "independent_leakage_components": n_components,
        "n_folds": folds,
        "split_rule": "GroupKFold over connected components; every exact target sequence, individual peptide sequence, and unordered peptide pair is confined to one fold",
        "split_leakage_checks": {
            "all_folds_have_zero_target_overlap": all(x["target_sequence_overlap"] == 0 for x in split_rows),
            "all_folds_have_zero_exact_pair_overlap": all(x["exact_pair_signature_overlap"] == 0 for x in split_rows),
            "all_folds_have_zero_peptide_sequence_overlap": all(x["peptide_sequence_overlap"] == 0 for x in split_rows),
            "all_rows_have_one_oof_prediction": True,
        },
        "model": {"name": "Ridge", "alpha": RIDGE_ALPHA, "feature_count": int(X.shape[1]),
                  "feature_definition": "directional amino-acid substitution, peptide length, normalized position, and flanking residues; sequence-only; no target/assay identifiers",
                  "hyperparameter_tuning": "none; alpha fixed before CV"},
        "metrics": {
            "zero_change_baseline": metric_dict(y, pred_zero),
            "training_fold_mean_baseline": metric_dict(y, pred_fold_mean),
            "sequence_only_ridge": metric_dict(y, pred_ridge),
        },
        "paired_bootstrap_ridge_minus_fold_mean": _cluster_bootstrap_differences(
            out, "delta_pKd", "prediction_sequence_only_ridge", "prediction_training_fold_mean",
            "leakage_component_id", bootstrap_replicates, seed),
        "bootstrap_seed": seed,
        "limitations": [
            "Retrospective benchmark on deposited affinity data; not prospective validation.",
            "Sequence-only comparator is separate from the production native c_score and is not integrated into candidate generation.",
            "Holding target sequences and exact peptide pairs together prevents those forms of split leakage, but does not establish generalization to HERV-K Env.",
            "Pair orientation is lexicographic; delta_pKd is pKd_b minus pKd_a and positive values mean peptide_b has stronger measured affinity.",
            "Experimental annotations may differ in construct, assay conditions, and uncertainty even after exact target-sequence grouping.",
        ],
        "result_directory": str(run_dir.resolve()),
    }
    (run_dir / "metrics.json").write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    report = [
        "# Leakage-resistant sequence-only mutation-pair validation",
        "",
        f"- Input file: `{Path(input_csv).resolve()}`",
        f"- Collapsed target/assay/peptide records: **{len(measurements)}**",
        f"- Single-substitution pairs: **{len(pairs)}**",
        f"- Unique target sequences: **{pairs['target_sequence'].nunique()}**",
        f"- Leakage components: **{n_components}**; folds: **{folds}**",
        "- Split rule: connected components link all rows sharing a target sequence, an individual peptide sequence, or an exact unordered peptide pair.",
        "- Each fold was checked for zero shared target sequences, individual peptide sequences, and exact pair signatures.",
        f"- Ridge alpha: **{RIDGE_ALPHA}**, fixed without tuning on held-out folds.",
        "",
        "## Out-of-fold metrics",
        "",
        "| Model | n | MAE (pKd units) | RMSE | Spearman rho | Informative direction coverage | Direction accuracy among covered | Direction score (ties half) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for label, key in [("Zero-change baseline", "zero_change_baseline"),
                       ("Training-fold mean", "training_fold_mean_baseline"),
                       ("Sequence-only Ridge", "sequence_only_ridge")]:
        m = metrics["metrics"][key]
        def fmt(x: Any) -> str:
            return "NA" if x is None else f"{x:.4f}"
        report.append(f"| {label} | {m['n']} | {fmt(m['mae'])} | {fmt(m['rmse'])} | {fmt(m['spearman_rho'])} | {fmt(m['informative_direction_coverage'])} | {fmt(m['informative_direction_accuracy'])} | {fmt(m['informative_direction_accuracy_ties_half'])} |")
    report += [
        "",
        "## Interpretation guardrails",
        "",
        "This is a separate sequence-only research comparator. It does not validate the native `c_score`, prove the adaptive optimizer improves, or validate the target-aware structural scorer. A positive result supports only the reported held-out statistical relationship under this split and feature definition.",
        "",
        "Do not report target-binding, biological efficacy, or HERV-K-specific accuracy from this benchmark. A null/weak result is a valid scientific outcome.",
        "",
    ]
    (run_dir / "REPORT.md").write_text("\n".join(report), encoding="utf-8")
    print(f"RUN_DIR={run_dir.resolve()}")
    print(json.dumps(metrics["metrics"], indent=2))
    print(f"Leakage audit: {metrics['split_leakage_checks']}")
    return run_dir


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Filtered exact-Kd CSV produced by ppikb_kd_benchmark.py")
    parser.add_argument("--out-root", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--bootstrap-replicates", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20261026)
    args = parser.parse_args(argv)
    if args.folds < 3:
        parser.error("--folds must be at least 3")
    if args.bootstrap_replicates < 0:
        parser.error("--bootstrap-replicates must be nonnegative")
    try:
        run_validation(args.input, args.out_root, args.folds, args.bootstrap_replicates, args.seed)
    except Exception as exc:
        print(f"FAIL: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
