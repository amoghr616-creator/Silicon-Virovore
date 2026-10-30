from __future__ import annotations

import csv
import json
import math
import re
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import rankdata, spearmanr
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "prototype3" / "results"
SOURCE = RESULTS / "ppikb_kd_benchmark" / "filtered_kd_measurements.csv"
STAMP = datetime.now().strftime("%Y%m%dT%H%M%S")
OUT = RESULTS / f"p3_ppikb_grouped_comparison_{STAMP}"
OUT.mkdir(parents=True, exist_ok=False)

AA = "ACDEFGHIKLMNPQRSTVWY"
MAX_LEN = 20
ALPHAS = [0.1, 1.0, 10.0, 100.0]
RNG = np.random.default_rng(271027)


def read_rows(path):
    with path.open(newline="", errors="replace") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    rows = list(rows)
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def clean_peptide(x):
    s = re.sub(r"\s+", "", str(x or "").upper())
    return s if s and re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", s) else ""


def clean_protein(x):
    return re.sub(r"\s+", "", str(x or "").upper())


def number(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (ValueError, TypeError):
        return None


def is_true(x):
    return str(x or "").strip().lower() in {"true", "1", "yes"}


def is_false(x):
    return str(x or "").strip().lower() in {"false", "0", "no"}


def rank_corr(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    if len(a) < 3 or len(np.unique(a)) < 2 or len(np.unique(b)) < 2:
        return np.nan
    return float(spearmanr(a, b).statistic)


def boot_ci(values, repetitions=5000):
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if not len(x):
        return None, None
    means = np.array([
        np.mean(RNG.choice(x, size=len(x), replace=True))
        for _ in range(repetitions)
    ])
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def paired_signflip_p(differences, repetitions=20000):
    d = np.asarray(differences, dtype=float)
    d = d[np.isfinite(d)]
    if len(d) < 2:
        return None
    observed = abs(float(np.mean(d)))
    signs = RNG.choice([-1.0, 1.0], size=(repetitions, len(d)))
    null = np.abs(np.mean(signs * d, axis=1))
    return float((1 + np.sum(null >= observed)) / (repetitions + 1))


def load_measurements():
    if not SOURCE.exists():
        raise SystemExit(f"Required existing dataset not found: {SOURCE}")

    raw = read_rows(SOURCE)
    if not raw:
        raise SystemExit("Input CSV is empty.")

    required = {
        "_peptide", "_protein", "_assay", "_p_affinity",
        "_native_score", "_molar", "_affinity_type",
    }
    missing = required - set(raw[0])
    if missing:
        raise SystemExit(f"ABORT: expected verified PPIKB columns missing: {sorted(missing)}")

    rejected = defaultdict(int)
    buckets = defaultdict(list)

    for row in raw:
        peptide = clean_peptide(row.get("_peptide"))
        protein = clean_protein(row.get("_protein"))
        assay = str(row.get("_assay", "")).strip()
        p_affinity = number(row.get("_p_affinity"))
        molar = number(row.get("_molar"))
        score = number(row.get("_native_score"))
        affinity_type = str(row.get("_affinity_type", "")).strip().upper()

        if not peptide or not (9 <= len(peptide) <= 20):
            rejected["noncanonical_or_length"] += 1
            continue
        if not protein or not assay:
            rejected["missing_target_or_assay"] += 1
            continue
        if affinity_type != "KD":
            rejected["not_KD"] += 1
            continue
        if molar is None or molar <= 0 or p_affinity is None:
            rejected["invalid_affinity"] += 1
            continue
        if not is_true(row.get("_eligible")):
            rejected["not_eligible"] += 1
            continue
        if not is_true(row.get("_target_valid")):
            rejected["invalid_target"] += 1
            continue
        if not is_true(row.get("_canonical")):
            rejected["noncanonical_flag"] += 1
            continue
        if not is_true(row.get("_length_ok")):
            rejected["length_flag"] += 1
            continue
        if not is_true(row.get("_linear")) or not is_true(row.get("_unmodified")):
            rejected["modified_or_non_linear"] += 1
            continue
        if not is_false(row.get("_censored")):
            rejected["censored"] += 1
            continue
        if score is None:
            rejected["missing_native_score"] += 1
            continue

        # Higher pKd means stronger measured affinity.
        # Context is exact protein sequence plus assay.
        context = protein + "\x1f" + assay
        buckets[(context, peptide)].append({
            "context": context,
            "protein": protein,
            "assay": assay,
            "sequence": peptide,
            "length": len(peptide),
            "pKd": p_affinity,
            "native_score": score,
            "source_id": row.get("ID", ""),
            "n_source_rows": 1,
        })

    data = []
    for (context, peptide), items in buckets.items():
        data.append({
            "context": context,
            "protein": items[0]["protein"],
            "assay": items[0]["assay"],
            "sequence": peptide,
            "length": len(peptide),
            "pKd": float(np.median([x["pKd"] for x in items])),
            "native_score": float(np.median([x["native_score"] for x in items])),
            "source_rows_collapsed": len(items),
            "source_ids": ";".join(sorted({x["source_id"] for x in items if x["source_id"]})),
        })

    info = {
        "input_rows": len(raw),
        "aggregated_context_sequence_rows": len(data),
        "unique_peptide_sequences": len({r["sequence"] for r in data}),
        "unique_target_assay_contexts": len({r["context"] for r in data}),
        "rejected_row_counts": dict(rejected),
        "source": str(SOURCE.relative_to(ROOT)),
        "affinity_measure": "uncensored measured Kd converted to pKd; higher pKd means stronger affinity",
        "aggregation": "median per exact protein-sequence + assay + canonical peptide sequence",
    }
    return data, info


def components(rows):
    """Keep every exact target-assay context and peptide sequence in one fold."""
    parent = list(range(len(rows)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        x, y = find(a), find(b)
        if x != y:
            parent[y] = x

    first_context, first_sequence = {}, {}
    for i, r in enumerate(rows):
        if r["context"] in first_context:
            union(i, first_context[r["context"]])
        else:
            first_context[r["context"]] = i
        if r["sequence"] in first_sequence:
            union(i, first_sequence[r["sequence"]])
        else:
            first_sequence[r["sequence"]] = i

    roots = [find(i) for i in range(len(rows))]
    remap = {root: j for j, root in enumerate(sorted(set(roots), key=str))}
    return np.array([remap[x] for x in roots], dtype=int)


def make_feature_vector(r, assays):
    seq = r["sequence"]
    protein = r["protein"]
    out = []

    # Position-specific peptide amino-acid encoding.
    for i in range(MAX_LEN):
        residue = seq[i] if i < len(seq) else "-"
        out.extend([1.0 if residue == aa else 0.0 for aa in AA])

    # Amino-acid composition.
    for aa in AA:
        out.append(seq.count(aa) / len(seq))

    # Dipeptide composition.
    di = defaultdict(int)
    for i in range(len(seq) - 1):
        di[seq[i:i + 2]] += 1
    denom = max(1, len(seq) - 1)
    for a in AA:
        for b in AA:
            out.append(di[a + b] / denom)

    # Target composition, not an exact target identity shortcut.
    valid_target = [x for x in protein if x in AA]
    for aa in AA:
        out.append(valid_target.count(aa) / max(1, len(valid_target)))
    out.append(len(valid_target) / 2000.0)

    # Assay one-hot features.
    out.extend([1.0 if r["assay"] == assay else 0.0 for assay in assays])

    # Existing sequence-only score as one predictor feature.
    out.append(r["native_score"])
    out.append(len(seq) / MAX_LEN)
    return np.asarray(out, dtype=float)


def fill_scale_fit_predict(X_train, y_train, X_test, alpha):
    med = np.nanmedian(X_train, axis=0)
    med[~np.isfinite(med)] = 0.0
    A = np.where(np.isfinite(X_train), X_train, med)
    B = np.where(np.isfinite(X_test), X_test, med)

    mean = A.mean(axis=0)
    sd = A.std(axis=0)
    sd[sd < 1e-12] = 1.0

    A = (A - mean) / sd
    B = (B - mean) / sd
    model = Ridge(alpha=alpha)
    model.fit(A, y_train)
    return model.predict(B)


def choose_alpha(X, y, groups):
    unique = len(np.unique(groups))
    if unique < 3:
        return 10.0
    n = min(3, unique)
    splitter = GroupKFold(n_splits=n)
    loss = {}
    for alpha in ALPHAS:
        errors = []
        for tr, va in splitter.split(X, y, groups):
            pred = fill_scale_fit_predict(X[tr], y[tr], X[va], alpha)
            errors.extend(np.abs(y[va] - pred).tolist())
        loss[alpha] = float(np.mean(errors)) if errors else float("inf")
    return min(loss, key=loss.get)


def context_metrics(rows, policies):
    by_context = defaultdict(list)
    for r in rows:
        by_context[r["context"]].append(r)

    selection_rows = []
    performance_rows = []

    for context, group in by_context.items():
        if len(group) < 3:
            continue

        actual = np.array([r["pKd"] for r in group])
        ranks = rankdata(actual, method="average")
        pct = (ranks - 1) / max(1, len(group) - 1)

        native = np.array([r["native_score"] for r in group])
        perf = {
            "context": context,
            "assay": group[0]["assay"],
            "n_distinct_peptides": len(group),
            "native_spearman": rank_corr(native, actual),
        }

        if all("ridge_prediction" in r for r in group):
            ridge = np.array([r["ridge_prediction"] for r in group])
            mean_pred = np.array([r["training_fold_mean_prediction"] for r in group])
            perf.update({
                "ridge_spearman": rank_corr(ridge, actual),
                "ridge_mae": float(np.mean(np.abs(actual - ridge))),
                "ridge_rmse": float(np.sqrt(np.mean((actual - ridge) ** 2))),
                "training_fold_mean_mae": float(np.mean(np.abs(actual - mean_pred))),
                "training_fold_mean_rmse": float(np.sqrt(np.mean((actual - mean_pred) ** 2))),
            })

        performance_rows.append(perf)

        if len(group) < 5:
            continue

        k = max(1, int(math.ceil(0.20 * len(group))))
        random_expected = float(np.mean(pct))

        available = [
            ("native_score_policy", native),
        ]
        if all("ridge_prediction" in r for r in group):
            available.append((
                "grouped_ridge_policy",
                np.array([r["ridge_prediction"] for r in group])
            ))

        for name, policy_score in available:
            order = np.argsort(-policy_score, kind="stable")[:k]
            selected = float(np.mean(pct[order]))
            selection_rows.append({
                "context": context,
                "assay": group[0]["assay"],
                "policy": name,
                "n_peptides": len(group),
                "selection_budget_k": k,
                "mean_selected_pKd_percentile": selected,
                "random_selection_expected_percentile": random_expected,
                "lift_over_random_expectation": selected - random_expected,
            })

    return performance_rows, selection_rows


def main():
    data, load_info = load_measurements()
    write_csv(OUT / "eligible_aggregated_measurements.csv", data)

    # Primary native score analysis at the context level.
    by_context = defaultdict(list)
    for r in data:
        by_context[r["context"]].append(r)

    corr_rows = []
    for context, group in by_context.items():
        if len(group) < 3:
            continue
        rho = rank_corr(
            [r["native_score"] for r in group],
            [r["pKd"] for r in group],
        )
        corr_rows.append({
            "context": context,
            "assay": group[0]["assay"],
            "distinct_peptides": len(group),
            "spearman_native_score_vs_pKd": rho,
        })

    write_csv(OUT / "native_score_context_correlations.csv", corr_rows)

    native_means = {}
    for minimum in (3, 5):
        vals = [
            r["spearman_native_score_vs_pKd"]
            for r in corr_rows
            if r["distinct_peptides"] >= minimum
            and np.isfinite(r["spearman_native_score_vs_pKd"])
        ]
        low, high = boot_ci(vals)
        native_means[str(minimum)] = {
            "minimum_distinct_peptides": minimum,
            "contexts": len(vals),
            "mean_spearman": float(np.mean(vals)) if vals else None,
            "median_spearman": float(np.median(vals)) if vals else None,
            "positive_contexts": int(sum(x > 0 for x in vals)),
            "negative_contexts": int(sum(x < 0 for x in vals)),
            "bootstrap_95ci_low": low,
            "bootstrap_95ci_high": high,
        }

    # The strict predictor validation uses connected components.
    labels = components(data)
    component_count = len(np.unique(labels))
    ridge_info = {
        "status": "not_estimable",
        "connected_components": component_count,
        "split_rule": "Exact peptide sequences and exact target-assay contexts cannot cross outer folds.",
        "features": "peptide sequence, composition, dipeptide composition, target composition, assay, native score",
        "target": "measured pKd from uncensored Kd records only",
    }

    predictions = []
    selection_rows = []
    performance_rows = []

    if component_count >= 5:
        assays = sorted({r["assay"] for r in data})
        X = np.vstack([make_feature_vector(r, assays) for r in data])
        y = np.array([r["pKd"] for r in data], dtype=float)
        nfolds = min(5, component_count)
        outer = GroupKFold(n_splits=nfolds)

        for fold, (train, test) in enumerate(outer.split(X, y, labels), 1):
            alpha = choose_alpha(X[train], y[train], labels[train])
            ridge_pred = fill_scale_fit_predict(
                X[train], y[train], X[test], alpha
            )
            fold_mean = float(np.mean(y[train]))

            for index, pred in zip(test, ridge_pred):
                row = dict(data[index])
                row.update({
                    "outer_fold": fold,
                    "leakage_component": int(labels[index]),
                    "ridge_prediction": float(pred),
                    "training_fold_mean_prediction": fold_mean,
                    "selected_alpha": alpha,
                })
                predictions.append(row)

        write_csv(OUT / "grouped_oof_predictions.csv", predictions)
        performance_rows, selection_rows = context_metrics(
            predictions, ["native_score_policy", "grouped_ridge_policy"]
        )
        write_csv(OUT / "heldout_context_performance.csv", performance_rows)
        write_csv(OUT / "heldout_top20_selection.csv", selection_rows)

        mae_diffs = [
            r["ridge_mae"] - r["training_fold_mean_mae"]
            for r in performance_rows
            if "ridge_mae" in r
        ]
        low, high = boot_ci(mae_diffs)
        ridge_rhos = [
            r["ridge_spearman"] for r in performance_rows
            if np.isfinite(r.get("ridge_spearman", np.nan))
        ]
        low_rho, high_rho = boot_ci(ridge_rhos)

        ridge_info.update({
            "status": "completed",
            "outer_folds": nfolds,
            "contexts_with_oof_predictions": len(performance_rows),
            "mean_paired_context_MAE_difference_ridge_minus_mean": (
                float(np.mean(mae_diffs)) if mae_diffs else None
            ),
            "bootstrap_95ci_MAE_difference_low": low,
            "bootstrap_95ci_MAE_difference_high": high,
            "mean_heldout_context_ridge_spearman": (
                float(np.mean(ridge_rhos)) if ridge_rhos else None
            ),
            "bootstrap_95ci_ridge_spearman_low": low_rho,
            "bootstrap_95ci_ridge_spearman_high": high_rho,
            "interpretation": (
                "Negative MAE difference favors Ridge. A correlation estimate is included "
                "only for contexts with enough prediction and outcome variation."
            ),
        })

    write_csv(OUT / "top20_selection_summary.csv", selection_rows)

    selection_summary = []
    for policy in ("native_score_policy", "grouped_ridge_policy"):
        vals = [
            r["lift_over_random_expectation"]
            for r in selection_rows if r["policy"] == policy
        ]
        selected = [
            r["mean_selected_pKd_percentile"]
            for r in selection_rows if r["policy"] == policy
        ]
        low, high = boot_ci(vals)
        selection_summary.append({
            "policy": policy,
            "contexts": len(vals),
            "mean_selected_affinity_percentile": (
                float(np.mean(selected)) if selected else None
            ),
            "mean_lift_over_uniform_expectation": (
                float(np.mean(vals)) if vals else None
            ),
            "bootstrap_95ci_lift_low": low,
            "bootstrap_95ci_lift_high": high,
            "two_sided_sign_flip_p_descriptive": paired_signflip_p(vals),
        })
    write_csv(OUT / "top20_selection_summary_by_policy.csv", selection_summary)

    # Graph 1: target-assay group native-score correlation distribution.
    vals3 = [
        r["spearman_native_score_vs_pKd"] for r in corr_rows
        if r["distinct_peptides"] >= 3
        and np.isfinite(r["spearman_native_score_vs_pKd"])
    ]
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.hist(vals3, bins=25)
    ax.axvline(0, linestyle="--", label="No association")
    ax.axvline(np.mean(vals3), linestyle=":", label=f"Mean = {np.mean(vals3):.3f}")
    ax.set_title("Native C-score versus measured pKd within target-assay groups")
    ax.set_xlabel("Within-group Spearman correlation")
    ax.set_ylabel("Number of target-assay groups")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "01_native_score_group_correlation_distribution.png", dpi=180)
    plt.close(fig)

    # Graph 2: native score and Ridge selection policy compared with uniform.
    plot_rows = []
    for policy in ("native_score_policy", "grouped_ridge_policy"):
        subset = [r for r in selection_rows if r["policy"] == policy]
        if subset:
            lifts = [r["lift_over_random_expectation"] for r in subset]
            low, high = boot_ci(lifts)
            plot_rows.append((policy, float(np.mean(lifts)), low, high))
    if plot_rows:
        fig, ax = plt.subplots(figsize=(8, 4.8))
        names = [r[0].replace("_", "\n") for r in plot_rows]
        means = [r[1] for r in plot_rows]
        lo = [r[1] - r[2] for r in plot_rows]
        hi = [r[3] - r[1] for r in plot_rows]
        ax.bar(names, means, yerr=[lo, hi], capsize=5)
        ax.axhline(0, linestyle="--", label="Uniform-selection expectation")
        ax.set_title("Top-20% selection lift on held-out experimental peptide groups")
        ax.set_ylabel("Affinity-percentile lift over uniform selection")
        ax.legend()
        fig.tight_layout()
        fig.savefig(OUT / "02_top20_selection_lift.png", dpi=180)
        plt.close(fig)

    # Graph 3: Ridge versus fold mean MAE, aggregated by held-out context.
    if ridge_info["status"] == "completed":
        valid = [
            r for r in performance_rows
            if "ridge_mae" in r and "training_fold_mean_mae" in r
        ]
        if valid:
            a = np.array([r["ridge_mae"] for r in valid])
            b = np.array([r["training_fold_mean_mae"] for r in valid])
            fig, ax = plt.subplots(figsize=(7, 5))
            ax.scatter(b, a, alpha=0.65)
            low_bound = min(float(np.min(a)), float(np.min(b)))
            high_bound = max(float(np.max(a)), float(np.max(b)))
            ax.plot([low_bound, high_bound], [low_bound, high_bound], linestyle="--")
            ax.set_xlabel("Training-fold-mean MAE (pKd)")
            ax.set_ylabel("Grouped Ridge MAE (pKd)")
            ax.set_title("Held-out context prediction error")
            fig.tight_layout()
            fig.savefig(OUT / "03_ridge_vs_fold_mean_mae.png", dpi=180)
            plt.close(fig)

    # Preserve existing exact-21-aa and receptor-context results.
    prior_files = {
        "exact21_metrics": RESULTS / "p3_21aa_affinity_audit" / "native_c_score_21aa" / "metrics.json",
        "receptor_context_summary": RESULTS / "pepbi_target_mutation_invariance" / "summary.json",
        "arise_position_audit": RESULTS / "arise_experimental_position_audit" / "metrics.json",
        "pepbi_summary": RESULTS / "pepbi_corrected_benchmark" / "benchmark_summary.json",
    }
    prior = {}
    for key, path in prior_files.items():
        if path.exists():
            try:
                prior[key] = json.loads(path.read_text())
            except json.JSONDecodeError:
                prior[key] = {"status": "unreadable_json", "path": str(path)}
    (OUT / "prior_validation_metrics.json").write_text(json.dumps(prior, indent=2))

    receptor = prior.get("receptor_context_summary", {})
    if receptor:
        names = ["Native score changes", "Experimental ΔG changes\n>0.5 kcal/mol"]
        values = [
            receptor.get("pairs_with_native_score_difference_over_1e-12", 0),
            receptor.get("pairs_with_experimental_delta_g_difference_over_0_5_kcal_mol", 0),
        ]
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.bar(names, values)
        ax.set_title("Existing receptor-context invariance control")
        ax.set_ylabel("Matched target-variant pairs")
        fig.tight_layout()
        fig.savefig(OUT / "04_receptor_context_invariance.png", dpi=180)
        plt.close(fig)

    summary = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "input_summary": load_info,
        "native_score_group_correlations": native_means,
        "grouped_ridge_validation": ridge_info,
        "top20_selection_summary": selection_summary,
        "prior_evidence": {
            "exact21_available": "exact21_metrics" in prior,
            "receptor_context_available": "receptor_context_summary" in prior,
            "arise_position_audit_available": "arise_position_audit" in prior,
            "pepbi_summary_available": "pepbi_summary" in prior,
        },
        "interpretation_guardrails": [
            "PPIKB is retrospective and has been examined in previous project analyses; this is not a pristine prospective test.",
            "Native-score correlations concern the score, not the optimizer's complete candidate-generation policy.",
            "Grouped Ridge measures a candidate alternative scoring approach and does not overwrite the original adaptive optimizer.",
            "No result here proves that any unmeasured candidate binds HERV-K Env.",
            "An inconclusive strict grouped validation is reported as inconclusive, not repaired by leaking target or peptide information."
        ],
    }
    (OUT / "P3_SUMMARY.json").write_text(json.dumps(summary, indent=2))

    # Main report.
    lines = [
        "# P3 PPIKB grouped validation and selection comparison",
        "",
        f"Generated: {summary['generated']}",
        "",
        "## Experimental design",
        "",
        "This is a retrospective, controlled comparison of score-guided selection policies on measured PPIKB peptide data. The optimizer and docking were not run. Kd values are converted to pKd; higher pKd means stronger affinity. Duplicate records are collapsed by exact protein sequence, assay, and peptide sequence.",
        "",
        "Native-score selection, grouped Ridge prediction, and the uniform-selection reference use the same top-20% selection budget in each eligible held-out context. For Ridge, exact peptide sequences and target-assay contexts are kept within the same connected component before outer fold assignment; model preprocessing and alpha selection occur inside training folds.",
        "",
        "## Data audit",
        "",
        f"- Input rows: {load_info['input_rows']}",
        f"- Aggregated target-assay-peptide rows: {load_info['aggregated_context_sequence_rows']}",
        f"- Distinct peptide sequences: {load_info['unique_peptide_sequences']}",
        f"- Target-assay contexts: {load_info['unique_target_assay_contexts']}",
        "",
        "## Test 3: Native score against measured affinity",
        "",
    ]
    for key, metrics in native_means.items():
        lines += [
            f"### Groups with at least {metrics['minimum_distinct_peptides']} peptides",
            "",
            f"- Groups: {metrics['contexts']}",
            f"- Mean within-group Spearman: {metrics['mean_spearman']}",
            f"- Median: {metrics['median_spearman']}",
            f"- Positive/negative groups: {metrics['positive_contexts']}/{metrics['negative_contexts']}",
            f"- Group-bootstrap 95% CI: [{metrics['bootstrap_95ci_low']}, {metrics['bootstrap_95ci_high']}]",
            "",
        ]

    lines += [
        "## Tests 6 and 7: Grouped prediction and measured-peptide selection",
        "",
        f"Grouped Ridge status: **{ridge_info['status']}**.",
        "",
    ]
    if ridge_info["status"] == "completed":
        lines += [
            f"- Leakage-separated connected components: {ridge_info['connected_components']}.",
            f"- Outer folds: {ridge_info['outer_folds']}.",
            f"- Contexts with held-out predictions: {ridge_info['contexts_with_oof_predictions']}.",
            f"- Ridge-minus-fold-mean paired context MAE difference: {ridge_info['mean_paired_context_MAE_difference_ridge_minus_mean']}.",
            f"- Group-bootstrap 95% interval for MAE difference: [{ridge_info['bootstrap_95ci_MAE_difference_low']}, {ridge_info['bootstrap_95ci_MAE_difference_high']}].",
            f"- Mean held-out context Ridge Spearman: {ridge_info['mean_heldout_context_ridge_spearman']}.",
            "",
            "Negative Ridge-minus-mean MAE favors Ridge. Interpret this as an exploratory predictor comparison, not proof of a target-specific binding advantage.",
            "",
        ]
    else:
        lines += [
            "Too few leakage-separated connected components were available for a defensible grouped predictor comparison. The native-score and selection-policy audits remain available; the Ridge result must not be claimed.",
            "",
        ]

    lines += [
        "### Top-20% selection",
        "",
        "| Policy | Contexts | Mean selected percentile | Mean lift over uniform expectation | Group-bootstrap 95% CI |",
        "|---|---:|---:|---:|---:|",
    ]
    for r in selection_summary:
        lines.append(
            f"| {r['policy']} | {r['contexts']} | "
            f"{r['mean_selected_affinity_percentile']} | "
            f"{r['mean_lift_over_uniform_expectation']} | "
            f"[{r['bootstrap_95ci_lift_low']}, {r['bootstrap_95ci_lift_high']}] |"
        )

    lines += [
        "",
        "A positive lift means selected peptides had higher measured affinity percentiles than the uniform-selection expectation. A confidence interval containing zero is inconclusive. These are already-measured peptide landscapes, not evidence about newly generated HERV-K Env candidates.",
        "",
        "## Tests 4, 5, and 9: supporting controls",
        "",
        "The exact-21-aa PPIKB analysis remains a separate, sparse subgroup. The saved ARISE mutation-position audit is exploratory and previously examined; it is not independent validation of a newly tuned policy. The receptor-context audit shows whether a sequence-only score responds to target changes. Its result is a limitation check, not a policy victory metric.",
        "",
        "See `prior_validation_metrics.json` for the recorded summaries.",
        "",
        "## Overall interpretation",
        "",
        "The current native-score benchmark is unfavorable: the mean within-context association with pKd is weakly negative. This justifies testing a separately trained, leakage-controlled affinity predictor as a potential selection component. It does not justify retroactively altering the original P2 claims or relabeling the new predictor as the original adaptive algorithm.",
        "",
        "A positive held-out retrospective selection result would support prioritization among measured peptide variants in the evaluated datasets. It would not prove that new candidates bind HERV-K Env. That claim requires valid target-specific experimental evidence.",
        "",
        "## Output inventory",
        "",
    ]
    for p in sorted(OUT.iterdir()):
        if p.is_file():
            lines.append(f"- `{p.name}`")
    (OUT / "P3_REPORT.md").write_text("\n".join(lines) + "\n")

    data_log = [
        f"# Logbook — Data and results ({datetime.now().strftime('%Y-%m-%d')})",
        "",
        "**Experiment:** PPIKB grouped affinity and selection-policy validation.",
        "",
        "**Objective:** Measure whether the current native C-score ranks experimental peptide affinity within target-assay groups and compare it with a grouped Ridge predictor and a uniform-selection reference.",
        "",
        "**Design:** Uncensored Kd records for canonical, unmodified, linear peptides were filtered using the stored eligibility flags. Records were aggregated by exact protein sequence, assay, and peptide sequence using median pKd and score. Higher pKd indicates stronger measured affinity. The selection metric was the mean measured affinity percentile of the top 20% selected within each eligible context.",
        "",
        f"**Dataset counts:** {load_info['input_rows']} input rows; {load_info['aggregated_context_sequence_rows']} aggregated context-sequence rows; {load_info['unique_peptide_sequences']} unique peptide sequences; {load_info['unique_target_assay_contexts']} exact target-assay contexts.",
        "",
        "**Native score results:**",
    ]
    for key, m in native_means.items():
        data_log.append(
            f"- Groups with at least {m['minimum_distinct_peptides']} peptides: "
            f"n={m['contexts']}; mean Spearman={m['mean_spearman']}; "
            f"95% group-bootstrap CI=[{m['bootstrap_95ci_low']}, {m['bootstrap_95ci_high']}]; "
            f"positive/negative={m['positive_contexts']}/{m['negative_contexts']}."
        )
    data_log += [
        "",
        f"**Grouped Ridge status:** {ridge_info['status']}.",
        "",
        f"**Top-selection summary:** {json.dumps(selection_summary, ensure_ascii=False)}",
        "",
        "**Graphs:** `01_native_score_group_correlation_distribution.png`, `02_top20_selection_lift.png`, `03_ridge_vs_fold_mean_mae.png` when estimable, and `04_receptor_context_invariance.png` when prior audit data are available.",
        "",
        f"**Output directory:** `{OUT.relative_to(ROOT)}`",
    ]
    (OUT / "P3_LOGBOOK_DATA_ENTRY.md").write_text("\n".join(data_log) + "\n")

    analysis_log = [
        f"# Logbook — Analysis and interpretation ({datetime.now().strftime('%Y-%m-%d')})",
        "",
        "The native score's within-group correlation with pKd assesses the score, not the full adaptive optimizer. Weak or negative group-level association means a higher internal score has not been shown to prioritize peptides with stronger measured affinity.",
        "",
        "The equal-budget top-20% replay measures selection performance on already-measured peptides. Positive lift over uniform selection suggests that a policy finds stronger measured variants in those contexts; a confidence interval crossing zero is inconclusive. These retrospective records cannot establish binding of new candidates to HERV-K Env.",
        "",
        "The Ridge model is evaluated separately using grouped held-out predictions. A genuine improvement requires lower held-out error than the training-fold-mean benchmark and no leakage of the same peptide sequence or target-assay context across folds. If strict grouping leaves too few independent components, the model result is not estimable.",
        "",
        "This PPIKB dataset has been examined in prior project work, so the results must be described as retrospective exploratory validation rather than an untouched prospective confirmation. The original P2 optimizer comparison remains unchanged.",
        "",
        "The next model change should be based on this evidence and evaluated in a new, pre-specified comparison—not tuned repeatedly against the same hold-out results.",
    ]
    (OUT / "P3_LOGBOOK_ANALYSIS_ENTRY.md").write_text("\n".join(analysis_log) + "\n")

    print("\nPPIKB grouped comparison completed.")
    print("Optimizer/docking runs: NONE")
    print("Output directory:", OUT.relative_to(ROOT))
    print("\nDataset:", json.dumps(load_info, indent=2))
    print("\nNative-score within-group results:")
    print(json.dumps(native_means, indent=2))
    print("\nGrouped Ridge:")
    print(json.dumps(ridge_info, indent=2))
    print("\nTop-20 selection summary:")
    print(json.dumps(selection_summary, indent=2))
    print("\nMain report:", OUT.relative_to(ROOT) / "P3_REPORT.md")
    print("Logbook data:", OUT.relative_to(ROOT) / "P3_LOGBOOK_DATA_ENTRY.md")
    print("Logbook analysis:", OUT.relative_to(ROOT) / "P3_LOGBOOK_ANALYSIS_ENTRY.md")


if __name__ == "__main__":
    main()
