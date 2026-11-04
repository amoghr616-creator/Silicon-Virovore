from pathlib import Path
from itertools import product
import sqlite3
import json
import csv
import math

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from scipy.stats import binomtest

ROOT = Path("prototype4/results/model_a_validation_20261030")
DB = ROOT / "p4_model_a_working.sqlite"
OUT = ROOT / "p4_nested_threshold_directional_20261101"

SEED = 20261101
THRESHOLD = 0.5
C_VALUES = [0.001, 0.01, 0.1, 1.0, 10.0, 100.0]
CLASS_WEIGHTS = [None, "balanced"]
PROBABILITY_THRESHOLDS = np.linspace(0.20, 0.80, 13)
INNER_FOLDS = 4
BOOTSTRAPS = 20000


def write_csv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=columns,
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def make_model(C, class_weight):
    return make_pipeline(
        StandardScaler(),
        LogisticRegression(
            C=C,
            penalty="l2",
            solver="liblinear",
            class_weight=class_weight,
            max_iter=5000,
            random_state=SEED,
        )
    )


def fit_probability_model(X_train, y_train, X_valid, C, class_weight):
    unique = np.unique(y_train)

    if len(unique) < 2:
        return np.full(
            len(X_valid),
            float(np.mean(y_train))
        ), "single_class_training_fallback"

    model = make_model(C, class_weight)
    model.fit(X_train, y_train)

    return model.predict_proba(X_valid)[:, 1], "logistic_regression"


def macro_group_accuracy(y, pred, groups):
    values = []

    for group in sorted(set(groups)):
        idx = np.asarray(groups) == group

        if np.any(idx):
            values.append(
                float(np.mean(
                    np.asarray(y)[idx] == np.asarray(pred)[idx]
                ))
            )

    return float(np.mean(values)) if values else float("nan")


def balanced_accuracy(y, pred):
    y = np.asarray(y, dtype=int)
    pred = np.asarray(pred, dtype=int)

    recalls = []

    for label in (0, 1):
        mask = y == label
        if not np.any(mask):
            return float("nan")
        recalls.append(float(np.mean(pred[mask] == label)))

    return float(np.mean(recalls))


def classification_metrics(y, pred, groups):
    y = np.asarray(y, dtype=int)
    pred = np.asarray(pred, dtype=int)

    tn = int(np.sum((y == 0) & (pred == 0)))
    fp = int(np.sum((y == 0) & (pred == 1)))
    fn = int(np.sum((y == 1) & (pred == 0)))
    tp = int(np.sum((y == 1) & (pred == 1)))

    return {
        "n": int(len(y)),
        "accuracy": float(np.mean(y == pred)),
        "balanced_accuracy": balanced_accuracy(y, pred),
        "macro_IF30_accuracy": macro_group_accuracy(
            y, pred, groups
        ),
        "negative_recall": tn / (tn + fp) if tn + fp else None,
        "positive_recall": tp / (tp + fn) if tp + fn else None,
        "tn": tn,
        "fp": fp,
        "fn": fn,
        "tp": tp,
    }


def inner_oof_probabilities(
    X,
    y,
    groups,
    feature_mask,
    C,
    class_weight
):
    unique_groups = sorted(set(groups))
    n_splits = min(INNER_FOLDS, len(unique_groups))

    if n_splits < 2:
        raise RuntimeError("Insufficient inner groups.")

    cv = GroupKFold(n_splits=n_splits)
    probabilities = np.full(len(y), np.nan)

    for train_idx, valid_idx in cv.split(X, y, groups):
        prob, _ = fit_probability_model(
            X[train_idx][:, feature_mask],
            y[train_idx],
            X[valid_idx][:, feature_mask],
            C,
            class_weight
        )
        probabilities[valid_idx] = prob

    if not np.isfinite(probabilities).all():
        raise RuntimeError("Incomplete inner OOF probabilities.")

    return probabilities


def choose_model_and_threshold(X, y, groups, masks):
    """
    Select feature set, C, class weight, and probability threshold
    using only inner group-held-out predictions.

    Primary objective: macro-averaged accuracy across IF30 groups.
    Secondary objective: pooled balanced accuracy.
    Tie breaks favor threshold 0.5, fewer features, and stronger
    regularization.
    """
    best = None

    for feature_set, C, class_weight in product(
        masks.keys(), C_VALUES, CLASS_WEIGHTS
    ):
        mask = masks[feature_set]

        probs = inner_oof_probabilities(
            X, y, groups, mask, C, class_weight
        )

        for threshold in PROBABILITY_THRESHOLDS:
            pred = (probs >= threshold).astype(int)

            macro_acc = macro_group_accuracy(
                y, pred, groups
            )
            pooled_ba = balanced_accuracy(y, pred)

            if not math.isfinite(macro_acc):
                continue

            if not math.isfinite(pooled_ba):
                continue

            key = (
                macro_acc,
                pooled_ba,
                -abs(float(threshold) - 0.5),
                -int(mask.sum()),
                -float(C),
                class_weight is None,
            )

            if best is None or key > best["key"]:
                best = {
                    "key": key,
                    "macro_accuracy": macro_acc,
                    "balanced_accuracy": pooled_ba,
                    "feature_set": feature_set,
                    "feature_count": int(mask.sum()),
                    "C": float(C),
                    "class_weight": class_weight,
                    "threshold": float(threshold),
                    "mask": mask,
                }

    if best is None:
        raise RuntimeError("No valid inner model configuration.")

    return best


def exact_group_swap(group_differences):
    """
    Exact paired sign-flip test on family-level accuracy differences.
    """
    differences = np.asarray(group_differences, dtype=float)
    observed = float(np.sum(differences))

    null = np.array([
        np.dot(np.asarray(signs), differences)
        for signs in product((-1.0, 1.0), repeat=len(differences))
    ])

    p = float(np.mean(
        np.abs(null) >= abs(observed) - 1e-12
    ))

    return {
        "test": "exact paired IF30-group sign-flip test",
        "n_groups": len(differences),
        "observed_sum_group_differences": observed,
        "two_sided_p": p,
    }


def family_bootstrap(group_rows, model_key, comparator_key):
    rng = np.random.default_rng(SEED)
    n_groups = len(group_rows)

    rates = np.array([
        (
            r[model_key] - r[comparator_key]
        ) / r["n"]
        for r in group_rows
    ])

    count_differences = np.array([
        r[model_key] - r[comparator_key]
        for r in group_rows
    ])

    sizes = np.array([r["n"] for r in group_rows])

    macro_draws = np.empty(BOOTSTRAPS)
    pooled_draws = np.empty(BOOTSTRAPS)

    for i in range(BOOTSTRAPS):
        idx = rng.integers(
            0, n_groups, size=n_groups
        )

        macro_draws[i] = np.mean(rates[idx])
        pooled_draws[i] = (
            np.sum(count_differences[idx])
            / np.sum(sizes[idx])
        )

    return {
        "replicates": BOOTSTRAPS,
        "seed": SEED,
        "macro_delta_95_ci": np.quantile(
            macro_draws, [0.025, 0.975]
        ).tolist(),
        "pooled_delta_95_ci": np.quantile(
            pooled_draws, [0.025, 0.975]
        ).tolist(),
        "note": "Resampling unit is the IF30 group."
    }


def mcnemar(y, pred_a, pred_b):
    correct_a = np.asarray(pred_a) == np.asarray(y)
    correct_b = np.asarray(pred_b) == np.asarray(y)

    a_only = int(np.sum(correct_a & ~correct_b))
    b_only = int(np.sum(~correct_a & correct_b))
    discordant = a_only + b_only

    p = (
        float(binomtest(
            min(a_only, b_only),
            n=discordant,
            p=0.5
        ).pvalue)
        if discordant else 1.0
    )

    return {
        "a_only_correct": a_only,
        "b_only_correct": b_only,
        "discordant_pairs": discordant,
        "exact_two_sided_p": p,
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    if not DB.exists():
        raise FileNotFoundError(DB)

    con = sqlite3.connect(
        f"file:{DB.resolve()}?mode=ro",
        uri=True
    )
    con.row_factory = sqlite3.Row

    records = con.execute("""
        SELECT
            f.feature_uid,
            f.cycle_uid,
            f.family_id,
            f.complex_id,
            f.feature_names_json,
            f.feature_values_json,
            f.experimental_epistasis_kcal_mol,
            p.heldout_group,
            p.split_unit,
            p.train_groups_json,
            p.y_pred,
            p.baseline_pred
        FROM cycle_features f
        JOIN loco_predictions p
          ON f.feature_uid = p.feature_uid
        ORDER BY f.family_id, f.feature_uid
    """).fetchall()

    con.close()

    if len(records) != 242:
        raise RuntimeError(
            f"Expected 242 rows, received {len(records)}"
        )

    feature_names = json.loads(
        records[0]["feature_names_json"]
    )

    X_all = []

    for r in records:
        if json.loads(r["feature_names_json"]) != feature_names:
            raise RuntimeError("Feature names differ between rows.")

        values = np.asarray(
            json.loads(r["feature_values_json"]),
            dtype=float
        )

        if len(values) != len(feature_names):
            raise RuntimeError("Feature vector length mismatch.")

        X_all.append(values)

    X_all = np.vstack(X_all)

    if not np.isfinite(X_all).all():
        raise RuntimeError("Non-finite feature values.")

    all_groups = sorted({
        r["family_id"] for r in records
    })

    # Verify that the original outer predictions are group isolated.
    for r in records:
        if r["heldout_group"] != r["family_id"]:
            raise RuntimeError(
                f"Held-out mismatch: {r['feature_uid']}"
            )

        if r["split_unit"] != "family":
            raise RuntimeError(
                f"Unexpected split unit: {r['feature_uid']}"
            )

        expected = set(all_groups) - {r["heldout_group"]}
        train_groups = set(json.loads(r["train_groups_json"]))

        if train_groups != expected:
            raise RuntimeError(
                f"Training-group isolation failure: "
                f"{r['feature_uid']}"
            )

    # Predefined informative subset.
    selected = [
        i for i, r in enumerate(records)
        if abs(float(r["experimental_epistasis_kcal_mol"]))
        >= THRESHOLD
    ]

    selected_records = [records[i] for i in selected]
    X = X_all[selected]

    y = np.array([
        int(
            float(r["experimental_epistasis_kcal_mol"]) > 0
        )
        for r in selected_records
    ], dtype=int)

    groups = np.array([
        r["family_id"] for r in selected_records
    ])

    existing_ridge_sign = np.array([
        int(float(r["y_pred"]) > 0)
        for r in selected_records
    ])

    outer_groups = sorted(set(groups))

    if len(outer_groups) != 16:
        raise RuntimeError(
            f"Expected 16 informative groups, got {len(outer_groups)}"
        )

    names = np.asarray(feature_names)

    masks = {
        "all_descriptors": np.ones(
            len(names), dtype=bool
        ),
        "physicochemical_and_context": np.array([
            not name.startswith("sub_")
            for name in names
        ]),
        "physicochemical_only": np.array([
            name.startswith(("sum_", "absdiff_", "product_"))
            for name in names
        ]),
    }

    nested_pred = np.full(len(y), -1, dtype=int)
    nested_prob = np.full(len(y), np.nan)
    chosen_threshold = np.full(len(y), np.nan)

    fold_rows = []

    print("PASS: 242 original family-held-out records")
    print("PASS: existing training-group isolation")
    print("Informative records:", len(y))
    print("Informative IF30 groups:", len(outer_groups))
    print("Positive / negative:", int(y.sum()), int(len(y) - y.sum()))

    for fold_number, group in enumerate(outer_groups, start=1):
        test_idx = np.where(groups == group)[0]
        train_idx = np.where(groups != group)[0]

        X_train = X[train_idx]
        y_train = y[train_idx]
        groups_train = groups[train_idx]
        X_test = X[test_idx]

        if len(np.unique(y_train)) < 2:
            majority = int(y_train[0])
            nested_pred[test_idx] = majority
            nested_prob[test_idx] = float(majority)
            chosen_threshold[test_idx] = 0.5

            fold_rows.append({
                "heldout_group": group,
                "n_train": len(train_idx),
                "n_test": len(test_idx),
                "status": "single_class_training_fallback",
                "selected_feature_set": "",
                "selected_threshold": 0.5,
            })
            continue

        best = choose_model_and_threshold(
            X_train,
            y_train,
            groups_train,
            masks
        )

        probability, status = fit_probability_model(
            X_train[:, best["mask"]],
            y_train,
            X_test[:, best["mask"]],
            best["C"],
            best["class_weight"]
        )

        pred = (
            probability >= best["threshold"]
        ).astype(int)

        nested_pred[test_idx] = pred
        nested_prob[test_idx] = probability
        chosen_threshold[test_idx] = best["threshold"]

        fold_rows.append({
            "heldout_group": group,
            "n_train": len(train_idx),
            "n_test": len(test_idx),
            "inner_macro_group_accuracy": best["macro_accuracy"],
            "inner_pooled_balanced_accuracy": best["balanced_accuracy"],
            "selected_feature_set": best["feature_set"],
            "selected_feature_count": best["feature_count"],
            "selected_C": best["C"],
            "selected_class_weight": (
                best["class_weight"]
                if best["class_weight"] is not None
                else "none"
            ),
            "selected_threshold": best["threshold"],
            "status": status,
        })

        print(
            f"[{fold_number:02d}/{len(outer_groups):02d}] "
            f"{group}: {best['feature_set']}, "
            f"C={best['C']}, "
            f"threshold={best['threshold']:.2f}, "
            f"inner_macro={best['macro_accuracy']:.3f}"
        )

    if np.any(nested_pred < 0):
        raise RuntimeError("Missing outer predictions.")

    if not np.isfinite(nested_prob).all():
        raise RuntimeError("Invalid outer probabilities.")

    # Recreate outer training-fold majority predictions.
    majority_pred = np.full(len(y), -1, dtype=int)

    for group in outer_groups:
        test_idx = np.where(groups == group)[0]
        train_idx = np.where(groups != group)[0]

        majority = int(np.mean(y[train_idx]) >= 0.5)
        majority_pred[test_idx] = majority

    if np.any(majority_pred < 0):
        raise RuntimeError("Missing majority baseline predictions.")

    methods = {
        "nested_threshold_classifier": nested_pred,
        "existing_ridge_sign": existing_ridge_sign,
        "outer_training_fold_majority": majority_pred,
    }

    method_metrics = {
        name: classification_metrics(y, pred, groups)
        for name, pred in methods.items()
    }

    group_rows = []

    for group in outer_groups:
        idx = np.where(groups == group)[0]

        row = {
            "IF30_group": group,
            "n": len(idx),
            "positive_n": int(y[idx].sum()),
            "negative_n": int(len(idx) - y[idx].sum()),
        }

        for name, pred in methods.items():
            row[name + "_correct"] = int(
                np.sum(pred[idx] == y[idx])
            )
            row[name + "_accuracy"] = float(
                np.mean(pred[idx] == y[idx])
            )

        group_rows.append(row)

    # Family-equal comparisons.
    group_deltas_vs_ridge = [
        r["nested_threshold_classifier_accuracy"]
        - r["existing_ridge_sign_accuracy"]
        for r in group_rows
    ]

    group_deltas_vs_majority = [
        r["nested_threshold_classifier_accuracy"]
        - r["outer_training_fold_majority_accuracy"]
        for r in group_rows
    ]

    summary = {
        "experiment": (
            "Nested IF30-group-held-out directional classifier "
            "with inner family-balanced threshold selection"
        ),
        "n_cycles": len(y),
        "n_informative_groups": len(outer_groups),
        "informative_threshold_abs_epsilon": THRESHOLD,
        "label_definition": (
            "1 when experimental epistasis > 0; 0 otherwise"
        ),
        "primary_inner_selection_metric": (
            "Macro-averaged accuracy across inner held-out IF30 groups"
        ),
        "methods": method_metrics,
        "classifier_vs_ridge": {
            "exact_family_block_test": exact_group_swap(
                group_deltas_vs_ridge
            ),
            "family_cluster_bootstrap": family_bootstrap(
                group_rows,
                "nested_threshold_classifier_correct",
                "existing_ridge_sign_correct"
            ),
            "paired_mcnemar": mcnemar(
                y, nested_pred, existing_ridge_sign
            ),
        },
        "classifier_vs_majority": {
            "exact_family_block_test": exact_group_swap(
                group_deltas_vs_majority
            ),
            "family_cluster_bootstrap": family_bootstrap(
                group_rows,
                "nested_threshold_classifier_correct",
                "outer_training_fold_majority_correct"
            ),
            "paired_mcnemar": mcnemar(
                y, nested_pred, majority_pred
            ),
        },
        "thresholds_selected": {
            str(r["heldout_group"]): r.get("selected_threshold")
            for r in fold_rows
        },
        "selection_method": (
            "Feature set, regularization, class weight, and threshold "
            "chosen exclusively from inner group-held-out predictions "
            "inside each outer training partition."
        ),
        "validation_scope": (
            "IF30 sequence-similarity groups, not CATH H-level superfamilies."
        ),
        "biological_scope": (
            "SKEMPI-derived PPI mutation-cycle epistasis; "
            "not HERV-K Env peptide-binding validation."
        ),
        "status": "exploratory; requires independent replication",
    }

    (OUT / "thresholded_directional_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8"
    )

    write_csv(
        OUT / "thresholded_directional_fold_choices.csv",
        fold_rows,
        list(fold_rows[0].keys())
    )

    prediction_rows = []

    for i, r in enumerate(selected_records):
        prediction_rows.append({
            "feature_uid": r["feature_uid"],
            "cycle_uid": r["cycle_uid"],
            "IF30_group": r["family_id"],
            "experimental_epistasis": r[
                "experimental_epistasis_kcal_mol"
            ],
            "true_sign": int(y[i]),
            "nested_prediction": int(nested_pred[i]),
            "nested_probability_positive": float(nested_prob[i]),
            "selected_threshold": float(chosen_threshold[i]),
            "existing_ridge_sign": int(existing_ridge_sign[i]),
            "outer_majority_prediction": int(majority_pred[i]),
        })

    write_csv(
        OUT / "thresholded_directional_oof_predictions.csv",
        prediction_rows,
        list(prediction_rows[0].keys())
    )

    write_csv(
        OUT / "thresholded_directional_family_metrics.csv",
        group_rows,
        list(group_rows[0].keys())
    )

    try:
        import matplotlib.pyplot as plt

        names_for_plot = [
            "Outer majority",
            "Existing Ridge sign",
            "Nested threshold classifier",
        ]

        keys = [
            "outer_training_fold_majority",
            "existing_ridge_sign",
            "nested_threshold_classifier",
        ]

        acc = [
            method_metrics[k]["accuracy"] * 100 for k in keys
        ]
        ba = [
            method_metrics[k]["balanced_accuracy"] * 100 for k in keys
        ]

        x = np.arange(len(keys))
        width = 0.36

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.bar(x - width / 2, acc, width, label="Accuracy")
        ax.bar(x + width / 2, ba, width, label="Balanced accuracy")
        ax.set_xticks(x)
        ax.set_xticklabels(names_for_plot)
        ax.set_ylim(0, 100)
        ax.set_ylabel("Score (%)")
        ax.set_title(
            "P4: Nested threshold-selected directional prediction"
        )
        ax.legend()
        ax.grid(axis="y", alpha=0.25)
        fig.tight_layout()
        fig.savefig(
            OUT / "p4_thresholded_directional_comparison.png",
            dpi=220
        )
        plt.close(fig)

        ordered = sorted(
            group_rows,
            key=lambda r: (
                r["nested_threshold_classifier_accuracy"]
                - r["existing_ridge_sign_accuracy"]
            )
        )

        fig, ax = plt.subplots(figsize=(9, 7))
        labels = [
            f"{r['IF30_group']} (n={r['n']})"
            for r in ordered
        ]
        deltas = [
            100 * (
                r["nested_threshold_classifier_accuracy"]
                - r["existing_ridge_sign_accuracy"]
            )
            for r in ordered
        ]

        ax.barh(labels, deltas)
        ax.axvline(0, linestyle="--")
        ax.set_xlabel(
            "Classifier minus Ridge sign accuracy (percentage points)"
        )
        ax.set_title(
            "P4: Threshold-selected classifier effect by held-out IF30 group"
        )
        ax.grid(axis="x", alpha=0.25)
        fig.tight_layout()
        fig.savefig(
            OUT / "p4_thresholded_directional_family_effects.png",
            dpi=220
        )
        plt.close(fig)

        print("\nFigures created:")
        print(OUT / "p4_thresholded_directional_comparison.png")
        print(OUT / "p4_thresholded_directional_family_effects.png")

    except ImportError:
        print("\nMatplotlib unavailable; tables and JSON were saved.")

    print("\n=== THRESHOLD-SELECTED RESULTS ===")
    print(json.dumps(summary, indent=2))

    print("\n=== OUTPUT DIRECTORY ===")
    print(OUT)


if __name__ == "__main__":
    main()
