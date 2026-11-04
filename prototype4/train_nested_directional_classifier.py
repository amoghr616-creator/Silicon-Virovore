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
OUT = ROOT / "p4_nested_directional_classifier_20261101"

SEED = 20261101
THRESHOLD = 0.5
C_VALUES = [0.001, 0.01, 0.1, 1.0, 10.0, 100.0]
CLASS_WEIGHTS = [None, "balanced"]
INNER_FOLDS = 4
BOOTSTRAPS = 20000


def write_csv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


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


def classification_metrics(y, pred):
    y = np.asarray(y, dtype=int)
    pred = np.asarray(pred, dtype=int)

    tn = int(np.sum((y == 0) & (pred == 0)))
    fp = int(np.sum((y == 0) & (pred == 1)))
    fn = int(np.sum((y == 1) & (pred == 0)))
    tp = int(np.sum((y == 1) & (pred == 1)))

    return {
        "n": len(y),
        "accuracy": float(np.mean(y == pred)),
        "balanced_accuracy": balanced_accuracy(y, pred),
        "negative_recall": tn / (tn + fp) if tn + fp else None,
        "positive_recall": tp / (tp + fn) if tp + fn else None,
        "tn": tn,
        "fp": fp,
        "fn": fn,
        "tp": tp,
    }


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


def safe_fit_predict(X_train, y_train, X_valid, C, class_weight):
    """
    Fit the classifier. If a training partition has only one class,
    use its training-partition majority class instead.
    """
    unique = np.unique(y_train)

    if len(unique) < 2:
        label = int(unique[0])
        return (
            np.full(len(X_valid), label, dtype=int),
            np.full(len(X_valid), float(label)),
            "single_class_training_fallback"
        )

    model = make_model(C, class_weight)
    model.fit(X_train, y_train)

    pred = model.predict(X_valid).astype(int)
    probability = model.predict_proba(X_valid)[:, 1]

    return pred, probability, "logistic_regression"


def exact_family_swap(group_deltas):
    """
    Exact paired sign-flip test on equally weighted family-level
    accuracy differences.
    """
    group_deltas = np.asarray(group_deltas, dtype=float)
    observed = float(np.sum(group_deltas))

    null_statistics = []

    for signs in product((-1.0, 1.0), repeat=len(group_deltas)):
        null_statistics.append(
            float(np.dot(np.asarray(signs), group_deltas))
        )

    null_statistics = np.asarray(null_statistics)

    p = float(np.mean(
        np.abs(null_statistics) >= abs(observed) - 1e-12
    ))

    return {
        "test": "exact IF30-group paired sign-flip test",
        "groups": len(group_deltas),
        "two_sided_p": p,
        "observed_sum_of_group_differences": observed,
    }


def family_bootstrap(group_rows, model_correct_key, baseline_correct_key):
    rng = np.random.default_rng(SEED)

    n_groups = len(group_rows)

    delta_rates = np.array([
        (r[model_correct_key] - r[baseline_correct_key])
        / r["n"]
        for r in group_rows
    ])

    delta_counts = np.array([
        r[model_correct_key] - r[baseline_correct_key]
        for r in group_rows
    ])

    group_sizes = np.array([
        r["n"] for r in group_rows
    ])

    macro_draws = np.empty(BOOTSTRAPS)
    pooled_draws = np.empty(BOOTSTRAPS)

    for i in range(BOOTSTRAPS):
        idx = rng.integers(0, n_groups, size=n_groups)

        macro_draws[i] = np.mean(delta_rates[idx])

        pooled_draws[i] = (
            np.sum(delta_counts[idx])
            / np.sum(group_sizes[idx])
        )

    return {
        "replicates": BOOTSTRAPS,
        "seed": SEED,
        "macro_difference_95_ci": np.quantile(
            macro_draws, [0.025, 0.975]
        ).tolist(),
        "pooled_difference_95_ci": np.quantile(
            pooled_draws, [0.025, 0.975]
        ).tolist(),
        "note": (
            "Bootstrap resamples whole IF30 groups, not individual cycles."
        ),
    }


def paired_mcnemar(y, pred_a, pred_b):
    a_correct = np.asarray(pred_a) == np.asarray(y)
    b_correct = np.asarray(pred_b) == np.asarray(y)

    a_only = int(np.sum(a_correct & ~b_correct))
    b_only = int(np.sum(~a_correct & b_correct))
    discordant = a_only + b_only

    p = (
        float(binomtest(
            min(a_only, b_only),
            n=discordant,
            p=0.5,
            alternative="two-sided"
        ).pvalue)
        if discordant else 1.0
    )

    return {
        "model_a_only_correct": a_only,
        "model_b_only_correct": b_only,
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
            p.y_true,
            p.y_pred,
            p.baseline_pred
        FROM cycle_features f
        JOIN loco_predictions p
          ON f.feature_uid = p.feature_uid
        ORDER BY f.family_id, f.feature_uid
    """).fetchall()

    con.close()

    if not records:
        raise RuntimeError("No joined cycle features and predictions found.")

    names = json.loads(records[0]["feature_names_json"])
    X_rows = []

    for r in records:
        current_names = json.loads(r["feature_names_json"])
        values = json.loads(r["feature_values_json"])

        if current_names != names:
            raise RuntimeError(
                f"Feature schema mismatch: {r['feature_uid']}"
            )

        if len(values) != len(names):
            raise RuntimeError(
                f"Feature count mismatch: {r['feature_uid']}"
            )

        X_rows.append(np.asarray(values, dtype=float))

    X_all = np.vstack(X_rows)

    if not np.isfinite(X_all).all():
        raise RuntimeError("Non-finite descriptor values found.")

    # Verify the existing family-held-out predictions.
    all_groups = sorted({
        r["family_id"] for r in records
    })

    integrity_errors = []

    for r in records:
        if r["family_id"] != r["heldout_group"]:
            integrity_errors.append(
                f"Family mismatch: {r['feature_uid']}"
            )

        if r["split_unit"] != "family":
            integrity_errors.append(
                f"Unexpected split unit: {r['feature_uid']}"
            )

        if not math.isclose(
            float(r["y_true"]),
            float(r["experimental_epistasis_kcal_mol"]),
            rel_tol=0.0,
            abs_tol=1e-8
        ):
            integrity_errors.append(
                f"Outcome mismatch: {r['feature_uid']}"
            )

        train_groups = set(json.loads(r["train_groups_json"]))
        expected = set(all_groups) - {r["heldout_group"]}

        if train_groups != expected:
            integrity_errors.append(
                f"Existing split mismatch: {r['feature_uid']}"
            )

    if integrity_errors:
        print("INTEGRITY FAILURES:")
        for error in integrity_errors[:30]:
            print(error)
        raise SystemExit("Stop: existing OOF integrity failed.")

    print("PASS: existing IF30 fold integrity")
    print("Feature dimensions:", len(names))
    print("All cycle rows:", len(records))

    # Define the informative endpoint before fitting models.
    informative_indices = [
        i for i, r in enumerate(records)
        if abs(float(r["experimental_epistasis_kcal_mol"])) >= THRESHOLD
    ]

    if not informative_indices:
        raise RuntimeError("No informative cycles.")

    X = X_all[informative_indices]
    selected_records = [records[i] for i in informative_indices]

    y = np.array([
        int(float(r["experimental_epistasis_kcal_mol"]) > 0)
        for r in selected_records
    ])

    groups = np.array([
        r["family_id"] for r in selected_records
    ])

    feature_names = np.asarray(names)

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

    for name, mask in masks.items():
        if not np.any(mask):
            raise RuntimeError(f"Empty feature group: {name}")

    print("Informative cycles:", len(y))
    print("Positive class:", int(np.sum(y == 1)))
    print("Negative class:", int(np.sum(y == 0)))
    print("Informative IF30 groups:", len(set(groups)))

    # Original regression predictions are retained only as a comparator.
    existing_ridge_sign = np.array([
        int(float(r["y_pred"]) > 0)
        for r in selected_records
    ])

    existing_mean_sign = np.array([
        int(float(r["baseline_pred"]) > 0)
        for r in selected_records
    ])

    nested_predictions = np.full(len(y), -1, dtype=int)
    nested_probabilities = np.full(len(y), np.nan)
    majority_predictions = np.full(len(y), -1, dtype=int)

    fold_rows = []
    parameter_grid = list(product(
        masks.keys(),
        C_VALUES,
        CLASS_WEIGHTS
    ))

    outer_groups = sorted(set(groups))

    for outer_number, outer_group in enumerate(outer_groups, start=1):
        test_idx = np.where(groups == outer_group)[0]
        train_idx = np.where(groups != outer_group)[0]

        X_train = X[train_idx]
        y_train = y[train_idx]
        groups_train = groups[train_idx]
        X_test = X[test_idx]

        if len(np.unique(y_train)) < 2:
            majority = int(y_train[0])
            nested_predictions[test_idx] = majority
            nested_probabilities[test_idx] = float(majority)
            majority_predictions[test_idx] = majority

            fold_rows.append({
                "heldout_group": outer_group,
                "n_train": len(train_idx),
                "n_test": len(test_idx),
                "inner_score": "",
                "feature_set": "single_class_fallback",
                "selected_C": "",
                "class_weight": "",
                "status": "single_class_outer_training",
            })
            continue

        majority = int(
            np.mean(y_train) >= 0.5
        )
        majority_predictions[test_idx] = majority

        n_inner = min(
            INNER_FOLDS,
            len(np.unique(groups_train))
        )

        if n_inner < 2:
            raise RuntimeError(
                f"Insufficient inner groups for {outer_group}"
            )

        inner_cv = GroupKFold(n_splits=n_inner)
        inner_splits = list(
            inner_cv.split(X_train, y_train, groups_train)
        )

        best = None

        # All feature-set and hyperparameter decisions are made
        # using only inner group-held-out predictions.
        for feature_set, C, class_weight in parameter_grid:
            mask = masks[feature_set]
            inner_oof = np.full(len(train_idx), -1, dtype=int)

            failed = False

            for inner_train, inner_valid in inner_splits:
                pred, _, _ = safe_fit_predict(
                    X_train[inner_train][:, mask],
                    y_train[inner_train],
                    X_train[inner_valid][:, mask],
                    C,
                    class_weight
                )

                inner_oof[inner_valid] = pred

            if np.any(inner_oof < 0):
                failed = True

            if failed:
                continue

            score = balanced_accuracy(y_train, inner_oof)

            # Prefer fewer features and stronger regularization when
            # inner scores tie.
            key = (
                -score,
                int(np.sum(mask)),
                C,
                class_weight is not None
            )

            if best is None or key < best["key"]:
                best = {
                    "key": key,
                    "score": float(score),
                    "feature_set": feature_set,
                    "C": C,
                    "class_weight": class_weight,
                    "mask": mask,
                }

        if best is None:
            raise RuntimeError(
                f"No valid inner model selected for {outer_group}"
            )

        pred, prob, fit_status = safe_fit_predict(
            X_train[:, best["mask"]],
            y_train,
            X_test[:, best["mask"]],
            best["C"],
            best["class_weight"]
        )

        nested_predictions[test_idx] = pred
        nested_probabilities[test_idx] = prob

        fold_rows.append({
            "heldout_group": outer_group,
            "n_train": len(train_idx),
            "n_test": len(test_idx),
            "inner_balanced_accuracy": best["score"],
            "feature_set": best["feature_set"],
            "feature_count": int(np.sum(best["mask"])),
            "selected_C": best["C"],
            "class_weight": (
                best["class_weight"]
                if best["class_weight"] is not None
                else "none"
            ),
            "status": fit_status,
        })

        print(
            f"[{outer_number:02d}/{len(outer_groups):02d}] "
            f"{outer_group}: "
            f"{best['feature_set']}, "
            f"C={best['C']}, "
            f"inner_BA={best['score']:.3f}"
        )

    if np.any(nested_predictions < 0):
        raise RuntimeError("Some outer predictions were not generated.")

    if not np.isfinite(nested_probabilities).all():
        raise RuntimeError("Some outer probabilities are invalid.")

    metrics = {
        "nested_directional_classifier":
            classification_metrics(y, nested_predictions),
        "outer_training_fold_majority":
            classification_metrics(y, majority_predictions),
        "existing_ridge_regression_sign":
            classification_metrics(y, existing_ridge_sign),
        "existing_training_fold_mean_sign":
            classification_metrics(y, existing_mean_sign),
    }

    family_rows = []

    for group in outer_groups:
        idx = np.where(groups == group)[0]

        family_rows.append({
            "IF30_group": group,
            "n": len(idx),
            "positive_cycles": int(np.sum(y[idx] == 1)),
            "negative_cycles": int(np.sum(y[idx] == 0)),
            "classifier_correct": int(
                np.sum(nested_predictions[idx] == y[idx])
            ),
            "majority_correct": int(
                np.sum(majority_predictions[idx] == y[idx])
            ),
            "ridge_sign_correct": int(
                np.sum(existing_ridge_sign[idx] == y[idx])
            ),
            "classifier_accuracy": float(
                np.mean(nested_predictions[idx] == y[idx])
            ),
            "majority_accuracy": float(
                np.mean(majority_predictions[idx] == y[idx])
            ),
            "ridge_sign_accuracy": float(
                np.mean(existing_ridge_sign[idx] == y[idx])
            ),
        })

    # Family-level comparison against the majority baseline.
    family_deltas = [
        r["classifier_accuracy"] - r["majority_accuracy"]
        for r in family_rows
    ]

    family_deltas_ridge = [
        r["classifier_accuracy"] - r["ridge_sign_accuracy"]
        for r in family_rows
    ]

    summary = {
        "study": (
            "Exploratory nested IF30-group-held-out classification "
            "of experimental epistasis sign"
        ),
        "endpoint": (
            "Positive vs negative experimental epistasis; "
            f"only cycles with |epsilon| >= {THRESHOLD} kcal/mol"
        ),
        "n_informative_cycles": len(y),
        "n_informative_groups": len(outer_groups),
        "n_total_groups_in_source": len(all_groups),
        "feature_names_count": len(names),
        "feature_set_sizes": {
            key: int(mask.sum())
            for key, mask in masks.items()
        },
        "metrics": metrics,
        "family_block_test_classifier_vs_majority":
            exact_family_swap(family_deltas),
        "family_block_test_classifier_vs_existing_ridge_sign":
            exact_family_swap(family_deltas_ridge),
        "family_cluster_bootstrap_classifier_vs_majority":
            family_bootstrap(
                family_rows,
                "classifier_correct",
                "majority_correct"
            ),
        "family_cluster_bootstrap_classifier_vs_ridge_sign":
            family_bootstrap(
                family_rows,
                "classifier_correct",
                "ridge_sign_correct"
            ),
        "mcnemar_classifier_vs_majority":
            paired_mcnemar(
                y,
                nested_predictions,
                majority_predictions
            ),
        "mcnemar_classifier_vs_existing_ridge_sign":
            paired_mcnemar(
                y,
                nested_predictions,
                existing_ridge_sign
            ),
        "feature_and_hyperparameter_selection":
            "Nested GroupKFold within each outer training partition",
        "outer_split":
            "Leave-one-IF30-group-out",
        "cath_H_superfamily_LOFO": False,
        "biological_scope": (
            "SKEMPI-derived protein-protein interaction mutation cycles; "
            "not HERV-K Env peptide-binding validation"
        ),
        "interpretation_warning": (
            "Exploratory model development on the same benchmark. "
            "Requires independent replication before confirmatory claims."
        )
    }

    (OUT / "nested_directional_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8"
    )

    write_csv(
        OUT / "nested_directional_fold_choices.csv",
        fold_rows,
        [
            "heldout_group", "n_train", "n_test",
            "inner_balanced_accuracy", "feature_set",
            "feature_count", "selected_C", "class_weight",
            "status"
        ]
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
            "nested_classifier_prediction": int(
                nested_predictions[i]
            ),
            "nested_classifier_probability_positive": float(
                nested_probabilities[i]
            ),
            "training_fold_majority_prediction": int(
                majority_predictions[i]
            ),
            "existing_ridge_sign_prediction": int(
                existing_ridge_sign[i]
            ),
            "existing_mean_sign_prediction": int(
                existing_mean_sign[i]
            ),
        })

    write_csv(
        OUT / "nested_directional_oof_predictions.csv",
        prediction_rows,
        list(prediction_rows[0].keys())
    )

    write_csv(
        OUT / "nested_directional_family_metrics.csv",
        family_rows,
        list(family_rows[0].keys())
    )

    try:
        import matplotlib.pyplot as plt

        model_names = [
            "Training-fold\nmajority",
            "Existing Ridge\nsign",
            "Nested sign\nclassifier"
        ]

        balanced_values = [
            metrics["outer_training_fold_majority"]["balanced_accuracy"],
            metrics["existing_ridge_regression_sign"]["balanced_accuracy"],
            metrics["nested_directional_classifier"]["balanced_accuracy"],
        ]

        accuracy_values = [
            metrics["outer_training_fold_majority"]["accuracy"],
            metrics["existing_ridge_regression_sign"]["accuracy"],
            metrics["nested_directional_classifier"]["accuracy"],
        ]

        fig, ax = plt.subplots(figsize=(8, 5))
        x = np.arange(len(model_names))
        width = 0.36

        ax.bar(
            x - width / 2,
            np.array(accuracy_values) * 100,
            width,
            label="Accuracy"
        )
        ax.bar(
            x + width / 2,
            np.array(balanced_values) * 100,
            width,
            label="Balanced accuracy"
        )

        ax.set_xticks(x)
        ax.set_xticklabels(model_names)
        ax.set_ylabel("Score (%)")
        ax.set_ylim(0, 100)
        ax.set_title(
            "P4: Nested family-held-out epistasis-sign prediction"
        )
        ax.legend()
        ax.grid(axis="y", alpha=0.25)
        fig.tight_layout()
        fig.savefig(
            OUT / "p4_nested_directional_model_comparison.png",
            dpi=220
        )
        plt.close(fig)

        ordered = sorted(
            family_rows,
            key=lambda r: (
                r["classifier_accuracy"] - r["majority_accuracy"]
            )
        )

        fig, ax = plt.subplots(figsize=(9, 7))
        labels = [
            f"{r['IF30_group']} (n={r['n']})"
            for r in ordered
        ]
        deltas = [
            100 * (
                r["classifier_accuracy"] - r["majority_accuracy"]
            )
            for r in ordered
        ]

        ax.barh(labels, deltas)
        ax.axvline(0, linestyle="--")
        ax.set_xlabel(
            "Classifier minus majority accuracy (percentage points)"
        )
        ax.set_title(
            "P4: Directional classifier performance by held-out IF30 group"
        )
        ax.grid(axis="x", alpha=0.25)
        fig.tight_layout()
        fig.savefig(
            OUT / "p4_nested_directional_family_effects.png",
            dpi=220
        )
        plt.close(fig)

        print("\nFigures created:")
        print(OUT / "p4_nested_directional_model_comparison.png")
        print(OUT / "p4_nested_directional_family_effects.png")

    except ImportError:
        print("\nMatplotlib unavailable; numerical outputs were still saved.")

    print("\n=== NESTED DIRECTIONAL RESULTS ===")
    print(json.dumps(summary, indent=2))

    print("\n=== OUTPUT DIRECTORY ===")
    print(OUT)


if __name__ == "__main__":
    main()
