from pathlib import Path
import importlib.util
import sqlite3
import hashlib
import json
import csv

import numpy as np

PROJECT = Path.cwd()
P4_SOURCE = PROJECT / "prototype4/p4.py"
DB = PROJECT / (
    "prototype4/results/model_a_validation_20261030/"
    "p4_model_a_working.sqlite"
)

OUT = PROJECT / (
    "prototype4/results/model_a_validation_20261030/"
    "p4_split_sensitivity_20261101"
)

EXPECTED_SOURCE_SHA256 = (
    "4ccb1f37b4be33b23636eab45584021b1df491178f7ac23b62f2609b7cca1181"
)

N_REPETITIONS = 100
TEST_FRACTION = 0.20
BOOTSTRAP_MODELS = 30
RIDGE_ALPHA = 10.0
SEED = 20261101


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"No rows to write: {path}")

    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=list(rows[0].keys()),
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def load_original_p4():
    digest = sha256(P4_SOURCE)

    if digest != EXPECTED_SOURCE_SHA256:
        raise RuntimeError(
            "p4.py SHA-256 differs from the validated original source. "
            f"Found {digest}"
        )

    spec = importlib.util.spec_from_file_location(
        "original_p4_model",
        P4_SOURCE
    )

    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load original p4.py")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    return module, digest


def load_data():
    con = sqlite3.connect(
        f"file:{DB.resolve()}?mode=ro",
        uri=True
    )
    con.row_factory = sqlite3.Row

    rows = con.execute("""
        SELECT
            f.feature_uid,
            f.cycle_uid,
            f.complex_id,
            f.family_id,
            f.feature_names_json,
            f.feature_values_json,
            f.experimental_epistasis_kcal_mol,
            p.heldout_group,
            p.split_unit,
            p.y_true,
            p.y_pred,
            p.baseline_pred
        FROM cycle_features AS f
        JOIN loco_predictions AS p
          ON p.feature_uid = f.feature_uid
        WHERE p.split_unit = 'family'
        ORDER BY f.family_id, f.feature_uid
    """).fetchall()

    con.close()

    if len(rows) != 242:
        raise RuntimeError(
            f"Expected 242 joined IF30 records; found {len(rows)}"
        )

    names = json.loads(rows[0]["feature_names_json"])
    X = []

    for r in rows:
        if r["family_id"] != r["heldout_group"]:
            raise RuntimeError(
                f"Family mismatch in existing OOF record: {r['feature_uid']}"
            )

        if json.loads(r["feature_names_json"]) != names:
            raise RuntimeError("Feature schema varies across records.")

        values = np.asarray(
            json.loads(r["feature_values_json"]),
            dtype=float
        )

        if len(values) != len(names):
            raise RuntimeError("Descriptor-vector length mismatch.")

        if not np.isfinite(values).all():
            raise RuntimeError("Non-finite descriptor value.")

        if not np.isclose(
            float(r["experimental_epistasis_kcal_mol"]),
            float(r["y_true"]),
            atol=1e-8,
            rtol=0
        ):
            raise RuntimeError(
                f"Experimental endpoint mismatch: {r['feature_uid']}"
            )

        X.append(values)

    return rows, names, np.vstack(X)


def within_family_split(groups, rng):
    """
    Hold out approximately 20% of cycles within each IF30 group.

    Every group with >=2 cycles retains at least one training cycle.
    Singleton groups remain in training and are not tested.
    """
    train_parts = []
    test_parts = []

    for group in sorted(set(groups)):
        idx = np.flatnonzero(groups == group)
        shuffled = rng.permutation(idx)

        if len(idx) < 2:
            train_parts.extend(shuffled.tolist())
            continue

        n_test = max(
            1,
            int(np.floor(TEST_FRACTION * len(idx)))
        )
        n_test = min(n_test, len(idx) - 1)

        test_parts.extend(shuffled[:n_test].tolist())
        train_parts.extend(shuffled[n_test:].tolist())

    train_idx = np.asarray(train_parts, dtype=int)
    test_idx = np.asarray(test_parts, dtype=int)

    rng.shuffle(train_idx)
    rng.shuffle(test_idx)

    if len(set(train_idx) & set(test_idx)):
        raise AssertionError("Train/test cycle overlap.")

    train_groups = set(groups[train_idx])
    test_groups = set(groups[test_idx])

    return train_idx, test_idx, train_groups, test_groups


def summarize_group_mae(y, pred, groups):
    values = []

    for group in sorted(set(groups)):
        idx = groups == group

        if np.any(idx):
            values.append(float(np.mean(
                np.abs(y[idx] - pred[idx])
            )))

    return {
        "macro_mae": float(np.mean(values)),
        "groups_with_test_cycles": len(values),
    }


def metric_values(p4, y, pred, groups):
    metrics = p4.metric_summary(y, pred)
    group_metrics = summarize_group_mae(y, pred, groups)

    return {
        "n_test": int(len(y)),
        "mae": float(metrics["mae_kcal_mol"]),
        "rmse": float(metrics["rmse_kcal_mol"]),
        "spearman_rho": (
            float(metrics["spearman_rho"])
            if metrics["spearman_rho"] is not None
            else None
        ),
        "direction_accuracy": (
            float(metrics["direction_accuracy_informative"])
            if metrics["direction_accuracy_informative"] is not None
            else None
        ),
        "n_informative": int(
            metrics["n_informative_abs_epsilon_ge_0_5"]
        ),
        "macro_group_mae": group_metrics["macro_mae"],
        "groups_with_test_cycles":
            group_metrics["groups_with_test_cycles"],
    }


def distribution(values):
    arr = np.asarray(
        [v for v in values if v is not None and np.isfinite(v)],
        dtype=float
    )

    if not len(arr):
        return None

    return {
        "n_splits": int(len(arr)),
        "mean": float(np.mean(arr)),
        "median": float(np.median(arr)),
        "sd": float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0,
        "q025": float(np.quantile(arr, 0.025)),
        "q975": float(np.quantile(arr, 0.975)),
        "interpretation": (
            "Empirical split-to-split distribution, not a confidence interval."
        ),
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    if not DB.exists():
        raise FileNotFoundError(DB)

    p4, source_hash = load_original_p4()
    records, feature_names, X = load_data()

    y = np.asarray([
        float(r["experimental_epistasis_kcal_mol"])
        for r in records
    ])

    groups = np.asarray([
        str(r["family_id"]) for r in records
    ])

    feature_ids = [r["feature_uid"] for r in records]
    if len(feature_ids) != len(set(feature_ids)):
        raise RuntimeError("Duplicate feature IDs.")

    print("PASS: exact original p4.py source hash")
    print("PASS: 242 source and OOF rows joined")
    print("PASS: 410-dimensional descriptors:", len(feature_names))
    print("PASS: existing family-held-out endpoints match")
    print("IF30 groups:", len(set(groups)))

    # Recompute existing family-held-out reference metrics.
    lofo_pred = np.asarray([
        float(r["y_pred"]) for r in records
    ])
    lofo_base = np.asarray([
        float(r["baseline_pred"]) for r in records
    ])

    lofo_model = metric_values(p4, y, lofo_pred, groups)
    lofo_baseline = metric_values(p4, y, lofo_base, groups)

    lofo_group_rows = []

    for group in sorted(set(groups)):
        idx = groups == group

        lofo_group_rows.append({
            "IF30_group": group,
            "n": int(idx.sum()),
            "model_mae": float(np.mean(
                np.abs(y[idx] - lofo_pred[idx])
            )),
            "baseline_mae": float(np.mean(
                np.abs(y[idx] - lofo_base[idx])
            )),
        })

    write_csv(
        OUT / "existing_if30_lofo_reference_by_group.csv",
        lofo_group_rows
    )

    random_rows = []
    prediction_rows = []

    for repetition in range(N_REPETITIONS):
        split_rng = np.random.default_rng(SEED + repetition)

        train_idx, test_idx, train_groups, test_groups = (
            within_family_split(groups, split_rng)
        )

        train_group_array = groups[train_idx]

        # Every tested group must also be represented in training.
        overlap_groups = train_groups & test_groups

        if overlap_groups != test_groups:
            raise RuntimeError(
                "This is not a within-family random split."
            )

        if len(np.unique(train_group_array)) < 3:
            raise RuntimeError(
                "Insufficient training groups for nested calibration."
            )

        # The exact original P4 prediction/calibration function is used.
        # No hyperparameter or feature selection is performed.
        pred, raw_u, calibrated, status, calibration_n = (
            p4.nested_fold_predictions(
                X[train_idx],
                y[train_idx],
                train_group_array,
                X[test_idx],
                n_boot=BOOTSTRAP_MODELS,
                seed=SEED + 100003 + repetition * 7919,
                alpha=RIDGE_ALPHA
            )
        )

        baseline = np.full(
            len(test_idx),
            float(np.mean(y[train_idx]))
        )

        y_test = y[test_idx]
        test_family_ids = groups[test_idx]

        model_metrics = metric_values(
            p4, y_test, pred, test_family_ids
        )
        baseline_metrics = metric_values(
            p4, y_test, baseline, test_family_ids
        )

        random_rows.append({
            "repetition": repetition + 1,
            "seed": SEED + repetition,
            "n_train": int(len(train_idx)),
            "n_test": int(len(test_idx)),
            "n_train_groups": int(len(train_groups)),
            "n_test_groups": int(len(test_groups)),
            "n_overlapping_groups": int(len(overlap_groups)),
            "test_cycle_family_overlap_fraction": float(
                sum(
                    groups[i] in train_groups
                    for i in test_idx
                ) / len(test_idx)
            ),
            "calibration_status": status,
            "inner_calibration_n": int(calibration_n),
            **{
                "model_" + k: v
                for k, v in model_metrics.items()
            },
            **{
                "baseline_" + k: v
                for k, v in baseline_metrics.items()
            },
            "model_minus_baseline_mae": (
                model_metrics["mae"] - baseline_metrics["mae"]
            ),
            "model_minus_baseline_macro_group_mae": (
                model_metrics["macro_group_mae"]
                - baseline_metrics["macro_group_mae"]
            ),
        })

        for local_i, global_i in enumerate(test_idx):
            prediction_rows.append({
                "repetition": repetition + 1,
                "feature_uid": records[global_i]["feature_uid"],
                "cycle_uid": records[global_i]["cycle_uid"],
                "complex_id": records[global_i]["complex_id"],
                "IF30_group": groups[global_i],
                "y_true": float(y_test[local_i]),
                "y_pred": float(pred[local_i]),
                "baseline_pred": float(baseline[local_i]),
                "raw_uncertainty": float(raw_u[local_i]),
                "calibrated_abs_error": float(calibrated[local_i]),
            })

        if (repetition + 1) % 10 == 0:
            print(
                f"Completed {repetition + 1}/{N_REPETITIONS} "
                f"random splits"
            )

    write_csv(
        OUT / "within_family_random_split_metrics.csv",
        random_rows
    )

    write_csv(
        OUT / "within_family_random_split_predictions.csv",
        prediction_rows
    )

    keys = [
        "model_mae",
        "model_rmse",
        "model_spearman_rho",
        "model_direction_accuracy",
        "model_macro_group_mae",
        "baseline_mae",
        "baseline_rmse",
        "baseline_spearman_rho",
        "baseline_direction_accuracy",
        "baseline_macro_group_mae",
        "model_minus_baseline_mae",
        "model_minus_baseline_macro_group_mae",
    ]

    random_distributions = {
        key: distribution([r[key] for r in random_rows])
        for key in keys
    }

    # Descriptive random-versus-LOFO differences.
    lofo_values = {
        "model_mae": lofo_model["mae"],
        "model_rmse": lofo_model["rmse"],
        "model_spearman_rho": lofo_model["spearman_rho"],
        "model_direction_accuracy": lofo_model["direction_accuracy"],
        "model_macro_group_mae": lofo_model["macro_group_mae"],
        "baseline_mae": lofo_baseline["mae"],
        "baseline_rmse": lofo_baseline["rmse"],
        "baseline_spearman_rho": lofo_baseline["spearman_rho"],
        "baseline_direction_accuracy": lofo_baseline["direction_accuracy"],
        "baseline_macro_group_mae": lofo_baseline["macro_group_mae"],
    }

    median_gaps = {}

    for key, lofo_value in lofo_values.items():
        random_key = key

        values = [
            r[random_key] for r in random_rows
            if r[random_key] is not None
        ]

        median_gaps[key] = {
            "random_median": (
                float(np.median(values)) if values else None
            ),
            "if30_lofo_reference": lofo_value,
            "random_median_minus_lofo": (
                float(np.median(values) - lofo_value)
                if values and lofo_value is not None else None
            ),
        }

    summary = {
        "experiment": (
            "Repeated within-IF30-group random cycle split "
            "versus existing leave-one-IF30-group-out evaluation"
        ),
        "source_p4_sha256": source_hash,
        "n_cycles": len(y),
        "n_IF30_groups": int(len(set(groups))),
        "descriptor_count": int(X.shape[1]),
        "random_split_repetitions": N_REPETITIONS,
        "random_split_test_fraction_target": TEST_FRACTION,
        "random_split_design": (
            "Within each group, approximately 20% of cycles are randomly "
            "held out; at least one cycle per group remains in training. "
            "Singleton groups remain in training."
        ),
        "model": (
            "Original p4.nested_fold_predictions; group-bootstrap Ridge "
            f"alpha={RIDGE_ALPHA}; {BOOTSTRAP_MODELS} bootstrap fits per model"
        ),
        "random_split": {
            "metric_distributions": random_distributions,
        },
        "existing_if30_lofo": {
            "ridge": lofo_model,
            "training_fold_mean": lofo_baseline,
            "n_test_cycles": len(y),
            "n_groups": int(len(set(groups))),
        },
        "random_median_minus_if30_lofo": median_gaps,
        "interpretation_limits": [
            "Random split tests interpolation within represented IF30 groups.",
            "IF30 LOFO tests transfer to an entirely held-out IF30 group.",
            "Different test-set sizes and compositions limit direct equivalence.",
            "Repeated random splits reuse the same dataset and are not independent replications.",
            "A performance gap indicates split-strategy sensitivity; it does not by itself prove leakage.",
            "This SKEMPI-derived PPI-cycle benchmark is not HERV-K Env peptide-binding validation."
        ],
    }

    (OUT / "split_sensitivity_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8"
    )

    # Separate figures for each metric.
    try:
        import matplotlib.pyplot as plt

        plots = [
            ("model_mae", "MAE (kcal/mol)", "Lower is better"),
            ("model_rmse", "RMSE (kcal/mol)", "Lower is better"),
            ("model_spearman_rho", "Spearman correlation", "Higher is better"),
            ("model_direction_accuracy", "Directional accuracy", "Higher is better"),
            ("model_macro_group_mae", "Macro IF30-group MAE (kcal/mol)", "Lower is better"),
        ]

        for key, ylabel, note in plots:
            values = [
                r[key] for r in random_rows
                if r[key] is not None and np.isfinite(r[key])
            ]

            if not values:
                continue

            lofo_key = key
            lofo_value = lofo_values[lofo_key]

            fig, ax = plt.subplots(figsize=(7, 5))
            ax.boxplot(
                values,
                tick_labels=["Within-family random split"],
                showmeans=True
            )

            if lofo_value is not None:
                ax.axhline(
                    lofo_value,
                    linestyle="--",
                    label=f"IF30 LOFO = {lofo_value:.3f}"
                )
                ax.legend()

            ax.set_ylabel(ylabel)
            ax.set_title(
                "P4 split sensitivity: " + ylabel + "\n" + note
            )
            ax.grid(axis="y", alpha=0.25)
            fig.tight_layout()

            filename = "p4_split_sensitivity_" + key + ".png"
            fig.savefig(OUT / filename, dpi=220)
            plt.close(fig)

        print("\nSplit-sensitivity figures saved.")

    except ImportError:
        print("Matplotlib unavailable; CSV and JSON still saved.")

    print("\n=== WITHIN-FAMILY RANDOM SPLIT MEDIANS ===")
    print(json.dumps(random_distributions, indent=2))

    print("\n=== EXISTING IF30 LOFO REFERENCE ===")
    print(json.dumps({
        "ridge": lofo_model,
        "training_fold_mean": lofo_baseline,
    }, indent=2))

    print("\n=== RANDOM MEDIAN MINUS LOFO ===")
    print(json.dumps(median_gaps, indent=2))

    print("\nOUTPUT:")
    print(OUT)


if __name__ == "__main__":
    main()
