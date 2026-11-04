from pathlib import Path
from collections import defaultdict
import importlib.util
import sqlite3
import hashlib
import json
import csv

import numpy as np
from itertools import product
from scipy.stats import binomtest

ROOT = Path("prototype4/results/model_a_validation_20261030")
DB = ROOT / "p4_model_a_working.sqlite"
OUT = ROOT / "p4_equal_training_size_20261101"

EXPECTED_HASH = (
    "4ccb1f37b4be33b23636eab45584021b1df491178f7ac23b62f2609b7cca1181"
)

TRAIN_N = 150
REPEATS = 10
TEST_FRACTION = 0.20
BOOTSTRAP_MODELS = 30
ALPHA = 10.0
SEED = 20261102
BOOTSTRAPS = 20000


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise RuntimeError(f"No rows to write: {path}")

    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=list(rows[0].keys()),
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def exact_sign_flip(values):
    values = np.asarray(values, dtype=float)
    observed = float(np.sum(values))

    null = np.fromiter(
        (
            float(np.dot(np.asarray(signs), values))
            for signs in product((-1.0, 1.0), repeat=len(values))
        ),
        dtype=float
    )

    return {
        "test": "exact paired IF30-family sign-flip permutation",
        "n_families": len(values),
        "observed_sum": observed,
        "two_sided_p": float(np.mean(
            np.abs(null) >= abs(observed) - 1e-12
        ))
    }


def cluster_bootstrap(values):
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(SEED)

    draws = np.empty(BOOTSTRAPS)

    for i in range(BOOTSTRAPS):
        idx = rng.integers(0, len(values), size=len(values))
        draws[i] = np.mean(values[idx])

    return {
        "replicates": BOOTSTRAPS,
        "seed": SEED,
        "mean": float(np.mean(values)),
        "ci95": np.quantile(draws, [0.025, 0.975]).tolist(),
        "resampling_unit": "IF30 family"
    }


def sign_test(values):
    values = np.asarray(values)
    positive = int(np.sum(values > 0))
    negative = int(np.sum(values < 0))
    ties = int(np.sum(values == 0))

    p = (
        float(binomtest(
            min(positive, negative),
            n=positive + negative,
            p=0.5
        ).pvalue)
        if positive + negative else 1.0
    )

    return {
        "positive": positive,
        "negative": negative,
        "ties": ties,
        "two_sided_p": p
    }


def metric_gaps(family_records):
    output = []

    for family, cycles in sorted(family_records.items()):
        random_mae = []
        excluded_mae = []
        random_mse = []
        excluded_mse = []
        random_base_mae = []
        excluded_base_mae = []
        random_direction = []
        excluded_direction = []

        for cycle in cycles.values():
            y = cycle["y_true"]

            random_mae.append(np.mean(cycle["random_abs"]))
            excluded_mae.append(np.mean(cycle["excluded_abs"]))

            random_mse.append(np.mean(cycle["random_sq"]))
            excluded_mse.append(np.mean(cycle["excluded_sq"]))

            random_base_mae.append(
                np.mean(cycle["random_base_abs"])
            )
            excluded_base_mae.append(
                np.mean(cycle["excluded_base_abs"])
            )

            if abs(y) >= 0.5:
                random_direction.append(
                    np.mean(cycle["random_direction"])
                )
                excluded_direction.append(
                    np.mean(cycle["excluded_direction"])
                )

        r_mae = float(np.mean(random_mae))
        e_mae = float(np.mean(excluded_mae))

        r_rmse = float(np.sqrt(np.mean(random_mse)))
        e_rmse = float(np.sqrt(np.mean(excluded_mse)))

        r_base = float(np.mean(random_base_mae))
        e_base = float(np.mean(excluded_base_mae))

        mae_gap = e_mae - r_mae
        base_gap = e_base - r_base

        rmse_base_random = float(np.sqrt(
            np.mean([
                np.mean(c["random_base_sq"])
                for c in cycles.values()
            ])
        ))
        rmse_base_excluded = float(np.sqrt(
            np.mean([
                np.mean(c["excluded_base_sq"])
                for c in cycles.values()
            ])
        ))

        rmse_gap = e_rmse - r_rmse
        baseline_rmse_gap = rmse_base_excluded - rmse_base_random

        output.append({
            "IF30_group": family,
            "n_test_cycles": len(cycles),
            "random_Ridge_MAE": r_mae,
            "excluded_Ridge_MAE": e_mae,
            "Ridge_MAE_gap": mae_gap,
            "random_baseline_MAE": r_base,
            "excluded_baseline_MAE": e_base,
            "baseline_MAE_gap": base_gap,
            "excess_MAE_gap": mae_gap - base_gap,
            "random_Ridge_RMSE": r_rmse,
            "excluded_Ridge_RMSE": e_rmse,
            "Ridge_RMSE_gap": rmse_gap,
            "baseline_RMSE_gap": baseline_rmse_gap,
            "excess_RMSE_gap": rmse_gap - baseline_rmse_gap,
            "random_direction_accuracy": (
                float(np.mean(random_direction))
                if random_direction else None
            ),
            "excluded_direction_accuracy": (
                float(np.mean(excluded_direction))
                if excluded_direction else None
            )
        })

    return output


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    p4_path = Path("prototype4/p4.py")
    source_hash = sha256(p4_path)

    if source_hash != EXPECTED_HASH:
        raise RuntimeError(
            "Original p4.py SHA-256 changed. Stop and inspect."
        )

    spec = importlib.util.spec_from_file_location(
        "original_p4_equal_size",
        p4_path
    )

    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load original p4.py")

    p4 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(p4)

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
            f.experimental_epistasis_kcal_mol
        FROM cycle_features f
        ORDER BY f.family_id, f.feature_uid
    """).fetchall()

    con.close()

    if len(records) != 242:
        raise RuntimeError(
            f"Expected 242 cycles; found {len(records)}"
        )

    names = json.loads(records[0]["feature_names_json"])
    X = []

    for r in records:
        if json.loads(r["feature_names_json"]) != names:
            raise RuntimeError("Feature schema mismatch.")

        values = np.asarray(
            json.loads(r["feature_values_json"]),
            dtype=float
        )

        if len(values) != len(names):
            raise RuntimeError("Descriptor dimension mismatch.")

        X.append(values)

    X = np.vstack(X)

    if not np.isfinite(X).all():
        raise RuntimeError("Non-finite feature values.")

    y = np.asarray([
        float(r["experimental_epistasis_kcal_mol"])
        for r in records
    ])

    groups = np.asarray([
        str(r["family_id"]) for r in records
    ])

    ids = np.asarray([
        str(r["feature_uid"]) for r in records
    ])

    unique_groups = sorted(set(groups))

    if len(unique_groups) != 19:
        raise RuntimeError("Expected 19 IF30 groups.")

    largest = max(
        int(np.sum(groups == group))
        for group in unique_groups
    )

    available_min = len(records) - largest

    if TRAIN_N > available_min:
        raise RuntimeError(
            f"TRAIN_N={TRAIN_N} exceeds smallest leave-one-family-out "
            f"training pool ({available_min})."
        )

    # Only groups with at least two cycles can have a within-family split.
    test_groups = [
        group for group in unique_groups
        if np.sum(groups == group) >= 2
    ]

    if len(test_groups) != 14:
        raise RuntimeError(
            f"Expected 14 testable groups, got {len(test_groups)}"
        )

    cycle_records = defaultdict(dict)
    split_manifest = []
    rng = np.random.default_rng(SEED)

    for repeat in range(REPEATS):
        for group_index, group in enumerate(test_groups):
            family_idx = np.flatnonzero(groups == group)

            shuffled_family = rng.permutation(family_idx)

            n_test = max(
                1,
                int(np.floor(TEST_FRACTION * len(family_idx)))
            )
            n_test = min(n_test, len(family_idx) - 1)

            test_idx = shuffled_family[:n_test]
            same_family_train = shuffled_family[n_test:]

            outside_idx = np.flatnonzero(groups != group)

            if len(outside_idx) < TRAIN_N:
                raise RuntimeError(
                    f"Insufficient outside-family training data: {group}"
                )

            # Draw one shared external pool. The two arms differ only
            # in whether same-family examples replace external examples.
            outside_pool = rng.choice(
                outside_idx,
                size=TRAIN_N,
                replace=False
            )

            n_same = len(same_family_train)
            n_outside_random = TRAIN_N - n_same

            if n_outside_random < 1:
                raise RuntimeError(
                    f"Same-family training subset too large: {group}"
                )

            random_train = np.concatenate([
                same_family_train,
                outside_pool[:n_outside_random]
            ])

            excluded_train = outside_pool.copy()

            if len(random_train) != TRAIN_N:
                raise AssertionError("Random arm train-size mismatch.")

            if len(excluded_train) != TRAIN_N:
                raise AssertionError("Excluded arm train-size mismatch.")

            if set(test_idx) & set(random_train):
                raise AssertionError("Test leakage into random arm.")

            if set(test_idx) & set(excluded_train):
                raise AssertionError("Test leakage into excluded arm.")

            if np.any(groups[excluded_train] == group):
                raise AssertionError(
                    "Held-out family leaked into excluded arm."
                )

            if not np.any(groups[random_train] == group):
                raise AssertionError(
                    "Within-family arm lacks same-family training examples."
                )

            seed_base = (
                SEED
                + repeat * 100000
                + group_index * 1009
            )

            random_pred, _, _, _ = p4.fit_ridge_ensemble(
                X[random_train],
                y[random_train],
                groups[random_train],
                X[test_idx],
                n_boot=BOOTSTRAP_MODELS,
                seed=seed_base + 1,
                alpha=ALPHA
            )

            excluded_pred, _, _, _ = p4.fit_ridge_ensemble(
                X[excluded_train],
                y[excluded_train],
                groups[excluded_train],
                X[test_idx],
                n_boot=BOOTSTRAP_MODELS,
                seed=seed_base + 2,
                alpha=ALPHA
            )

            random_baseline = float(np.mean(y[random_train]))
            excluded_baseline = float(np.mean(y[excluded_train]))

            split_manifest.append({
                "repeat": repeat + 1,
                "IF30_group": group,
                "n_test": len(test_idx),
                "random_train_n": len(random_train),
                "excluded_train_n": len(excluded_train),
                "same_family_train_n": n_same,
                "outside_family_train_random_n": n_outside_random,
                "outside_family_train_excluded_n": TRAIN_N,
                "random_training_groups": len(
                    set(groups[random_train])
                ),
                "excluded_training_groups": len(
                    set(groups[excluded_train])
                )
            })

            for j, idx in enumerate(test_idx):
                fid = ids[idx]

                if fid not in cycle_records[group]:
                    cycle_records[group][fid] = {
                        "y_true": float(y[idx]),
                        "random_abs": [],
                        "excluded_abs": [],
                        "random_sq": [],
                        "excluded_sq": [],
                        "random_base_abs": [],
                        "excluded_base_abs": [],
                        "random_base_sq": [],
                        "excluded_base_sq": [],
                        "random_direction": [],
                        "excluded_direction": []
                    }

                r = cycle_records[group][fid]
                yr = float(y[idx])
                rp = float(random_pred[j])
                ep = float(excluded_pred[j])

                r["random_abs"].append(abs(yr - rp))
                r["excluded_abs"].append(abs(yr - ep))

                r["random_sq"].append((yr - rp) ** 2)
                r["excluded_sq"].append((yr - ep) ** 2)

                r["random_base_abs"].append(
                    abs(yr - random_baseline)
                )
                r["excluded_base_abs"].append(
                    abs(yr - excluded_baseline)
                )

                r["random_base_sq"].append(
                    (yr - random_baseline) ** 2
                )
                r["excluded_base_sq"].append(
                    (yr - excluded_baseline) ** 2
                )

                if abs(yr) >= 0.5:
                    r["random_direction"].append(
                        int(np.sign(yr) == np.sign(rp))
                    )
                    r["excluded_direction"].append(
                        int(np.sign(yr) == np.sign(ep))
                    )

        print(
            f"Completed matched-size repetition "
            f"{repeat + 1}/{REPEATS}"
        )

    family_rows = metric_gaps(cycle_records)

    if len(family_rows) != 14:
        raise RuntimeError("Expected 14 family summary rows.")

    def analyze(key):
        values = np.asarray([
            row[key] for row in family_rows
        ], dtype=float)

        return {
            "mean_family_gap": float(np.mean(values)),
            "families_positive": int(np.sum(values > 0)),
            "families_negative": int(np.sum(values < 0)),
            "families_tied": int(np.sum(values == 0)),
            "exact_sign_flip": exact_sign_flip(values),
            "sign_test": sign_test(values),
            "cluster_bootstrap": cluster_bootstrap(values)
        }

    summary = {
        "experiment": (
            "Training-size-matched comparison of within-family "
            "interpolation versus held-out-family prediction"
        ),
        "source_p4_sha256": source_hash,
        "total_cycles": len(records),
        "IF30_groups_total": len(unique_groups),
        "testable_IF30_groups": len(test_groups),
        "matched_training_size_per_arm": TRAIN_N,
        "repeats_per_family": REPEATS,
        "bootstrap_models_per_fit": BOOTSTRAP_MODELS,
        "ridge_alpha": ALPHA,
        "test_fraction_within_family": TEST_FRACTION,
        "n_family_cycle_test_observations": int(sum(
            row["n_test"] for row in split_manifest
        )),
        "primary_Ridge_MAE_gap": analyze("Ridge_MAE_gap"),
        "secondary_Ridge_RMSE_gap": analyze("Ridge_RMSE_gap"),
        "baseline_MAE_gap": analyze("baseline_MAE_gap"),
        "baseline_RMSE_gap": analyze("baseline_RMSE_gap"),
        "excess_MAE_gap": analyze("excess_MAE_gap"),
        "excess_RMSE_gap": analyze("excess_RMSE_gap"),
        "direction_accuracy": {
            "mean_random": float(np.mean([
                r["random_direction_accuracy"]
                for r in family_rows
                if r["random_direction_accuracy"] is not None
            ])),
            "mean_family_excluded": float(np.mean([
                r["excluded_direction_accuracy"]
                for r in family_rows
                if r["excluded_direction_accuracy"] is not None
            ]))
        },
        "design": (
            "For each testable IF30 group, the same held-out cycles "
            "are evaluated by both arms. Both arms use exactly 150 "
            "training cycles. The within-family arm includes all "
            "non-test cycles from the held-out family and fills its "
            "remaining training slots from a shared external pool. "
            "The family-excluded arm uses 150 cycles from outside "
            "the held-out family."
        ),
        "limitations": [
            "Only 14 non-singleton IF30 groups can be tested.",
            "The 10 repeats reuse the same underlying cycles.",
            "The two arms intentionally differ in family composition.",
            "The result concerns IF30 sequence-similarity groups, not CATH H superfamilies.",
            "This is a SKEMPI-derived PPI epistasis benchmark, not HERV-K Env binding validation."
        ]
    }

    write_csv(
        OUT / "equal_size_family_gap_by_group.csv",
        family_rows
    )

    write_csv(
        OUT / "equal_size_split_manifest.csv",
        split_manifest
    )

    (OUT / "equal_size_family_gap_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8"
    )

    try:
        import matplotlib.pyplot as plt

        ordered = sorted(
            family_rows,
            key=lambda r: r["Ridge_MAE_gap"]
        )

        fig, ax = plt.subplots(figsize=(9, 7))

        ax.barh(
            [
                f"{r['IF30_group']} (n={r['n_test_cycles']})"
                for r in ordered
            ],
            [r["Ridge_MAE_gap"] for r in ordered]
        )

        ax.axvline(0, linestyle="--")
        ax.set_xlabel(
            "Family-excluded MAE minus within-family MAE (kcal/mol)"
        )
        ax.set_title(
            "P4: Matched-training-size family generalization gap"
        )
        ax.grid(axis="x", alpha=0.25)
        fig.tight_layout()

        fig.savefig(
            OUT / "p4_equal_size_family_gap.png",
            dpi=220
        )
        plt.close(fig)

        print("Figure:", OUT / "p4_equal_size_family_gap.png")

    except ImportError:
        print("Matplotlib unavailable; summary data still exported.")

    print("\n=== EQUAL TRAINING-SIZE RESULTS ===")
    print(json.dumps(summary, indent=2))
    print("\nOUTPUT:", OUT)


if __name__ == "__main__":
    main()
