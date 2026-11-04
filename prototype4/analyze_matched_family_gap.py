from pathlib import Path
from collections import defaultdict
from itertools import product
import sqlite3
import hashlib
import json
import csv

import numpy as np

ROOT = Path(
    "prototype4/results/model_a_validation_20261030"
)

DB = ROOT / "p4_model_a_working.sqlite"

RANDOM = (
    ROOT / "p4_split_sensitivity_20261101/"
    "within_family_random_split_predictions.csv"
)

OUT = (
    ROOT / "p4_split_sensitivity_20261101/"
    "matched_family_gap_20261101"
)

EXPECTED_P4_HASH = (
    "4ccb1f37b4be33b23636eab45584021b1df491178f7ac23b62f2609b7cca1181"
)

BOOTSTRAPS = 20000
SEED = 20261101
INFORMATIVE_THRESHOLD = 0.5


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_csv(path):
    if not path.exists():
        raise FileNotFoundError(path)

    with path.open(
        "r", encoding="utf-8", newline=""
    ) as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)

    if not rows:
        raise RuntimeError(f"No rows to write: {path}")

    with path.open(
        "w", encoding="utf-8", newline=""
    ) as f:
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
            for signs in product(
                (-1.0, 1.0),
                repeat=len(values)
            )
        ),
        dtype=float
    )

    return {
        "test": "exact paired IF30-family sign-flip permutation",
        "n_families": len(values),
        "observed_sum": observed,
        "two_sided_p": float(np.mean(
            np.abs(null) >= abs(observed) - 1e-12
        )),
    }


def family_bootstrap(values, rng):
    values = np.asarray(values, dtype=float)
    n = len(values)

    draws = np.empty(BOOTSTRAPS)

    for i in range(BOOTSTRAPS):
        idx = rng.integers(0, n, size=n)
        draws[i] = np.mean(values[idx])

    return {
        "replicates": BOOTSTRAPS,
        "seed": SEED,
        "mean_difference": float(np.mean(values)),
        "ci95": np.quantile(
            draws, [0.025, 0.975]
        ).tolist(),
        "resampling_unit": "IF30 family",
        "interpretation": (
            "Conditional on the saved random-split predictions."
        )
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    if not DB.exists() or not RANDOM.exists():
        raise FileNotFoundError(
            f"Missing input:\n{DB}\n{RANDOM}"
        )

    # Verify the exact original model source.
    p4_source = Path("prototype4/p4.py")
    p4_hash = sha256(p4_source)

    if p4_hash != EXPECTED_P4_HASH:
        raise RuntimeError(
            "Original p4.py source hash has changed."
        )

    random_rows = read_csv(RANDOM)

    con = sqlite3.connect(
        f"file:{DB.resolve()}?mode=ro",
        uri=True
    )
    con.row_factory = sqlite3.Row

    lofo_rows = con.execute("""
        SELECT
            f.feature_uid,
            f.cycle_uid,
            f.family_id,
            f.complex_id,
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

    if len(lofo_rows) != 242:
        raise RuntimeError(
            f"Expected 242 LOFO records; got {len(lofo_rows)}"
        )

    lofo = {}

    for row in lofo_rows:
        r = dict(row)
        fid = r["feature_uid"]

        if r["family_id"] != r["heldout_group"]:
            raise RuntimeError(
                f"LOFO family mismatch: {fid}"
            )

        if fid in lofo:
            raise RuntimeError(
                f"Duplicate LOFO feature: {fid}"
            )

        if not np.isclose(
            float(r["y_true"]),
            float(r["experimental_epistasis_kcal_mol"]),
            atol=1e-8,
            rtol=0
        ):
            raise RuntimeError(
                f"LOFO target mismatch: {fid}"
            )

        lofo[fid] = r

    # Build repeated within-family predictions.
    random_by_family = defaultdict(lambda: defaultdict(list))
    split_ids = set()

    for r in random_rows:
        fid = r["feature_uid"]
        group = r["IF30_group"]

        if fid not in lofo:
            raise RuntimeError(
                f"Random prediction not found in LOFO source: {fid}"
            )

        if group != lofo[fid]["family_id"]:
            raise RuntimeError(
                f"Random prediction family mismatch: {fid}"
            )

        if not np.isclose(
            float(r["y_true"]),
            float(lofo[fid]["y_true"]),
            atol=1e-8,
            rtol=0
        ):
            raise RuntimeError(
                f"Random target mismatch: {fid}"
            )

        repetition = int(r["repetition"])
        split_ids.add(repetition)

        random_by_family[group][fid].append({
            "repetition": repetition,
            "y_true": float(r["y_true"]),
            "y_pred": float(r["y_pred"]),
            "baseline_pred": float(r["baseline_pred"]),
        })

    # Only families with random held-out observations are eligible.
    random_families = sorted(random_by_family)

    if len(random_families) < 5:
        raise RuntimeError(
            "Too few matched families for the planned analysis."
        )

    family_rows = []

    for group in random_families:
        cycle_predictions = random_by_family[group]

        lofo_ids = {
            fid for fid, record in lofo.items()
            if record["family_id"] == group
        }

        random_ids = set(cycle_predictions)

        if random_ids != lofo_ids:
            missing = sorted(lofo_ids - random_ids)
            extra = sorted(random_ids - lofo_ids)

            raise RuntimeError(
                f"Incomplete random-cycle coverage for {group}; "
                f"missing={missing[:5]}, extra={extra[:5]}"
            )

        per_cycle = []

        for fid in sorted(lofo_ids):
            repetitions = cycle_predictions[fid]
            observed = float(lofo[fid]["y_true"])
            lofo_prediction = float(lofo[fid]["y_pred"])
            lofo_baseline = float(lofo[fid]["baseline_pred"])

            abs_errors = []
            squared_errors = []
            baseline_abs_errors = []
            baseline_squared_errors = []
            direction_correct = []
            baseline_direction_correct = []

            for prediction in repetitions:
                yp = prediction["y_pred"]
                bp = prediction["baseline_pred"]

                abs_errors.append(abs(observed - yp))
                squared_errors.append((observed - yp) ** 2)

                baseline_abs_errors.append(
                    abs(observed - bp)
                )
                baseline_squared_errors.append(
                    (observed - bp) ** 2
                )

                if abs(observed) >= INFORMATIVE_THRESHOLD:
                    direction_correct.append(
                        int(np.sign(observed) == np.sign(yp))
                    )
                    baseline_direction_correct.append(
                        int(np.sign(observed) == np.sign(bp))
                    )

            informative = (
                abs(observed) >= INFORMATIVE_THRESHOLD
            )

            per_cycle.append({
                "feature_uid": fid,
                "n_random_predictions": len(repetitions),
                "y_true": observed,
                "random_abs_error_mean": float(
                    np.mean(abs_errors)
                ),
                "random_squared_error_mean": float(
                    np.mean(squared_errors)
                ),
                "random_baseline_abs_error_mean": float(
                    np.mean(baseline_abs_errors)
                ),
                "random_baseline_squared_error_mean": float(
                    np.mean(baseline_squared_errors)
                ),
                "random_direction_accuracy": (
                    float(np.mean(direction_correct))
                    if informative else None
                ),
                "random_baseline_direction_accuracy": (
                    float(np.mean(baseline_direction_correct))
                    if informative else None
                ),
                "lofo_abs_error": abs(
                    observed - lofo_prediction
                ),
                "lofo_squared_error": (
                    observed - lofo_prediction
                ) ** 2,
                "lofo_baseline_abs_error": abs(
                    observed - lofo_baseline
                ),
                "lofo_baseline_squared_error": (
                    observed - lofo_baseline
                ) ** 2,
                "lofo_direction_correct": (
                    int(
                        np.sign(observed)
                        == np.sign(lofo_prediction)
                    )
                    if informative else None
                ),
                "lofo_baseline_direction_correct": (
                    int(
                        np.sign(observed)
                        == np.sign(lofo_baseline)
                    )
                    if informative else None
                ),
            })

        random_mae = np.mean([
            r["random_abs_error_mean"] for r in per_cycle
        ])
        lofo_mae = np.mean([
            r["lofo_abs_error"] for r in per_cycle
        ])

        random_rmse = np.sqrt(np.mean([
            r["random_squared_error_mean"] for r in per_cycle
        ]))
        lofo_rmse = np.sqrt(np.mean([
            r["lofo_squared_error"] for r in per_cycle
        ]))

        random_base_mae = np.mean([
            r["random_baseline_abs_error_mean"]
            for r in per_cycle
        ])
        lofo_base_mae = np.mean([
            r["lofo_baseline_abs_error"]
            for r in per_cycle
        ])

        random_base_rmse = np.sqrt(np.mean([
            r["random_baseline_squared_error_mean"]
            for r in per_cycle
        ]))
        lofo_base_rmse = np.sqrt(np.mean([
            r["lofo_baseline_squared_error"]
            for r in per_cycle
        ]))

        informative_cycles = [
            r for r in per_cycle
            if r["random_direction_accuracy"] is not None
        ]

        if informative_cycles:
            random_direction = float(np.mean([
                r["random_direction_accuracy"]
                for r in informative_cycles
            ]))
            lofo_direction = float(np.mean([
                r["lofo_direction_correct"]
                for r in informative_cycles
            ]))
            random_base_direction = float(np.mean([
                r["random_baseline_direction_accuracy"]
                for r in informative_cycles
            ]))
            lofo_base_direction = float(np.mean([
                r["lofo_baseline_direction_correct"]
                for r in informative_cycles
            ]))
        else:
            random_direction = None
            lofo_direction = None
            random_base_direction = None
            lofo_base_direction = None

        repeat_counts = [
            r["n_random_predictions"] for r in per_cycle
        ]

        family_rows.append({
            "IF30_group": group,
            "n_cycles": len(per_cycle),
            "n_informative_cycles": len(informative_cycles),
            "random_prediction_rows": sum(repeat_counts),
            "min_repeats_per_cycle": min(repeat_counts),
            "median_repeats_per_cycle": float(
                np.median(repeat_counts)
            ),
            "random_Ridge_MAE": float(random_mae),
            "LOFO_Ridge_MAE": float(lofo_mae),
            "MAE_generalization_gap": float(
                lofo_mae - random_mae
            ),
            "random_Ridge_RMSE": float(random_rmse),
            "LOFO_Ridge_RMSE": float(lofo_rmse),
            "RMSE_generalization_gap": float(
                lofo_rmse - random_rmse
            ),
            "random_Ridge_direction_accuracy": random_direction,
            "LOFO_Ridge_direction_accuracy": lofo_direction,
            "direction_accuracy_drop": (
                random_direction - lofo_direction
                if random_direction is not None else None
            ),
            "random_baseline_MAE": float(random_base_mae),
            "LOFO_baseline_MAE": float(lofo_base_mae),
            "baseline_MAE_gap": float(
                lofo_base_mae - random_base_mae
            ),
            "random_baseline_RMSE": float(random_base_rmse),
            "LOFO_baseline_RMSE": float(lofo_base_rmse),
            "baseline_RMSE_gap": float(
                lofo_base_rmse - random_base_rmse
            ),
            "random_baseline_direction_accuracy":
                random_base_direction,
            "LOFO_baseline_direction_accuracy":
                lofo_base_direction,
        })

    if not family_rows:
        raise RuntimeError("No matched family metrics created.")

    rng = np.random.default_rng(SEED)

    def compare_metric(gap_key):
        usable = [
            r for r in family_rows
            if r[gap_key] is not None
            and np.isfinite(r[gap_key])
        ]

        values = np.array([
            r[gap_key] for r in usable
        ], dtype=float)

        return {
            "n_families": len(values),
            "macro_gap": float(np.mean(values)),
            "families_with_positive_gap": int(
                np.sum(values > 0)
            ),
            "families_with_negative_gap": int(
                np.sum(values < 0)
            ),
            "families_with_zero_gap": int(
                np.sum(values == 0)
            ),
            "exact_family_sign_flip": exact_sign_flip(values),
            "family_cluster_bootstrap": family_bootstrap(
                values, rng
            ),
        }

    # Equal-family weighting is primary.
    comparisons = {
        "ridge_MAE_LOFO_minus_random": compare_metric(
            "MAE_generalization_gap"
        ),
        "ridge_RMSE_LOFO_minus_random": compare_metric(
            "RMSE_generalization_gap"
        ),
        "ridge_direction_random_minus_LOFO": compare_metric(
            "direction_accuracy_drop"
        ),
        "baseline_MAE_LOFO_minus_random": compare_metric(
            "baseline_MAE_gap"
        ),
        "baseline_RMSE_LOFO_minus_random": compare_metric(
            "baseline_RMSE_gap"
        ),
    }

    summary = {
        "experiment": (
            "Matched-family comparison of within-family random "
            "interpolation versus IF30 leave-one-family-out prediction"
        ),
        "model": (
            "Original p4.py group-bootstrap Ridge ensemble; "
            "alpha=10; 30 models per fit"
        ),
        "original_p4_sha256": p4_hash,
        "random_split_repetitions": len(split_ids),
        "random_split_design": (
            "Repeated 20% within-family holdout; remaining cycles "
            "from each testable family stay in training."
        ),
        "n_matched_IF30_families": len(family_rows),
        "n_matched_cycles": int(sum(
            r["n_cycles"] for r in family_rows
        )),
        "family_level_comparisons": comparisons,
        "interpretation": (
            "Positive MAE/RMSE gaps indicate higher error when an "
            "entire IF30 group is held out. Positive direction gap "
            "indicates lower direction accuracy under IF30 LOFO."
        ),
        "limitations": [
            "Random split repetitions reuse the same underlying cycles.",
            "Family bootstrap intervals are conditional on saved random predictions.",
            "This evaluates interpolation within IF30 groups versus transfer to held-out IF30 groups.",
            "It does not establish that every random split protocol causes leakage.",
            "IF30 groups are sequence-similarity groups, not CATH H-level superfamilies.",
            "This SKEMPI-derived PPI epistasis benchmark is not HERV-K Env peptide-binding validation."
        ]
    }

    write_csv(
        OUT / "matched_family_gap_by_group.csv",
        family_rows
    )

    with (
        OUT / "matched_family_gap_summary.json"
    ).open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("\n=== MATCHED-FAMILY COMPARISON ===")
    print(json.dumps(summary, indent=2))

    try:
        import matplotlib.pyplot as plt

        # Figure 1: paired family MAE.
        x = np.array([
            r["random_Ridge_MAE"] for r in family_rows
        ])
        y = np.array([
            r["LOFO_Ridge_MAE"] for r in family_rows
        ])

        upper = max(float(np.max(x)), float(np.max(y))) * 1.1

        fig, ax = plt.subplots(figsize=(7, 6))
        ax.scatter(x, y)

        ax.plot(
            [0, upper], [0, upper],
            linestyle="--"
        )

        ax.set_xlim(0, upper)
        ax.set_ylim(0, upper)
        ax.set_xlabel("Within-family random MAE (kcal/mol)")
        ax.set_ylabel("IF30-held-out MAE (kcal/mol)")
        ax.set_title(
            "P4: Matched-family prediction error comparison"
        )
        ax.grid(alpha=0.25)
        fig.tight_layout()
        fig.savefig(
            OUT / "p4_matched_family_mae_scatter.png",
            dpi=220
        )
        plt.close(fig)

        # Figure 2: per-family MAE generalization gap.
        ordered = sorted(
            family_rows,
            key=lambda r: r["MAE_generalization_gap"]
        )

        fig, ax = plt.subplots(figsize=(9, 7))
        labels = [
            f"{r['IF30_group']} (n={r['n_cycles']})"
            for r in ordered
        ]
        values = [
            r["MAE_generalization_gap"] for r in ordered
        ]

        ax.barh(labels, values)
        ax.axvline(0, linestyle="--")
        ax.set_xlabel(
            "LOFO MAE minus random-split MAE (kcal/mol)"
        )
        ax.set_title(
            "P4: Family-specific generalization gap"
        )
        ax.grid(axis="x", alpha=0.25)
        fig.tight_layout()
        fig.savefig(
            OUT / "p4_family_generalization_gap.png",
            dpi=220
        )
        plt.close(fig)

        print("\nFigures:")
        print(OUT / "p4_matched_family_mae_scatter.png")
        print(OUT / "p4_family_generalization_gap.png")

    except ImportError:
        print("Matplotlib unavailable; summary tables were saved.")

    print("\nOUTPUT DIRECTORY:")
    print(OUT)


if __name__ == "__main__":
    main()
