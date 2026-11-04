from pathlib import Path
from itertools import product
import csv
import json

import numpy as np
from scipy.stats import binomtest

ROOT = Path(
    "prototype4/results/model_a_validation_20261030/"
    "p4_split_sensitivity_20261101/"
    "matched_family_gap_20261101"
)
SOURCE = ROOT / "matched_family_gap_by_group.csv"
OUT = ROOT / "excess_generalization_gap_20261101"

N_BOOT = 20000
SEED = 20261101


def read_csv(path):
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def exact_sign_flip(values):
    values = np.asarray(values, dtype=float)
    observed = float(np.sum(values))

    null = np.fromiter(
        (
            float(np.dot(np.asarray(signs), values))
            for signs in product(
                (-1.0, 1.0), repeat=len(values)
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


def cluster_bootstrap(values):
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(SEED)
    n = len(values)

    draws = np.empty(N_BOOT)

    for i in range(N_BOOT):
        idx = rng.integers(0, n, size=n)
        draws[i] = np.mean(values[idx])

    return {
        "replicates": N_BOOT,
        "seed": SEED,
        "mean": float(np.mean(values)),
        "ci95": np.quantile(draws, [0.025, 0.975]).tolist(),
        "resampling_unit": "IF30 family",
        "note": "Conditional on the existing saved split predictions."
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    if not SOURCE.exists():
        raise FileNotFoundError(SOURCE)

    rows = read_csv(SOURCE)

    required = {
        "IF30_group",
        "n_cycles",
        "random_Ridge_MAE",
        "LOFO_Ridge_MAE",
        "random_baseline_MAE",
        "LOFO_baseline_MAE",
        "random_Ridge_RMSE",
        "LOFO_Ridge_RMSE",
        "random_baseline_RMSE",
        "LOFO_baseline_RMSE",
    }

    missing = required - set(rows[0])

    if missing:
        raise RuntimeError(
            f"Missing columns: {sorted(missing)}"
        )

    family_results = []

    for r in rows:
        ridge_mae_gap = (
            float(r["LOFO_Ridge_MAE"])
            - float(r["random_Ridge_MAE"])
        )

        baseline_mae_gap = (
            float(r["LOFO_baseline_MAE"])
            - float(r["random_baseline_MAE"])
        )

        ridge_rmse_gap = (
            float(r["LOFO_Ridge_RMSE"])
            - float(r["random_Ridge_RMSE"])
        )

        baseline_rmse_gap = (
            float(r["LOFO_baseline_RMSE"])
            - float(r["random_baseline_RMSE"])
        )

        family_results.append({
            "IF30_group": r["IF30_group"],
            "n_cycles": int(r["n_cycles"]),
            "ridge_MAE_gap": ridge_mae_gap,
            "baseline_MAE_gap": baseline_mae_gap,
            "excess_MAE_generalization_gap": (
                ridge_mae_gap - baseline_mae_gap
            ),
            "ridge_RMSE_gap": ridge_rmse_gap,
            "baseline_RMSE_gap": baseline_rmse_gap,
            "excess_RMSE_generalization_gap": (
                ridge_rmse_gap - baseline_rmse_gap
            ),
        })

    if len(family_results) != 14:
        raise RuntimeError(
            f"Expected 14 matched families; got {len(family_results)}"
        )

    def analyze(metric):
        values = np.array([
            r[metric] for r in family_results
        ], dtype=float)

        positive = int(np.sum(values > 0))
        negative = int(np.sum(values < 0))
        ties = int(np.sum(values == 0))

        sign_p = (
            float(binomtest(
                min(positive, negative),
                n=positive + negative,
                p=0.5
            ).pvalue)
            if positive + negative else 1.0
        )

        return {
            "mean_excess_gap": float(np.mean(values)),
            "positive_families": positive,
            "negative_families": negative,
            "tied_families": ties,
            "exact_family_sign_flip": exact_sign_flip(values),
            "family_sign_test_p": sign_p,
            "family_cluster_bootstrap": cluster_bootstrap(values),
        }

    summary = {
        "experiment": (
            "Difference-in-differences comparison of Ridge and "
            "training-fold-mean generalization gaps"
        ),
        "n_matched_IF30_families": len(family_results),
        "n_matched_cycles": int(sum(
            r["n_cycles"] for r in family_results
        )),
        "primary_MAE_analysis": analyze(
            "excess_MAE_generalization_gap"
        ),
        "secondary_RMSE_analysis": analyze(
            "excess_RMSE_generalization_gap"
        ),
        "definition": (
            "(Ridge LOFO error - Ridge random error) - "
            "(baseline LOFO error - baseline random error)"
        ),
        "interpretation": (
            "Positive excess gap indicates that Ridge's increase in "
            "error under family-held-out evaluation exceeds the "
            "corresponding increase for the training-fold-mean baseline."
        ),
        "limitations": [
            "Exploratory analysis on an existing dataset.",
            "The random split and LOFO training sets differ in size.",
            "Baseline adjustment does not eliminate every training-set difference.",
            "Bootstrap uncertainty is conditional on saved random predictions.",
            "IF30 grouping is not CATH H-level superfamily grouping.",
            "This is a SKEMPI-derived PPI-cycle benchmark, not HERV-K Env affinity validation."
        ],
    }

    with (OUT / "excess_generalization_gap_summary.json").open(
        "w", encoding="utf-8"
    ) as f:
        json.dump(summary, f, indent=2)

    with (OUT / "excess_generalization_gap_by_family.csv").open(
        "w", newline="", encoding="utf-8"
    ) as f:
        writer = csv.DictWriter(
            f,
            fieldnames=list(family_results[0].keys())
        )
        writer.writeheader()
        writer.writerows(family_results)

    try:
        import matplotlib.pyplot as plt

        ordered = sorted(
            family_results,
            key=lambda r: r["excess_MAE_generalization_gap"]
        )

        fig, ax = plt.subplots(figsize=(9, 7))

        labels = [
            f"{r['IF30_group']} (n={r['n_cycles']})"
            for r in ordered
        ]

        values = [
            r["excess_MAE_generalization_gap"]
            for r in ordered
        ]

        ax.barh(labels, values)
        ax.axvline(0, linestyle="--")
        ax.set_xlabel(
            "Ridge excess generalization gap (kcal/mol)"
        )
        ax.set_title(
            "P4: Ridge-specific family generalization penalty"
        )
        ax.grid(axis="x", alpha=0.25)
        fig.tight_layout()
        fig.savefig(
            OUT / "p4_excess_generalization_gap.png",
            dpi=220
        )
        plt.close(fig)

        print("Figure:", OUT / "p4_excess_generalization_gap.png")

    except ImportError:
        print("Matplotlib unavailable; data exports completed.")

    print("\n=== EXCESS GENERALIZATION GAP ===")
    print(json.dumps(summary, indent=2))
    print("\nOUTPUT:", OUT)


if __name__ == "__main__":
    main()
