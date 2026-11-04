from pathlib import Path
from itertools import product
import csv
import json

import numpy as np
from scipy.stats import binomtest

ROOT = Path(
    "prototype4/results/model_a_validation_20261030/"
    "p4_nested_directional_classifier_20261101"
)

FAMILY_FILE = ROOT / "nested_directional_family_metrics.csv"
PRED_FILE = ROOT / "nested_directional_oof_predictions.csv"
OUT = ROOT / "ridge_family_advantage_20261101"

BOOTSTRAPS = 20000
SEED = 20261101


def read_csv(path):
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def exact_sign_flip(differences):
    differences = np.asarray(differences, dtype=float)
    observed = float(np.sum(differences))

    null = np.array([
        np.dot(np.asarray(signs), differences)
        for signs in product((-1.0, 1.0), repeat=len(differences))
    ])

    return {
        "test": "exact paired IF30-family sign-flip permutation",
        "groups": len(differences),
        "observed_sum_family_differences": observed,
        "two_sided_p": float(np.mean(
            np.abs(null) >= abs(observed) - 1e-12
        )),
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    family_rows = read_csv(FAMILY_FILE)
    prediction_rows = read_csv(PRED_FILE)

    required = {
        "IF30_group", "n",
        "ridge_sign_correct", "majority_correct",
        "ridge_sign_accuracy", "majority_accuracy",
    }

    missing = required - set(family_rows[0])
    if missing:
        raise RuntimeError(
            f"Missing family CSV columns: {sorted(missing)}"
        )

    groups = []

    for r in family_rows:
        n = int(r["n"])
        ridge_correct = int(r["ridge_sign_correct"])
        majority_correct = int(r["majority_correct"])

        if n <= 0:
            raise RuntimeError("Invalid family sample count.")

        groups.append({
            "IF30_group": r["IF30_group"],
            "n": n,
            "ridge_correct": ridge_correct,
            "majority_correct": majority_correct,
            "ridge_accuracy": ridge_correct / n,
            "majority_accuracy": majority_correct / n,
            "delta": (
                ridge_correct - majority_correct
            ) / n,
            "net_correct": ridge_correct - majority_correct,
        })

    if len(groups) != 16:
        raise RuntimeError(
            f"Expected 16 informative groups; found {len(groups)}"
        )

    # Primary family-level statistics.
    deltas = np.array([r["delta"] for r in groups])
    sizes = np.array([r["n"] for r in groups])
    net_counts = np.array([r["net_correct"] for r in groups])

    macro_delta = float(np.mean(deltas))
    pooled_delta = float(np.sum(net_counts) / np.sum(sizes))

    wins = sum(x > 0 for x in deltas)
    losses = sum(x < 0 for x in deltas)
    ties = sum(x == 0 for x in deltas)

    sign_test = binomtest(
        min(wins, losses),
        n=wins + losses,
        p=0.5,
        alternative="two-sided"
    )

    block_test = exact_sign_flip(deltas)

    # Family-cluster bootstrap.
    rng = np.random.default_rng(SEED)
    n_groups = len(groups)

    macro_boot = np.empty(BOOTSTRAPS)
    pooled_boot = np.empty(BOOTSTRAPS)

    for i in range(BOOTSTRAPS):
        idx = rng.integers(0, n_groups, size=n_groups)

        macro_boot[i] = np.mean(deltas[idx])
        pooled_boot[i] = (
            np.sum(net_counts[idx]) / np.sum(sizes[idx])
        )

    macro_ci = np.quantile(macro_boot, [0.025, 0.975])
    pooled_ci = np.quantile(pooled_boot, [0.025, 0.975])

    # Leave-one-family-out influence analysis.
    loo_macro = [
        float(np.mean(np.delete(deltas, i)))
        for i in range(n_groups)
    ]

    # Observation-level paired McNemar comparison.
    ridge_only = 0
    majority_only = 0
    ridge_total = 0
    majority_total = 0

    for r in prediction_rows:
        y = int(r["true_sign"])
        ridge = int(r["existing_ridge_sign_prediction"])
        majority = int(r["training_fold_majority_prediction"])

        ridge_correct = ridge == y
        majority_correct = majority == y

        ridge_total += int(ridge_correct)
        majority_total += int(majority_correct)

        ridge_only += int(ridge_correct and not majority_correct)
        majority_only += int(majority_correct and not ridge_correct)

    discordant = ridge_only + majority_only

    mcnemar_p = (
        float(binomtest(
            min(ridge_only, majority_only),
            n=discordant,
            p=0.5
        ).pvalue)
        if discordant else 1.0
    )

    summary = {
        "comparison": (
            "Existing Ridge regression sign versus "
            "outer training-fold majority classifier"
        ),
        "validation_unit": "IF30 family",
        "n_informative_groups": n_groups,
        "n_informative_cycles": int(np.sum(sizes)),
        "ridge_pooled_accuracy": ridge_total / len(prediction_rows),
        "majority_pooled_accuracy": (
            majority_total / len(prediction_rows)
        ),
        "ridge_family_macro_accuracy": float(np.mean([
            r["ridge_accuracy"] for r in groups
        ])),
        "majority_family_macro_accuracy": float(np.mean([
            r["majority_accuracy"] for r in groups
        ])),
        "macro_difference_ridge_minus_majority": macro_delta,
        "pooled_difference_ridge_minus_majority": pooled_delta,
        "family_wins": wins,
        "family_losses": losses,
        "family_ties": ties,
        "family_sign_test_p": float(sign_test.pvalue),
        "family_block_permutation": block_test,
        "family_cluster_bootstrap": {
            "replicates": BOOTSTRAPS,
            "seed": SEED,
            "macro_difference_95_ci": macro_ci.tolist(),
            "pooled_difference_95_ci": pooled_ci.tolist(),
        },
        "leave_one_family_out_macro_difference_range": [
            min(loo_macro), max(loo_macro)
        ],
        "observation_level_paired_comparison": {
            "ridge_only_correct": ridge_only,
            "majority_only_correct": majority_only,
            "discordant_pairs": discordant,
            "mcnemar_exact_p": mcnemar_p,
        },
        "scope_warning": (
            "Exploratory SKEMPI-derived PPI epistasis sign evaluation; "
            "not HERV-K Env peptide-binding validation. "
            "Only 16 of 19 IF30 groups have informative cycles."
        ),
    }

    with (OUT / "ridge_family_advantage_summary.json").open(
        "w", encoding="utf-8"
    ) as f:
        json.dump(
            summary,
            f,
            indent=2,
            default=lambda value: (
                value.item()
                if hasattr(value, "item")
                else str(value)
            )
        )

    with (OUT / "ridge_family_advantage_by_group.csv").open(
        "w", newline="", encoding="utf-8"
    ) as f:
        writer = csv.DictWriter(
            f,
            fieldnames=list(groups[0].keys())
        )
        writer.writeheader()
        writer.writerows(groups)

    try:
        import matplotlib.pyplot as plt

        ordered = sorted(
            groups, key=lambda r: r["delta"]
        )

        fig, ax = plt.subplots(figsize=(9, 7))
        labels = [
            f"{r['IF30_group']} (n={r['n']})"
            for r in ordered
        ]
        values = [100 * r["delta"] for r in ordered]

        ax.barh(labels, values)
        ax.axvline(0, linestyle="--")
        ax.set_xlabel(
            "Ridge minus majority accuracy (percentage points)"
        )
        ax.set_title(
            "P4: Existing Ridge directional advantage by held-out IF30 group"
        )
        ax.grid(axis="x", alpha=0.25)
        fig.tight_layout()
        fig.savefig(
            OUT / "p4_ridge_vs_majority_by_family.png",
            dpi=220
        )
        plt.close(fig)

        print("Figure:", OUT / "p4_ridge_vs_majority_by_family.png")

    except ImportError:
        print("Matplotlib unavailable; CSV and JSON still saved.")

    print("\n=== RIDGE FAMILY ADVANTAGE ===")
    print(json.dumps(
        summary,
        indent=2,
        default=lambda value: (
            value.item() if hasattr(value, "item") else str(value)
        )
    ))
    print("\nOUTPUT:", OUT)


if __name__ == "__main__":
    main()
