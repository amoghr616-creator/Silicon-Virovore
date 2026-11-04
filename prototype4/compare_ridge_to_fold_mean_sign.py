from pathlib import Path
from itertools import product
import csv
import json

import numpy as np
from scipy.stats import binomtest

SOURCE = Path(
    "prototype4/results/model_a_validation_20261030/"
    "p4_nested_directional_classifier_20261101/"
    "nested_directional_oof_predictions.csv"
)

OUT = Path(
    "prototype4/results/model_a_validation_20261030/"
    "p4_nested_directional_classifier_20261101/"
    "ridge_vs_fold_mean_sign_20261101"
)

BOOTSTRAPS = 20000
SEED = 20261101


def read_csv(path):
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=list(rows[0].keys()),
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def exact_family_sign_flip(differences):
    differences = np.asarray(differences, dtype=float)
    observed = float(np.sum(differences))

    null = np.array([
        np.dot(np.asarray(signs), differences)
        for signs in product((-1.0, 1.0), repeat=len(differences))
    ])

    return {
        "test": "exact paired IF30-family sign-flip test",
        "n_families": len(differences),
        "observed_sum": observed,
        "two_sided_p": float(np.mean(
            np.abs(null) >= abs(observed) - 1e-12
        )),
    }


def family_bootstrap(group_rows, key_a, key_b):
    rng = np.random.default_rng(SEED)
    n_groups = len(group_rows)

    deltas = np.array([
        r[key_a] - r[key_b]
        for r in group_rows
    ])

    sizes = np.array([r["n"] for r in group_rows])
    count_differences = np.rint(deltas * sizes).astype(int)

    macro = np.empty(BOOTSTRAPS)
    pooled = np.empty(BOOTSTRAPS)

    for i in range(BOOTSTRAPS):
        idx = rng.integers(0, n_groups, size=n_groups)
        macro[i] = np.mean(deltas[idx])
        pooled[i] = (
            np.sum(count_differences[idx])
            / np.sum(sizes[idx])
        )

    return {
        "replicates": BOOTSTRAPS,
        "seed": SEED,
        "macro_delta_95_ci": np.quantile(
            macro, [0.025, 0.975]
        ).tolist(),
        "pooled_delta_95_ci": np.quantile(
            pooled, [0.025, 0.975]
        ).tolist(),
        "resampling_unit": "IF30 family",
    }


def mcnemar(y, a, b):
    a_correct = np.asarray(a) == np.asarray(y)
    b_correct = np.asarray(b) == np.asarray(y)

    a_only = int(np.sum(a_correct & ~b_correct))
    b_only = int(np.sum(~a_correct & b_correct))
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
        "exact_p": p,
        "warning": "Cycle-level test; does not account for family clustering."
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    if not SOURCE.exists():
        raise FileNotFoundError(SOURCE)

    rows = read_csv(SOURCE)

    required = {
        "IF30_group",
        "true_sign",
        "existing_ridge_sign_prediction",
        "existing_mean_sign_prediction",
        "training_fold_majority_prediction",
    }

    missing = required - set(rows[0])
    if missing:
        raise RuntimeError(f"Missing columns: {sorted(missing)}")

    groups = sorted({r["IF30_group"] for r in rows})

    if len(groups) != 16:
        raise RuntimeError(
            f"Expected 16 informative IF30 groups, found {len(groups)}"
        )

    # Verify unique cycle-level records.
    feature_ids = [r["feature_uid"] for r in rows]
    if len(feature_ids) != len(set(feature_ids)):
        raise RuntimeError("Duplicate feature IDs found.")

    family_rows = []

    for group in groups:
        subset = [r for r in rows if r["IF30_group"] == group]

        y = np.array([int(r["true_sign"]) for r in subset])
        ridge = np.array([
            int(r["existing_ridge_sign_prediction"])
            for r in subset
        ])
        fold_mean = np.array([
            int(r["existing_mean_sign_prediction"])
            for r in subset
        ])
        majority = np.array([
            int(r["training_fold_majority_prediction"])
            for r in subset
        ])

        family_rows.append({
            "IF30_group": group,
            "n": len(subset),
            "ridge_correct": int(np.sum(ridge == y)),
            "fold_mean_sign_correct": int(np.sum(fold_mean == y)),
            "majority_class_correct": int(np.sum(majority == y)),
            "ridge_accuracy": float(np.mean(ridge == y)),
            "fold_mean_sign_accuracy": float(np.mean(fold_mean == y)),
            "majority_class_accuracy": float(np.mean(majority == y)),
            "ridge_minus_fold_mean": float(
                np.mean(ridge == y) - np.mean(fold_mean == y)
            ),
            "ridge_minus_majority": float(
                np.mean(ridge == y) - np.mean(majority == y)
            ),
        })

    # Primary comparison: Ridge sign versus original regression baseline sign.
    deltas = np.array([
        r["ridge_minus_fold_mean"] for r in family_rows
    ])

    wins = int(np.sum(deltas > 0))
    losses = int(np.sum(deltas < 0))
    ties = int(np.sum(deltas == 0))

    sign_test = (
        float(binomtest(
            min(wins, losses),
            n=wins + losses,
            p=0.5
        ).pvalue)
        if wins + losses else 1.0
    )

    total_n = sum(r["n"] for r in family_rows)
    ridge_correct = sum(r["ridge_correct"] for r in family_rows)
    mean_correct = sum(
        r["fold_mean_sign_correct"] for r in family_rows
    )

    y = [int(r["true_sign"]) for r in rows]
    ridge = [
        int(r["existing_ridge_sign_prediction"])
        for r in rows
    ]
    fold_mean = [
        int(r["existing_mean_sign_prediction"])
        for r in rows
    ]

    summary = {
        "primary_comparison": (
            "Existing Ridge regression sign versus "
            "training-fold-mean regression prediction sign"
        ),
        "validation_unit": "leave-one-IF30-group-out",
        "n_informative_cycles": total_n,
        "n_informative_groups": len(groups),
        "ridge_pooled_accuracy": ridge_correct / total_n,
        "fold_mean_sign_pooled_accuracy": mean_correct / total_n,
        "pooled_difference": (
            ridge_correct - mean_correct
        ) / total_n,
        "ridge_macro_IF30_accuracy": float(np.mean([
            r["ridge_accuracy"] for r in family_rows
        ])),
        "fold_mean_sign_macro_IF30_accuracy": float(np.mean([
            r["fold_mean_sign_accuracy"] for r in family_rows
        ])),
        "macro_difference": float(np.mean(deltas)),
        "family_wins": wins,
        "family_losses": losses,
        "family_ties": ties,
        "family_sign_test_p": sign_test,
        "family_block_permutation": exact_family_sign_flip(deltas),
        "family_cluster_bootstrap": family_bootstrap(
            family_rows,
            "ridge_accuracy",
            "fold_mean_sign_accuracy"
        ),
        "cycle_level_mcnemar_descriptive": mcnemar(
            y, ridge, fold_mean
        ),
        "secondary_comparison": {
            "ridge_macro_accuracy": float(np.mean([
                r["ridge_accuracy"] for r in family_rows
            ])),
            "majority_macro_accuracy": float(np.mean([
                r["majority_class_accuracy"] for r in family_rows
            ])),
            "macro_difference": float(np.mean([
                r["ridge_minus_majority"] for r in family_rows
            ])),
        },
        "scope": (
            "SKEMPI-derived PPI mutation-cycle epistasis sign; "
            "not HERV-K Env peptide-binding validation."
        ),
        "interpretation": (
            "Exploratory family-held-out comparison. "
            "A positive point estimate is not evidence of a robust "
            "advantage unless family-level uncertainty supports it."
        ),
    }

    with (OUT / "ridge_vs_fold_mean_sign_summary.json").open(
        "w", encoding="utf-8"
    ) as f:
        json.dump(summary, f, indent=2)

    write_csv(
        OUT / "ridge_vs_fold_mean_sign_by_family.csv",
        family_rows
    )

    try:
        import matplotlib.pyplot as plt

        ordered = sorted(
            family_rows,
            key=lambda r: r["ridge_minus_fold_mean"]
        )

        fig, ax = plt.subplots(figsize=(9, 7))
        labels = [
            f"{r['IF30_group']} (n={r['n']})"
            for r in ordered
        ]
        values = [
            100 * r["ridge_minus_fold_mean"]
            for r in ordered
        ]

        ax.barh(labels, values)
        ax.axvline(0, linestyle="--")
        ax.set_xlabel(
            "Ridge minus fold-mean sign accuracy "
            "(percentage points)"
        )
        ax.set_title(
            "P4: Ridge versus regression-baseline direction by IF30 group"
        )
        ax.grid(axis="x", alpha=0.25)
        fig.tight_layout()
        fig.savefig(
            OUT / "p4_ridge_vs_fold_mean_sign_by_family.png",
            dpi=220
        )
        plt.close(fig)

        print("Figure:", OUT / "p4_ridge_vs_fold_mean_sign_by_family.png")

    except ImportError:
        print("Matplotlib unavailable; data exports completed.")

    print("\n=== RIDGE VS FOLD-MEAN SIGN ===")
    print(json.dumps(summary, indent=2))
    print("\nOUTPUT:", OUT)


if __name__ == "__main__":
    main()
