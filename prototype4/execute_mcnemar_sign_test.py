
# prototype4/execute_mcnemar_sign_test.py
#
# Paired directional comparison:
#   Ridge LOFO predictions vs training-fold-mean predictions.
#
# Reports:
#   1. Paired contingency counts
#   2. Exact McNemar test
#   3. Continuity-corrected McNemar statistic
#   4. Pooled accuracy difference
#   5. Family-macro difference and family-cluster bootstrap CI
#
# This does not make the cycles statistically independent.
# Interpret the family bootstrap cautiously because there are only
# 16 informative families.

import sqlite3
import math
from pathlib import Path
from collections import defaultdict

import numpy as np
from scipy.stats import binomtest, chi2

DB_PATH = Path(
    "prototype4/results/model_a_validation_20261030/"
    "p4_model_a_working.sqlite"
)

SEED = 20261102
BOOTSTRAP_REPS = 20000
THRESHOLD = 0.5


def main():
    if not DB_PATH.exists():
        raise FileNotFoundError(DB_PATH)

    with sqlite3.connect(DB_PATH) as conn:
        columns = {
            row[1]
            for row in conn.execute(
                "PRAGMA table_info(loco_predictions)"
            )
        }

        required = {
            "y_true",
            "y_pred",
            "baseline_pred",
            "heldout_group",
        }

        missing = required - columns
        if missing:
            raise RuntimeError(
                f"Missing required columns: {sorted(missing)}"
            )

        rows = conn.execute(
            """
            SELECT y_true, y_pred, baseline_pred, heldout_group
            FROM loco_predictions
            WHERE ABS(y_true) >= ?
            """,
            (THRESHOLD,),
        ).fetchall()

    if not rows:
        raise RuntimeError("No informative validation cycles found.")

    y_true = np.asarray([float(r[0]) for r in rows])
    y_ridge = np.asarray([float(r[1]) for r in rows])
    y_base = np.asarray([float(r[2]) for r in rows])
    families = np.asarray([str(r[3]) for r in rows])

    if not (
        np.isfinite(y_true).all()
        and np.isfinite(y_ridge).all()
        and np.isfinite(y_base).all()
    ):
        raise RuntimeError("Non-finite prediction values detected.")

    if np.any(y_true == 0):
        raise RuntimeError(
            "Unexpected zero-valued truth in informative subset."
        )

    ridge_correct = np.sign(y_ridge) == np.sign(y_true)
    base_correct = np.sign(y_base) == np.sign(y_true)

    # Discordant pairs
    b = int(np.sum(ridge_correct & ~base_correct))
    c = int(np.sum(~ridge_correct & base_correct))

    n = len(rows)
    ridge_accuracy = float(ridge_correct.mean())
    base_accuracy = float(base_correct.mean())

    print("========== PAIRED DIRECTIONAL RESULTS ==========")
    print(f"Informative cycles: {n}")
    print(f"Informative families: {len(set(families))}")
    print(f"Ridge correct: {int(ridge_correct.sum())}/{n}")
    print(f"Baseline correct: {int(base_correct.sum())}/{n}")
    print(f"Pooled accuracy difference: "
          f"{ridge_accuracy - base_accuracy:+.4%}")

    print("\n========== MCNEMAR CONTINGENCY ==========")
    print(f"Ridge only correct (b): {b}")
    print(f"Baseline only correct (c): {c}")
    print(f"Both correct: "
          f"{int(np.sum(ridge_correct & base_correct))}")
    print(f"Both incorrect: "
          f"{int(np.sum(~ridge_correct & ~base_correct))}")

    if b + c == 0:
        exact_p = 1.0
        statistic = 0.0
        asymptotic_p = 1.0
    else:
        # Exact two-sided McNemar test
        exact_p = float(
            binomtest(
                b,
                n=b + c,
                p=0.5,
                alternative="two-sided",
            ).pvalue
        )

        # Continuity-corrected McNemar chi-square
        numerator = max(abs(b - c) - 1, 0) ** 2
        statistic = numerator / (b + c)
        asymptotic_p = float(chi2.sf(statistic, df=1))

    print(f"Exact McNemar p-value: {exact_p:.6f}")
    print(f"Corrected chi-square: {statistic:.6f}")
    print(f"Corrected asymptotic p-value: {asymptotic_p:.6f}")

    # Family-macro difference:
    # Calculate each family's Ridge accuracy minus baseline accuracy.
    per_family = defaultdict(list)

    for family, r_ok, b_ok in zip(
        families, ridge_correct, base_correct
    ):
        per_family[family].append(
            (float(r_ok), float(b_ok))
        )

    family_differences = []

    for values in per_family.values():
        values = np.asarray(values)
        family_differences.append(
            float(values[:, 0].mean() - values[:, 1].mean())
        )

    family_differences = np.asarray(family_differences)

    print("\n========== FAMILY-MACRO ANALYSIS ==========")
    print(f"Families: {len(family_differences)}")
    print(
        "Mean family accuracy difference: "
        f"{family_differences.mean():+.4%}"
    )

    # Resample whole families, not individual cycles.
    rng = np.random.default_rng(SEED)
    bootstrap = np.empty(BOOTSTRAP_REPS)

    for i in range(BOOTSTRAP_REPS):
        sample = rng.choice(
            family_differences,
            size=len(family_differences),
            replace=True,
        )
        bootstrap[i] = sample.mean()

    lower, upper = np.quantile(
        bootstrap, [0.025, 0.975]
    )

    print(
        "Family-cluster bootstrap 95% CI: "
        f"[{lower:+.4%}, {upper:+.4%}]"
    )


    # Bootstrap pooled accuracy difference by resampling
    # entire structural families, preserving their cycles.

    unique_families = np.unique(families)

    family_indices = {
        family: np.flatnonzero(families == family)
        for family in unique_families
    }

    pooled_bootstrap = np.empty(BOOTSTRAP_REPS)
    pooled_rng = np.random.default_rng(SEED + 1)

    for i in range(BOOTSTRAP_REPS):
        sampled_families = pooled_rng.choice(
            unique_families,
            size=len(unique_families),
            replace=True,
        )

        sampled_indices = np.concatenate([
            family_indices[family]
            for family in sampled_families
        ])

        pooled_bootstrap[i] = (
            ridge_correct[sampled_indices].mean()
            - base_correct[sampled_indices].mean()
        )

    pooled_lower, pooled_upper = np.quantile(
        pooled_bootstrap,
        [0.025, 0.975],
    )

    print("\n========== POOLED FAMILY-CLUSTER BOOTSTRAP ==========")
    print(
        f"Pooled difference: "
        f"{ridge_correct.mean() - base_correct.mean():+.4%}"
    )
    print(
        f"Family-cluster bootstrap 95% CI: "
        f"[{pooled_lower:+.4%}, {pooled_upper:+.4%}]"
    )

    print("\nInterpretation:")
    print(
        "The exact McNemar test is a paired-cycle diagnostic. "
        "The family bootstrap addresses family grouping, but its "
        "precision is limited by the small number of families."
    )


if __name__ == "__main__":
    main()
