from pathlib import Path
from collections import Counter
import sqlite3
import json
import math
import csv

import numpy as np
from scipy.stats import binomtest

ROOT = Path("prototype4/results/model_a_validation_20261030")
DB = ROOT / "p4_model_a_working.sqlite"
METRICS = ROOT / "model_a_outputs/p4_validation_metrics.json"
OUT = ROOT / "p4_directional_analysis_20261101"

EXPECTED_N = 242
EXPECTED_GROUPS = 19
THRESHOLD = 0.5
TOL = 1e-8
BOOTSTRAPS = 20000
SEED = 20261101


def write_csv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def exact_group_swap_test(group_differences):
    """
    Exact paired family-block swap test.

    Each group's Ridge and baseline labels are exchanged as a block.
    The statistic is the total net number of correct directional
    predictions across informative cycles.
    """
    distribution = Counter({0: 1})

    for difference in group_differences:
        updated = Counter()

        for score, count in distribution.items():
            updated[score + difference] += count
            updated[score - difference] += count

        distribution = updated

    observed = sum(group_differences)
    extreme = sum(
        count for score, count in distribution.items()
        if abs(score) >= abs(observed)
    )

    total = 2 ** len(group_differences)

    return {
        "observed_net_correct": int(observed),
        "extreme_permutations": int(extreme),
        "total_permutations": int(total),
        "two_sided_p": extreme / total,
        "test": "exact paired IF30-family block-swap permutation",
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    if not DB.exists() or not METRICS.exists():
        raise FileNotFoundError(
            f"Required file missing:\n{DB}\n{METRICS}"
        )

    metrics = json.loads(METRICS.read_text(encoding="utf-8"))

    con = sqlite3.connect(
        f"file:{DB.resolve()}?mode=ro",
        uri=True
    )
    con.row_factory = sqlite3.Row

    result = con.execute("""
        SELECT
            p.prediction_uid,
            p.feature_uid,
            p.heldout_group,
            p.split_unit,
            p.train_groups_json,
            p.fold_index,
            p.y_true,
            p.y_pred,
            p.baseline_pred,
            p.calibrated_abs_error,
            f.cycle_uid,
            f.family_id,
            f.experimental_epistasis_kcal_mol
        FROM loco_predictions p
        JOIN cycle_features f
          ON p.feature_uid = f.feature_uid
        ORDER BY p.heldout_group, p.feature_uid
    """).fetchall()

    rows = [dict(r) for r in result]
    con.close()

    if len(rows) != EXPECTED_N:
        raise SystemExit(
            f"STOP: expected {EXPECTED_N} joined rows, got {len(rows)}"
        )

    groups = sorted({r["family_id"] for r in rows})

    if len(groups) != EXPECTED_GROUPS:
        raise SystemExit(
            f"STOP: expected {EXPECTED_GROUPS} groups, got {len(groups)}"
        )

    problems = []

    for r in rows:
        if not r["family_id"]:
            problems.append("Missing family ID")
            continue

        if r["split_unit"] != "family":
            problems.append(
                f"Unexpected split unit: {r['split_unit']}"
            )

        if r["heldout_group"] != r["family_id"]:
            problems.append(
                f"Held-out/family mismatch: {r['feature_uid']}"
            )

        if not math.isclose(
            float(r["y_true"]),
            float(r["experimental_epistasis_kcal_mol"]),
            rel_tol=0.0,
            abs_tol=TOL
        ):
            problems.append(
                f"Target mismatch: {r['feature_uid']}"
            )

        try:
            train_groups = set(json.loads(r["train_groups_json"]))
        except Exception:
            problems.append(
                f"Invalid training group list: {r['feature_uid']}"
            )
            continue

        expected_train = set(groups) - {r["heldout_group"]}

        if train_groups != expected_train:
            problems.append(
                f"Training-group isolation mismatch: "
                f"{r['heldout_group']}"
            )

        for field in ("y_true", "y_pred", "baseline_pred"):
            if not math.isfinite(float(r[field])):
                problems.append(
                    f"Non-finite value {field}: {r['feature_uid']}"
                )

    if problems:
        print("=== VALIDATION FAILURES ===")
        for problem in problems[:30]:
            print(problem)
        raise SystemExit(
            f"STOP: {len(problems)} validation problems."
        )

    print("PASS: 242 joined predictions")
    print("PASS: 19 distinct IF30 groups")
    print("PASS: held-out group matches feature family")
    print("PASS: training groups exclude held-out group")
    print("PASS: all targets match source cycle outcomes")

    # Restrict directional analysis to the pre-defined informative subset.
    informative = [
        r for r in rows
        if abs(float(r["y_true"])) >= THRESHOLD
    ]

    if not informative:
        raise SystemExit("No informative cycles at the specified threshold.")

    for r in informative:
        true_sign = np.sign(float(r["y_true"]))
        ridge_sign = np.sign(float(r["y_pred"]))
        base_sign = np.sign(float(r["baseline_pred"]))

        r["ridge_correct"] = int(ridge_sign == true_sign)
        r["baseline_correct"] = int(base_sign == true_sign)
        r["direction_delta"] = (
            r["ridge_correct"] - r["baseline_correct"]
        )

    n = len(informative)
    ridge_correct = sum(r["ridge_correct"] for r in informative)
    baseline_correct = sum(r["baseline_correct"] for r in informative)

    ridge_accuracy = ridge_correct / n
    baseline_accuracy = baseline_correct / n
    pooled_delta = ridge_accuracy - baseline_accuracy

    # Paired, observation-level exact McNemar test.
    ridge_only = sum(
        r["ridge_correct"] == 1 and r["baseline_correct"] == 0
        for r in informative
    )
    baseline_only = sum(
        r["ridge_correct"] == 0 and r["baseline_correct"] == 1
        for r in informative
    )

    discordant = ridge_only + baseline_only

    mcnemar_p = (
        binomtest(
            min(ridge_only, baseline_only),
            n=discordant,
            p=0.5,
            alternative="two-sided"
        ).pvalue
        if discordant else 1.0
    )

    # Aggregate paired outcomes within each held-out family.
    family_rows = []

    for group in groups:
        subset = [
            r for r in informative
            if r["heldout_group"] == group
        ]

        if not subset:
            continue

        n_group = len(subset)
        r_correct = sum(x["ridge_correct"] for x in subset)
        b_correct = sum(x["baseline_correct"] for x in subset)

        family_rows.append({
            "IF30_group": group,
            "n_informative": n_group,
            "ridge_correct": r_correct,
            "baseline_correct": b_correct,
            "ridge_accuracy": r_correct / n_group,
            "baseline_accuracy": b_correct / n_group,
            "delta_percentage_points":
                100 * (r_correct - b_correct) / n_group,
            "net_correct_difference": r_correct - b_correct,
        })

    family_differences = [
        row["net_correct_difference"]
        for row in family_rows
    ]

    family_swap = exact_group_swap_test(family_differences)

    # Family-cluster bootstrap. Resample whole held-out groups,
    # never individual cycles independently.
    rng = np.random.default_rng(SEED)
    n_groups_with_data = len(family_rows)

    pooled_boot = np.empty(BOOTSTRAPS)
    macro_boot = np.empty(BOOTSTRAPS)

    for b in range(BOOTSTRAPS):
        selected = rng.integers(
            0, n_groups_with_data,
            size=n_groups_with_data
        )

        sampled = [family_rows[i] for i in selected]

        total_n = sum(x["n_informative"] for x in sampled)
        total_ridge = sum(x["ridge_correct"] for x in sampled)
        total_baseline = sum(x["baseline_correct"] for x in sampled)

        pooled_boot[b] = (
            total_ridge - total_baseline
        ) / total_n

        macro_boot[b] = np.mean([
            x["ridge_accuracy"] - x["baseline_accuracy"]
            for x in sampled
        ])

    pooled_ci = np.quantile(pooled_boot, [0.025, 0.975])
    macro_ci = np.quantile(macro_boot, [0.025, 0.975])

    # Summarize results.
    summary = {
        "analysis": "P4 directional prediction evaluation",
        "split_unit": "leave-one-IF30-group-out",
        "not_cath_H_superfamily_LOFO": True,
        "informative_threshold_abs_epsilon_kcal_mol": THRESHOLD,
        "n_total_cycles": len(rows),
        "n_informative_cycles": n,
        "n_groups_total": len(groups),
        "n_groups_with_informative_cycles": n_groups_with_data,
        "ridge_correct": ridge_correct,
        "baseline_correct": baseline_correct,
        "ridge_direction_accuracy": ridge_accuracy,
        "baseline_direction_accuracy": baseline_accuracy,
        "absolute_improvement_percentage_points": 100 * pooled_delta,
        "ridge_only_correct": ridge_only,
        "baseline_only_correct": baseline_only,
        "observation_level_mcnemar_exact_p": float(mcnemar_p),
        "family_block_swap_test": family_swap,
        "family_cluster_bootstrap": {
            "replicates": BOOTSTRAPS,
            "seed": SEED,
            "pooled_delta_95_ci": pooled_ci.tolist(),
            "macro_delta_95_ci": macro_ci.tolist(),
            "note": (
                "Percentile intervals from resampling whole IF30 groups; "
                "only 19 groups are available."
            ),
        },
        "family_balanced_delta_percentage_points":
            100 * float(np.mean([
                x["ridge_accuracy"] - x["baseline_accuracy"]
                for x in family_rows
            ])),
        "stored_model_mae": metrics[
            "sequence_mutation_descriptor_model"
        ]["mae_kcal_mol"],
        "stored_baseline_mae": metrics[
            "baseline_training_fold_mean"
        ]["mae_kcal_mol"],
        "stored_risk_coverage": metrics["risk_coverage"],
        "stored_uncertainty_calibration": metrics[
            "uncertainty_calibration"
        ],
    }

    (OUT / "p4_directional_statistics.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8"
    )

    write_csv(
        OUT / "p4_directional_by_if30.csv",
        family_rows,
        [
            "IF30_group",
            "n_informative",
            "ridge_correct",
            "baseline_correct",
            "ridge_accuracy",
            "baseline_accuracy",
            "delta_percentage_points",
            "net_correct_difference",
        ]
    )

    # Export informative cycle-level paired outcomes.
    write_csv(
        OUT / "p4_directional_paired_predictions.csv",
        informative,
        [
            "prediction_uid",
            "feature_uid",
            "cycle_uid",
            "heldout_group",
            "y_true",
            "y_pred",
            "baseline_pred",
            "ridge_correct",
            "baseline_correct",
            "direction_delta",
        ]
    )

    try:
        import matplotlib.pyplot as plt

        # Figure 1: directional accuracy by model.
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.bar(
            ["Training-fold mean", "Descriptor Ridge"],
            [100 * baseline_accuracy, 100 * ridge_accuracy]
        )
        ax.set_ylabel("Directional accuracy (%)")
        ax.set_ylim(0, 100)
        ax.set_title(
            f"P4: Directional prediction on {n} informative cycles\n"
            "Leave-one-IF30-group-out evaluation"
        )
        ax.grid(axis="y", alpha=0.25)
        fig.tight_layout()
        fig.savefig(
            OUT / "p4_directional_accuracy.png",
            dpi=220
        )
        plt.close(fig)

        # Figure 2: family-specific paired differences.
        ordered = sorted(
            family_rows,
            key=lambda x: x["delta_percentage_points"]
        )

        fig, ax = plt.subplots(figsize=(9, 7))
        labels = [
            f"{x['IF30_group']} (n={x['n_informative']})"
            for x in ordered
        ]
        values = [
            x["delta_percentage_points"] for x in ordered
        ]

        ax.barh(labels, values)
        ax.axvline(0, linestyle="--")
        ax.set_xlabel(
            "Ridge minus baseline directional accuracy (percentage points)"
        )
        ax.set_title(
            "P4: Directional advantage varies across held-out IF30 groups"
        )
        ax.grid(axis="x", alpha=0.25)
        fig.tight_layout()
        fig.savefig(
            OUT / "p4_directional_family_differences.png",
            dpi=220
        )
        plt.close(fig)

        print("\nFigures:")
        print(OUT / "p4_directional_accuracy.png")
        print(OUT / "p4_directional_family_differences.png")

    except ImportError:
        print("\nMatplotlib unavailable; numerical exports completed.")

    print("\n=== DIRECTIONAL RESULTS ===")
    print(json.dumps(summary, indent=2))

    print("\n=== OUTPUT ===")
    print(OUT)


if __name__ == "__main__":
    main()
