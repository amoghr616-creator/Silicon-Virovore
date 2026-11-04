from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import sqlite3
import sys
from pathlib import Path

SEEDS = list(range(631, 641))
ARMS = ["P1", "P2", "P3", "P4"]

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "unified_ablation_20261102"
DB = OUT / "unified_ablation_matrix.sqlite"

FIELDS = [
    "config_hash",
    "initial_pool_hash",
    "actual_evaluations",
    "mean_common_score",
    "best_common_score",
    "mutation_localization",
    "unique_sequences",
    "pairwise_identity",
    "gate_abstention_count",
    "runtime_seconds",
    "run_status",
]

NUMERIC_FIELDS = [
    "actual_evaluations",
    "mean_common_score",
    "best_common_score",
    "mutation_localization",
    "unique_sequences",
    "pairwise_identity",
    "gate_abstention_count",
    "runtime_seconds",
]


def connect():
    OUT.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(DB)
    conn.row_factory = sqlite3.Row
    return conn


def init_db():
    with connect() as conn:
        conn.execute("""
            CREATE TABLE IF NOT EXISTS unified_ablation_matrix (
                seed INTEGER NOT NULL,
                arm TEXT NOT NULL,
                config_hash TEXT,
                initial_pool_hash TEXT,
                actual_evaluations INTEGER,
                mean_common_score REAL,
                best_common_score REAL,
                mutation_localization REAL,
                unique_sequences INTEGER,
                pairwise_identity REAL,
                gate_abstention_count INTEGER,
                runtime_seconds REAL,
                run_status TEXT NOT NULL
                    CHECK(run_status IN ('PASS','INCOMPLETE','FAILED')),
                PRIMARY KEY(seed, arm)
            )
        """)

        for seed in SEEDS:
            for arm in ARMS:
                conn.execute("""
                    INSERT OR IGNORE INTO unified_ablation_matrix
                    (seed, arm, run_status)
                    VALUES (?, ?, 'INCOMPLETE')
                """, (seed, arm))

    print(f"Initialized: {DB}")
    print("40 records created or preserved.")
    print("All unexecuted records remain INCOMPLETE.")


def valid_hash(value):
    if not isinstance(value, str) or len(value) != 64:
        return False
    try:
        int(value, 16)
        return True
    except ValueError:
        return False


def validate_metrics(seed, arm, data):
    status = data.get("run_status")

    if status not in {"PASS", "INCOMPLETE", "FAILED"}:
        raise ValueError("Invalid run_status.")

    if status != "PASS":
        return

    for field in ("config_hash", "initial_pool_hash"):
        if not valid_hash(data.get(field)):
            raise ValueError(f"{field} must be a valid SHA-256 hash.")

    missing = [
        field for field in NUMERIC_FIELDS
        if data.get(field) is None
    ]
    if missing:
        raise ValueError(f"PASS record is missing: {missing}")

    for field in NUMERIC_FIELDS:
        value = data[field]
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"{field} must be numeric.")
        if not math.isfinite(float(value)):
            raise ValueError(f"{field} must be finite.")

    if data["actual_evaluations"] <= 0:
        raise ValueError("actual_evaluations must be positive.")

    if data["unique_sequences"] < 1:
        raise ValueError("unique_sequences must be positive.")

    if not 0 <= data["mutation_localization"] <= 1:
        raise ValueError("mutation_localization must be between 0 and 1.")

    if not 0 <= data["pairwise_identity"] <= 1:
        raise ValueError("pairwise_identity must be between 0 and 1.")

    if data["gate_abstention_count"] < 0:
        raise ValueError("gate_abstention_count cannot be negative.")

    if data["runtime_seconds"] < 0:
        raise ValueError("runtime_seconds cannot be negative.")

    if arm != "P4" and data["gate_abstention_count"] != 0:
        raise ValueError("P1-P3 must have zero P4 gate abstentions.")


def record_result(seed, arm, metrics_path):
    if seed not in SEEDS or arm not in ARMS:
        raise ValueError("Seed or arm is outside the registered design.")

    data = json.loads(metrics_path.read_text())
    validate_metrics(seed, arm, data)

    values = [data.get(field) for field in FIELDS]

    with connect() as conn:
        conn.execute(
            f"""
            UPDATE unified_ablation_matrix
            SET {", ".join(field + " = ?" for field in FIELDS)}
            WHERE seed = ? AND arm = ?
            """,
            values + [seed, arm],
        )

    print(f"Recorded seed={seed}, arm={arm}, status={data['run_status']}")


def exact_sign_flip(differences):
    n = len(differences)
    observed = abs(sum(differences) / n)
    extreme = 0
    total = 2 ** n

    for signs in itertools.product((-1, 1), repeat=n):
        permuted = abs(
            sum(sign * value for sign, value in zip(signs, differences)) / n
        )
        if permuted >= observed - 1e-12:
            extreme += 1

    return extreme / total


def holm_bonferroni(p_values):
    ordered = sorted(p_values.items(), key=lambda item: item[1])
    adjusted = {}
    running_max = 0.0
    count = len(ordered)

    for rank, (name, p_value) in enumerate(ordered):
        candidate = min(1.0, (count - rank) * p_value)
        running_max = max(running_max, candidate)
        adjusted[name] = running_max

    return adjusted


def analyze(budget):
    with connect() as conn:
        rows = conn.execute("""
            SELECT * FROM unified_ablation_matrix
            ORDER BY seed, arm
        """).fetchall()

    if len(rows) != 40:
        raise RuntimeError(f"Expected 40 rows; found {len(rows)}.")

    unfinished = [
        (row["seed"], row["arm"], row["run_status"])
        for row in rows
        if row["run_status"] != "PASS"
    ]

    if unfinished:
        print("ANALYSIS BLOCKED: actual results are incomplete.")
        for seed, arm, status in unfinished:
            print(f"  seed={seed}, arm={arm}, status={status}")
        print("No missing results will be imputed.")
        return 2

    data = {
        (row["seed"], row["arm"]): dict(row)
        for row in rows
    }

    for seed in SEEDS:
        block = [data[(seed, arm)] for arm in ARMS]

        pool_hashes = {row["initial_pool_hash"] for row in block}
        evaluation_counts = {row["actual_evaluations"] for row in block}

        if len(pool_hashes) != 1:
            raise RuntimeError(
                f"Initial-pool mismatch in seed block {seed}."
            )

        if evaluation_counts != {budget}:
            raise RuntimeError(
                f"Budget mismatch in seed block {seed}: "
                f"{evaluation_counts}; expected {budget}."
            )

    contrasts = {
        "Delta_21": ("P1", "P2"),
        "Delta_32": ("P2", "P3"),
        "Delta_43": ("P3", "P4"),
    }

    results = {}
    p_values = {}
    paired_rows = []

    for name, (control, treatment) in contrasts.items():
        differences = []

        for seed in SEEDS:
            difference = (
                data[(seed, treatment)]["mean_common_score"]
                - data[(seed, control)]["mean_common_score"]
            )
            differences.append(difference)

            paired_rows.append({
                "contrast": name,
                "seed": seed,
                "control_arm": control,
                "treatment_arm": treatment,
                "paired_difference": difference,
            })

        p_value = exact_sign_flip(differences)
        p_values[name] = p_value

        results[name] = {
            "control": control,
            "treatment": treatment,
            "n_seed_pairs": len(differences),
            "mean_paired_difference": sum(differences) / len(differences),
            "median_paired_difference": sorted(differences)[
                len(differences) // 2
            ],
            "positive_pairs": sum(value > 0 for value in differences),
            "negative_pairs": sum(value < 0 for value in differences),
            "ties": sum(value == 0 for value in differences),
            "exact_two_sided_sign_flip_p": p_value,
        }

    adjusted = holm_bonferroni(p_values)

    for name in results:
        results[name]["holm_adjusted_p"] = adjusted[name]

    report = {
        "design": "Unified four-arm matched-seed cumulative ablation",
        "seeds": SEEDS,
        "arms": ARMS,
        "budget_per_run": budget,
        "actual_passed_runs": len(rows),
        "independent_seed_blocks": len(SEEDS),
        "primary_endpoint": "mean_common_score",
        "statistical_method": "Exact two-sided paired sign-flip permutation",
        "multiple_testing": "Holm-Bonferroni across three contrasts",
        "contrasts": results,
        "imputation": "None",
    }

    OUT.mkdir(parents=True, exist_ok=True)

    (OUT / "paired_contrast_report.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )

    with (OUT / "paired_seed_differences.csv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "contrast", "seed", "control_arm",
                "treatment_arm", "paired_difference",
            ],
        )
        writer.writeheader()
        writer.writerows(paired_rows)

    print(json.dumps(report, indent=2))
    print(f"\nSaved analysis to {OUT}")
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["init", "record", "analyze"])
    parser.add_argument("--seed", type=int)
    parser.add_argument("--arm", choices=ARMS)
    parser.add_argument("--metrics-json", type=Path)
    parser.add_argument("--budget", type=int, default=500)
    args = parser.parse_args()

    if args.action == "init":
        init_db()

    elif args.action == "record":
        if args.seed is None or args.arm is None or args.metrics_json is None:
            parser.error("record requires --seed, --arm, and --metrics-json")
        if not args.metrics_json.is_file():
            parser.error("Metrics JSON file does not exist.")
        record_result(args.seed, args.arm, args.metrics_json)

    elif args.action == "analyze":
        return analyze(args.budget)

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, RuntimeError, sqlite3.Error) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
