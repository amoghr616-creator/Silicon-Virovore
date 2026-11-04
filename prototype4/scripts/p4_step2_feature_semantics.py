#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import sqlite3

from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path.cwd()

BASE = ROOT / "prototype4" / "results" / "skempi_audit_20261029"
DB_PATH = BASE / "p4_evidence.sqlite"
ELIGIBLE_PATH = BASE / "eligible_cycles.csv"
FAMILY_PATH = BASE / "family_map.csv"

OUTPUT_DIR = (
    ROOT / "prototype4" / "results" / "p4_step2_feature_semantics"
)

# kcal mol^-1 K^-1
R_KCAL = 0.00198720425864083

EPSILON_TOLERANCE = 1e-6

SEARCH_TERMS = re.compile(
    r"feature_names_json|feature_values_json|"
    r"sub_[A-Z]{2}|feature_names|feature_values",
    re.IGNORECASE,
)

SKIP_PARTS = {
    ".git",
    ".venv",
    "venv",
    "__pycache__",
    "node_modules",
    "results",
    "data",
}

SOURCE_EXTENSIONS = {
    ".py", ".c", ".h", ".cpp", ".hpp", ".sql"
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)

    return digest.hexdigest()


def read_csv(path: Path) -> list[dict]:
    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def parse_json(value):
    try:
        return json.loads(value)
    except (TypeError, json.JSONDecodeError):
        return None


def as_float(value):
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError):
        return None

    if not math.isfinite(result):
        return None

    return result


def safe_log_ratio(kd0, kdi, kdj, kdij, temperature):
    values = [kd0, kdi, kdj, kdij, temperature]

    if any(value is None for value in values):
        return None

    if any(value <= 0 for value in values):
        return None

    # epsilon = RT ln(Kd_Pij * Kd_P0 / (Kd_Pi * Kd_Pj))
    ratio = (kdij * kd0) / (kdi * kdj)

    if not math.isfinite(ratio) or ratio <= 0:
        return None

    return R_KCAL * temperature * math.log(ratio)


def summarize_numeric(values):
    values = [
        float(value)
        for value in values
        if value is not None and math.isfinite(float(value))
    ]

    if not values:
        return {
            "count": 0,
            "min": None,
            "median": None,
            "max": None,
            "mean": None,
        }

    ordered = sorted(values)
    n = len(ordered)

    if n % 2:
        median = ordered[n // 2]
    else:
        median = (
            ordered[n // 2 - 1] + ordered[n // 2]
        ) / 2

    return {
        "count": n,
        "min": min(ordered),
        "median": median,
        "max": max(ordered),
        "mean": sum(ordered) / n,
    }


def search_feature_source():
    hits = []

    search_roots = [
        ROOT / "prototype4",
        ROOT / "src",
    ]

    for search_root in search_roots:
        if not search_root.exists():
            continue

        for path in search_root.rglob("*"):
            if not path.is_file():
                continue

            if path.suffix.lower() not in SOURCE_EXTENSIONS:
                continue

            if any(
                part in SKIP_PARTS
                for part in path.parts
            ):
                continue

            try:
                lines = path.read_text(
                    encoding="utf-8",
                    errors="replace",
                ).splitlines()
            except OSError:
                continue

            for line_number, line in enumerate(
                lines,
                start=1,
            ):
                if SEARCH_TERMS.search(line):
                    hits.append({
                        "path": str(path.relative_to(ROOT)),
                        "line": line_number,
                        "text": line.strip()[:350],
                    })

                    if len(hits) >= 150:
                        return hits

    return hits


def main() -> int:
    for path in [DB_PATH, ELIGIBLE_PATH, FAMILY_PATH]:
        if not path.is_file():
            raise FileNotFoundError(path)

    eligible = read_csv(ELIGIBLE_PATH)
    family_rows = read_csv(FAMILY_PATH)

    eligible_by_id = {
        row["cycle_uid"]: row
        for row in eligible
    }

    family_by_complex = {
        row["complex_id"]: row
        for row in family_rows
    }

    cycle_to_family = {}

    for cycle in eligible:
        family_record = family_by_complex.get(
            cycle["complex_id"]
        )

        if family_record:
            cycle_to_family[cycle["cycle_uid"]] = (
                family_record.get("family_id", "").strip()
            )

    uri = (
        f"file:{DB_PATH.resolve()}?mode=ro"
    )

    connection = sqlite3.connect(uri, uri=True)
    connection.row_factory = sqlite3.Row

    try:
        feature_rows = connection.execute(
            """
            SELECT
                feature_uid,
                cycle_uid,
                family_id,
                complex_id,
                feature_names_json,
                feature_values_json,
                experimental_epistasis_kcal_mol
            FROM cycle_features
            ORDER BY cycle_uid
            """
        ).fetchall()

        node_rows = connection.execute(
            """
            SELECT
                cycle_uid,
                node_state,
                delta_g_kcal_mol
            FROM cycle_nodes
            """
        ).fetchall()

    finally:
        connection.close()

    # ---------------------------------------------------------
    # A. Feature vector schema and numeric integrity
    # ---------------------------------------------------------

    feature_names_reference = None
    feature_schema_mismatches = []
    feature_value_length_mismatches = []
    invalid_feature_json = []
    non_numeric_features = []
    nonfinite_features = []

    feature_vectors = {}
    nonzero_counts = []
    value_distribution = Counter()
    active_feature_counts = Counter()

    feature_row_examples = []

    for row in feature_rows:
        cycle_id = row["cycle_uid"]

        names = parse_json(row["feature_names_json"])
        values = parse_json(row["feature_values_json"])

        if not isinstance(names, list) or not isinstance(values, list):
            invalid_feature_json.append(cycle_id)
            continue

        if len(names) != len(values):
            feature_value_length_mismatches.append({
                "cycle_uid": cycle_id,
                "names_length": len(names),
                "values_length": len(values),
            })
            continue

        if feature_names_reference is None:
            feature_names_reference = list(names)
        elif names != feature_names_reference:
            feature_schema_mismatches.append(cycle_id)

        numeric_values = []
        invalid_indices = []

        for index, value in enumerate(values):
            parsed = as_float(value)

            if parsed is None:
                invalid_indices.append(index)
                continue

            numeric_values.append(parsed)

        if invalid_indices:
            non_numeric_features.append({
                "cycle_uid": cycle_id,
                "invalid_indices": invalid_indices[:20],
                "invalid_count": len(invalid_indices),
            })
            continue

        if len(numeric_values) != len(values):
            continue

        nonzero_indices = [
            index
            for index, value in enumerate(numeric_values)
            if abs(value) > 1e-12
        ]

        nonzero_names = [
            str(names[index])
            for index in nonzero_indices
        ]

        for index in nonzero_indices:
            active_feature_counts[str(names[index])] += 1

        value_distribution.update(
            "zero" if abs(value) <= 1e-12 else "nonzero"
            for value in numeric_values
        )

        feature_vectors[cycle_id] = tuple(numeric_values)
        nonzero_counts.append(len(nonzero_indices))

        if len(feature_row_examples) < 12:
            cycle = eligible_by_id.get(cycle_id, {})

            feature_row_examples.append({
                "cycle_uid": cycle_id,
                "mutation_i": cycle.get("mutation_i"),
                "mutation_j": cycle.get("mutation_j"),
                "mutation_topology": None,
                "nonzero_feature_count": len(nonzero_indices),
                "active_features": nonzero_names[:20],
            })

    feature_name_counts = Counter(
        str(name)
        for name in (feature_names_reference or [])
    )

    feature_prefix_counts = Counter()

    for name in (feature_names_reference or []):
        name = str(name)

        if "_" in name:
            prefix = name.split("_", 1)[0]
        else:
            prefix = "(no_prefix)"

        feature_prefix_counts[prefix] += 1

    duplicate_feature_names = [
        name
        for name, count in feature_name_counts.items()
        if count > 1
    ]

    vector_frequencies = Counter(
        feature_vectors.values()
    )

    # ---------------------------------------------------------
    # B. Reconstruct experimental epistasis from the four Kd
    # measurements and compare with the stored target values.
    # ---------------------------------------------------------

    epsilon_comparisons = []
    invalid_kd_cycles = []

    for cycle_id, row in eligible_by_id.items():
        kd0 = as_float(row.get("kd_p0_m"))
        kdi = as_float(row.get("kd_pi_m"))
        kdj = as_float(row.get("kd_pj_m"))
        kdij = as_float(row.get("kd_pij_m"))
        temperature = as_float(row.get("temperature_k"))

        epsilon_csv = as_float(
            row.get("epsilon_exp_kcal_mol")
        )

        epsilon_recomputed = safe_log_ratio(
            kd0, kdi, kdj, kdij, temperature
        )

        if (
            epsilon_recomputed is None
            or epsilon_csv is None
        ):
            invalid_kd_cycles.append(cycle_id)
            continue

        epsilon_comparisons.append({
            "cycle_uid": cycle_id,
            "epsilon_csv": epsilon_csv,
            "epsilon_recomputed": epsilon_recomputed,
            "absolute_difference": abs(
                epsilon_csv - epsilon_recomputed
            ),
        })

    node_delta_g = defaultdict(dict)

    for row in node_rows:
        state = str(row["node_state"]).strip().upper()

        node_delta_g[row["cycle_uid"]][state] = as_float(
            row["delta_g_kcal_mol"]
        )

    node_epsilon_comparisons = []
    invalid_node_epsilon_cycles = []

    for cycle_id, nodes in node_delta_g.items():
        expected = {"P0", "PI", "PJ", "PIJ"}

        if set(nodes) != expected:
            invalid_node_epsilon_cycles.append(cycle_id)
            continue

        if any(nodes[state] is None for state in expected):
            invalid_node_epsilon_cycles.append(cycle_id)
            continue

        epsilon_nodes = (
            nodes["PIJ"]
            - nodes["PI"]
            - nodes["PJ"]
            + nodes["P0"]
        )

        node_epsilon_comparisons.append({
            "cycle_uid": cycle_id,
            "epsilon_nodes": epsilon_nodes,
        })

    feature_target_comparisons = []

    feature_by_cycle = {
        row["cycle_uid"]: row
        for row in feature_rows
    }

    for cycle_id, comparison in {
        row["cycle_uid"]: row
        for row in epsilon_comparisons
    }.items():
        feature_row = feature_by_cycle.get(cycle_id)

        if feature_row is None:
            continue

        stored_target = as_float(
            feature_row["experimental_epistasis_kcal_mol"]
        )

        if stored_target is not None:
            feature_target_comparisons.append({
                "cycle_uid": cycle_id,
                "epsilon_csv": comparison["epsilon_csv"],
                "feature_target": stored_target,
                "absolute_difference": abs(
                    comparison["epsilon_csv"] - stored_target
                ),
            })

    # ---------------------------------------------------------
    # C. Assay and family overlap accounting
    # ---------------------------------------------------------

    context_family = defaultdict(set)
    reference_family = defaultdict(set)

    for cycle in eligible:
        family = cycle_to_family.get(cycle["cycle_uid"])

        if not family:
            continue

        context_family[cycle["context_uid"]].add(family)
        reference_family[cycle["reference"]].add(family)

    contexts_shared_across_families = {
        key: sorted(value)
        for key, value in context_family.items()
        if len(value) > 1
    }

    references_shared_across_families = {
        key: sorted(value)
        for key, value in reference_family.items()
        if len(value) > 1
    }

    # ---------------------------------------------------------
    # D. Search source code for feature-generation logic
    # ---------------------------------------------------------

    source_hits = search_feature_source()

    # ---------------------------------------------------------
    # E. Report
    # ---------------------------------------------------------

    csv_epsilon_differences = [
        row["absolute_difference"]
        for row in epsilon_comparisons
    ]

    stored_target_differences = [
        row["absolute_difference"]
        for row in feature_target_comparisons
    ]

    report = {
        "report_type": "P4_STEP2_FEATURE_SEMANTICS_AUDIT",
        "created_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "input_hashes": {
            str(path.relative_to(ROOT)): sha256_file(path)
            for path in [DB_PATH, ELIGIBLE_PATH, FAMILY_PATH]
        },
        "feature_schema": {
            "feature_rows": len(feature_rows),
            "valid_feature_vectors": len(feature_vectors),
            "feature_dimension": (
                len(feature_names_reference)
                if feature_names_reference is not None
                else 0
            ),
            "feature_names": feature_names_reference or [],
            "feature_prefix_counts": dict(feature_prefix_counts),
            "duplicate_feature_names": duplicate_feature_names,
            "rows_with_schema_mismatch": (
                feature_schema_mismatches
            ),
            "rows_with_invalid_json": invalid_feature_json,
            "rows_with_length_mismatch": (
                feature_value_length_mismatches
            ),
            "rows_with_non_numeric_values": non_numeric_features,
            "nonzero_feature_count_per_row": summarize_numeric(
                nonzero_counts
            ),
            "unique_feature_vectors": len(vector_frequencies),
            "most_common_vector_counts": [
                {
                    "frequency": count,
                    "vector_preview": list(vector[:20]),
                }
                for vector, count in vector_frequencies.most_common(5)
            ],
            "most_frequently_active_features": [
                {
                    "feature": name,
                    "active_rows": count,
                }
                for name, count in active_feature_counts.most_common(30)
            ],
            "examples": feature_row_examples,
        },
        "thermodynamic_reconstruction": {
            "gas_constant_kcal_mol_K": R_KCAL,
            "formula": (
                "RT*ln(Kd_Pij*Kd_P0/(Kd_Pi*Kd_Pj))"
            ),
            "valid_csv_recomputations": len(
                epsilon_comparisons
            ),
            "invalid_or_nonpositive_kd_cycles": (
                invalid_kd_cycles
            ),
            "csv_vs_recomputed_epsilon_difference": (
                summarize_numeric(csv_epsilon_differences)
            ),
            "cycles_over_tolerance": [
                row for row in epsilon_comparisons
                if row["absolute_difference"] > EPSILON_TOLERANCE
            ],
            "valid_node_epistasis_reconstructions": len(
                node_epsilon_comparisons
            ),
            "invalid_node_epistasis_cycles": (
                invalid_node_epsilon_cycles
            ),
            "csv_vs_database_feature_target_difference": (
                summarize_numeric(stored_target_differences)
            ),
            "feature_targets_over_tolerance": [
                row for row in feature_target_comparisons
                if row["absolute_difference"] > EPSILON_TOLERANCE
            ],
        },
        "assay_overlap": {
            "contexts_shared_across_families": (
                contexts_shared_across_families
            ),
            "references_shared_across_families": (
                references_shared_across_families
            ),
        },
        "feature_generator_source_search": {
            "hits": source_hits,
            "hit_count": len(source_hits),
        },
        "database_modified": False,
        "model_fitted": False,
        "synthetic_data_created": False,
    }

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    output_path = (
        OUTPUT_DIR / "p4_step2_feature_semantics.json"
    )

    output_path.write_text(
        json.dumps(report, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    print("\n======= P4 FEATURE SEMANTICS AUDIT =======")

    print("\n--- FEATURE MATRIX ---")
    print(
        "Feature dimension: "
        f"{report['feature_schema']['feature_dimension']}"
    )
    print(
        "Valid vectors: "
        f"{report['feature_schema']['valid_feature_vectors']}"
    )
    print(
        "Feature prefixes: "
        f"{dict(feature_prefix_counts)}"
    )
    print(
        "Nonzero features per row: "
        f"{report['feature_schema']['nonzero_feature_count_per_row']}"
    )
    print(
        "Unique vectors: "
        f"{report['feature_schema']['unique_feature_vectors']}"
    )
    print(
        "Most active features: "
        f"{report['feature_schema']['most_frequently_active_features'][:15]}"
    )

    print("\nFirst active-feature examples:")

    for example in feature_row_examples[:8]:
        print(
            f"  {example['cycle_uid']} | "
            f"{example['mutation_i']} + "
            f"{example['mutation_j']} | "
            f"{example['active_features']}"
        )

    print("\n--- THERMODYNAMIC RECONSTRUCTION ---")
    print(
        "Valid Kd recomputations: "
        f"{len(epsilon_comparisons)}"
    )
    print(
        "Invalid Kd cycles: "
        f"{len(invalid_kd_cycles)}"
    )
    print(
        "CSV epsilon vs recomputed epsilon: "
        f"{summarize_numeric(csv_epsilon_differences)}"
    )
    print(
        "Cycles exceeding tolerance: "
        f"{len(report['thermodynamic_reconstruction']['cycles_over_tolerance'])}"
    )
    print(
        "CSV epsilon vs database feature target: "
        f"{summarize_numeric(stored_target_differences)}"
    )
    print(
        "Feature target discrepancies over tolerance: "
        f"{len(report['thermodynamic_reconstruction']['feature_targets_over_tolerance'])}"
    )

    print("\n--- ASSAY OVERLAP ---")
    print(
        "Contexts shared across families: "
        f"{len(contexts_shared_across_families)}"
    )
    print(
        "References shared across families: "
        f"{len(references_shared_across_families)}"
    )

    print("\n--- FEATURE GENERATOR SOURCE ---")

    for hit in source_hits[:35]:
        print(
            f"{hit['path']}:{hit['line']}: "
            f"{hit['text']}"
        )

    print("\n==========================================")
    print(f"Full JSON report: {output_path}")
    print("Database modified: NO")
    print("Model fitted: NO")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
