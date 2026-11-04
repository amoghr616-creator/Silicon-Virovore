#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import json
import sqlite3

from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path.cwd()

BASE = ROOT / "prototype4" / "results" / "skempi_audit_20261029"

DB_PATH = BASE / "p4_evidence.sqlite"
ELIGIBLE_PATH = BASE / "eligible_cycles.csv"
FAMILY_PATH = BASE / "family_map.csv"
MUTATION_PATH = BASE / "mutation_chain_audit.csv"

OUTPUT_DIR = (
    ROOT / "prototype4" / "results" / "p4_step2_model_a_audit"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)

    return digest.hexdigest()


def read_csv(path: Path) -> list[dict]:
    if not path.is_file():
        raise FileNotFoundError(path)

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def parse_json(value):
    if value is None or value == "":
        return None

    try:
        return json.loads(value)
    except (TypeError, json.JSONDecodeError):
        return "__INVALID_JSON__"


def truthy(value) -> bool:
    return str(value).strip().lower() in {
        "true", "1", "yes", "y"
    }


def summarize_json(value, max_items=12):
    parsed = parse_json(value)

    if parsed is None:
        return {"type": "missing"}

    if parsed == "__INVALID_JSON__":
        return {"type": "invalid_json"}

    result = {
        "type": type(parsed).__name__,
    }

    if isinstance(parsed, dict):
        result["length"] = len(parsed)
        result["keys"] = list(parsed.keys())[:max_items]

    elif isinstance(parsed, list):
        result["length"] = len(parsed)
        result["preview"] = parsed[:max_items]

    else:
        result["value"] = parsed

    return result


def inspect_feature_alignment(rows):
    invalid_names = 0
    invalid_values = 0
    missing_names = 0
    missing_values = 0
    aligned_lists = 0
    aligned_dicts = 0
    mismatched_pairs = 0
    unknown_shapes = 0

    examples = []

    for row in rows:
        names = parse_json(row["feature_names_json"])
        values = parse_json(row["feature_values_json"])

        if names is None:
            missing_names += 1

        if values is None:
            missing_values += 1

        if names == "__INVALID_JSON__":
            invalid_names += 1

        if values == "__INVALID_JSON__":
            invalid_values += 1

        if (
            isinstance(names, list)
            and isinstance(values, list)
        ):
            if len(names) == len(values):
                aligned_lists += 1
            else:
                mismatched_pairs += 1

        elif (
            isinstance(names, dict)
            and isinstance(values, dict)
        ):
            if set(names) == set(values):
                aligned_dicts += 1
            else:
                mismatched_pairs += 1

        elif names is not None and values is not None:
            unknown_shapes += 1

        if len(examples) < 3:
            examples.append({
                "feature_uid": row["feature_uid"],
                "feature_names": summarize_json(
                    row["feature_names_json"]
                ),
                "feature_values": summarize_json(
                    row["feature_values_json"]
                ),
            })

    return {
        "rows": len(rows),
        "missing_names": missing_names,
        "missing_values": missing_values,
        "invalid_names_json": invalid_names,
        "invalid_values_json": invalid_values,
        "aligned_list_pairs": aligned_lists,
        "aligned_dictionary_pairs": aligned_dicts,
        "mismatched_pairs": mismatched_pairs,
        "unknown_shapes": unknown_shapes,
        "examples": examples,
    }


def inspect_unique_parent_key(connection, table, column):
    table_q = '"' + table.replace('"', '""') + '"'

    indexes = connection.execute(
        f"PRAGMA index_list({table_q})"
    ).fetchall()

    inspected = []

    for index in indexes:
        index_name = index["name"]
        index_q = '"' + index_name.replace('"', '""') + '"'

        columns = [
            item["name"]
            for item in connection.execute(
                f"PRAGMA index_info({index_q})"
            ).fetchall()
        ]

        unique = bool(index["unique"])

        inspected.append({
            "index_name": index_name,
            "unique": unique,
            "columns": columns,
        })

        if unique and columns == [column]:
            return {
                "has_unique_single_column_parent_key": True,
                "indexes": inspected,
            }

    return {
        "has_unique_single_column_parent_key": False,
        "indexes": inspected,
    }


def main() -> int:
    required_paths = [
        DB_PATH,
        ELIGIBLE_PATH,
        FAMILY_PATH,
        MUTATION_PATH,
    ]

    missing = [
        str(path)
        for path in required_paths
        if not path.is_file()
    ]

    if missing:
        raise FileNotFoundError(
            "Required audit inputs are missing:\n"
            + "\n".join(missing)
        )

    eligible = read_csv(ELIGIBLE_PATH)
    family_rows = read_csv(FAMILY_PATH)
    mutation_rows = read_csv(MUTATION_PATH)

    eligible_by_id = {
        row["cycle_uid"]: row
        for row in eligible
    }

    family_by_complex = {
        row["complex_id"]: row
        for row in family_rows
    }

    duplicate_family_complexes = (
        len(family_rows) - len(family_by_complex)
    )

    mapped_cycle_families = {}
    unmapped_complexes = set()

    for cycle in eligible:
        complex_id = cycle["complex_id"]
        family_record = family_by_complex.get(complex_id)

        if family_record is None:
            unmapped_complexes.add(complex_id)
            continue

        mapped_cycle_families[cycle["cycle_uid"]] = (
            family_record.get("family_id", "").strip()
        )

    family_cycle_counts = Counter(
        value
        for value in mapped_cycle_families.values()
        if value
    )

    mutation_by_id = {
        row["cycle_uid"]: row
        for row in mutation_rows
    }

    mutation_topology_counts = Counter()
    mutation_partner_pairs = Counter()
    both_mapped_count = 0
    mutation_cycles_missing = []

    for cycle_id in eligible_by_id:
        row = mutation_by_id.get(cycle_id)

        if row is None:
            mutation_cycles_missing.append(cycle_id)
            continue

        if truthy(row.get("both_mutations_mapped")):
            both_mapped_count += 1

        partner_i = (
            row.get("mutation_i_partner") or "MISSING"
        ).strip()

        partner_j = (
            row.get("mutation_j_partner") or "MISSING"
        ).strip()

        mutation_partner_pairs[
            (partner_i, partner_j)
        ] += 1

        is_cross = truthy(
            row.get("cross_partner_cycle")
        )

        if is_cross:
            mutation_topology_counts["cross_partner"] += 1
        elif partner_i != "MISSING" and partner_j != "MISSING":
            if partner_i == partner_j:
                mutation_topology_counts[
                    f"same_partner:{partner_i}"
                ] += 1
            else:
                mutation_topology_counts[
                    "partner_disagreement_without_cross_flag"
                ] += 1
        else:
            mutation_topology_counts["unclassified"] += 1

    uri = (
        f"file:{DB_PATH.resolve()}?mode=ro"
    )

    connection = sqlite3.connect(uri, uri=True)
    connection.row_factory = sqlite3.Row

    try:
        connection.execute("PRAGMA foreign_keys = ON")

        features = connection.execute(
            """
            SELECT
                feature_uid,
                cycle_uid,
                family_id,
                complex_id,
                node_uid_p0,
                node_uid_pi,
                node_uid_pj,
                node_uid_pij,
                feature_names_json,
                feature_values_json,
                experimental_epistasis_kcal_mol
            FROM cycle_features
            ORDER BY cycle_uid
            """
        ).fetchall()

        nodes = connection.execute(
            """
            SELECT
                node_uid,
                cycle_uid,
                node_state,
                source_row_uid,
                complex_id,
                context_uid,
                reference,
                temperature_k,
                measurement_method,
                kd_m,
                delta_g_kcal_mol
            FROM cycle_nodes
            ORDER BY cycle_uid, node_state
            """
        ).fetchall()

        raw_status = connection.execute(
            """
            SELECT import_status, COUNT(*) AS n
            FROM skempi_raw
            GROUP BY import_status
            ORDER BY import_status
            """
        ).fetchall()

        raw_hashes = connection.execute(
            """
            SELECT COUNT(DISTINCT source_sha256)
            FROM skempi_raw
            """
        ).fetchone()[0]

        raw_missing_source = connection.execute(
            """
            SELECT COUNT(*)
            FROM skempi_raw
            WHERE source_sha256 IS NULL
               OR source_row_number IS NULL
            """
        ).fetchone()[0]

        db_family_counts = Counter(
            (row["family_id"] or "").strip()
            for row in features
        )

        feature_cycle_counts = Counter(
            row["cycle_uid"]
            for row in features
        )

        family_null_rows = sum(
            1 for row in features
            if not (row["family_id"] or "").strip()
        )

        family_disagreements = []

        for row in features:
            expected = mapped_cycle_families.get(
                row["cycle_uid"]
            )
            observed = (row["family_id"] or "").strip()

            if expected is not None and observed:
                if expected != observed:
                    family_disagreements.append({
                        "cycle_uid": row["cycle_uid"],
                        "database_family": observed,
                        "family_map_family": expected,
                    })

        node_states_by_cycle = defaultdict(list)

        for row in nodes:
            node_states_by_cycle[row["cycle_uid"]].append(
                str(row["node_state"]).strip().upper()
            )

        expected_states = {"P0", "PI", "PJ", "PIJ"}

        incomplete_cycles = []

        for cycle_id in eligible_by_id:
            states = node_states_by_cycle.get(cycle_id, [])
            state_set = set(states)

            if (
                len(states) != 4
                or state_set != expected_states
            ):
                incomplete_cycles.append({
                    "cycle_uid": cycle_id,
                    "node_count": len(states),
                    "states": states,
                })

        eligible_ids = set(eligible_by_id)
        feature_ids = set(feature_cycle_counts)
        node_ids = set(node_states_by_cycle)

        feature_uid_duplicates = (
            len(features)
            - len({row["feature_uid"] for row in features})
        )

        cycle_feature_duplicates = {
            key: count
            for key, count in feature_cycle_counts.items()
            if count != 1
        }

        feature_summary = inspect_feature_alignment(
            features
        )

        parent_key_audit = inspect_unique_parent_key(
            connection,
            "cycle_features",
            "cycle_uid",
        )

        structure_fk = connection.execute(
            """
            PRAGMA foreign_key_list(structure_runs)
            """
        ).fetchall()

        foreign_key_violations = [
            list(row)
            for row in connection.execute(
                "PRAGMA foreign_key_check"
            ).fetchall()
        ]

        integrity = [
            row[0]
            for row in connection.execute(
                "PRAGMA quick_check"
            ).fetchall()
        ]

    finally:
        connection.close()

    family_map_summary = {
        "complexes_in_family_map": len(family_rows),
        "duplicate_complex_rows": duplicate_family_complexes,
        "distinct_family_groups": len(family_cycle_counts),
        "cycles_with_mapped_complex": len(
            mapped_cycle_families
        ),
        "cycles_with_unmapped_complex": len(
            unmapped_complexes
        ),
        "cycles_by_family": dict(
            sorted(family_cycle_counts.items())
        ),
        "database_family_null_rows": family_null_rows,
        "database_family_counts_including_blank": dict(
            db_family_counts
        ),
        "database_vs_csv_family_disagreements": (
            family_disagreements
        ),
    }

    node_summary = {
        "database_node_rows": len(nodes),
        "expected_node_rows": 4 * len(eligible),
        "eligible_cycles_missing_feature_rows": sorted(
            eligible_ids - feature_ids
        ),
        "feature_cycles_not_in_eligible_csv": sorted(
            feature_ids - eligible_ids
        ),
        "eligible_cycles_missing_node_rows": sorted(
            eligible_ids - node_ids
        ),
        "node_cycles_not_in_eligible_csv": sorted(
            node_ids - eligible_ids
        ),
        "incomplete_or_invalid_four_state_cycles": (
            incomplete_cycles
        ),
        "cycle_feature_duplicate_counts": (
            cycle_feature_duplicates
        ),
        "duplicate_feature_uid_rows": feature_uid_duplicates,
    }

    mutation_summary = {
        "mutation_audit_rows": len(mutation_rows),
        "eligible_cycles_missing_mutation_audit": (
            mutation_cycles_missing
        ),
        "both_mutations_mapped": both_mapped_count,
        "topology_counts": dict(mutation_topology_counts),
        "partner_pair_counts": {
            f"{i} -> {j}": count
            for (i, j), count in mutation_partner_pairs.items()
        },
    }

    context_summary = {
        "distinct_contexts": len({
            row["context_uid"] for row in eligible
        }),
        "distinct_references": len({
            row["reference"] for row in eligible
        }),
        "distinct_methods": len({
            row["method"] for row in eligible
        }),
        "distinct_temperatures_K": sorted({
            row["temperature_k"] for row in eligible
        }),
        "method_counts": dict(Counter(
            row["method"] for row in eligible
        )),
        "reference_counts": dict(Counter(
            row["reference"] for row in eligible
        )),
    }

    report = {
        "report_type": "P4_STEP2_MODEL_A_SOURCE_AUDIT",
        "created_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "inputs": {
            str(path.relative_to(ROOT)): sha256_file(path)
            for path in required_paths
        },
        "source_inventory": {
            "eligible_cycles_csv_rows": len(eligible),
            "family_map_csv_rows": len(family_rows),
            "mutation_audit_csv_rows": len(mutation_rows),
        },
        "database_inventory": {
            "cycle_feature_rows": len(features),
            "cycle_node_rows": len(nodes),
            "raw_import_status_counts": {
                str(row["import_status"]): row["n"]
                for row in raw_status
            },
            "distinct_raw_source_hashes": raw_hashes,
            "raw_rows_missing_source_provenance": (
                raw_missing_source
            ),
        },
        "family_mapping": family_map_summary,
        "four_state_node_integrity": node_summary,
        "mutation_topology": mutation_summary,
        "assay_context": context_summary,
        "feature_payload_alignment": feature_summary,
        "structure_runs_schema": {
            "foreign_keys": [
                dict(row) for row in structure_fk
            ],
            "cycle_uid_parent_key_audit": (
                parent_key_audit
            ),
        },
        "sqlite_integrity": {
            "foreign_key_violations": (
                foreign_key_violations
            ),
            "quick_check": integrity,
        },
        "database_modified": False,
        "synthetic_data_created": False,
    }

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    output_path = (
        OUTPUT_DIR / "p4_step2_model_a_source_audit.json"
    )

    output_path.write_text(
        json.dumps(report, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    print("\n========== MODEL A SOURCE AUDIT ==========")

    print("\n--- CYCLE AND NODE INTEGRITY ---")
    print(f"Eligible cycles: {len(eligible)}")
    print(f"Feature rows: {len(features)}")
    print(f"Node rows: {len(nodes)}")
    print(
        "Eligible cycles missing feature rows: "
        f"{len(node_summary['eligible_cycles_missing_feature_rows'])}"
    )
    print(
        "Eligible cycles missing node rows: "
        f"{len(node_summary['eligible_cycles_missing_node_rows'])}"
    )
    print(
        "Incomplete four-state cycles: "
        f"{len(incomplete_cycles)}"
    )
    print(
        "Duplicate feature rows per cycle: "
        f"{len(cycle_feature_duplicates)}"
    )

    print("\n--- FAMILY MAPPING ---")
    print(
        "Mapped cycles: "
        f"{len(mapped_cycle_families)} / {len(eligible)}"
    )
    print(
        "Distinct mapped groups: "
        f"{len(family_cycle_counts)}"
    )
    print(
        "Database family_id NULL rows: "
        f"{family_null_rows}"
    )
    print(
        "Database/CSV family disagreements: "
        f"{len(family_disagreements)}"
    )

    print("\nCycles per family:")

    for family_id, count in sorted(
        family_cycle_counts.items(),
        key=lambda item: (-item[1], item[0]),
    ):
        print(f"  {family_id}: {count}")

    print("\n--- MUTATION TOPOLOGY ---")
    print(
        f"Both mutations mapped: {both_mapped_count} "
        f"/ {len(eligible)}"
    )
    print(
        "Topology counts: "
        f"{dict(mutation_topology_counts)}"
    )
    print(
        "Missing mutation audit rows: "
        f"{len(mutation_cycles_missing)}"
    )

    print("\n--- ASSAY CONTEXT ---")
    print(
        f"Contexts: {context_summary['distinct_contexts']}"
    )
    print(
        f"References: {context_summary['distinct_references']}"
    )
    print(
        f"Methods: {context_summary['distinct_methods']}"
    )
    print(
        "Temperatures (K): "
        f"{context_summary['distinct_temperatures_K']}"
    )

    print("\n--- FEATURE JSON ---")
    print(
        json.dumps(
            {
                key: value
                for key, value in feature_summary.items()
                if key != "examples"
            },
            indent=2,
        )
    )

    print("\nFeature examples:")
    print(
        json.dumps(
            feature_summary["examples"],
            indent=2,
        )
    )

    print("\n--- STRUCTURE RUNS FOREIGN KEY ---")
    print(
        "Unique cycle_uid parent key exists: "
        f"{parent_key_audit['has_unique_single_column_parent_key']}"
    )

    print("\n--- SQLITE ---")
    print(f"Quick check: {integrity}")
    print(
        f"Foreign-key violations: "
        f"{len(foreign_key_violations)}"
    )

    print("\n==========================================")
    print(f"Report: {output_path}")
    print("Database modified: NO")
    print("Synthetic data created: NO")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
