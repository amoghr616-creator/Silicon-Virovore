#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import json
import platform
import sqlite3
import sys

from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote


ROOT = Path.cwd()
P4 = ROOT / "prototype4"

SKIP_DIRS = {
    ".git",
    ".venv",
    "__pycache__",
    "node_modules",
}

DATA_EXTENSIONS = {
    ".csv",
    ".tsv",
    ".sqlite",
    ".sqlite3",
    ".db",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)

    return digest.hexdigest()


def quote_identifier(value: str) -> str:
    return '"' + value.replace('"', '""') + '"'


def detect_delimiter(path: Path) -> str:
    with path.open(
        "r",
        encoding="utf-8-sig",
        errors="replace",
        newline="",
    ) as handle:
        sample = handle.read(16384)

    if not sample.strip():
        return "," if path.suffix.lower() == ".csv" else "\t"

    try:
        dialect = csv.Sniffer().sniff(
            sample,
            delimiters=",;\t|",
        )
        return dialect.delimiter

    except csv.Error:
        first_line = sample.splitlines()[0]

        counts = {
            delimiter: first_line.count(delimiter)
            for delimiter in (",", ";", "\t", "|")
        }

        selected = max(counts, key=counts.get)

        if counts[selected] == 0:
            raise ValueError(
                f"Could not determine delimiter for {path}"
            )

        return selected


def inspect_csv(path: Path) -> dict:
    delimiter = detect_delimiter(path)

    with path.open(
        "r",
        encoding="utf-8-sig",
        errors="replace",
        newline="",
    ) as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header = next(reader, [])
        row_count = sum(1 for _ in reader)

    normalized_headers = {
        str(column).strip().lower(): str(column)
        for column in header
    }

    group_candidates = [
        "sequence_similarity_group",
        "family_id",
        "primary_family",
        "family_group",
        "group_id",
    ]

    cycle_candidates = [
        "cycle_uid",
        "cycle_id",
        "cycle_uid",
        "cycleid",
    ]

    group_column = next(
        (
            normalized_headers[name]
            for name in group_candidates
            if name in normalized_headers
        ),
        None,
    )

    cycle_column = next(
        (
            normalized_headers[name]
            for name in cycle_candidates
            if name in normalized_headers
        ),
        None,
    )

    group_counts = {}
    missing_group_count = None
    duplicate_cycle_count = None

    if group_column or cycle_column:
        groups = set()
        cycle_ids = set()

        missing_groups = 0
        duplicate_cycles = 0

        with path.open(
            "r",
            encoding="utf-8-sig",
            errors="replace",
            newline="",
        ) as handle:
            reader = csv.DictReader(
                handle,
                delimiter=delimiter,
            )

            for row in reader:
                if group_column:
                    value = (row.get(group_column) or "").strip()

                    if value:
                        groups.add(value)
                    else:
                        missing_groups += 1

                if cycle_column:
                    value = (row.get(cycle_column) or "").strip()

                    if value:
                        if value in cycle_ids:
                            duplicate_cycles += 1
                        cycle_ids.add(value)

        if group_column:
            group_counts = {
                "group_column": group_column,
                "distinct_groups": len(groups),
                "rows_missing_group": missing_groups,
            }

        if cycle_column:
            duplicate_cycle_count = {
                "cycle_column": cycle_column,
                "distinct_cycle_ids": len(cycle_ids),
                "duplicate_cycle_rows": duplicate_cycles,
            }

    return {
        "path": str(path.relative_to(ROOT)),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "delimiter": repr(delimiter),
        "header": header,
        "data_rows": row_count,
        "group_summary": group_counts,
        "cycle_id_summary": duplicate_cycle_count,
    }


def inspect_sqlite(path: Path) -> dict:
    uri = f"file:{quote(str(path.resolve()))}?mode=ro"

    result = {
        "path": str(path.relative_to(ROOT)),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "tables": {},
        "foreign_key_violations": [],
        "integrity_check": None,
        "inspection_error": None,
    }

    try:
        connection = sqlite3.connect(uri, uri=True)
        connection.row_factory = sqlite3.Row

        try:
            table_rows = connection.execute(
                """
                SELECT name
                FROM sqlite_master
                WHERE type = 'table'
                  AND name NOT LIKE 'sqlite_%'
                ORDER BY name
                """
            ).fetchall()

            for table_row in table_rows:
                table = table_row["name"]
                quoted = quote_identifier(table)

                columns = [
                    dict(row)
                    for row in connection.execute(
                        f"PRAGMA table_info({quoted})"
                    )
                ]

                foreign_keys = [
                    dict(row)
                    for row in connection.execute(
                        f"PRAGMA foreign_key_list({quoted})"
                    )
                ]

                row_count = connection.execute(
                    f"SELECT COUNT(*) FROM {quoted}"
                ).fetchone()[0]

                result["tables"][table] = {
                    "row_count": row_count,
                    "columns": columns,
                    "foreign_keys": foreign_keys,
                }

            result["foreign_key_violations"] = [
                list(row)
                for row in connection.execute(
                    "PRAGMA foreign_key_check"
                ).fetchall()
            ]

            result["integrity_check"] = [
                row[0]
                for row in connection.execute(
                    "PRAGMA quick_check"
                ).fetchall()
            ]

        finally:
            connection.close()

    except Exception as exc:
        result["inspection_error"] = (
            f"{type(exc).__name__}: {exc}"
        )

    return result


def main() -> int:
    if not P4.is_dir():
        print(f"[ERROR] Could not find: {P4}")
        print("Run this script from the Silicon-Virovore repository root.")
        return 2

    files = []

    for path in P4.rglob("*"):
        if not path.is_file():
            continue

        if any(part in SKIP_DIRS for part in path.parts):
            continue

        if path.suffix.lower() not in DATA_EXTENSIONS:
            continue

        files.append(path)

    files.sort()

    csv_results = []
    database_results = []
    failures = []

    for path in files:
        try:
            if path.suffix.lower() in {".csv", ".tsv"}:
                csv_results.append(inspect_csv(path))
            else:
                database_results.append(inspect_sqlite(path))

        except Exception as exc:
            failures.append({
                "path": str(path.relative_to(ROOT)),
                "error": f"{type(exc).__name__}: {exc}",
            })

    report = {
        "report_type": "P4_STEP2_READ_ONLY_PREFLIGHT",
        "created_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "repository": str(ROOT),
        "python_version": sys.version,
        "platform": platform.platform(),
        "sqlite_version": sqlite3.sqlite_version,
        "discovered_file_count": len(files),
        "csv_files": csv_results,
        "sqlite_databases": database_results,
        "inspection_failures": failures,
        "database_mutations_performed": False,
        "synthetic_data_created": False,
    }

    output_dir = P4 / "results" / "p4_step2_preflight"
    output_dir.mkdir(parents=True, exist_ok=True)

    output_path = output_dir / "p4_step2_preflight.json"

    output_path.write_text(
        json.dumps(report, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    print("\n========== P4 STEP 2 PREFLIGHT ==========")

    print(f"Repository: {ROOT}")
    print(f"Discovered files: {len(files)}")
    print(f"CSV/TSV files: {len(csv_results)}")
    print(f"SQLite databases: {len(database_results)}")
    print(f"Inspection failures: {len(failures)}")

    print("\n--- CSV INVENTORY ---")

    for item in csv_results:
        print(
            f"\n{item['path']}\n"
            f"  Rows: {item['data_rows']}\n"
            f"  Delimiter: {item['delimiter']}\n"
            f"  Columns: {item['header']}"
        )

        if item["group_summary"]:
            print(f"  Groups: {item['group_summary']}")

        if item["cycle_id_summary"]:
            print(f"  Cycle IDs: {item['cycle_id_summary']}")

    print("\n--- SQLITE INVENTORY ---")

    for database in database_results:
        print(f"\nDATABASE: {database['path']}")

        if database["inspection_error"]:
            print(f"  ERROR: {database['inspection_error']}")
            continue

        for table, details in database["tables"].items():
            columns = [
                column["name"]
                for column in details["columns"]
            ]

            print(f"\n  TABLE: {table}")
            print(f"    Rows: {details['row_count']}")
            print(f"    Columns: {columns}")
            print(
                "    Foreign keys: "
                f"{details['foreign_keys']}"
            )

        print(
            "  Foreign-key violations: "
            f"{len(database['foreign_key_violations'])}"
        )

        print(
            "  Integrity check: "
            f"{database['integrity_check']}"
        )

    if failures:
        print("\n--- INSPECTION FAILURES ---")

        for failure in failures:
            print(failure)

    print("\n=========================================")
    print(f"JSON report: {output_path}")

    print(
        "\nNo SQLite tables were modified. "
        "No optimizer parameters were changed."
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
