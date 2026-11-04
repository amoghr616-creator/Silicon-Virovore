from pathlib import Path
import sqlite3
import json
import hashlib
from collections import Counter

ROOT = Path("prototype4")

DATABASES = [
    ROOT / "results/model_a_validation_20261030/p4_model_a_working.sqlite",
    ROOT / "data/p4_evidence.sqlite",
]

METRICS_JSON = (
    ROOT / "results/model_a_validation_20261030/"
    "model_a_outputs/p4_validation_metrics.json"
)


def section(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def load_json(value):
    if not value:
        return None
    try:
        return json.loads(value)
    except (TypeError, json.JSONDecodeError):
        return value


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(1024 * 1024), b""
        ):
            digest.update(block)
    return digest.hexdigest()


def inspect_metrics_file():
    section("1. MODEL A METRICS JSON")

    if not METRICS_JSON.exists():
        print("NOT FOUND:", METRICS_JSON)
        return

    raw = METRICS_JSON.read_text(encoding="utf-8")

    print("Path:", METRICS_JSON)
    print("SHA256:", sha256(METRICS_JSON))

    try:
        obj = json.loads(raw)
        print(json.dumps(obj, indent=2, sort_keys=True))
    except json.JSONDecodeError as exc:
        print("JSON parse error:", exc)
        print(raw[:10000])


def inspect_database(path):
    section(f"2. DATABASE: {path}")

    if not path.exists():
        print("NOT FOUND")
        return

    print("Size:", path.stat().st_size, "bytes")
    print("SHA256:", sha256(path))

    uri = f"file:{path.resolve()}?mode=ro"
    con = sqlite3.connect(uri, uri=True)

    tables = [
        row[0]
        for row in con.execute("""
            SELECT name
            FROM sqlite_master
            WHERE type='table'
            ORDER BY name
        """)
    ]

    print("\nTables:", tables)

    if "validation_manifest" in tables:
        section("VALIDATION MANIFEST")

        columns = [
            row[1]
            for row in con.execute(
                'PRAGMA table_info("validation_manifest")'
            )
        ]

        print("Columns:", columns)

        rows = con.execute(
            'SELECT * FROM "validation_manifest"'
        ).fetchall()

        print("Manifest rows:", len(rows))

        for index, row in enumerate(rows, start=1):
            record = dict(zip(columns, row))

            print(f"\nManifest #{index}")

            for key in (
                "manifest_id",
                "manifest_uid",
                "created_utc",
                "created_at",
                "input_dataset_sha256",
                "script_sha256",
                "split_unit",
                "family_count",
                "cycle_count",
                "integrity_note",
            ):
                if key in record:
                    print(f"{key}: {record[key]}")

            for key in ("config_json", "metrics_json", "software_json"):
                if key in record:
                    print(f"\n{key}:")
                    print(json.dumps(
                        load_json(record[key]),
                        indent=2,
                        sort_keys=True
                    ))

    if "cycle_features" in tables:
        section("CYCLE FEATURE FAMILY DISTRIBUTION")

        columns = [
            row[1]
            for row in con.execute(
                'PRAGMA table_info("cycle_features")'
            )
        ]

        print("Columns:", columns)

        if "family_id" in columns:
            total = con.execute(
                'SELECT COUNT(*) FROM cycle_features'
            ).fetchone()[0]

            null_count = con.execute("""
                SELECT COUNT(*)
                FROM cycle_features
                WHERE family_id IS NULL
                   OR TRIM(family_id) = ''
            """).fetchone()[0]

            print("Total feature/cycle rows:", total)
            print("Missing family IDs:", null_count)

            counts = con.execute("""
                SELECT family_id, COUNT(*)
                FROM cycle_features
                GROUP BY family_id
                ORDER BY COUNT(*) DESC, family_id
            """).fetchall()

            print("Distinct family IDs:", len(counts))

            for family, count in counts:
                print(f"  {family}: {count}")

            unique_cycles = (
                "cycle_uid"
                if "cycle_uid" in columns
                else "cycle_uid"
                if "cycle_uid" in columns
                else None
            )

            if unique_cycles:
                unique_count = con.execute(
                    'SELECT COUNT(DISTINCT cycle_uid) '
                    'FROM cycle_features'
                ).fetchone()[0]
                print("Distinct cycle UIDs:", unique_count)

    section("PREDICTION / FOLD / VALIDATION TABLES")

    candidate_tables = [
        table for table in tables
        if any(
            term in table.lower()
            for term in (
                "predict", "loco", "fold",
                "validation", "prediction", "result"
            )
        )
    ]

    if not candidate_tables:
        print("No separately named prediction/fold tables found.")

    for table in candidate_tables:
        escaped = table.replace('"', '""')

        columns = [
            row[1]
            for row in con.execute(
                f'PRAGMA table_info("{escaped}")'
            )
        ]

        print(f"\nTABLE: {table}")
        print("Columns:", columns)

        rows = con.execute(
            f'SELECT * FROM "{escaped}" LIMIT 3'
        ).fetchall()

        for row in rows:
            print("  ", dict(zip(columns, row)))

    con.close()


def main():
    inspect_metrics_file()

    for path in DATABASES:
        inspect_database(path)

    section("4. INTERPRETATION CHECKLIST")
    print("1. Confirm split_unit is actually family.")
    print("2. Confirm all 242 cycles have usable family IDs.")
    print("3. Confirm there are 19 distinct held-out groups.")
    print("4. Recover every baseline and model metric.")
    print("5. Determine whether per-cycle out-of-fold predictions exist.")
    print("6. Determine whether the family labels are IF30 or CATH H-level.")
    print("7. Do not call IF30 validation CATH-superfamily LOFO.")

    section("AUDIT COMPLETE")
    print("Read-only inspection; no files or databases modified.")


if __name__ == "__main__":
    main()
