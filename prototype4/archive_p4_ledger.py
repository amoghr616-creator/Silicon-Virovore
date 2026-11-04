import sqlite3
import hashlib
import json
import os
import stat
from datetime import datetime, timezone
from pathlib import Path

DB_PATH = Path(
    "prototype4/results/model_a_validation_20261030/"
    "p4_model_a_working.sqlite"
)

ARCHIVE_DIR = DB_PATH.parent / "evidence_archives"


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    if not DB_PATH.is_file():
        raise FileNotFoundError(f"Database not found: {DB_PATH}")

    ARCHIVE_DIR.mkdir(parents=True, exist_ok=True)

    stamp = datetime.now(timezone.utc).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    snapshot = ARCHIVE_DIR / f"p4_ledger_{stamp}.sqlite"
    manifest = ARCHIVE_DIR / f"p4_ledger_{stamp}.json"

    if snapshot.exists() or manifest.exists():
        raise FileExistsError("Archive destination already exists.")

    # Create a consistent SQLite backup using the backup API.
    source = sqlite3.connect(DB_PATH)
    destination = sqlite3.connect(snapshot)

    try:
        source.backup(destination)
    except Exception:
        destination.close()
        source.close()
        snapshot.unlink(missing_ok=True)
        raise
    else:
        destination.close()
        source.close()

    # Reopen the snapshot read-only and verify its contents.
    read_only_uri = snapshot.resolve().as_uri() + "?mode=ro"
    conn = sqlite3.connect(read_only_uri, uri=True)

    try:
        integrity = conn.execute(
            "PRAGMA integrity_check"
        ).fetchone()[0]

        if integrity != "ok":
            raise RuntimeError(
                f"Snapshot integrity check failed: {integrity}"
            )

        tables = [
            row[0]
            for row in conn.execute(
                """
                SELECT name
                FROM sqlite_master
                WHERE type = 'table'
                  AND name NOT LIKE 'sqlite_%'
                ORDER BY name
                """
            )
        ]

        table_counts = {}

        for table in tables:
            safe_name = '"' + table.replace('"', '""') + '"'
            table_counts[table] = conn.execute(
                f"SELECT COUNT(*) FROM {safe_name}"
            ).fetchone()[0]

        family_count = None
        cycle_count = table_counts.get("cycle_features")

        if "cycle_features" in tables:
            feature_columns = {
                row[1]
                for row in conn.execute(
                    "PRAGMA table_info(cycle_features)"
                )
            }

            if "family_id" in feature_columns:
                family_count = conn.execute(
                    """
                    SELECT COUNT(DISTINCT family_id)
                    FROM cycle_features
                    """
                ).fetchone()[0]

    finally:
        conn.close()

    digest = sha256_file(snapshot)

    # Read-only filesystem permission for the archived database.
    os.chmod(
        snapshot,
        stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH
    )

    record = {
        "created_utc": stamp,
        "source_database": str(DB_PATH),
        "snapshot_database": str(snapshot),
        "snapshot_sha256": digest,
        "snapshot_size_bytes": snapshot.stat().st_size,
        "integrity_check": integrity,
        "table_counts": table_counts,
        "cycle_features_count": cycle_count,
        "distinct_family_count": family_count,
        "original_database_modified": False,
        "purpose": (
            "Read-only evidence archive created using SQLite backup API."
        ),
    }

    manifest.write_text(json.dumps(record, indent=2) + "\n")

    print("========== P4 LEDGER ARCHIVE ==========")
    print(json.dumps(record, indent=2))
    print("\nArchive and manifest created.")
    print("Original working database was not modified.")


if __name__ == "__main__":
    main()
