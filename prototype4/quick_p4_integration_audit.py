from __future__ import annotations

import ast
import json
import os
import sqlite3
from pathlib import Path

ROOT = Path.cwd()

DATABASES = [
    ROOT / "prototype4/data/p4_evidence.sqlite",
    ROOT / "prototype4/data/p4_closed_loop_registry.sqlite",
    ROOT / "prototype4/results/model_a_validation_20261030/p4_model_a_working.sqlite",
]

SCAN_DIRS = [
    ROOT / "prototype4",
    ROOT / "src",
    ROOT / "tests",
]

SKIP_DIRS = {
    ".venv", "venv", ".git", "__pycache__", ".pytest_cache",
    "results", "node_modules", "data",
}

TERMS = (
    "target", "pdb", "chain", "residue", "site",
    "structure", "feature", "score", "rank", "optim",
    "dock", "reliab", "entropy", "fingerprint",
    "converg", "epistasis", "cycle", "mutant",
    "skempi", "evsi", "budget", "ood",
)


def inspect_database(path):
    print("\n" + "=" * 75)
    print("DATABASE:", path.relative_to(ROOT))

    if not path.exists():
        print("MISSING")
        return

    print("Size:", path.stat().st_size, "bytes")

    if path.stat().st_size == 0:
        print("EMPTY FILE")
        return

    try:
        uri = f"file:{path.resolve().as_posix()}?mode=ro"
        with sqlite3.connect(uri, uri=True, timeout=3) as con:
            objects = con.execute("""
                SELECT type, name, sql
                FROM sqlite_master
                WHERE name NOT LIKE 'sqlite_%'
                ORDER BY type, name
            """).fetchall()

            for obj_type, name, sql in objects:
                print(f"\n{obj_type.upper()}: {name}")
                if sql:
                    print("  SQL:", " ".join(sql.split())[:500])

                if obj_type == "table":
                    quoted = '"' + name.replace('"', '""') + '"'
                    columns = con.execute(
                        f"PRAGMA table_info({quoted})"
                    ).fetchall()

                    print("  Columns:", [
                        f"{c[1]} ({c[2]})" for c in columns
                    ])

                    try:
                        count = con.execute(
                            f"SELECT COUNT(*) FROM {quoted}"
                        ).fetchone()[0]
                        print("  Rows:", count)
                    except sqlite3.Error as exc:
                        print("  Count error:", exc)

    except sqlite3.Error as exc:
        print("DATABASE ERROR:", exc)


def inspect_source():
    print("\n" + "=" * 75)
    print("PYTHON SOURCE INVENTORY")

    seen = set()
    count = 0

    for base in SCAN_DIRS:
        print(f"\nDIRECTORY: {base.relative_to(ROOT)}")
        if not base.exists():
            print("  DOES NOT EXIST")
            continue

        for current, dirs, filenames in os.walk(base):
            dirs[:] = [
                d for d in dirs
                if d not in SKIP_DIRS
            ]

            for filename in filenames:
                if not filename.endswith(".py"):
                    continue

                path = Path(current) / filename
                if path in seen:
                    continue
                seen.add(path)
                count += 1

                try:
                    source = path.read_text(errors="replace")
                    tree = ast.parse(source)
                except (OSError, SyntaxError):
                    continue

                definitions = []
                for node in ast.walk(tree):
                    if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
                        if any(term in node.name.lower() for term in TERMS):
                            definitions.append(
                                f"{node.name}@{node.lineno}"
                            )

                if definitions:
                    print(
                        f"  {path.relative_to(ROOT)}: "
                        + ", ".join(definitions[:20])
                    )

    print("\nTotal Python files inspected:", count)


def inspect_preregistration():
    path = ROOT / "prototype4/configs/p4_preregistration.json"
    print("\n" + "=" * 75)
    print("P4 PREREGISTRATION")

    if not path.exists():
        print("MISSING")
        return

    try:
        data = json.loads(path.read_text())
        print(json.dumps(data, indent=2)[:10000])
    except Exception as exc:
        print("READ ERROR:", exc)


def main():
    inspect_source()

    for path in DATABASES:
        inspect_database(path)

    inspect_preregistration()

    print("\n========== AUDIT COMPLETE ==========")


if __name__ == "__main__":
    main()
