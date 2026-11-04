from pathlib import Path
import ast
import csv
import sqlite3

ROOT = Path("prototype4")
RESULTS = ROOT / "results"
AUDIT = RESULTS / "cath_boundary_domain_map.tsv"
ASSIGNER = ROOT / "assign_skempi_cath_superfamilies.py"

KEYWORDS = (
    "cath", "boundary", "domain", "interface",
    "residue", "family", "superfamily"
)


def section(title):
    print("\n" + "=" * 76)
    print(title)
    print("=" * 76)


def inspect_audit():
    section("1. VERIFIED BOUNDARY AUDIT")

    if not AUDIT.exists():
        print("Missing:", AUDIT)
        return

    with AUDIT.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    print("Rows:", len(rows))

    matched = [
        row for row in rows
        if row["classification_status"] == "CLASSIFIED"
    ]
    unresolved = [
        row for row in rows
        if row["classification_status"] != "CLASSIFIED"
    ]

    print("Exact classification matches:", len(matched))
    print("Unresolved classifications:", len(unresolved))

    print("\nDomain-level classifications:")
    for row in rows:
        print(
            f"  {row['domain_id']:10s} "
            f"CATH={row['cath_code'] or 'UNASSIGNED':16s} "
            f"S35={row['s35_cluster'] or 'UNASSIGNED':10s} "
            f"segments={row['segment_count']} "
            f"status={row['classification_status']}"
        )

    print("\nChain-level CATH H-level diversity:")
    by_chain = {}

    for row in rows:
        by_chain.setdefault(row["chain_id"], set())

        if row["cath_code"]:
            parts = row["cath_code"].split(".")
            if len(parts) == 4:
                by_chain[row["chain_id"]].add(
                    ".".join(parts)
                )

    for chain, families in sorted(by_chain.items()):
        print(
            f"  {chain}: "
            + (
                ", ".join(sorted(families))
                if families else "NO CLASSIFIED DOMAINS"
            )
        )

    print(
        "\nNOTE: CATH H-level IDs and S35 IDs are separate "
        "classification fields."
    )


def inspect_assigner():
    section("2. PRODUCTION ASSIGNMENT SCRIPT")

    if not ASSIGNER.exists():
        print("Not found:", ASSIGNER)
        return

    source = ASSIGNER.read_text(
        encoding="utf-8",
        errors="replace"
    )
    lines = source.splitlines()

    print("Path:", ASSIGNER)
    print("Total lines:", len(lines))

    try:
        tree = ast.parse(source)
    except SyntaxError as exc:
        print("AST parsing failed:", exc)
        return

    nodes = []

    for node in ast.walk(tree):
        if isinstance(
            node,
            (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)
        ):
            name = node.name.lower()

            if any(keyword in name for keyword in KEYWORDS):
                nodes.append(node)

    nodes.sort(key=lambda node: node.lineno)

    print("\nRelevant function/class definitions:")

    for node in nodes:
        print(
            f"  {node.lineno:4}-{node.end_lineno:4} "
            f"{node.name}"
        )

    # Print compact source windows around relevant definitions
    # and important CATH input/output references.
    windows = []

    for node in nodes:
        windows.append((
            max(1, node.lineno - 2),
            min(len(lines), node.end_lineno + 2)
        ))

    for i, line in enumerate(lines, start=1):
        low = line.lower()

        if (
            "cath-domain-boundaries" in low
            or "cath-domain-list" in low
            or "boundary_without_classification" in low
            or "assigned_h" in low
            or "no_cath_h" in low
            or "interface_residue_records" in low
        ):
            windows.append((
                max(1, i - 3),
                min(len(lines), i + 5)
            ))

    windows.sort()

    merged = []

    for start, end in windows:
        if merged and start <= merged[-1][1] + 2:
            merged[-1] = (
                merged[-1][0],
                max(merged[-1][1], end)
            )
        else:
            merged.append((start, end))

    print("\nRelevant source excerpts:")

    printed = 0

    for start, end in merged:
        if printed >= 220:
            print("  ... output capped at 220 source lines")
            break

        print(f"\n--- lines {start}-{end} ---")

        for number in range(start, end + 1):
            print(f"{number:4}: {lines[number - 1]}")

        printed += end - start + 1


def inspect_result_files():
    section("3. EXISTING CATH / FAMILY OUTPUT FILES")

    if not RESULTS.exists():
        print("No results directory:", RESULTS)
        return

    terms = (
        "cath", "family", "assignment",
        "interface", "skempi", "validation"
    )

    found = []

    for path in RESULTS.rglob("*"):
        if not path.is_file():
            continue

        if any(term in path.name.lower() for term in terms):
            found.append(path)

    for path in sorted(found)[:50]:
        print(
            f"{path} | {path.stat().st_size:,} bytes"
        )

        if path.suffix.lower() in (".tsv", ".csv"):
            try:
                with path.open(
                    encoding="utf-8",
                    errors="replace",
                    newline=""
                ) as handle:
                    header = next(csv.reader(
                        handle,
                        delimiter="\t" if path.suffix == ".tsv"
                        else ","
                    ), [])

                print("    columns:", header[:30])
            except Exception as exc:
                print("    header read failed:", exc)

    if not found:
        print("No matching output filenames found.")

    if len(found) > 50:
        print("Output capped at 50 files.")


def inspect_databases():
    section("4. SQLITE DATABASE SCHEMAS")

    candidates = []

    for path in ROOT.rglob("*"):
        if (
            path.is_file()
            and path.suffix.lower() in (
                ".db", ".sqlite", ".sqlite3"
            )
            and ".git" not in path.parts
        ):
            candidates.append(path)

    if not candidates:
        print("No SQLite databases found under", ROOT)
        return

    for path in sorted(candidates)[:8]:
        print(f"\nDATABASE: {path}")

        try:
            uri = f"file:{path.resolve()}?mode=ro"
            con = sqlite3.connect(uri, uri=True)

            tables = con.execute("""
                SELECT name
                FROM sqlite_master
                WHERE type='table'
                ORDER BY name
            """).fetchall()

            for (table,) in tables:
                low = table.lower()

                if not any(
                    word in low for word in (
                        "cath", "family", "interface",
                        "cycle", "validation", "skempi",
                        "structure"
                    )
                ):
                    continue

                print(f"\n  TABLE: {table}")

                columns = con.execute(
                    f'PRAGMA table_info("{table}")'
                ).fetchall()

                names = [column[1] for column in columns]
                print("  COLUMNS:", names)

                sample = con.execute(
                    f'SELECT * FROM "{table}" LIMIT 1'
                ).fetchone()

                if sample is not None:
                    print("  SAMPLE:")
                    for name, value in zip(names, sample):
                        print(
                            f"    {name}: "
                            f"{str(value)[:180]}"
                        )

            con.close()

        except Exception as exc:
            print("  DATABASE INSPECTION ERROR:", exc)

    if len(candidates) > 8:
        print("\nDatabase listing capped at 8 files.")


def main():
    inspect_audit()
    inspect_assigner()
    inspect_result_files()
    inspect_databases()

    section("5. END OF READ-ONLY AUDIT")
    print("No production files or databases were modified.")


if __name__ == "__main__":
    main()
