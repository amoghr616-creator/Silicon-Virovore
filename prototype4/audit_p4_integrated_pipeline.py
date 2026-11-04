from __future__ import annotations

import ast
import re
import sqlite3
from collections import defaultdict
from pathlib import Path

ROOT = Path.cwd()
OUT = ROOT / "prototype4/results/p4_integration_audit_20261101"
OUT.mkdir(parents=True, exist_ok=True)

EXCLUDED = {
    ".git", ".venv", "venv", "__pycache__",
    ".pytest_cache", "results", "node_modules",
}

SOURCE_SUFFIXES = {
    ".py", ".toml", ".yaml", ".yml", ".json",
    ".cfg", ".ini", ".md",
}

CATEGORIES = {
    "TARGET CONFIGURATION": [
        r"target_config", r"TargetConfig", r"pdb_path",
        r"pdb_id", r"receptor_chain", r"target_chain",
        r"target_residues", r"target_site", r"residue_map",
        r"PDBParser", r"MMCIFParser",
    ],
    "TARGET FEATURE GENERATION": [
        r"target_specific", r"target.score", r"anchor",
        r"hotspot", r"interface", r"contact_pairs",
        r"target_window", r"structure_features",
    ],
    "OPTIMIZER INTEGRATION": [
        r"CandidateRanker", r"fitness", r"objective",
        r"optimizer", r"optimization", r"candidate_score",
        r"docking_score_is_target_specific",
        r"strongest_anchor_delta",
    ],
    "RELIABILITY": [
        r"Shannon", r"entropy", r"fingerprint",
        r"convergence", r"reproducib", r"uncertainty",
        r"calibration", r"out.of.distribution",
        r"\bOOD\b", r"failure_boundary",
    ],
    "EPISTASIS": [
        r"epistasis", r"epsilon", r"mutant_cycle",
        r"double_mutant", r"cycle_features",
        r"cycle_nodes", r"SKEMPI", r"delta_delta_g",
    ],
    "BUDGET AND EVSI": [
        r"\bEVSI\b", r"expected.value",
        r"budget", r"cost_effective",
        r"evaluation_budget", r"acquisition",
    ],
}

FUNCTION_TERMS = re.compile(
    r"target|pdb|chain|residue|site|structure|"
    r"feature|score|rank|optim|dock|reliab|"
    r"entropy|fingerprint|converg|epistasis|"
    r"cycle|mutant|skempi|evsi|budget|ood",
    re.I,
)


def excluded(path: Path) -> bool:
    return any(part in EXCLUDED for part in path.parts)


def source_files():
    files = []
    for base in (ROOT, ROOT / "src", ROOT / "prototype4", ROOT / "tests"):
        if not base.exists():
            continue
        for path in base.rglob("*"):
            if (
                path.is_file()
                and path.suffix.lower() in SOURCE_SUFFIXES
                and not excluded(path)
                and path.stat().st_size < 5_000_000
            ):
                files.append(path)
    return sorted(set(files))


def scan_source(files):
    hits = defaultdict(list)
    symbols = defaultdict(list)
    test_files = []

    for path in files:
        rel = str(path.relative_to(ROOT))

        try:
            lines = path.read_text(errors="replace").splitlines()
        except OSError:
            continue

        full_text = "\n".join(lines)

        for category, patterns in CATEGORIES.items():
            for line_no, line in enumerate(lines, 1):
                if any(re.search(p, line, re.I) for p in patterns):
                    hits[category].append(
                        (rel, line_no, line.strip()[:220])
                    )

        if path.suffix == ".py":
            try:
                tree = ast.parse(full_text)
            except SyntaxError as exc:
                symbols["PARSE ERRORS"].append(
                    (rel, exc.lineno or 0, str(exc))
                )
                continue

            for node in ast.walk(tree):
                if isinstance(
                    node,
                    (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef),
                ):
                    if FUNCTION_TERMS.search(node.name):
                        args = ""
                        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                            args = ", ".join(
                                a.arg for a in node.args.args
                            )
                        symbols["RELEVANT SYMBOLS"].append(
                            (rel, node.lineno, f"{node.name}({args})")
                        )

            if path.name.startswith("test_") or "tests" in path.parts:
                if re.search(
                    r"target|switch|epistasis|cycle|reliab|"
                    r"fingerprint|converg|pdb",
                    full_text,
                    re.I,
                ):
                    test_files.append(rel)

    return hits, symbols, sorted(set(test_files))


def find_configs():
    found = []
    for path in ROOT.rglob("*"):
        if not path.is_file() or excluded(path):
            continue
        if path.suffix.lower() not in {
            ".yaml", ".yml", ".toml", ".cfg", ".ini", ".json"
        }:
            continue
        if re.search(
            r"target|haddock|p4|optim|reliab|epistasis",
            path.name,
            re.I,
        ):
            found.append(str(path.relative_to(ROOT)))
    return sorted(found)


def inspect_databases():
    candidates = []
    for base in (
        ROOT / "prototype4" / "data",
        ROOT / "prototype4" / "results",
    ):
        if not base.exists():
            continue
        for path in base.rglob("*"):
            if (
                path.is_file()
                and path.suffix.lower() in {".sqlite", ".sqlite3", ".db"}
                and path.stat().st_size > 0
                and not excluded(path)
            ):
                candidates.append(path)

    report = []

    for path in sorted(candidates):
        try:
            uri = f"file:{path.resolve()}?mode=ro"
            with sqlite3.connect(uri, uri=True, timeout=3) as con:
                tables = con.execute(
                    """
                    SELECT name FROM sqlite_master
                    WHERE type='table' AND name NOT LIKE 'sqlite_%'
                    ORDER BY name
                    """
                ).fetchall()

                table_info = []
                for (name,) in tables:
                    columns = [
                        row[1]
                        for row in con.execute(
                            f'PRAGMA table_info("{name}")'
                        ).fetchall()
                    ]

                    count = None
                    if re.search(
                        r"cycle|skempi|manifest|structure|target|prediction",
                        name,
                        re.I,
                    ):
                        try:
                            count = con.execute(
                                f'SELECT COUNT(*) FROM "{name}"'
                            ).fetchone()[0]
                        except sqlite3.Error:
                            pass

                    table_info.append({
                        "name": name,
                        "columns": columns,
                        "rows": count,
                    })

            report.append({
                "path": str(path.relative_to(ROOT)),
                "size_bytes": path.stat().st_size,
                "tables": table_info,
            })

        except Exception as exc:
            report.append({
                "path": str(path.relative_to(ROOT)),
                "error": str(exc),
            })

    return report


def main():
    files = source_files()
    hits, symbols, test_files = scan_source(files)
    configs = find_configs()
    databases = inspect_databases()

    sections = []

    def add(title, lines):
        sections.append(f"\n{'=' * 78}\n{title}\n{'=' * 78}")
        sections.extend(lines)

    add("SOURCE INVENTORY", [
        f"Source/configuration files scanned: {len(files)}",
        *[f"  {x}" for x in files[:100]],
    ])

    for category, records in hits.items():
        add(f"{category} ({len(records)} matches)", [
            f"{path}:{line}: {text}"
            for path, line, text in records[:45]
        ] or ["No static source matches found."])

    add("RELEVANT CLASSES AND FUNCTIONS", [
        f"{path}:{line}: {symbol}"
        for path, line, symbol in symbols["RELEVANT SYMBOLS"][:160]
    ] or ["No matching Python definitions found."])

    add("PARSE ERRORS", [
        f"{path}:{line}: {message}"
        for path, line, message in symbols["PARSE ERRORS"]
    ] or ["None detected."])

    add("RELEVANT TEST FILES", [
        *test_files
    ] or ["No obviously named target/epistasis/reliability tests found."])

    add("RELEVANT CONFIGURATION FILES", [
        *configs
    ] or ["No matching configuration files found."])

    db_lines = []
    for db in databases:
        db_lines.append(f"DATABASE: {db['path']}")
        db_lines.append(f"  size: {db.get('size_bytes')}")
        if "error" in db:
            db_lines.append(f"  error: {db['error']}")
            continue
        for table in db["tables"]:
            db_lines.append(
                f"  {table['name']}: rows={table['rows']}, "
                f"columns={table['columns']}"
            )

    add("SQLITE DATABASE SCHEMAS", db_lines or ["No SQLite databases found."])

    report = "\n".join(sections)
    output = OUT / "integration_audit.txt"
    output.write_text(report)

    print(report)
    print("\nREPORT:", output)


if __name__ == "__main__":
    main()
