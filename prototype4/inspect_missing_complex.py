from __future__ import annotations

import ast
import json
import re
from pathlib import Path

ROOT = Path.cwd()

SOURCE_FILES = [
    ROOT / "src/p4_target_integration.py",
    ROOT / "run_pipeline.py",
]

REPORTS = [
    ROOT / (
        "results/p4_live_gate_smoke_20261102/"
        "p4_live_gate_smoke_20261102/p4_target_evidence/"
        "generation_0001.json"
    ),
    ROOT / (
        "results/p4_live_gate_smoke_20261102/"
        "p4_live_gate_smoke_20261102/p4_integrated_assessment/"
        "generation_0001.json"
    ),
]

TERMS = re.compile(
    r"missing_complex|complex_path|complex_pdb|complex_file|"
    r"complex_results|target_evidence|target_status|"
    r"target_geometry|sequence_hash|CandidateRanker|"
    r"p4_target_integration|HADDOCK|generation",
    re.IGNORECASE,
)


def print_source_structure(path: Path):
    print("\n" + "=" * 90)
    print("SOURCE STRUCTURE:", path.relative_to(ROOT))

    if not path.is_file():
        print("FILE MISSING")
        return

    source = path.read_text(encoding="utf-8")
    lines = source.splitlines()

    try:
        tree = ast.parse(source)
    except SyntaxError as exc:
        print("AST parse error:", exc)
        tree = None

    if tree:
        print("\nCLASSES AND FUNCTIONS")
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef):
                print(f"  class {node.name} (line {node.lineno})")
            elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                args = [
                    arg.arg for arg in node.args.args
                ]
                print(
                    f"  def {node.name}({', '.join(args)}) "
                    f"(line {node.lineno})"
                )

    # Print merged source windows around relevant code.
    hits = [
        i for i, line in enumerate(lines)
        if TERMS.search(line)
    ]

    intervals = []
    for index in hits:
        start = max(0, index - 7)
        end = min(len(lines), index + 11)

        if intervals and start <= intervals[-1][1] + 2:
            intervals[-1] = (intervals[-1][0], max(intervals[-1][1], end))
        else:
            intervals.append((start, end))

    print("\nRELEVANT SOURCE CONTEXT")
    if not intervals:
        print("  No matching terms found.")
        return

    for start, end in intervals[:22]:
        print(f"\n  --- Lines {start + 1}-{end} ---")
        for i in range(start, end):
            print(f"  {i + 1:5d}: {lines[i]}")

    if len(intervals) > 22:
        print(f"\n  ... {len(intervals) - 22} additional source blocks omitted.")


def candidate_records(data):
    candidates = data.get("candidates", [])
    if isinstance(candidates, dict):
        candidates = candidates.get("items", [])
    return candidates if isinstance(candidates, list) else []


def relevant_nested_fields(obj, prefix="candidate"):
    """Collect path, status, error, complex, and evidence-related fields."""
    found = []

    if isinstance(obj, dict):
        for key, value in obj.items():
            path = f"{prefix}.{key}"
            k = key.lower()

            relevant = any(term in k for term in (
                "path", "file", "complex", "status", "reason",
                "error", "source", "evidence", "structure",
            ))

            if relevant and not isinstance(value, (dict, list)):
                found.append((path, value))
            elif isinstance(value, (dict, list)):
                found.extend(relevant_nested_fields(value, path))

    elif isinstance(obj, list):
        for i, value in enumerate(obj):
            found.extend(relevant_nested_fields(value, f"{prefix}[{i}]"))

    return found


def inspect_report(path: Path):
    print("\n" + "=" * 90)
    print("GENERATION REPORT:", path.relative_to(ROOT))

    if not path.is_file():
        print("REPORT MISSING")
        return

    data = json.loads(path.read_text(encoding="utf-8"))
    print("Generation:", data.get("generation"))
    print("Run ID:", data.get("run_id"))
    print("Candidate count:", data.get("candidate_count"))
    print("Target evidence count:", data.get("target_evidence_count"))
    print("Target evidence coverage:", data.get("target_evidence_coverage"))
    print("Target selection reason:", data.get("reason"))
    print("Manifest:", data.get("manifest", data.get("target_manifest")))

    candidates = candidate_records(data)
    print("\nCandidate records:", len(candidates))

    status_counts = {}
    for candidate in candidates:
        status = (
            candidate.get("target_status")
            or candidate.get("status")
            or "<NO STATUS>"
        )
        status_counts[status] = status_counts.get(status, 0) + 1

    print("Candidate status counts:", status_counts)

    for i, candidate in enumerate(candidates):
        if not isinstance(candidate, dict):
            print(f"\nCandidate {i}: non-dictionary record")
            continue

        print(f"\n--- Candidate {i + 1} ---")

        for key in (
            "candidate_id", "id", "name", "sequence_hash",
            "sequence", "target_status", "status",
            "target_structure_id", "target_chain",
            "target_geometry_score", "target_gate_valid",
            "target_eligible", "complex_count",
        ):
            if key in candidate:
                value = candidate[key]
                if key == "sequence" and isinstance(value, str):
                    value = f"{value[:30]}... (length={len(value)})"
                print(f"  {key}: {value!r}")

        print("  All candidate keys:", sorted(candidate.keys()))

        fields = relevant_nested_fields(candidate)
        for field, value in fields[:45]:
            if isinstance(value, str) and len(value) > 220:
                value = value[:220] + "..."
            print(f"  {field}: {value!r}")

        if len(fields) > 45:
            print(f"  ... {len(fields) - 45} additional nested fields omitted.")


def main():
    print("MISSING-COMPLEX DIAGNOSTIC")
    print("Project root:", ROOT)

    for path in SOURCE_FILES:
        print_source_structure(path)

    for path in REPORTS:
        inspect_report(path)

    print("\n" + "=" * 90)
    print("Additional project-wide missing_complex references")

    for directory in (ROOT / "src", ROOT / "prototype4"):
        if not directory.exists():
            continue

        for path in directory.rglob("*.py"):
            if any(x in path.parts for x in (".venv", "__pycache__", "results")):
                continue

            try:
                lines = path.read_text(encoding="utf-8").splitlines()
            except (OSError, UnicodeDecodeError):
                continue

            for i, line in enumerate(lines):
                if "missing_complex" in line.lower():
                    print(f"\n{path.relative_to(ROOT)}:{i + 1}")
                    for j in range(max(0, i - 5), min(len(lines), i + 8)):
                        print(f"  {j + 1:5d}: {lines[j]}")


if __name__ == "__main__":
    main()
