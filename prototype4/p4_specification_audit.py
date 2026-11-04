from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
import sqlite3
from collections import Counter
from pathlib import Path

ROOT = Path.cwd()
P4 = ROOT / "prototype4"
RESULTS = P4 / "results"
OUT = RESULTS / "p4_specification_audit_20261101"
OUT.mkdir(parents=True, exist_ok=True)

SKIP_DIRS = {
    ".git", "__pycache__", ".pytest_cache", "node_modules",
    "hybridock-pep", "site-packages", "build", "dist",
}

EXPECTED_TABLES = {
    "skempi_raw",
    "cycle_nodes",
    "structure_runs",
    "cycle_features",
    "loco_predictions",
    "validation_manifest",
}

TEXT_EXT = {
    ".py", ".json", ".jsonl", ".csv", ".md", ".txt",
    ".sql", ".yaml", ".yml", ".toml",
}

TRAJECTORY_EXT = {
    ".xtc", ".dcd", ".trr", ".nc", ".netcdf",
    ".mdcrd", ".lammpstrj", ".traj",
}


def walk_files(base: Path):
    if not base.exists():
        return

    for current, dirs, files in os.walk(base):
        dirs[:] = [
            d for d in dirs
            if d not in SKIP_DIRS
            and not d.startswith(".venv")
            and d not in {"venv", "env"}
        ]

        for filename in files:
            yield Path(current) / filename


def load_json(path: Path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def load_csv(path: Path):
    try:
        with path.open(
            encoding="utf-8-sig", newline=""
        ) as handle:
            return list(csv.DictReader(handle))
    except Exception:
        return []


def sha256(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def find_key_values(obj, wanted, prefix=""):
    found = []

    if isinstance(obj, dict):
        for key, value in obj.items():
            location = f"{prefix}.{key}" if prefix else key

            if key.lower() in wanted:
                found.append((location, value))

            found.extend(find_key_values(value, wanted, location))

    elif isinstance(obj, list):
        for index, value in enumerate(obj[:100]):
            found.extend(
                find_key_values(value, wanted, f"{prefix}[{index}]")
            )

    return found


def numeric(value):
    try:
        number = float(value)
        return number if number == number and abs(number) != float("inf") else None
    except (TypeError, ValueError):
        return None


def record(label, passed, detail):
    status = "PASS" if passed else "BLOCKED"
    checks.append({
        "requirement": label,
        "status": status,
        "detail": str(detail),
    })
    print(f"{status:8s} {label}: {detail}")


checks = []
all_files = list(walk_files(ROOT))
p4_files = [
    p for p in all_files
    if "prototype4" in p.relative_to(ROOT).parts
]
source_files = [
    p for p in all_files
    if p.suffix == ".py"
    and (
        "src" in p.relative_to(ROOT).parts
        or "prototype4" in p.relative_to(ROOT).parts
        or p.parent == ROOT
    )
]

print("\n========== P4 SPECIFICATION AUDIT ==========")
print("Project:", ROOT)
print("Project files indexed:", len(all_files))


# ============================================================
# 1. LOFO FEATURE MODEL AND UNCERTAINTY CALIBRATION
# ============================================================

MODEL_DIR = RESULTS / "model_a_validation_20261030"
MODEL_OUT = MODEL_DIR / "model_a_outputs"

LOCO = load_csv(MODEL_OUT / "loco_oof_predictions.csv")
FOLDS = load_csv(MODEL_OUT / "outer_fold_manifest.csv")
RISK = load_csv(MODEL_OUT / "risk_coverage.csv")
CALIBRATION = load_csv(
    MODEL_OUT / "uncertainty_calibration_bins.csv"
)
VALIDATION = load_json(MODEL_OUT / "p4_validation_metrics.json")

GROUPS = {
    row.get("heldout_group")
    for row in FOLDS
    if row.get("heldout_group")
}

print("\n========== LOFO MODEL ==========")
print("LOFO prediction rows:", len(LOCO))
print("Outer fold rows:", len(FOLDS))
print("Unique held-out groups:", len(GROUPS))
print("Risk-coverage rows:", len(RISK))
print("Uncertainty calibration bins:", len(CALIBRATION))

feature_dimension_evidence = []

wanted_feature_keys = {
    "feature_names",
    "feature_count",
    "n_features",
    "n_features_in",
    "feature_dimension",
    "descriptor_count",
    "n_descriptors",
}

for path in p4_files:
    if (
        path.suffix == ".json"
        and path.stat().st_size < 5_000_000
        and any(k in path.name.lower() for k in (
            "validation", "feature", "model", "manifest"
        ))
    ):
        obj = load_json(path)
        if obj is None:
            continue

        for key, value in find_key_values(obj, wanted_feature_keys):
            if isinstance(value, list):
                feature_dimension_evidence.append(
                    (str(path.relative_to(ROOT)), key, len(value))
                )
            elif isinstance(value, (int, float)):
                feature_dimension_evidence.append(
                    (str(path.relative_to(ROOT)), key, value)
                )

feature_csv_candidates = [
    p for p in p4_files
    if p.suffix == ".csv"
    and "feature" in p.name.lower()
    and p.stat().st_size < 30_000_000
]

feature_csv_info = []

for path in feature_csv_candidates:
    try:
        with path.open(
            encoding="utf-8-sig", newline=""
        ) as handle:
            reader = csv.DictReader(handle)
            columns = reader.fieldnames or []
            sample = []
            for row in reader:
                sample.append(row)
                if len(sample) >= 30:
                    break

        metadata_names = {
            "feature_uid", "cycle_uid", "node_uid",
            "complex_id", "context_uid", "protein1",
            "protein2", "family_id", "target",
            "y_true", "y_pred", "epsilon_exp_kcal_mol",
        }

        numeric_columns = []
        for col in columns:
            if col.lower() in metadata_names:
                continue
            values = [numeric(row.get(col)) for row in sample]
            if values and all(v is not None for v in values):
                numeric_columns.append(col)

        feature_csv_info.append({
            "path": str(path.relative_to(ROOT)),
            "columns": len(columns),
            "sample_numeric_columns": len(numeric_columns),
            "sample_rows": len(sample),
        })
    except Exception:
        pass

print("Feature dimension evidence:", feature_dimension_evidence[:12])
print("Feature CSV inspection:", feature_csv_info[:12])

validation_summary = {}
if isinstance(VALIDATION, dict):
    validation_summary = {
        key: VALIDATION.get(key)
        for key in (
            "n_cycles",
            "n_independent_groups",
            "split_unit",
            "model_minus_baseline_mae",
            "group_bootstrap_ci_model_minus_baseline_mae",
            "uncertainty_calibration",
            "risk_coverage",
            "validation_strength",
            "interpretation",
        )
        if key in VALIDATION
    }

print("Validation summary:")
print(json.dumps(validation_summary, indent=2, default=str)[:6000])

better_than_random_mean = 0
better_than_random_lower_ci = 0
risk_comparable = 0

for row in RISK:
    selected = numeric(row.get("uncertainty_selected_mae"))
    random_mean = numeric(row.get("random_retention_mae_mean"))
    random_low = numeric(row.get("random_retention_mae_ci95_low"))

    if selected is not None and random_mean is not None:
        risk_comparable += 1
        if selected < random_mean:
            better_than_random_mean += 1

        if random_low is not None and selected < random_low:
            better_than_random_lower_ci += 1

print(
    "Risk-coverage rows where uncertainty selection beats random mean:",
    f"{better_than_random_mean}/{risk_comparable}",
)
print(
    "Rows below random-retention 95% lower bound:",
    f"{better_than_random_lower_ci}/{risk_comparable}",
)

has_lofo = bool(LOCO) and len(GROUPS) >= 2
has_risk_curve = len(RISK) >= 3 and risk_comparable >= 3

record(
    "LOFO held-out validation",
    has_lofo,
    f"{len(LOCO)} predictions across {len(GROUPS)} held-out groups",
)
record(
    "Risk-coverage analysis",
    has_risk_curve,
    f"{len(RISK)} coverage levels; uncertainty MAE beats random mean at {better_than_random_mean}/{risk_comparable}",
)

feature_410_found = any(
    value == 410
    for _, _, value in feature_dimension_evidence
) or any(
    row["sample_numeric_columns"] == 410
    for row in feature_csv_info
)

record(
    "410-dimensional feature vector verified",
    feature_410_found,
    (
        "Exact 410 feature dimension found"
        if feature_410_found
        else "No exact 410-dimensional feature declaration verified"
    ),
)

has_calibration_bins = len(CALIBRATION) >= 3
record(
    "Uncertainty calibration diagnostics",
    has_calibration_bins,
    f"{len(CALIBRATION)} calibration bins found",
)


# ============================================================
# 2. POSE / TRAJECTORY ENSEMBLE ENTROPY
# ============================================================

fingerprints_path = (
    RESULTS
    / "haddock_run_reproducibility_20261101"
    / "pose_contact_fingerprints.csv"
)

pairwise_path = (
    RESULTS
    / "haddock_run_reproducibility_20261101"
    / "pairwise_interface_similarity.csv"
)

entropy_path = (
    RESULTS
    / "haddock_pose_audit_20261031"
    / "p4_target_scored_poses"
    / "contact_entropy_summary.csv"
)

fingerprints = load_csv(fingerprints_path)
pairwise = load_csv(pairwise_path)
entropy_rows = load_csv(entropy_path)

print("\n========== ENSEMBLE ENTROPY ==========")
print("Fingerprint records:", len(fingerprints))
print("Pairwise similarity records:", len(pairwise))
print("Entropy summary records:", len(entropy_rows))

contact_counts = Counter()
zero_contact_rows = 0

for row in fingerprints:
    run = row.get("run", "unknown")
    count = numeric(row.get("contact_pair_count"))
    if count is not None:
        contact_counts[run] += 1
        if count == 0:
            zero_contact_rows += 1

print("Fingerprint rows by run:", dict(contact_counts))
print("Rows with zero detected contacts:", zero_contact_rows)

cross_run_values = []

for row in pairwise:
    comparison = row.get("comparison", "").lower()
    similarity = numeric(row.get("jaccard_similarity"))

    if (
        similarity is not None
        and ("pilot_vs_repeat" in comparison
             or "pilot vs repeat" in comparison)
    ):
        cross_run_values.append(similarity)

print("Pilot-repeat pairwise values:", len(cross_run_values))
if cross_run_values:
    print(
        "Pilot-repeat mean Jaccard:",
        sum(cross_run_values) / len(cross_run_values),
    )
    print(
        "Pilot-repeat positive Jaccard:",
        sum(x > 0 for x in cross_run_values),
    )

trajectory_files = [
    p for p in all_files
    if p.suffix.lower() in TRAJECTORY_EXT
    and any(term in str(p).lower() for term in (
        "prototype4", "haddock", "cabs", "trajectory", "md/"
    ))
]

print("Trajectory-format files:", len(trajectory_files))
for path in trajectory_files[:20]:
    print(" ", path.relative_to(ROOT))

has_pose_entropy = bool(fingerprints) and bool(entropy_rows)
has_trajectory_files = bool(trajectory_files)

record(
    "Pose/contact entropy engine",
    has_pose_entropy,
    f"{len(fingerprints)} fingerprints and {len(entropy_rows)} entropy summaries",
)
record(
    "Post-equilibration trajectory input",
    has_trajectory_files,
    f"{len(trajectory_files)} trajectory-format files found",
)

record(
    "Nonzero pilot-repeat interface agreement",
    bool(cross_run_values) and any(x > 0 for x in cross_run_values),
    (
        f"{sum(x > 0 for x in cross_run_values)}/{len(cross_run_values)} "
        "pilot-repeat Jaccard comparisons are positive"
        if cross_run_values
        else "No pilot-repeat Jaccard comparisons found"
    ),
)


# ============================================================
# 3. EVSI / DECISION CONTROLLER / RESOURCE ALLOCATION
# ============================================================

production_sources = [
    p for p in source_files
    if "tests" not in p.relative_to(ROOT).parts
    and "test_" not in p.name
]

evsi_defs = []
evsi_calls = []
prune_defs = []
prune_calls = []
gate_defs = []
budget_references = []

for path in production_sources:
    try:
        lines = path.read_text(
            encoding="utf-8", errors="replace"
        ).splitlines()
    except OSError:
        continue

    for number, line in enumerate(lines, start=1):
        stripped = line.strip()

        if re.search(r"def\s+estimate_evsi\s*\(", line):
            evsi_defs.append((path, number, stripped))
        elif re.search(r"estimate_evsi\s*\(", line):
            evsi_calls.append((path, number, stripped))

        if re.search(r"def\s+candidate_prune_decision\s*\(", line):
            prune_defs.append((path, number, stripped))
        elif re.search(r"candidate_prune_decision\s*\(", line):
            prune_calls.append((path, number, stripped))

        if re.search(r"def\s+operator_gate\s*\(", line):
            gate_defs.append((path, number, stripped))

        if re.search(
            r"cpu_hours|gpu_hours|compute_budget|runtime_budget|"
            r"budget_remaining|resource_budget|cost_per_hour",
            line,
            re.I,
        ):
            budget_references.append((path, number, stripped))

print("\n========== DECISION-THEORETIC CONTROL ==========")
print("EVSI definitions:", len(evsi_defs))
print("EVSI references/calls:", len(evsi_calls))
print("Pruning definitions:", len(prune_defs))
print("Pruning references/calls:", len(prune_calls))
print("Operator gate definitions:", len(gate_defs))
print("Budget references:", len(budget_references))

for label, entries in (
    ("EVSI", evsi_calls),
    ("PRUNING", prune_calls),
    ("BUDGET", budget_references),
):
    print(f"\n{label} sites:")
    for path, line, content in entries[:15]:
        print(f" {path.relative_to(ROOT)}:{line}: {content[:180]}")

integration_evsi_calls = [
    item for item in evsi_calls
    if item[0].name not in {"p4_closed_loop.py"}
]

record(
    "EVSI helper implemented",
    bool(evsi_defs),
    f"{len(evsi_defs)} definition(s)",
)
record(
    "EVSI called outside its defining module",
    bool(integration_evsi_calls),
    f"{len(integration_evsi_calls)} external call site(s)",
)
record(
    "Candidate pruning function implemented",
    bool(prune_defs),
    f"{len(prune_defs)} definition(s)",
)
record(
    "Explicit compute-budget accounting",
    bool(budget_references),
    f"{len(budget_references)} source references",
)


# ============================================================
# 4. RELATIONAL LEDGER / DATABASE INTEGRITY
# ============================================================

database_files = [
    p for p in all_files
    if p.suffix.lower() in {".db", ".sqlite", ".sqlite3", ".duckdb"}
]

database_results = []

print("\n========== DATABASE LEDGER ==========")
print("Database files discovered:", len(database_files))

for db_path in database_files:
    if db_path.suffix.lower() == ".duckdb":
        database_results.append({
            "path": str(db_path.relative_to(ROOT)),
            "status": "DuckDB found; SQLite audit not applicable",
        })
        continue

    try:
        connection = sqlite3.connect(
            f"file:{db_path.resolve()}?mode=ro",
            uri=True,
        )
        objects = connection.execute(
            "SELECT type, name, sql FROM sqlite_master "
            "WHERE type IN ('table', 'trigger')"
        ).fetchall()

        tables = {
            name for kind, name, _ in objects if kind == "table"
        }

        triggers = [
            {
                "name": name,
                "sql": sql or "",
            }
            for kind, name, sql in objects
            if kind == "trigger"
        ]

        present = sorted(EXPECTED_TABLES & tables)

        immutable_trigger_tables = []
        for trigger in triggers:
            sql = trigger["sql"].upper()
            if (
                "BEFORE UPDATE" in sql
                or "BEFORE DELETE" in sql
                or "BEFORE INSERT" in sql
            ):
                immutable_trigger_tables.append(trigger["name"])

        database_results.append({
            "path": str(db_path.relative_to(ROOT)),
            "tables": sorted(tables),
            "expected_tables_present": present,
            "expected_table_count": len(present),
            "triggers": triggers,
            "immutability_trigger_names": immutable_trigger_tables,
        })

        print("\nDatabase:", db_path.relative_to(ROOT))
        print("Tables:", sorted(tables))
        print("Expected tables present:", present)
        print("Triggers:", [x["name"] for x in triggers])

        connection.close()

    except Exception as exc:
        database_results.append({
            "path": str(db_path.relative_to(ROOT)),
            "error": f"{type(exc).__name__}: {exc}",
        })

max_tables = max(
    (
        item.get("expected_table_count", 0)
        for item in database_results
    ),
    default=0,
)

has_ledger = max_tables == len(EXPECTED_TABLES)

has_immutability = any(
    item.get("immutability_trigger_names")
    for item in database_results
)

record(
    "Six-table relational ledger",
    has_ledger,
    f"Maximum expected tables co-located in one database: {max_tables}/6",
)
record(
    "Database immutability triggers",
    has_immutability,
    (
        "Update/delete protection triggers detected"
        if has_immutability
        else "No append-only trigger protection detected"
    ),
)


# ============================================================
# 5. ENTROPY INFORMATION GAIN: OUT-OF-SAMPLE DELTA RMSE
# ============================================================

information_gain_patterns = re.compile(
    r"delta.?rmse|rmse.{0,40}entropy|entropy.{0,40}rmse|"
    r"single.?pose|with.?entropy|without.?entropy|"
    r"trajectory_entropy|entropy_added",
    re.I,
)

information_gain_hits = []

for path in p4_files:
    if path.suffix.lower() not in TEXT_EXT:
        continue

    if path.stat().st_size > 3_000_000:
        continue

    try:
        content = path.read_text(
            encoding="utf-8", errors="replace"
        )
    except OSError:
        continue

    matching = []
    for line_number, line in enumerate(content.splitlines(), start=1):
        if information_gain_patterns.search(line):
            matching.append((line_number, line.strip()[:200]))

    if matching:
        information_gain_hits.append({
            "path": str(path.relative_to(ROOT)),
            "matches": matching[:8],
        })

print("\n========== ENTROPY INFORMATION GAIN ==========")
print("Files with entropy/RMSE comparison references:",
      len(information_gain_hits))

for item in information_gain_hits[:20]:
    print("\n", item["path"])
    for line_number, line in item["matches"]:
        print(f"  {line_number}: {line}")

# Require an actual paired, out-of-sample metric comparison
# artifact, not merely code comments or a single entropy value.
validated_information_gain = any(
    any(
        ("rmse" in line.lower() and "entropy" in line.lower())
        for _, line in item["matches"]
    )
    and any(
        word in item["path"].lower()
        for word in ("validation", "benchmark", "comparison", "ablation")
    )
    for item in information_gain_hits
)

record(
    "Entropy-added out-of-sample RMSE comparison",
    validated_information_gain,
    f"{len(information_gain_hits)} candidate files; statistical validation still requires source inspection",
)


# ============================================================
# 6. STRUCTURAL SENSITIVITY CONTROLS
# ============================================================

sensitivity_keywords = (
    "sensitivity",
    "distal",
    "solvent_exposed",
    "solvent-exposed",
    "pocket_probe",
    "pocket-probe",
    "distal_control",
    "mutation_probe",
)

sensitivity_artifacts = []

for path in p4_files:
    lower = str(path).lower()

    if any(word in lower for word in sensitivity_keywords):
        sensitivity_artifacts.append(str(path.relative_to(ROOT)))

print("\n========== SENSITIVITY AUDIT ==========")
print("Sensitivity-named artifacts:", len(sensitivity_artifacts))

for path in sensitivity_artifacts[:40]:
    print(" ", path)

has_named_sensitivity_audit = any(
    any(word in path.lower() for word in (
        "distal", "solvent_exposed", "solvent-exposed",
        "pocket_probe", "pocket-probe", "mutation_probe",
    ))
    for path in sensitivity_artifacts
)

record(
    "Distal-versus-pocket structural sensitivity audit",
    has_named_sensitivity_audit,
    f"{len(sensitivity_artifacts)} sensitivity-related paths; dedicated control audit not confirmed",
)


# ============================================================
# 7. TESTS / EXISTING TARGET EVIDENCE
# ============================================================

test_count = sum(
    1 for p in all_files
    if p.suffix == ".py"
    and p.name.startswith("test_")
    and "tests" in p.relative_to(ROOT).parts
)

generation_reports = [
    p for p in all_files
    if p.name.startswith("generation_")
    and p.suffix == ".json"
    and "p4_target_evidence" in p.parts
]

print("\n========== IMPLEMENTATION COVERAGE ==========")
print("Test files found:", test_count)
print("Generation target-evidence reports:", len(generation_reports))

record(
    "P4 tests present",
    test_count > 0,
    f"{test_count} test files discovered",
)
record(
    "Generation-level target evidence reports",
    bool(generation_reports),
    f"{len(generation_reports)} reports discovered",
)


# ============================================================
# 8. FINAL REPORT
# ============================================================

passed = sum(x["status"] == "PASS" for x in checks)
blocked = sum(x["status"] == "BLOCKED" for x in checks)

audit = {
    "project": "AmphiMize P4 specification audit",
    "checks_passed": passed,
    "checks_blocked": blocked,
    "checks": checks,
    "lofo": {
        "predictions": len(LOCO),
        "outer_folds": len(FOLDS),
        "unique_groups": len(GROUPS),
        "risk_coverage_rows": len(RISK),
        "calibration_bins": len(CALIBRATION),
        "feature_dimension_evidence": feature_dimension_evidence,
        "feature_csv_info": feature_csv_info,
        "validation_summary": validation_summary,
        "risk_better_than_random_mean": better_than_random_mean,
        "risk_comparison_rows": risk_comparable,
    },
    "ensemble": {
        "fingerprint_rows": len(fingerprints),
        "pairwise_rows": len(pairwise),
        "entropy_summary_rows": len(entropy_rows),
        "zero_contact_rows": zero_contact_rows,
        "pilot_repeat_jaccard_values": cross_run_values,
        "trajectory_files": [
            str(x.relative_to(ROOT))
            for x in trajectory_files
        ],
    },
    "decision_control": {
        "evsi_definitions": len(evsi_defs),
        "evsi_calls": len(evsi_calls),
        "external_evsi_calls": len(integration_evsi_calls),
        "pruning_definitions": len(prune_defs),
        "pruning_calls": len(prune_calls),
        "budget_references": len(budget_references),
    },
    "databases": database_results,
    "information_gain_candidates": information_gain_hits,
    "sensitivity_artifacts": sensitivity_artifacts,
}

json_path = OUT / "p4_specification_audit.json"
json_path.write_text(
    json.dumps(audit, indent=2, default=str),
    encoding="utf-8",
)

md_path = OUT / "P4_SPECIFICATION_AUDIT.md"

lines = [
    "# P4 Specification Audit",
    "",
    f"Checks passed: {passed}",
    f"Checks blocked: {blocked}",
    "",
    "| Requirement | Status | Evidence |",
    "|---|---|---|",
]

for item in checks:
    detail = item["detail"].replace("|", "\\|")
    lines.append(
        f"| {item['requirement']} | {item['status']} | {detail} |"
    )

lines.extend([
    "",
    "## Interpretation",
    "",
    "A passing software implementation check is not equivalent to scientific validation.",
    "Pose-based contact entropy is not automatically trajectory ensemble entropy.",
    "A helper function is not an operational resource allocator unless production",
    "code calls it and records the resulting decisions.",
    "A risk-coverage table is not sufficient evidence of calibration without",
    "appropriate held-out evaluation and uncertainty analysis.",
    "A database is not immutable merely because it stores provenance.",
    "",
    "## Files",
    "",
    f"- JSON audit: `{json_path.relative_to(ROOT)}`",
    f"- Markdown audit: `{md_path.relative_to(ROOT)}`",
    "",
])

md_path.write_text("\n".join(lines), encoding="utf-8")

print("\n========== FINAL P4 AUDIT ==========")
print("Checks passed:", passed)
print("Checks blocked:", blocked)
print("JSON:", json_path)
print("Markdown:", md_path)
print("\nThis is an evidence audit, not a P4 completion certificate.")
