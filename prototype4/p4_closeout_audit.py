from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
P4 = ROOT / "prototype4"
RESULTS = P4 / "results"
OUT = RESULTS / "p4_closeout_20261101"
OUT.mkdir(parents=True, exist_ok=True)


def read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        return {"_error": f"{type(exc).__name__}: {exc}"}


def read_csv(path: Path) -> list[dict]:
    try:
        with path.open(
            "r", encoding="utf-8-sig", newline=""
        ) as handle:
            return list(csv.DictReader(handle))
    except Exception:
        return []


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_source(raw_path: str, manifest_path: Path) -> Path | None:
    if not raw_path:
        return None

    supplied = Path(raw_path).expanduser()

    candidates = (
        [supplied]
        if supplied.is_absolute()
        else [
            ROOT / supplied,
            manifest_path.parent / supplied,
            P4 / supplied,
            RESULTS / supplied,
        ]
    )

    for candidate in candidates:
        try:
            resolved = candidate.resolve()
            if resolved.is_file():
                return resolved
        except OSError:
            continue

    return None


def normalize_complexes(value: Any) -> list[dict]:
    if isinstance(value, list):
        return [x for x in value if isinstance(x, dict)]
    if isinstance(value, dict):
        return [x for x in value.values() if isinstance(x, dict)]
    return []


def find_latest_target_report() -> Path | None:
    candidates = []

    for path in (ROOT / "results").rglob("generation_0001.json"):
        if "p4_target_evidence" in path.parts:
            candidates.append(path)

    if not candidates:
        return None

    return max(candidates, key=lambda x: x.stat().st_mtime)


def short_metrics(obj: Any, keys: tuple[str, ...]) -> dict:
    if not isinstance(obj, dict):
        return {}

    result = {}
    for key in keys:
        if key in obj:
            result[key] = obj[key]
    return result


# ---------------------------------------------------------
# 1. Historical structural pose manifest and file integrity
# ---------------------------------------------------------

manifest_path = (
    RESULTS
    / "real_structure_integration_check_20261101"
    / "historical_pose_manifest.json"
)

manifest = read_json(manifest_path) if manifest_path.exists() else {}
complexes = normalize_complexes(manifest.get("complexes", []))

manifest_rows = []
independent_records = []
independent_ids = set()
independent_seeds = set()
protocol_hashes = set()

for index, item in enumerate(complexes, start=1):
    raw_path = str(item.get("path", ""))
    resolved = resolve_source(raw_path, manifest_path)

    expected_hash = str(item.get("sha256", "")).lower()
    actual_hash = sha256(resolved) if resolved else ""

    hash_valid = bool(
        expected_hash
        and actual_hash
        and expected_hash == actual_hash
    )

    run_id = item.get("independent_run_id")
    seed = item.get("random_seed")
    protocol = item.get("protocol_sha256")

    complete_independent_metadata = all(
        value is not None and str(value).strip()
        for value in (run_id, seed, protocol)
    )

    if complete_independent_metadata:
        independent_records.append(item)
        independent_ids.add(str(run_id))
        independent_seeds.add(str(seed))
        protocol_hashes.add(str(protocol))

    manifest_rows.append({
        "manifest_index": index,
        "relative_or_supplied_path": raw_path,
        "resolved_path": (
            str(resolved.relative_to(ROOT))
            if resolved and resolved.is_relative_to(ROOT)
            else str(resolved) if resolved else ""
        ),
        "file_exists": resolved is not None,
        "expected_sha256": expected_hash,
        "actual_sha256": actual_hash,
        "sha256_valid": hash_valid,
        "sequence": item.get("sequence", ""),
        "independent_run_id": run_id,
        "random_seed": seed,
        "protocol_sha256": protocol,
        "independent_metadata_complete": complete_independent_metadata,
    })

manifest_csv = OUT / "pose_manifest_integrity.csv"

with manifest_csv.open("w", newline="", encoding="utf-8") as handle:
    fields = list(manifest_rows[0].keys()) if manifest_rows else [
        "manifest_index", "relative_or_supplied_path",
        "file_exists", "sha256_valid",
    ]
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    writer.writerows(manifest_rows)

existing_pose_files = sum(x["file_exists"] for x in manifest_rows)
valid_pose_hashes = sum(x["sha256_valid"] for x in manifest_rows)

# ---------------------------------------------------------
# 2. HADDOCK pose and contact-entropy evidence
# ---------------------------------------------------------

pose_dir = (
    RESULTS
    / "haddock_pose_audit_20261031"
)

pose_comparison_path = (
    pose_dir / "p4_target_scored_poses"
    / "ten_pose_target_comparison.csv"
)

contact_entropy_path = (
    pose_dir / "p4_target_scored_poses"
    / "contact_entropy_summary.csv"
)

pose_comparison = read_csv(pose_comparison_path)
contact_entropy = read_csv(contact_entropy_path)

protomer_counts: dict[str, int] = {}

for row in pose_comparison:
    protomer = (
        row.get("source_protomer")
        or row.get("chain")
        or row.get("source_chain")
        or "unknown"
    )
    protomer_counts[protomer] = protomer_counts.get(protomer, 0) + 1

# ---------------------------------------------------------
# 3. HADDOCK run reproducibility evidence
# ---------------------------------------------------------

repro_path = (
    RESULTS
    / "haddock_run_reproducibility_20261101"
    / "reproducibility_summary.json"
)

repro = read_json(repro_path) if repro_path.exists() else {}

repro_metadata = repro.get("run_metadata", {})
repro_comparisons = repro.get("comparisons", {})
repro_best_matches = repro.get("best_cross_run_match_summary", {})

# ---------------------------------------------------------
# 4. SKEMPI independent benchmark and LOFO outputs
# ---------------------------------------------------------

skempi_dir = RESULTS / "skempi_audit_featurefix_20261030"

cycles_path = skempi_dir / "eligible_cycles.csv"
family_map_path = skempi_dir / "family_map_template.csv"

cycles = read_csv(cycles_path)
family_map = read_csv(family_map_path)

family_ids = {
    row.get("family_id")
    for row in family_map
    if row.get("family_id")
}

model_dir = RESULTS / "model_a_validation_20261030"

loco_dir = model_dir / "model_a_outputs"
loco_predictions_path = loco_dir / "loco_oof_predictions.csv"
outer_metrics_path = loco_dir / "outer_group_metrics.csv"
validation_metrics_path = loco_dir / "p4_validation_metrics.json"

loco_predictions = read_csv(loco_predictions_path)
outer_metrics = read_csv(outer_metrics_path)
validation_metrics = (
    read_json(validation_metrics_path)
    if validation_metrics_path.exists()
    else {}
)

directional_path = (
    model_dir
    / "p4_directional_analysis_20261101"
    / "p4_directional_statistics.json"
)

directional = (
    read_json(directional_path)
    if directional_path.exists()
    else {}
)

ridge_family_path = (
    model_dir
    / "p4_nested_directional_classifier_20261101"
    / "ridge_family_advantage_20261101"
    / "ridge_family_advantage_summary.json"
)

ridge_family = (
    read_json(ridge_family_path)
    if ridge_family_path.exists()
    else {}
)

# ---------------------------------------------------------
# 5. Inspect latest generation-level target-evidence report
# ---------------------------------------------------------

latest_target_path = find_latest_target_report()

latest_target = (
    read_json(latest_target_path)
    if latest_target_path
    else {}
)

target_summary = short_metrics(
    latest_target,
    (
        "generation",
        "candidate_count",
        "target_evidence_count",
        "target_evidence_coverage",
        "reliability_diagnostic_candidate_count",
        "reliability_feedback_active",
        "reliability_calibrated",
        "pareto_requested",
        "pareto_reordered_candidates",
        "reason",
    ),
)

# ---------------------------------------------------------
# 6. Evidence gates
# ---------------------------------------------------------

validation_artifacts_present = all([
    bool(cycles),
    bool(loco_predictions),
    bool(outer_metrics),
    bool(validation_metrics),
    bool(directional),
    bool(ridge_family),
])

structural_pose_artifacts_present = all([
    manifest_path.exists(),
    bool(pose_comparison),
    bool(contact_entropy),
    bool(repro),
])

independent_runs_verified = (
    len(independent_ids) >= 2
    and len(independent_seeds) >= 2
    and len(protocol_hashes) == 1
    and len(independent_records) >= 2
)

reliability_calibrated = (
    latest_target.get("reliability_calibrated") is True
)

reliability_feedback_active = (
    latest_target.get("reliability_feedback_active") is True
)

checks = {
    "manifest_exists": manifest_path.exists(),
    "all_manifest_pose_files_exist": (
        bool(manifest_rows)
        and existing_pose_files == len(manifest_rows)
    ),
    "all_manifest_sha256_valid": (
        bool(manifest_rows)
        and valid_pose_hashes == len(manifest_rows)
    ),
    "at_least_two_independent_runs": len(independent_ids) >= 2,
    "at_least_two_independent_seeds": len(independent_seeds) >= 2,
    "single_consistent_protocol_hash": len(protocol_hashes) == 1,
    "pose_comparison_available": bool(pose_comparison),
    "contact_entropy_available": bool(contact_entropy),
    "run_reproducibility_report_available": bool(repro),
    "skempi_eligible_cycles_available": bool(cycles),
    "family_map_available": bool(family_map),
    "loco_predictions_available": bool(loco_predictions),
    "outer_group_metrics_available": bool(outer_metrics),
    "directional_analysis_available": bool(directional),
    "family_directional_analysis_available": bool(ridge_family),
    "reliability_calibrated": reliability_calibrated,
    "reliability_feedback_active": reliability_feedback_active,
}

# Scientific completion intentionally requires evidence, not code presence.
scientific_closeout_ready = all([
    validation_artifacts_present,
    structural_pose_artifacts_present,
    independent_runs_verified,
    reliability_calibrated,
])

# ---------------------------------------------------------
# 7. Build a reproducibility inventory of key source files
# ---------------------------------------------------------

source_paths = [
    manifest_path,
    pose_comparison_path,
    contact_entropy_path,
    repro_path,
    cycles_path,
    family_map_path,
    loco_predictions_path,
    outer_metrics_path,
    validation_metrics_path,
    directional_path,
    ridge_family_path,
]

artifact_rows = []

for path in source_paths:
    if path.is_file():
        artifact_rows.append({
            "relative_path": str(path.relative_to(ROOT)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
        })

artifact_csv = OUT / "p4_source_artifact_manifest.csv"

with artifact_csv.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["relative_path", "size_bytes", "sha256"],
    )
    writer.writeheader()
    writer.writerows(artifact_rows)

# ---------------------------------------------------------
# 8. Machine-readable audit
# ---------------------------------------------------------

audit = {
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "project_root": str(ROOT),
    "manifest": {
        "path": str(manifest_path.relative_to(ROOT)),
        "complex_count": len(complexes),
        "pose_files_found": existing_pose_files,
        "valid_sha256_count": valid_pose_hashes,
        "independent_metadata_records": len(independent_records),
        "unique_independent_run_ids": len(independent_ids),
        "unique_random_seeds": len(independent_seeds),
        "unique_protocol_hashes": len(protocol_hashes),
    },
    "pose_audit": {
        "target_comparison_rows": len(pose_comparison),
        "poses_by_protomer": protomer_counts,
        "contact_entropy_rows": len(contact_entropy),
    },
    "reproducibility": {
        "report_present": bool(repro),
        "run_metadata": repro_metadata,
        "comparisons": repro_comparisons,
        "best_cross_run_matches": repro_best_matches,
    },
    "benchmark": {
        "skempi_eligible_cycles": len(cycles),
        "family_map_complexes": len(family_map),
        "unique_family_ids": len(family_ids),
        "loco_predictions": len(loco_predictions),
        "outer_group_metric_rows": len(outer_metrics),
        "validation_metrics": validation_metrics,
        "directional_statistics": directional,
        "ridge_family_advantage": ridge_family,
    },
    "latest_target_integration_report": {
        "path": (
            str(latest_target_path.relative_to(ROOT))
            if latest_target_path
            else None
        ),
        "summary": target_summary,
    },
    "checks": checks,
    "validation_artifacts_present": validation_artifacts_present,
    "structural_pose_artifacts_present": structural_pose_artifacts_present,
    "independent_runs_verified": independent_runs_verified,
    "scientific_closeout_ready": scientific_closeout_ready,
    "source_artifacts": artifact_rows,
}

audit_json_path = OUT / "p4_closeout_audit.json"
audit_json_path.write_text(
    json.dumps(audit, indent=2, default=str),
    encoding="utf-8",
)

# ---------------------------------------------------------
# 9. Human-readable closeout report
# ---------------------------------------------------------

def status(value: bool) -> str:
    return "PASS" if value else "BLOCKED"


lines = [
    "# P4 Evidence and Reliability Closeout Audit",
    "",
    f"Audit timestamp: {audit['created_utc']}",
    "",
    "## 1. Structural manifest integrity",
    "",
    f"- Manifest complexes: {len(complexes)}",
    f"- Pose files found: {existing_pose_files}/{len(complexes)}",
    f"- Valid SHA-256 checksums: {valid_pose_hashes}/{len(complexes)}",
    f"- Complete independent-run metadata records: {len(independent_records)}",
    f"- Unique independent run IDs: {len(independent_ids)}",
    f"- Unique random seeds: {len(independent_seeds)}",
    f"- Unique protocol hashes: {len(protocol_hashes)}",
    "",
    "## 2. Structural pose audit",
    "",
    f"- Target-scored pose rows: {len(pose_comparison)}",
    f"- Poses by source protomer: `{json.dumps(protomer_counts)}`",
    f"- Contact-entropy summary rows: {len(contact_entropy)}",
    f"- Reproducibility report present: {bool(repro)}",
    "",
    "## 3. Independent benchmark",
    "",
    f"- SKEMPI eligible cycles: {len(cycles)}",
    f"- Family-map complexes: {len(family_map)}",
    f"- Unique family identifiers: {len(family_ids)}",
    f"- LOFO prediction rows: {len(loco_predictions)}",
    f"- Outer-group metric rows: {len(outer_metrics)}",
    "",
    "### Validation metric summary",
    "",
    "```json",
    json.dumps(
        short_metrics(
            validation_metrics,
            (
                "n_cycles",
                "n_independent_groups",
                "split_unit",
                "model_minus_baseline_mae",
                "validation_strength",
                "interpretation",
            ),
        ),
        indent=2,
    ),
    "```",
    "",
    "### Directional analysis summary",
    "",
    "```json",
    json.dumps(
        short_metrics(
            ridge_family,
            (
                "n_informative_groups",
                "n_informative_cycles",
                "ridge_pooled_accuracy",
                "majority_pooled_accuracy",
                "ridge_family_macro_accuracy",
                "majority_family_macro_accuracy",
                "macro_difference_ridge_minus_majority",
                "family_block_permutation",
                "family_cluster_bootstrap",
                "scope_warning",
            ),
        ),
        indent=2,
    ),
    "```",
    "",
    "## 4. Latest target-specific integration report",
    "",
    f"Report: `{latest_target_path}`" if latest_target_path else "No generation report found.",
    "",
    "```json",
    json.dumps(target_summary, indent=2, default=str),
    "```",
    "",
    "## 5. Completion gates",
    "",
]

for name, value in checks.items():
    lines.append(f"- **{status(value)}** `{name}`")

lines.extend([
    "",
    "## 6. Scientific interpretation",
    "",
    "The SKEMPI family-aware benchmark and structural-pose diagnostics are "
    "separate evidence streams. The SKEMPI benchmark evaluates generalization "
    "on its own experimental dataset; it is not direct experimental validation "
    "of HERV-K Env candidate binding.",
    "",
    "Historical poses can support geometric and interface-consistency analysis. "
    "They cannot establish independent-run reliability unless run identifiers, "
    "random seeds, and protocol hashes demonstrate genuinely separate runs.",
    "",
    "A reliability score must not be described as calibrated unless the "
    "calibration procedure and independent evidence support that claim.",
    "",
    "## 7. Closeout decision",
    "",
    f"- Validation artifacts present: **{validation_artifacts_present}**",
    f"- Structural pose artifacts present: **{structural_pose_artifacts_present}**",
    f"- Independent runs verified: **{independent_runs_verified}**",
    f"- Reliability calibrated: **{reliability_calibrated}**",
    f"- Scientific closeout ready: **{scientific_closeout_ready}**",
    "",
    "Do not bypass the evidence gates or enable feedback/pruning solely to "
    "make the pipeline appear complete. Missing independent structural "
    "replicates require new, properly documented computational runs; missing "
    "experimental validation cannot be created by post-processing.",
    "",
    "## 8. Audit artifacts",
    "",
    f"- Pose manifest audit: `{manifest_csv.relative_to(ROOT)}`",
    f"- Source file checksums: `{artifact_csv.relative_to(ROOT)}`",
    f"- Machine-readable audit: `{audit_json_path.relative_to(ROOT)}`",
    "",
])

report_path = OUT / "P4_CLOSEOUT_REPORT.md"
report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

# ---------------------------------------------------------
# 10. Terminal summary
# ---------------------------------------------------------

print("\n========== P4 CLOSEOUT AUDIT ==========")
print("Historical complexes:", len(complexes))
print("Pose files found:", existing_pose_files)
print("Valid checksums:", valid_pose_hashes)
print("Complete independent metadata records:", len(independent_records))
print("Unique independent run IDs:", len(independent_ids))
print("Unique random seeds:", len(independent_seeds))
print("Unique protocol hashes:", len(protocol_hashes))
print("HADDOCK pose rows:", len(pose_comparison))
print("Contact entropy rows:", len(contact_entropy))
print("SKEMPI cycles:", len(cycles))
print("Family groups:", len(family_ids))
print("LOFO prediction rows:", len(loco_predictions))
print("Outer group metrics:", len(outer_metrics))

print("\n========== GATES ==========")
for name, value in checks.items():
    print(f"{status(value):8s} {name}")

print("\nScientific closeout ready:", scientific_closeout_ready)
print("\nCREATED:")
print(report_path)
print(artifact_csv)
print(manifest_csv)
print(audit_json_path)

if not scientific_closeout_ready:
    print(
        "\nP4 is not scientifically closed. "
        "Use the report to resolve evidence gaps; do not bypass gates."
    )
