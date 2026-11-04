from __future__ import annotations

import ast
import csv
import hashlib
import importlib.util
import json
import shutil
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path.cwd().resolve()
RESULTS = ROOT / "prototype4/results"
OUT = RESULTS / "p4_final_reconciliation_20261101"
OUT.mkdir(parents=True, exist_ok=True)

def read_csv(path, delimiter=","):
    if not path.is_file():
        return []
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))

def read_json(path):
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        print("JSON ERROR:", path, exc)
        return None

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def write_csv(path, rows):
    if not rows:
        path.write_text("No records found.\n", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=list(rows[0].keys()), extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)

def num(value):
    try:
        x = float(value)
        return x if x == x and abs(x) != float("inf") else None
    except (TypeError, ValueError):
        return None

checks = []
warnings = []

def check(name, passed, detail):
    checks.append({
        "requirement": name,
        "status": "PASS" if passed else "BLOCKED",
        "detail": str(detail),
    })
    print(
        f"{'PASS' if passed else 'BLOCKED':8s} "
        f"{name}: {detail}"
    )

print("\n========== 1. FAMILY MAP RECONCILIATION ==========")

canonical_path = RESULTS / "skempi_audit_20261029/family_map.csv"
recovered_path = RESULTS / "data_isolation_20261101/family_map_recovered.csv"
cycles_path = RESULTS / "skempi_audit_featurefix_20261030/eligible_cycles.csv"
fold_path = RESULTS / "model_a_validation_20261030/model_a_outputs/outer_fold_manifest.csv"
cath_path = RESULTS / "data_isolation_20261101/cath/partner_interface_cath_assignments.tsv"

canonical = read_csv(canonical_path)
recovered = read_csv(recovered_path)
cycles = read_csv(cycles_path)
folds = read_csv(fold_path)
cath = read_csv(cath_path, "\t")

canonical_map = {
    row["complex_id"]: row["family_id"]
    for row in canonical
    if row.get("complex_id") and row.get("family_id")
}

recovered_map = {
    row["complex_id"]: row["family_id"]
    for row in recovered
    if row.get("complex_id") and row.get("family_id")
}

canonical_families = set(canonical_map.values())
recovered_families = set(recovered_map.values())
lofo_groups = {
    row.get("heldout_group")
    for row in folds
    if row.get("heldout_group")
}

mapping_mismatches = [
    {
        "complex_id": complex_id,
        "canonical_family": canonical_map.get(complex_id),
        "recovered_family": recovered_map.get(complex_id),
    }
    for complex_id in sorted(set(canonical_map) | set(recovered_map))
    if canonical_map.get(complex_id) != recovered_map.get(complex_id)
]

cycle_family_counts = Counter()
unmapped_cycles = []

for row in cycles:
    complex_id = row.get("complex_id", "")
    family = canonical_map.get(complex_id)
    if family:
        cycle_family_counts[family] += 1
    else:
        unmapped_cycles.append(row.get("cycle_uid", complex_id))

family_cycle_rows = [
    {
        "family_id": family,
        "complex_count": sum(
            1 for cid, fid in canonical_map.items() if fid == family
        ),
        "eligible_cycle_count": cycle_count,
    }
    for family, cycle_count in sorted(cycle_family_counts.items())
]

write_csv(OUT / "family_cycle_reconciliation.csv", family_cycle_rows)
write_csv(OUT / "family_mapping_mismatches.csv", mapping_mismatches)

print("Canonical complexes:", len(canonical_map))
print("Recovered complexes:", len(recovered_map))
print("Canonical family IDs:", len(canonical_families))
print("Recovered family IDs:", len(recovered_families))
print("LOFO held-out groups:", len(lofo_groups))
print("Eligible SKEMPI cycles:", len(cycles))
print("Cycles mapped to a family:", sum(cycle_family_counts.values()))
print("Unmapped cycles:", len(unmapped_cycles))
print("Canonical/recovered mapping mismatches:", len(mapping_mismatches))
print("CATH partner assignment rows:", len(cath))

check(
    "Populated canonical family map",
    len(canonical_map) == 35 and len(canonical_families) == 19,
    f"{len(canonical_map)} complexes / {len(canonical_families)} families",
)
check(
    "Recovered map agrees with canonical map",
    not mapping_mismatches and canonical_map == recovered_map,
    f"{len(mapping_mismatches)} mismatches",
)
check(
    "LOFO group count agrees with family count",
    len(lofo_groups) == len(canonical_families) == 19,
    f"{len(lofo_groups)} LOFO groups / {len(canonical_families)} families",
)
check(
    "SKEMPI cycles have family assignments",
    len(cycles) > 0 and not unmapped_cycles,
    f"{len(cycle_family_counts)} families represented; {len(unmapped_cycles)} unmapped cycles",
)

cath_status = Counter(
    row.get("CATH_assignment_status", "missing")
    for row in cath
)
print("CATH assignment status:", dict(cath_status))


print("\n========== 2. LOFO / RISK-COVERAGE ==========")

model_dir = RESULTS / "model_a_validation_20261030/model_a_outputs"
predictions = read_csv(model_dir / "loco_oof_predictions.csv")
risk = read_csv(model_dir / "risk_coverage.csv")
validation = read_json(model_dir / "p4_validation_metrics.json")

risk_comparison = []
for row in risk:
    selected = num(row.get("uncertainty_selected_mae"))
    random_mean = num(row.get("random_retention_mae_mean"))
    random_low = num(row.get("random_retention_mae_ci95_low"))

    risk_comparison.append({
        "coverage": row.get("coverage"),
        "uncertainty_selected_mae": selected,
        "random_retention_mae_mean": random_mean,
        "random_retention_ci95_low": random_low,
        "beats_random_mean": (
            selected < random_mean
            if selected is not None and random_mean is not None
            else ""
        ),
        "beats_random_ci_lower": (
            selected < random_low
            if selected is not None and random_low is not None
            else ""
        ),
    })

write_csv(OUT / "risk_coverage_reconciled.csv", risk_comparison)

valid_risk = [
    x for x in risk_comparison
    if isinstance(x["beats_random_mean"], bool)
]
risk_wins = sum(x["beats_random_mean"] for x in valid_risk)

print("LOFO prediction rows:", len(predictions))
print("LOFO groups:", len(lofo_groups))
print("Risk coverage levels:", len(risk_comparison))
print("Levels beating random mean:", f"{risk_wins}/{len(valid_risk)}")
if isinstance(validation, dict):
    print("Model minus baseline MAE:",
          validation.get("model_minus_baseline_mae"))
    print("Risk-coverage summary:",
          json.dumps(validation.get("risk_coverage", {}), indent=2))

check(
    "LOFO prediction artifact",
    len(predictions) == 242,
    f"{len(predictions)} out-of-fold predictions",
)
check(
    "Risk-coverage comparison artifact",
    len(risk_comparison) == 10,
    f"{len(risk_comparison)} coverage levels",
)


print("\n========== 3. MD TRAJECTORY AND TOPOLOGY INVENTORY ==========")

md_dir = ROOT / "results/md"

trajectory_suffixes = {".trr", ".xtc", ".dcd", ".nc"}
topology_suffixes = {
    ".tpr", ".gro", ".pdb", ".prmtop", ".psf", ".top"
}

trajectories = []
topologies = []

if md_dir.is_dir():
    for path in md_dir.rglob("*"):
        if not path.is_file():
            continue

        suffix = path.suffix.lower()

        if suffix in trajectory_suffixes:
            trajectories.append(path.resolve())
        elif suffix in topology_suffixes:
            topologies.append(path.resolve())

mdtraj_available = importlib.util.find_spec("mdtraj") is not None
mdanalysis_available = importlib.util.find_spec("MDAnalysis") is not None
gmx = shutil.which("gmx")

md_rows = []
for path in sorted(trajectories + topologies):
    md_rows.append({
        "relative_path": str(path.relative_to(ROOT)),
        "extension": path.suffix.lower(),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
        "role": (
            "trajectory"
            if path.suffix.lower() in trajectory_suffixes
            else "topology_candidate"
        ),
    })

write_csv(OUT / "md_input_inventory.csv", md_rows)

print("Trajectory files:", len(trajectories))
for path in sorted(trajectories):
    print("  ", path.relative_to(ROOT))

print("Topology candidates:", len(topologies))
for path in sorted(topologies):
    print("  ", path.relative_to(ROOT))

print("MDTraj installed:", mdtraj_available)
print("MDAnalysis installed:", mdanalysis_available)
print("GROMACS executable:", gmx or "not found")

check(
    "Trajectory files identified",
    len(trajectories) > 0,
    f"{len(trajectories)} trajectory files",
)
check(
    "Trajectory processing prerequisites",
    bool(topologies) and (mdtraj_available or mdanalysis_available),
    f"{len(topologies)} topology candidates; MDTraj={mdtraj_available}; MDAnalysis={mdanalysis_available}",
)


print("\n========== 4. ENTROPY-ADDED RMSE BENCHMARK ==========")

entropy_result_files = []
entropy_metrics = []

for path in (ROOT / "prototype4").rglob("*.json"):
    if path.stat().st_size > 10_000_000:
        continue

    obj = read_json(path)
    if not isinstance(obj, dict):
        continue

    keys = set(obj.keys())

    delta_keys = [
        key for key in keys
        if "delta_rmse" in key.lower()
        or "entropy_augmented_oof_rmse" in key.lower()
    ]

    if delta_keys:
        record = {
            "path": str(path.relative_to(ROOT)),
            "sha256": sha256(path),
            "matching_keys": ",".join(delta_keys),
            "values": json.dumps(
                {key: obj[key] for key in delta_keys},
                default=str,
            ),
        }
        entropy_result_files.append(record)
        entropy_metrics.append(obj)

write_csv(OUT / "entropy_rmse_result_inventory.csv", entropy_result_files)

print("JSON artifacts containing entropy/RMSE result keys:",
      len(entropy_result_files))
for item in entropy_result_files:
    print(item["path"], item["matching_keys"], item["values"][:500])

check(
    "Entropy-added RMSE result artifact",
    bool(entropy_result_files),
    f"{len(entropy_result_files)} result artifacts found; interpretation requires validation-split audit",
)


print("\n========== 5. EVSI / PRUNING CALL GRAPH ==========")

closed_loop_path = ROOT / "prototype4/p4_closed_loop.py"
production_path = ROOT / "prototype4/p4.py"

evsi_info = []
prune_info = []
entropy_function_info = []

for path in (closed_loop_path, production_path):
    if not path.is_file():
        continue

    source = path.read_text(encoding="utf-8", errors="replace")
    tree = ast.parse(source)

    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            name = node.name.lower()

            if "entropy" in name or "rmse" in name:
                entropy_function_info.append({
                    "file": str(path.relative_to(ROOT)),
                    "function": node.name,
                    "line": node.lineno,
                })

            if name in {
                "estimate_evsi",
                "candidate_prune_decision",
            }:
                info = {
                    "file": str(path.relative_to(ROOT)),
                    "function": node.name,
                    "line": node.lineno,
                }
                if name == "estimate_evsi":
                    evsi_info.append(info)
                else:
                    prune_info.append(info)

        elif isinstance(node, ast.Call):
            call_name = None

            if isinstance(node.func, ast.Name):
                call_name = node.func.id
            elif isinstance(node.func, ast.Attribute):
                call_name = node.func.attr

            if call_name in {
                "estimate_evsi",
                "candidate_prune_decision",
            }:
                entry = {
                    "file": str(path.relative_to(ROOT)),
                    "called_function": call_name,
                    "line": node.lineno,
                }
                if call_name == "estimate_evsi":
                    evsi_info.append(entry)
                else:
                    prune_info.append(entry)

print("EVSI definitions/calls:", json.dumps(evsi_info, indent=2))
print("Pruning definitions/calls:", json.dumps(prune_info, indent=2))
print("Entropy/RMSE functions:", json.dumps(entropy_function_info, indent=2))

check(
    "EVSI implementation",
    any(x.get("function") == "estimate_evsi" for x in evsi_info),
    f"{len(evsi_info)} EVSI references found",
)
check(
    "Pruning callsite",
    any(x.get("called_function") == "candidate_prune_decision" for x in prune_info),
    f"{sum(x.get('called_function') == 'candidate_prune_decision' for x in prune_info)} pruning calls found",
)


print("\n========== 6. FINAL REPORT ==========")

audit = {
    "generated_utc": datetime.now(timezone.utc).isoformat(),
    "family_mapping": {
        "canonical_complexes": len(canonical_map),
        "canonical_families": len(canonical_families),
        "recovered_complexes": len(recovered_map),
        "recovered_families": len(recovered_families),
        "mapping_mismatches": mapping_mismatches,
        "lofo_groups": len(lofo_groups),
        "eligible_cycles": len(cycles),
        "mapped_cycles": sum(cycle_family_counts.values()),
        "unmapped_cycles": unmapped_cycles,
        "cath_assignment_status": dict(cath_status),
    },
    "risk_coverage": {
        "levels": len(risk_comparison),
        "beating_random_mean": risk_wins,
        "validation_metrics": validation,
    },
    "md_inputs": {
        "trajectory_files": [
            str(p.relative_to(ROOT)) for p in trajectories
        ],
        "topology_candidates": [
            str(p.relative_to(ROOT)) for p in topologies
        ],
        "mdtraj_available": mdtraj_available,
        "mdanalysis_available": mdanalysis_available,
        "gromacs_executable": gmx,
        "trajectory_entropy_calculated": False,
    },
    "entropy_rmse_results": entropy_result_files,
    "evsi_references": evsi_info,
    "pruning_references": prune_info,
    "entropy_functions": entropy_function_info,
    "checks": checks,
    "warnings": warnings,
}

json_out = OUT / "P4_FINAL_RECONCILIATION.json"
json_out.write_text(json.dumps(audit, indent=2, default=str))

md = [
    "# P4 Final Evidence Reconciliation",
    "",
    f"Generated: {audit['generated_utc']}",
    "",
    "## Family map",
    "",
    f"- Canonical complex records: {len(canonical_map)}",
    f"- Populated canonical families: {len(canonical_families)}",
    f"- Recovered map families: {len(recovered_families)}",
    f"- LOFO held-out groups: {len(lofo_groups)}",
    f"- Eligible cycles: {len(cycles)}",
    f"- Cycles mapped to families: {sum(cycle_family_counts.values())}",
    f"- Unmapped cycles: {len(unmapped_cycles)}",
    "",
    "## Risk coverage",
    "",
    f"- Coverage levels: {len(risk_comparison)}",
    f"- Levels where uncertainty selection beats random mean: {risk_wins}",
    "",
    "## Trajectory inventory",
    "",
    f"- Trajectory files: {len(trajectories)}",
    f"- Topology candidates: {len(topologies)}",
    f"- MDTraj: {mdtraj_available}",
    f"- MDAnalysis: {mdanalysis_available}",
    f"- GROMACS: {gmx or 'not found'}",
    "",
    "Trajectory entropy has not been calculated by this reconciliation.",
    "The existence of trajectory files does not prove they are valid",
    "protein-complex trajectories with the required chain mapping.",
    "",
    "## EVSI and pruning",
    "",
    f"- EVSI references: {len(evsi_info)}",
    f"- Pruning references: {len(prune_info)}",
    "",
    "## Entropy-added RMSE artifacts",
    "",
    f"- Result artifacts found: {len(entropy_result_files)}",
    "",
    "## Gate table",
    "",
    "| Check | Status | Evidence |",
    "|---|---|---|",
]

for item in checks:
    detail = item["detail"].replace("|", "\\|")
    md.append(
        f"| {item['requirement']} | {item['status']} | {detail} |"
    )

md.extend([
    "",
    "## Scientific conclusion",
    "",
    "The populated family map resolves the earlier empty-template audit error.",
    "The 410-feature, family-aware LOFO validation and risk-coverage results",
    "remain separate from direct HERV-K Env candidate-binding validation.",
    "No positive entropy-added RMSE benefit or fully operational EVSI-pruning",
    "loop is claimed unless supported by the corresponding held-out evidence.",
    "",
    "No candidate binding affinity, inhibition, or therapeutic efficacy",
    "claim follows from these computational analyses.",
    "",
])

md_out = OUT / "P4_FINAL_RECONCILIATION.md"
md_out.write_text("\n".join(md) + "\n", encoding="utf-8")

print("\n========== FINAL RESULT ==========")
print("PASS:", sum(x["status"] == "PASS" for x in checks))
print("BLOCKED:", sum(x["status"] == "BLOCKED" for x in checks))
print("JSON:", json_out)
print("REPORT:", md_out)
print("CSV directory:", OUT)
print("Original experimental data were not modified.")
