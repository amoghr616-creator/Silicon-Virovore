#!/usr/bin/env python3
"""Matched P1–P4-TA comparative experiment planner, runner, and analyzer.

P1: uniform sequence-search baseline.
P2: hotspot/entropy/epistasis pair selection with random residue proposals.
P3: matched P2-style operator with learned substitution proposals; structures
    are assessed blind and cannot reorder parents.
P4-TA: same P3 sequence operator plus generation-1 target-aware parent ordering.

Target geometry is computational evidence, not binding affinity. The PEPBI
family-held-out validation summary is reported as external-transfer evidence,
not as direct HERV-K Env validation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import random
import shutil
import shlex
import subprocess
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean, median, stdev
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.config import DEFAULT_SEED_SEQUENCE  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"
DEFAULT_HISTORICAL_MANIFEST = (
    ROOT / "prototype4/results/real_structure_integration_check_20261101/"
    "historical_pose_manifest.json"
)
DEFAULT_POOL_DIR = ROOT / "results/unified_ablation_20261102/initial_pools"
DEFAULT_EXTERNAL_SUMMARY = (
    ROOT / "prototype4/results/target_context_external_validation_fresh/summary.json"
)
DEFAULT_RECEPTOR = ROOT / "prototype4/data/haddock3/receptors/9sya_source_chain_A.pdb"
DEFAULT_HADDOCK_CONFIG = ROOT / "prototype4/configs/p4_reliability_seed1917.cfg"

ARMS: dict[str, dict[str, Any]] = {
    "P1": {
        "condition": "baseline",
        "mutation_guidance_mode": "uniform",
        "adaptive_enabled": False,
        "arise_enabled": False,
        "ale_enabled": False,
        "hotspot_enabled": False,
        "structural_assessment": False,
        "target_ordering": False,
        "description": "Uniform sequence-search baseline.",
    },
    "P2": {
        "condition": "adaptive",
        "mutation_guidance_mode": "hotspot_epistasis_entropy_exploratory",
        "adaptive_enabled": True,
        "arise_enabled": True,
        "ale_enabled": True,
        "hotspot_enabled": True,
        "structural_assessment": False,
        "target_ordering": False,
        "description": "Entropy/epistasis pair selection with random residue proposals.",
    },
    "P3": {
        "condition": "adaptive",
        "mutation_guidance_mode": "hotspot_entropy_substitution_epistasis_exploratory_v2_matched",
        "adaptive_enabled": True,
        "arise_enabled": True,
        "ale_enabled": True,
        "hotspot_enabled": True,
        "structural_assessment": True,
        "target_ordering": False,
        "description": "Learned substitution proposals; generation-1 structures assessed without parent reordering.",
    },
    "P4-TA": {
        "condition": "adaptive",
        "mutation_guidance_mode": "hotspot_entropy_substitution_epistasis_exploratory_v2_matched",
        "adaptive_enabled": True,
        "arise_enabled": True,
        "ale_enabled": True,
        "hotspot_enabled": True,
        "structural_assessment": True,
        "target_ordering": True,
        "description": "P3 sequence operator plus exploratory generation-1 target-aware parent ordering.",
    },
}


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode()
    ).hexdigest()


def parse_int_list(value: str, label: str) -> list[int]:
    values = [int(piece.strip()) for piece in value.split(",") if piece.strip()]
    if not values or len(values) != len(set(values)):
        raise ValueError(f"{label} must contain unique comma-separated integers")
    if any(item < 0 or item > 0xFFFFFFFF for item in values):
        raise ValueError(f"{label} values must be between 0 and 2**32 - 1")
    return values


def read_fasta(path: Path) -> list[str]:
    sequences: list[str] = []
    current: list[str] | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        text = line.strip()
        if not text:
            continue
        if text.startswith(">"):
            if current is not None:
                sequences.append("".join(current).upper())
            current = []
        else:
            if current is None:
                raise ValueError(f"FASTA sequence precedes first header: {path}")
            current.append(text)
    if current is not None:
        sequences.append("".join(current).upper())
    return sequences


def validate_pool(path: Path, *, count: int, sequence_length: int,
                  seed_sequence: str, historical_sequence: str) -> list[str]:
    sequences = read_fasta(path)
    if len(sequences) != count or len(set(sequences)) != count:
        raise ValueError(f"{path}: expected {count} unique FASTA sequences")
    if any(len(seq) != sequence_length or set(seq) - set(AA) for seq in sequences):
        raise ValueError(f"{path}: invalid sequence length or amino-acid symbol")
    required = {seed_sequence, historical_sequence}
    if not required.issubset(set(sequences)):
        raise ValueError(f"{path}: pool must include seed and historical reference sequences")
    return sequences


def validate_historical_manifest(path: Path, seed_sequence: str) -> tuple[str, dict[str, Any]]:
    if not path.is_file():
        raise FileNotFoundError(f"Historical structural manifest is missing: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema_version") != 1:
        raise ValueError("Historical structural manifest must use schema_version=1")
    target = payload.get("target") or {}
    if (target.get("target_structure_id") != "9SYA"
            or target.get("receptor_chain") != "A"
            or target.get("peptide_chain") != "B"):
        raise ValueError("Historical manifest must identify 9SYA receptor chain A and peptide chain B")
    hotspot_set = {
        (str(h.get("chain")), int(h.get("residue_number", -1)), str(h.get("residue_name", "")).upper())
        for h in target.get("hotspots", [])
    }
    expected_hotspots = {("A", 203, "ARG"), ("A", 205, "ARG")}
    if not expected_hotspots.issubset(hotspot_set):
        raise ValueError("Historical manifest must include 9SYA chain-A Arg203 and Arg205")
    complexes = payload.get("complexes")
    if not isinstance(complexes, list) or not complexes:
        raise ValueError("Historical manifest has no complex records")
    sequences = set()
    for index, record in enumerate(complexes):
        seq = str(record.get("sequence", "")).strip().upper()
        raw_path = record.get("path")
        expected_hash = str(record.get("sha256", "")).lower()
        if not seq or len(seq) != len(seed_sequence) or set(seq) - set(AA):
            raise ValueError(f"Historical manifest complex {index} has an invalid peptide sequence")
        if not raw_path or len(expected_hash) != 64:
            raise ValueError(f"Historical manifest complex {index} lacks path or SHA-256")
        complex_path = Path(raw_path).expanduser()
        if not complex_path.is_absolute():
            complex_path = path.parent / complex_path
        complex_path = complex_path.resolve()
        if not complex_path.is_file():
            raise FileNotFoundError(f"Historical complex file is missing: {complex_path}")
        if sha256_file(complex_path) != expected_hash:
            raise ValueError(f"Historical complex SHA-256 mismatch: {complex_path}")
        sequences.add(seq)
    if len(sequences) != 1:
        raise ValueError(f"Expected one unique historical peptide sequence; found {len(sequences)}")
    historical_sequence = next(iter(sequences))
    if historical_sequence == seed_sequence:
        raise ValueError("Historical reference sequence duplicates the configured seed")
    return historical_sequence, payload


def generate_matched_pool(seed: int, count: int, seed_sequence: str,
                          historical_sequence: str) -> list[str]:
    neighbours = []
    for position, original in enumerate(seed_sequence):
        for replacement in AA:
            if replacement == original:
                continue
            variant = seed_sequence[:position] + replacement + seed_sequence[position + 1:]
            if variant != historical_sequence:
                neighbours.append(variant)
    needed = count - 2
    if needed < 0 or needed > len(neighbours):
        raise ValueError("Pool size must allow seed + historical reference + unique single mutants")
    rng = random.Random(seed)
    return [seed_sequence, historical_sequence, *rng.sample(neighbours, needed)]


def write_fasta(path: Path, seed: int, sequences: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(
        f">matched_seed_{seed}_candidate_{index:03d}\n{sequence}\n"
        for index, sequence in enumerate(sequences, 1)
    ), encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows and not fields:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = fields or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _run_id(arm: str, seed: int) -> str:
    return f"{arm.replace('-', '_')}_seed{seed}"


def _pool_hash_record(plan: dict[str, Any], seed: int) -> dict[str, Any]:
    return next(row for row in plan["pool_records"] if int(row["seed"]) == seed)


def snapshot_historical_target_manifest(
    source_manifest: Path,
    payload: dict[str, Any],
    destination: Path,
) -> dict[str, Any]:
    """Freeze the registered target manifest and every referenced complex file."""
    destination.mkdir(parents=True, exist_ok=True)
    manifest_copy = destination / "historical_pose_manifest.json"
    shutil.copy2(source_manifest, manifest_copy)
    complex_rows = []
    for index, record in enumerate(payload.get("complexes", []), start=1):
        source = Path(str(record["path"])).expanduser()
        if not source.is_absolute():
            source = source_manifest.parent / source
        source = source.resolve()
        digest = sha256_file(source)
        snapshot = destination / "complexes" / f"{digest}.pdb"
        snapshot.parent.mkdir(parents=True, exist_ok=True)
        if not snapshot.exists():
            shutil.copy2(source, snapshot)
        complex_rows.append({
            "manifest_record_index": index,
            "sequence": str(record.get("sequence", "")).strip().upper(),
            "source_path": str(source),
            "source_sha256": digest,
            "snapshot_path": str(snapshot),
            "snapshot_sha256": sha256_file(snapshot),
            "independent_run_id": record.get("independent_run_id"),
            "random_seed": record.get("random_seed"),
            "protocol_sha256": record.get("protocol_sha256"),
            "provenance": record.get("provenance"),
        })
    snapshot_record = {
        "status": "snapshotted_and_hash_verified",
        "source_manifest_path": str(source_manifest),
        "source_manifest_sha256": sha256_file(source_manifest),
        "snapshot_manifest_path": str(manifest_copy),
        "snapshot_manifest_sha256": sha256_file(manifest_copy),
        "complex_records": complex_rows,
    }
    (destination / "snapshot_manifest.json").write_text(
        json.dumps(snapshot_record, indent=2, default=str) + "\n",
        encoding="utf-8",
    )
    snapshot_record["snapshot_record_path"] = str(destination / "snapshot_manifest.json")
    snapshot_record["snapshot_record_sha256"] = sha256_file(destination / "snapshot_manifest.json")
    return snapshot_record


def snapshot_external_validation(source_summary: Path, destination: Path) -> dict[str, Any]:
    """Freeze the existing family-held-out PEPBI report, folds, and artifacts."""
    destination.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, Any] = {
        "status": "not_registered_at_plan_time",
        "source_summary_path": str(source_summary),
        "summary_path": str(destination / "NO_REGISTERED_SUMMARY.json"),
        "summary_sha256": None,
        "source_dataset_path": None,
        "source_dataset_sha256": None,
        "copied_artifacts": [],
        "missing_artifacts": [],
        "rows": None,
        "unique_peptide_sequences": None,
        "independent_binding_groups": None,
        "binding_groups": [],
        "validation": None,
        "fold_count": None,
        "folds": [],
    }
    if not source_summary.is_file():
        sentinel = destination / "NO_REGISTERED_SUMMARY.json"
        sentinel.write_text(json.dumps({
            "status": "not_registered_at_plan_time",
            "reason": "No external family-held-out summary was present when the plan was created.",
        }, indent=2) + "\n", encoding="utf-8")
        manifest["summary_sha256"] = sha256_file(sentinel)
        (destination / "manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )
        return manifest

    try:
        summary = json.loads(source_summary.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        sentinel = destination / "NO_REGISTERED_SUMMARY.json"
        sentinel.write_text(json.dumps({
            "status": "unreadable_external_validation_summary",
            "reason": f"{type(exc).__name__}: {exc}",
        }, indent=2) + "\n", encoding="utf-8")
        manifest.update({
            "status": "unreadable_external_validation_summary",
            "summary_sha256": sha256_file(sentinel),
            "error": f"{type(exc).__name__}: {exc}",
        })
        (destination / "manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )
        return manifest
    summary_destination = destination / "summary.json"
    shutil.copy2(source_summary, summary_destination)
    manifest.update({
        "status": "summary_snapshotted",
        "summary_path": str(summary_destination),
        "summary_sha256": sha256_file(summary_destination),
        "rows": summary.get("rows"),
        "unique_peptide_sequences": summary.get("unique_peptide_sequences"),
        "independent_binding_groups": summary.get("independent_binding_groups"),
        "binding_groups": summary.get("binding_groups", []),
        "validation": summary.get("validation"),
        "fold_count": len(summary.get("folds", [])),
        "folds": summary.get("folds", []),
    })

    source_value = summary.get("source")
    if source_value:
        source_dataset = Path(source_value).expanduser()
        if not source_dataset.is_absolute():
            source_dataset = ROOT / source_dataset
        if source_dataset.is_file():
            actual_hash = sha256_file(source_dataset)
            expected_hash = str(summary.get("source_sha256") or "").lower()
            if expected_hash and actual_hash != expected_hash:
                raise ValueError(
                    "External-validation source CSV hash does not match its frozen summary: "
                    f"{source_dataset}"
                )
            target = destination / "source_dataset" / source_dataset.name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source_dataset, target)
            manifest["source_dataset_path"] = str(target)
            manifest["source_dataset_sha256"] = sha256_file(target)
        else:
            manifest["missing_artifacts"].append(str(source_dataset))

    for role, filename in (summary.get("artifacts") or {}).items():
        source_artifact = source_summary.parent / str(filename)
        if source_artifact.is_file():
            target = destination / "artifacts" / source_artifact.name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source_artifact, target)
            manifest["copied_artifacts"].append({
                "role": role, "path": str(target), "sha256": sha256_file(target)
            })
        else:
            manifest["missing_artifacts"].append(str(source_artifact))

    gate_path = source_summary.parent / "external_validation_gate.json"
    if gate_path.is_file():
        target = destination / "external_validation_gate.json"
        shutil.copy2(gate_path, target)
        manifest["copied_artifacts"].append({
            "role": "external_validation_gate", "path": str(target), "sha256": sha256_file(target)
        })

    (destination / "manifest.json").write_text(
        json.dumps(manifest, indent=2, default=str) + "\n", encoding="utf-8"
    )
    return manifest


def plan_experiment(args: argparse.Namespace) -> int:
    seeds = parse_int_list(args.seeds, "--seeds")
    dock_seeds = parse_int_list(args.dock_seeds, "--dock-seeds")
    if len(dock_seeds) < 2:
        raise ValueError("P4 repeatability assessment requires at least two distinct HADDOCK seeds")
    if args.population_size < 3 or args.generations < 2:
        raise ValueError("Use population_size >= 3 and generations >= 2 for this matched optimization pilot")
    if args.shortlist_size < 1 or args.refill_cap < 0 or args.haddock_ncores < 1:
        raise ValueError("Invalid structural shortlist, refill, or HADDOCK core setting")
    if args.shortlist_size + args.refill_cap > args.population_size:
        raise ValueError("Shortlist plus refill cap cannot exceed the population")

    out_dir = args.experiment_dir.expanduser().resolve()
    if out_dir.exists() and any(out_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty experiment directory: {out_dir}")

    historical_manifest = args.historical_manifest.expanduser().resolve()
    seed_sequence = args.seed_sequence.strip().upper()
    if not seed_sequence or set(seed_sequence) - set(AA):
        raise ValueError("--seed-sequence must contain canonical amino acids only")
    historical_sequence, historical_payload = validate_historical_manifest(
        historical_manifest, seed_sequence
    )
    receptor = DEFAULT_RECEPTOR
    haddock_config = DEFAULT_HADDOCK_CONFIG
    if not receptor.is_file():
        raise FileNotFoundError(f"Registered 9SYA chain-A receptor PDB is missing: {receptor}")
    if not haddock_config.is_file():
        raise FileNotFoundError(f"Registered HADDOCK template config is missing: {haddock_config}")
    manifest_hash = sha256_file(historical_manifest)
    historical_snapshot = snapshot_historical_target_manifest(
        historical_manifest,
        historical_payload,
        out_dir / "historical_target_manifest_snapshot",
    )

    pool_root = out_dir / "initial_pools"
    pool_root.mkdir(parents=True, exist_ok=True)
    pool_records: list[dict[str, Any]] = []
    source_pool_dir = args.initial_pool_dir.expanduser().resolve() if args.initial_pool_dir else None

    for seed in seeds:
        source_pool = source_pool_dir / f"seed_{seed}.fasta" if source_pool_dir else None
        destination = pool_root / f"seed_{seed}.fasta"
        source_label = "generated_from_registered_historical_manifest"
        if source_pool and source_pool.is_file():
            sequences = validate_pool(
                source_pool, count=args.population_size,
                sequence_length=len(seed_sequence), seed_sequence=seed_sequence,
                historical_sequence=historical_sequence,
            )
            shutil.copy2(source_pool, destination)
            source_label = str(source_pool)
        else:
            sequences = generate_matched_pool(
                seed, args.population_size, seed_sequence, historical_sequence
            )
            write_fasta(destination, seed, sequences)
            validate_pool(
                destination, count=args.population_size,
                sequence_length=len(seed_sequence), seed_sequence=seed_sequence,
                historical_sequence=historical_sequence,
            )
        pool_records.append({
            "seed": seed,
            "path": str(destination),
            "sha256": sha256_file(destination),
            "canonical_sequences_sha256": stable_hash(sequences),
            "record_count": len(sequences),
            "unique_sequence_count": len(set(sequences)),
            "seed_sequence_present": seed_sequence in sequences,
            "historical_reference_sequence_present": historical_sequence in sequences,
            "source": source_label,
        })

    out_dir.mkdir(parents=True, exist_ok=True)
    receptor_sha256 = sha256_file(receptor)
    haddock_config_sha256 = sha256_file(haddock_config)
    external_summary = args.external_validation_summary.expanduser().resolve()
    external_snapshot = snapshot_external_validation(
        external_summary, out_dir / "external_validation"
    )
    p4_gate_summary_path = external_snapshot["summary_path"]
    plan_data = {
        "schema_version": 2,
        "created_utc": now_utc(),
        "design": "four-arm matched-seed cumulative ablation with generation-1 structural assessment",
        "arms": ARMS,
        "seeds": seeds,
        "seed_sequence": seed_sequence,
        "population_size": args.population_size,
        "generations": args.generations,
        "planned_sequence_evaluations_per_run": args.population_size * args.generations,
        "planned_sequence_evaluations_total": len(seeds) * len(ARMS) * args.population_size * args.generations,
        "primary_common_endpoint": "native_fitness / C sequence-level fitness; not target binding affinity",
        "primary_comparison": "mean native_fitness in final generation",
        "secondary_common_endpoints": ["best native_fitness", "unique sequences", "mean pairwise identity", "mutation counts"],
        "p4_structural_generation": 1,
        "p3_p4_structural_budget_matched": True,
        "shortlist_size": args.shortlist_size,
        "refill_cap": args.refill_cap,
        "haddock_seeds": dock_seeds,
        "haddock_ncores": args.haddock_ncores,
        "planned_haddock_execution_slots_total": sum(1 for arm in ARMS.values() if arm["structural_assessment"]) * len(seeds) * args.shortlist_size * len(dock_seeds),
        "planned_esmfold_candidate_attempts_max": sum(1 for arm in ARMS.values() if arm["structural_assessment"]) * len(seeds) * (args.shortlist_size + args.refill_cap),
        "memory_policy": "ARISE historical-memory reads disabled; each run writes to a unique run-specific memory path",
        "historical_target_manifest": {
            "path": str(historical_manifest),
            "sha256": manifest_hash,
            "snapshot": historical_snapshot,
            "target": historical_payload["target"],
            "historical_reference_sequence": historical_sequence,
            "registered_complex_records": len(historical_payload["complexes"]),
        },
        "haddock_receptor": str(receptor),
        "haddock_receptor_sha256": receptor_sha256,
        "haddock_template_config": str(haddock_config),
        "haddock_template_config_sha256": haddock_config_sha256,
        "external_family_validation_summary": {
            "path": p4_gate_summary_path,
            "source_summary_path": str(external_summary),
            "exists_at_plan_time": external_snapshot.get("status") == "summary_snapshotted",
            "sha256": external_snapshot.get("summary_sha256"),
            "snapshot_manifest_path": str(out_dir / "external_validation" / "manifest.json"),
            "rows": external_snapshot.get("rows"),
            "unique_peptide_sequences": external_snapshot.get("unique_peptide_sequences"),
            "independent_binding_groups": external_snapshot.get("independent_binding_groups"),
            "fold_count": external_snapshot.get("fold_count"),
            "binding_groups": external_snapshot.get("binding_groups"),
            "scope_note": "Frozen external family-held-out peptide-protein benchmark; not direct HERV-K Env validation.",
        },
        "p4_gate_summary_path": p4_gate_summary_path,
        "pool_records": pool_records,
        "statistical_scope": "Exploratory paired analysis with two matched seeds; no strong inferential claim.",
        "no_time_limit": True,
    }
    (out_dir / "experiment_plan.json").write_text(json.dumps(plan_data, indent=2) + "\n", encoding="utf-8")

    run_rows = []
    for seed in seeds:
        pool_record = _pool_hash_record(plan_data, seed)
        for arm in ARMS:
            cfg = dict(ARMS[arm]) | {
                "seed": seed,
                "run_id": _run_id(arm, seed),
                "pool_sha256": pool_record["sha256"],
                "population_size": args.population_size,
                "generations": args.generations,
                "shortlist_size": args.shortlist_size if ARMS[arm]["structural_assessment"] else 0,
                "refill_cap": args.refill_cap if ARMS[arm]["structural_assessment"] else 0,
                "haddock_seeds": dock_seeds if ARMS[arm]["structural_assessment"] else [],
                "target_manifest_sha256": manifest_hash if ARMS[arm]["structural_assessment"] else None,
            }
            run_rows.append({
                "seed": seed,
                "arm": arm,
                "run_id": cfg["run_id"],
                "initial_pool_sha256": pool_record["sha256"],
                "config_sha256": stable_hash(cfg),
                "status": "PLANNED",
                "exit_code": "",
                "actual_sequence_evaluations": "",
                "scored_sequence_evaluations": "",
                "score_coverage": "",
            })
            run_manifest_path = out_dir / "run_manifests" / f"{cfg['run_id']}.json"
            run_manifest_path.parent.mkdir(parents=True, exist_ok=True)
            run_manifest_path.write_text(json.dumps({
                "schema_version": 1,
                "run_id": cfg["run_id"],
                "seed": seed,
                "arm": arm,
                "status": "PLANNED",
                "config": cfg,
                "initial_pool_path": pool_record["path"],
                "initial_pool_sha256": pool_record["sha256"],
                "historical_manifest_sha256": manifest_hash if ARMS[arm]["structural_assessment"] else None,
            }, indent=2) + "\n", encoding="utf-8")

    write_csv(out_dir / "run_ledger.csv", run_rows)
    (out_dir / "launch.sh").write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\n"
        f"cd {shlex.quote(str(ROOT))}\n"
        f"{shlex.quote(sys.executable)} scripts/p1_p4_comparative_experiment.py run "
        f"--experiment-dir {shlex.quote(str(out_dir))}\n",
        encoding="utf-8",
    )
    (out_dir / "launch.sh").chmod(0o755)
    print(json.dumps({
        "experiment_plan": str(out_dir / "experiment_plan.json"),
        "initial_pools": [row["path"] for row in pool_records],
        "planned_runs": len(run_rows),
        "planned_sequence_evaluations": plan_data["planned_sequence_evaluations_total"],
        "planned_haddock_execution_slots": plan_data["planned_haddock_execution_slots_total"],
        "max_esmfold_candidate_attempts": plan_data["planned_esmfold_candidate_attempts_max"],
        "launcher": str(out_dir / "launch.sh"),
        "launch_not_started": True,
    }, indent=2))
    return 0


def load_plan(experiment_dir: Path) -> dict[str, Any]:
    plan_path = experiment_dir / "experiment_plan.json"
    if not plan_path.is_file():
        raise FileNotFoundError(f"Missing experiment plan: {plan_path}")
    return json.loads(plan_path.read_text(encoding="utf-8"))


def update_ledger(experiment_dir: Path, seed: int, arm: str, **changes: Any) -> None:
    ledger_path = experiment_dir / "run_ledger.csv"
    with ledger_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        if int(row["seed"]) == seed and row["arm"] == arm:
            for key, value in changes.items():
                row[key] = value
            break
    write_csv(ledger_path, rows)


def _build_command(plan_data: dict[str, Any], experiment_dir: Path,
                   seed: int, arm: str) -> tuple[list[str], dict[str, str], Path, Path, dict[str, Any]]:
    config = ARMS[arm]
    run_id = _run_id(arm, seed)
    pool_record = _pool_hash_record(plan_data, seed)
    pool_path = Path(pool_record["path"]).resolve()
    run_root = experiment_dir / "runs"
    run_dir = run_root / run_id
    log_path = experiment_dir / "logs" / f"{run_id}.log"
    memory_path = experiment_dir / "memory" / f"{run_id}.json"
    command = [
        sys.executable, "run_pipeline.py",
        "--condition", str(config["condition"]),
        "--population-size", str(plan_data["population_size"]),
        "--generations", str(plan_data["generations"]),
        "--seed", str(seed),
        "--run-id", run_id,
        "--output-dir", str(run_root),
        "--seed-sequence", str(plan_data["seed_sequence"]),
        "--mutation-guidance-mode", str(config["mutation_guidance_mode"]),
        "--arise-memory-path", str(memory_path),
        "--no-arise-memory-enabled", "--no-arise-memory-use-history",
        "--no-enable-vina-tier2", "--no-cache-enabled", "--no-generate-plots",
        "--preserve-best",
    ]
    for flag, enabled in (
        ("adaptive-enabled", config["adaptive_enabled"]),
        ("arise-enabled", config["arise_enabled"]),
        ("ale-enabled", config["ale_enabled"]),
        ("hotspot-enabled", config["hotspot_enabled"]),
    ):
        command.append(f"--{flag}" if enabled else f"--no-{flag}")

    env = os.environ.copy()
    for key in (
        "SILICON_VIROVORE_P4_TARGET_MANIFEST",
        "SILICON_VIROVORE_P4_GENERATE_COMPLEXES",
        "SILICON_VIROVORE_P4_STRUCTURE_LIMIT",
        "SILICON_VIROVORE_P4_MAX_DOCK_CANDIDATES",
        "SILICON_VIROVORE_P4_MAX_STRUCTURE_REFILL_CANDIDATES",
        "SILICON_VIROVORE_P4_DOCK_SEEDS",
        "SILICON_VIROVORE_P4_DOCK_NCORES",
        "SILICON_VIROVORE_P4_TARGET_PARETO",
        "SILICON_VIROVORE_P4_TARGET_OPTIMIZATION",
        "SILICON_VIROVORE_P4_STRUCTURAL_GENERATION",
        "SILICON_VIROVORE_P4_EXTERNAL_VALIDATION_SUMMARY",
    ):
        env.pop(key, None)
    env["SILICON_VIROVORE_INITIAL_POOL_FASTA"] = str(pool_path)
    env["SILICON_VIROVORE_P4_STRUCTURAL_GENERATION"] = "1"
    env["SILICON_VIROVORE_P4_EXTERNAL_VALIDATION_SUMMARY"] = str(
        plan_data["p4_gate_summary_path"]
    )

    if config["structural_assessment"]:
        manifest = plan_data["historical_target_manifest"]["path"]
        env["SILICON_VIROVORE_P4_TARGET_MANIFEST"] = manifest
        env["SILICON_VIROVORE_P4_GENERATE_COMPLEXES"] = "1"
        env["SILICON_VIROVORE_P4_STRUCTURE_LIMIT"] = str(plan_data["shortlist_size"])
        env["SILICON_VIROVORE_P4_MAX_DOCK_CANDIDATES"] = str(plan_data["shortlist_size"])
        env["SILICON_VIROVORE_P4_MAX_STRUCTURE_REFILL_CANDIDATES"] = str(plan_data["refill_cap"])
        env["SILICON_VIROVORE_P4_DOCK_SEEDS"] = ",".join(map(str, plan_data["haddock_seeds"]))
        env["SILICON_VIROVORE_P4_DOCK_NCORES"] = str(plan_data["haddock_ncores"])
        env["SILICON_VIROVORE_P4_TARGET_PARETO"] = "1" if config["target_ordering"] else "0"
    else:
        env["SILICON_VIROVORE_P4_TARGET_PARETO"] = "0"

    if config["target_ordering"]:
        command.append("--p4-exploratory-target-optimization")

    run_config = {
        "arm": arm,
        "seed": seed,
        "command": command,
        "mutation_guidance_mode": config["mutation_guidance_mode"],
        "initial_pool_sha256": pool_record["sha256"],
        "historical_manifest_sha256": plan_data["historical_target_manifest"]["sha256"] if config["structural_assessment"] else None,
        "p4_structural_generation": 1 if config["structural_assessment"] else None,
        "target_ordering": config["target_ordering"],
        "memory_read_enabled": False,
        "memory_write_path": str(memory_path),
        "external_validation_summary_path": plan_data["p4_gate_summary_path"],
        "external_validation_summary_sha256": plan_data.get("external_family_validation_summary", {}).get("sha256"),
    }
    return command, env, log_path, run_dir, run_config


def run_experiment(args: argparse.Namespace) -> int:
    experiment_dir = args.experiment_dir.expanduser().resolve()
    plan_data = load_plan(experiment_dir)
    if (experiment_dir / "LAUNCH_STARTED").exists():
        raise FileExistsError("LAUNCH_STARTED exists; refusing to rerun or overwrite this experiment")

    for record in plan_data["pool_records"]:
        pool_path = Path(record["path"])
        if not pool_path.is_file() or sha256_file(pool_path) != record["sha256"]:
            raise RuntimeError(f"Initial pool missing or changed: {pool_path}")
    manifest_path = Path(plan_data["historical_target_manifest"]["path"])
    if not manifest_path.is_file() or sha256_file(manifest_path) != plan_data["historical_target_manifest"]["sha256"]:
        raise RuntimeError("Historical structural manifest changed after planning")
    external_record = plan_data["external_family_validation_summary"]
    if external_record.get("exists_at_plan_time"):
        frozen_summary = Path(external_record["path"])
        if not frozen_summary.is_file() or sha256_file(frozen_summary) != external_record.get("sha256"):
            raise RuntimeError("Frozen external-validation summary is missing or changed")
    for label, path_value, expected_hash in (
        ("9SYA chain-A receptor", plan_data.get("haddock_receptor"), plan_data.get("haddock_receptor_sha256")),
        ("HADDOCK template config", plan_data.get("haddock_template_config"), plan_data.get("haddock_template_config_sha256")),
    ):
        source = Path(path_value) if path_value else None
        if source is None or not source.is_file():
            raise FileNotFoundError(f"P4 structural preflight failed: {label} is missing: {source}")
        if expected_hash and sha256_file(source) != expected_hash:
            raise RuntimeError(f"P4 structural preflight failed: {label} changed after planning")

    snapshot = plan_data["historical_target_manifest"].get("snapshot") or {}
    snapshot_record_path = Path(snapshot.get("snapshot_record_path", ""))
    if not snapshot_record_path.is_file() or sha256_file(snapshot_record_path) != snapshot.get("snapshot_record_sha256"):
        raise RuntimeError("Historical target-manifest snapshot is missing or changed after planning")
    manifest_copy = Path(snapshot.get("snapshot_manifest_path", ""))
    if not manifest_copy.is_file() or sha256_file(manifest_copy) != snapshot.get("snapshot_manifest_sha256"):
        raise RuntimeError("Copied historical target manifest is missing or changed after planning")
    for complex_record in snapshot.get("complex_records", []):
        source = Path(complex_record["source_path"])
        copy = Path(complex_record["snapshot_path"])
        if not source.is_file() or sha256_file(source) != complex_record["source_sha256"]:
            raise RuntimeError(f"Registered historical complex changed after planning: {source}")
        if not copy.is_file() or sha256_file(copy) != complex_record["snapshot_sha256"]:
            raise RuntimeError(f"Historical complex snapshot is missing or changed: {copy}")

    patch_requirements = (
        '"--arise-memory-path"',
        'SILICON_VIROVORE_P4_STRUCTURAL_GENERATION',
        'p4_structural_generation_active',
    )
    pipeline_text = (ROOT / "run_pipeline.py").read_text(encoding="utf-8")
    missing = [item for item in patch_requirements if item not in pipeline_text]
    if missing:
        raise RuntimeError(
            "Pipeline support patch is not installed. Missing markers: " + ", ".join(missing)
        )

    (experiment_dir / "LAUNCH_STARTED").write_text(now_utc() + "\n", encoding="utf-8")
    (experiment_dir / "logs").mkdir(parents=True, exist_ok=True)
    (experiment_dir / "memory").mkdir(parents=True, exist_ok=True)
    (experiment_dir / "runs").mkdir(parents=True, exist_ok=True)

    failed = 0
    for seed_index, seed in enumerate(plan_data["seeds"]):
        ordered_arms = list(ARMS) if seed_index % 2 == 0 else list(reversed(ARMS))
        for arm in ordered_arms:
            command, env, log_path, run_dir, run_config = _build_command(
                plan_data, experiment_dir, int(seed), arm
            )
            if run_dir.exists():
                raise FileExistsError(f"Refusing to overwrite existing run directory: {run_dir}")
            run_manifest_path = experiment_dir / "run_manifests" / f"{_run_id(arm, int(seed))}.json"
            run_manifest = json.loads(run_manifest_path.read_text(encoding="utf-8"))
            run_manifest.update({
                "status": "RUNNING",
                "started_utc": now_utc(),
                "execution": run_config,
                "log_path": str(log_path),
                "output_directory": str(run_dir),
            })
            run_manifest_path.write_text(json.dumps(run_manifest, indent=2) + "\n", encoding="utf-8")
            update_ledger(experiment_dir, int(seed), arm, status="RUNNING")
            print(f"\n===== START {arm} seed={seed} =====", flush=True)
            start = time.perf_counter()
            log_path.parent.mkdir(parents=True, exist_ok=True)
            spawn_error = None
            interrupted = False
            process = None
            with log_path.open("w", encoding="utf-8") as log_handle:
                try:
                    process = subprocess.Popen(
                        command, cwd=ROOT, env=env,
                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                        text=True, bufsize=1,
                    )
                    assert process.stdout is not None
                    for line in process.stdout:
                        sys.stdout.write(line)
                        sys.stdout.flush()
                        log_handle.write(line)
                        log_handle.flush()
                    return_code = process.wait()
                except KeyboardInterrupt:
                    interrupted = True
                    spawn_error = "Interrupted by user; current run stopped and artifacts preserved"
                    if process is not None and process.poll() is None:
                        process.terminate()
                        try:
                            process.wait(timeout=8)
                        except subprocess.TimeoutExpired:
                            process.kill()
                            process.wait()
                    return_code = 130
                    log_handle.write(f"\nRUN_INTERRUPTED: {spawn_error}\n")
                except OSError as exc:
                    return_code = 127
                    spawn_error = f"{type(exc).__name__}: {exc}"
                    log_handle.write(f"LAUNCH_ERROR: {spawn_error}\n")
            elapsed = time.perf_counter() - start
            audit_path = run_dir / "experiment_audit.json"
            candidate_path = run_dir / "candidate_history.json"
            actual_evaluations = None
            scored_evaluations = 0
            generations_completed = 0
            postprocess_error = None
            if audit_path.is_file():
                try:
                    audit = json.loads(audit_path.read_text(encoding="utf-8"))
                    actual_evaluations = int(audit.get("execution", {}).get("total_candidates", 0))
                    generations_completed = len(audit.get("optimization", {}).get("generation_summaries", []))
                except (OSError, ValueError, TypeError, AttributeError) as exc:
                    postprocess_error = f"audit_parse_error={type(exc).__name__}: {exc}"
            if candidate_path.is_file():
                try:
                    history = json.loads(candidate_path.read_text(encoding="utf-8"))
                    if not isinstance(history, list):
                        raise ValueError("candidate_history.json must contain a JSON list")
                    candidate_history_count = len(history)
                    scored_evaluations = _count_scored_native_fitness(history)
                except (OSError, ValueError, TypeError) as exc:
                    candidate_history_count = 0
                    scored_evaluations = 0
                    postprocess_error = (
                        f"{postprocess_error + '; ' if postprocess_error else ''}"
                        f"candidate_history_parse_error={type(exc).__name__}: {exc}"
                    )
            else:
                candidate_history_count = 0
                scored_evaluations = 0
            expected_evaluations = int(plan_data["population_size"]) * int(plan_data["generations"])
            complete = (
                return_code == 0
                and audit_path.is_file()
                and candidate_path.is_file()
                and generations_completed == int(plan_data["generations"])
                and candidate_history_count == expected_evaluations
                and scored_evaluations == expected_evaluations
                and actual_evaluations == candidate_history_count
            )
            status = "INTERRUPTED" if interrupted else ("COMPLETED" if complete else ("FAILED" if return_code else "INCOMPLETE"))
            run_manifest.update({
                "status": status,
                "finished_utc": now_utc(),
                "exit_code": return_code,
                "runtime_seconds": round(elapsed, 3),
                "actual_sequence_evaluations_audit": actual_evaluations,
                "candidate_history_records": candidate_history_count,
                "scored_candidate_history_records": scored_evaluations,
                "native_fitness_score_coverage": (
                    scored_evaluations / candidate_history_count
                    if candidate_history_count else 0.0
                ),
                "generations_completed": generations_completed,
                "reason": None if complete else (
                    f"{spawn_error + '; ' if spawn_error else ''}"
                    f"{postprocess_error + '; ' if postprocess_error else ''}"
                    f"exit_code={return_code}; audit={audit_path.is_file()}; "
                    f"history_records={candidate_history_count}; "
                    f"scored_records={scored_evaluations}; "
                    f"generations={generations_completed}"
                ),
            })
            run_manifest_path.write_text(json.dumps(run_manifest, indent=2) + "\n", encoding="utf-8")
            update_ledger(
                experiment_dir, int(seed), arm,
                status=status, exit_code=str(return_code),
                actual_sequence_evaluations=str(candidate_history_count),
                scored_sequence_evaluations=str(scored_evaluations),
                score_coverage=(
                    f"{scored_evaluations / candidate_history_count:.6f}"
                    if candidate_history_count else "0"
                ),
            )
            if not complete:
                failed += 1
            print(f"===== {status}: {arm} seed={seed}; evaluations={candidate_history_count}; {elapsed:.1f}s =====", flush=True)
            if interrupted:
                break
        if 'interrupted' in locals() and interrupted:
            break

    print(f"\nExperiment runs finished: incomplete/failed={failed}; planned={len(plan_data['seeds']) * len(ARMS)}")
    print(f"Analyze with: {sys.executable} scripts/p1_p4_comparative_experiment.py analyze --experiment-dir {experiment_dir}")
    return 0 if failed == 0 else 2


def _safe_float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _count_scored_native_fitness(history: list[Any]) -> int:
    """Count candidate-history rows with a finite native-fitness score."""
    return sum(
        isinstance(row, dict)
        and _safe_float(row.get("native_fitness")) is not None
        for row in history
    )


def _load_run_data(experiment_dir: Path, plan_data: dict[str, Any],
                   seed: int, arm: str) -> dict[str, Any]:
    run_id = _run_id(arm, seed)
    run_dir = experiment_dir / "runs" / run_id
    manifest_path = experiment_dir / "run_manifests" / f"{run_id}.json"
    audit_path = run_dir / "experiment_audit.json"
    history_path = run_dir / "candidate_history.json"
    record: dict[str, Any] = {
        "seed": seed,
        "arm": arm,
        "run_id": run_id,
        "run_status": "MISSING",
        "actual_sequence_evaluations": 0,
        "native_fitness_score_count": 0,
        "pool_sha256": _pool_hash_record(plan_data, seed)["sha256"],
        "generations": [],
        "p4_generation1": None,
        "structure_refill": None,
        "structure_attempts": [],
        "target_candidates": [],
    }
    if manifest_path.is_file():
        run_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        record["run_status"] = run_manifest.get("status", "PLANNED")
        record["exit_code"] = run_manifest.get("exit_code")
        record["runtime_seconds"] = run_manifest.get("runtime_seconds")
    refill_path = run_dir / "p4_structure_refill" / "generation_0001.json"
    if refill_path.is_file():
        refill = json.loads(refill_path.read_text(encoding="utf-8"))
        record["structure_refill"] = refill
        record["structure_attempts"] = [
            {
                "seed": seed, "arm": arm, "run_id": run_id,
                "queue_rank": row.get("queue_rank"),
                "role": row.get("role"),
                "candidate_id": row.get("candidate_id"),
                "sequence_sha256": row.get("sequence_sha256"),
                "attempted": row.get("attempted"),
                "status": row.get("status"),
                "structure_available": row.get("structure_available"),
                "error": row.get("error"),
                "selected": row.get("selected"),
                "selection_rank": row.get("selection_rank"),
            }
            for row in refill.get("candidates", [])
        ]
    else:
        record["structure_refill"] = None
        record["structure_attempts"] = []

    target_path = run_dir / "p4_target_evidence" / "generation_0001.json"
    if target_path.is_file():
        report = json.loads(target_path.read_text(encoding="utf-8"))
        rows = report.get("candidates", [])
        selected = [row for row in rows if (row.get("structure_shortlist") or {}).get("selected") is True]
        def exact_sequence_evidence_valid(row: dict[str, Any]) -> bool:
            generation = row.get("complex_generation") or {}
            return (
                row.get("status") == "success"
                and generation.get("sequence_verified") is True
                and generation.get("execution_manifests_verified") is True
            )
        valid = [row for row in selected if exact_sequence_evidence_valid(row)]
        eligible = [row for row in selected if row.get("target_selection_eligible") is True]
        reasons = Counter(
            str(row.get("target_selection_reason") or row.get("status") or "unspecified")
            for row in selected if row.get("target_selection_eligible") is not True
        )
        for row in selected:
            reliability = row.get("reliability") or {}
            record["target_candidates"].append({
                "seed": seed,
                "arm": arm,
                "run_id": run_id,
                "candidate_id": row.get("candidate_id"),
                "sequence": row.get("sequence"),
                "status": row.get("status"),
                "target_complex_source": row.get("target_complex_source"),
                "registered_complex_count": row.get("registered_complex_count"),
                "complex_generation_status": (row.get("complex_generation") or {}).get("status"),
                "exact_sequence_verified": (row.get("complex_generation") or {}).get("sequence_verified") is True,
                "execution_manifests_verified": (row.get("complex_generation") or {}).get("execution_manifests_verified") is True,
                "independent_run_count": (row.get("complex_generation") or {}).get("independent_run_count"),
                "verified_pose_count": (row.get("complex_generation") or {}).get("verified_pose_count"),
                "complex_generation_json": json.dumps(row.get("complex_generation"), sort_keys=True, default=str),
                "errors_json": json.dumps(row.get("errors", []), sort_keys=True, default=str),
                "target_geometry_score_mean": row.get("target_geometry_score_mean"),
                "hotspot_proximity_score": row.get("hotspot_proximity_score"),
                "arg203_mean_distance_A": (row.get("hotspot_distances_A_mean") or {}).get("A:203:ARG"),
                "arg205_mean_distance_A": (row.get("hotspot_distances_A_mean") or {}).get("A:205:ARG"),
                "mean_cross_run_jaccard_similarity": reliability.get("mean_cross_run_jaccard_similarity"),
                "reliability_proxy_for_reporting_only": row.get("reliability_proxy_for_reporting_only"),
                "target_selection_eligible": row.get("target_selection_eligible"),
                "abstention_reason": (
                    None if row.get("target_selection_eligible") is True
                    else row.get("target_selection_reason") or row.get("status")
                ),
            })
        record["p4_generation1"] = {
            "shortlist_count": len(selected),
            "valid_target_evidence_count": len(valid),
            "exact_sequence_evidence_coverage": len(valid) / len(selected) if selected else 0.0,
            "selection_eligible_count": len(eligible),
            "selection_eligibility_coverage": len(eligible) / len(selected) if selected else 0.0,
            "abstention_count": len(selected) - len(eligible),
            "failed_structural_count": len(selected) - len(valid),
            "structure_attempt_count": sum(bool(row.get("attempted")) for row in record.get("structure_attempts", [])),
            "structure_attempt_failure_count": sum(
                bool(row.get("attempted")) and row.get("structure_available") is False
                for row in record.get("structure_attempts", [])
            ),
            "refill_attempted_count": (record.get("structure_refill") or {}).get("refill_attempted_count", 0),
            "unattempted_queue_count": (record.get("structure_refill") or {}).get("unattempted_queue_count", 0),
            "target_evidence_coverage": len(valid) / len(selected) if selected else 0.0,
            "selection_eligibility_coverage": len(eligible) / len(selected) if selected else 0.0,
            "pareto_reordered_candidates": bool(report.get("pareto_reordered_candidates")),
            "selection_reason": report.get("reason"),
            "abstention_reasons": dict(reasons),
            "report_path": str(target_path),
        }
    if not audit_path.is_file() or not history_path.is_file():
        return record

    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    history = json.loads(history_path.read_text(encoding="utf-8"))
    if not isinstance(history, list):
        history = []
    record["audit"] = audit
    record["actual_sequence_evaluations"] = len(history)
    record["native_fitness_score_count"] = sum(
        _safe_float(row.get("native_fitness")) is not None for row in history
    )
    record["audit_total_candidates"] = audit.get("execution", {}).get("total_candidates")
    record["pool_sha256_audit"] = audit.get("configuration", {}).get("initial_pool_sha256")
    record["pool_hash_matches_plan"] = record["pool_sha256_audit"] == record["pool_sha256"]
    record["generations_completed"] = len(audit.get("optimization", {}).get("generation_summaries", []))
    grouped: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for item in history:
        try:
            generation = int(item.get("generation"))
        except (TypeError, ValueError):
            continue
        grouped[generation].append(item)
    summary_by_generation = {
        int(item.get("generation")): item
        for item in audit.get("optimization", {}).get("generation_summaries", [])
        if item.get("generation") is not None
    }
    for generation in sorted(grouped):
        items = grouped[generation]
        scores = [x for x in (_safe_float(row.get("native_fitness")) for row in items) if x is not None]
        summary = summary_by_generation.get(generation, {})
        diversity = summary.get("diversity") or {}
        record["generations"].append({
            "generation": generation,
            "sequence_evaluations": len(items),
            "scored_candidate_count": len(scores),
            "native_fitness_score_coverage": len(scores) / len(items) if items else 0.0,
            "mean_native_fitness": mean(scores) if scores else None,
            "best_native_fitness": max(scores) if scores else None,
            "median_native_fitness": median(scores) if scores else None,
            "unique_sequences": summary.get("unique_sequences", len({row.get("sequence") for row in items})),
            "mean_pairwise_identity": diversity.get("mean_pairwise_identity", summary.get("mean_pairwise_identity")),
            "mutation_count": summary.get("mutation_count"),
            "pairwise_operator_diagnostics": summary.get("pairwise_operator_diagnostics"),
            "hotspot_state": summary.get("hotspot_state"),
            "p4_selected_winner": summary.get("p4_selected_winner"),
            "proposal_operator_active": summary.get("proposal_operator_active"),
        })

    return record


def _exact_two_sided_signflip(values: list[float]) -> float | None:
    if not values:
        return None
    import itertools
    observed = abs(mean(values))
    extreme = 0
    total = 2 ** len(values)
    for signs in itertools.product((-1, 1), repeat=len(values)):
        permuted = abs(mean([sign * value for sign, value in zip(signs, values)]))
        if permuted >= observed - 1e-12:
            extreme += 1
    return extreme / total


def _write_figures(out_dir: Path, run_data: list[dict[str, Any]],
                   plan_data: dict[str, Any]) -> list[str]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        (out_dir / "GRAPH_GENERATION_UNAVAILABLE.txt").write_text(
            f"Matplotlib could not be imported: {type(exc).__name__}: {exc}\n",
            encoding="utf-8",
        )
        return []

    figures_dir = out_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    generated: list[str] = []

    def save(fig, stem: str) -> None:
        for extension in ("png", "pdf", "svg"):
            destination = figures_dir / f"{stem}.{extension}"
            if extension == "png":
                fig.savefig(destination, dpi=300, bbox_inches="tight")
            else:
                fig.savefig(destination, bbox_inches="tight")
            generated.append(str(destination))
        plt.close(fig)

    arms = list(ARMS)
    colors = {"P1": "#555555", "P2": "#2878B5", "P3": "#2A9D8F", "P4-TA": "#E07A25"}

    # Figure 1: common sequence-level fitness trajectory.
    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    for arm in arms:
        by_gen: dict[int, list[float]] = defaultdict(list)
        for run in run_data:
            if run["arm"] != arm or run.get("run_status") != "COMPLETED":
                continue
            for item in run.get("generations", []):
                value = _safe_float(item.get("mean_native_fitness"))
                if value is not None:
                    by_gen[int(item["generation"])].append(value)
        if by_gen:
            xs = sorted(by_gen)
            ys = [mean(by_gen[x]) for x in xs]
            ax.plot(xs, ys, marker="o", linewidth=2.2, markersize=6, label=arm, color=colors[arm])
            if all(len(by_gen[x]) > 1 for x in xs):
                errors = [stdev(by_gen[x]) for x in xs]
                ax.errorbar(xs, ys, yerr=errors, fmt="none", capsize=4, color=colors[arm], alpha=0.8)
    ax.set_title("Matched P1–P4 sequence-level fitness trajectory")
    ax.set_xlabel("Generation")
    ax.set_ylabel("Mean native C fitness (common endpoint)")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, ncol=4)
    ax.text(0.01, -0.22, "Sequence-level computational endpoint; not binding affinity. Error bars show between-seed SD where available.", transform=ax.transAxes, fontsize=9)
    save(fig, "01_common_endpoint_trajectory")

    # Figure 2: population diversity trajectory.
    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    for arm in arms:
        by_gen = defaultdict(list)
        for run in run_data:
            if run["arm"] != arm or run.get("run_status") != "COMPLETED":
                continue
            for item in run.get("generations", []):
                val = _safe_float(item.get("unique_sequences"))
                if val is not None:
                    by_gen[int(item["generation"])].append(val)
        if by_gen:
            xs = sorted(by_gen)
            ys = [mean(by_gen[x]) for x in xs]
            ax.plot(xs, ys, marker="o", linewidth=2.2, markersize=6, label=arm, color=colors[arm])
    ax.set_title("Population diversity across optimization generations")
    ax.set_xlabel("Generation")
    ax.set_ylabel("Unique sequences per generation")
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, ncol=4)
    ax.text(0.01, -0.22, "Means include completed matched-seed runs only; missing runs are not imputed.", transform=ax.transAxes, fontsize=9)
    save(fig, "02_population_diversity")

    # Figure 3: matched P3/P4-TA generation-1 structural evidence and abstentions.
    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    evidence_by_arm: dict[str, list[float]] = {"P3": [], "P4-TA": []}
    p4_eligibility: list[float] = []
    p4_abstentions = 0
    for run in run_data:
        if run["arm"] not in evidence_by_arm:
            continue
        report = run.get("p4_generation1")
        if not report:
            continue
        evidence_by_arm[run["arm"]].append(float(report["target_evidence_coverage"]))
        if run["arm"] == "P4-TA":
            p4_eligibility.append(float(report["selection_eligibility_coverage"]))
            p4_abstentions += int(report["abstention_count"])

    labels = ["P3 exact-sequence\ngeometry coverage", "P4-TA exact-sequence\ngeometry coverage", "P4-TA diagnostic-present\neligibility"]
    values = [
        mean(evidence_by_arm["P3"]) if evidence_by_arm["P3"] else math.nan,
        mean(evidence_by_arm["P4-TA"]) if evidence_by_arm["P4-TA"] else math.nan,
        mean(p4_eligibility) if p4_eligibility else math.nan,
    ]
    bars = ax.bar(range(3), values, width=0.62)
    bars[0].set_color(colors["P3"])
    bars[1].set_color(colors["P4-TA"])
    bars[2].set_color("#6A8EAE")
    sample_counts = [
        len(evidence_by_arm["P3"]),
        len(evidence_by_arm["P4-TA"]),
        len(p4_eligibility),
    ]
    for index, value in enumerate(values):
        if math.isfinite(value):
            ax.text(index, value + 0.025, f"{value:.0%} (n={sample_counts[index]})",
                    ha="center", va="bottom", fontsize=9)
        else:
            ax.text(index, 0.5, "Not available", ha="center", va="center",
                    rotation=90, fontsize=9)
    ax.set_xticks(range(3), labels)
    ax.set_ylim(0, 1.15)
    ax.set_ylabel("Fraction of generation-1 structural shortlist")
    ax.set_title("P3/P4-TA exact-sequence evidence and diagnostic availability")
    ax.grid(axis="y", alpha=0.25)
    ax.text(
        0.01, -0.25,
        f"P4-TA candidate-level abstentions in available generation-1 reports: {p4_abstentions}. "
        "P3 is a blinded geometry-assessment control; it does not use target eligibility for parent ordering. "
        "P4 eligibility requires valid geometry and a repeatability diagnostic to be present; no minimum Jaccard threshold is applied. "
        "The diagnostic is not calibrated reliability.",
        transform=ax.transAxes, fontsize=9,
    )
    save(fig, "03_target_evidence_and_abstention")

    # Figure 4: frozen external family-held-out P3/P4 transfer result.
    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    ext = Path(plan_data["external_family_validation_summary"]["path"])
    if ext.is_file():
        try:
            external = json.loads(ext.read_text(encoding="utf-8"))
            effect = external.get("paired_primary_effect") or {}
            value = _safe_float(effect.get("macro_mae_benefit_kcal_mol"))
            ci = effect.get("group_bootstrap_95ci") or []
            low = _safe_float(ci[0]) if len(ci) >= 2 else None
            high = _safe_float(ci[1]) if len(ci) >= 2 else None
            if value is not None and low is not None and high is not None:
                ax.errorbar([value], [0], xerr=[[value-low], [high-value]], fmt="o", capsize=8, markersize=8, color=colors["P4-TA"])
                ax.axvline(0, color="#333333", linestyle="--", linewidth=1.2)
                ax.set_yticks([0], ["P3 macro-MAE − P4-TA macro-MAE"])
                ax.set_xlabel("Paired benefit (kcal/mol); positive favors P4-TA")
                ax.set_title("External family-held-out target-context validation")
                ax.text(0.01, -0.22, f"Independent binding groups: {effect.get('n_independent_groups', external.get('independent_binding_groups', 'n/a'))}; one-sided sign-flip p={effect.get('exact_one_sided_signflip_p', 'n/a')}. External benchmark only; not HERV-K Env validation.", transform=ax.transAxes, fontsize=9)
            else:
                ax.text(0.5, 0.5, "External summary lacks the predeclared paired effect or CI.", ha="center", va="center", transform=ax.transAxes)
                ax.set_title("External family-held-out validation — unavailable")
        except Exception as exc:
            ax.text(0.5, 0.5, f"Could not read external summary: {exc}", ha="center", va="center", transform=ax.transAxes)
            ax.set_title("External family-held-out validation — unavailable")
    else:
        ax.text(0.5, 0.5, "No registered external family-held-out summary found.\nNo validation values were imputed.", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("External family-held-out validation — not supplied")
    ax.grid(axis="x", alpha=0.25)
    save(fig, "04_external_family_validation")

    return generated


def analyze_experiment(args: argparse.Namespace) -> int:
    experiment_dir = args.experiment_dir.expanduser().resolve()
    plan_data = load_plan(experiment_dir)
    analysis_dir = experiment_dir / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    run_data = [
        _load_run_data(experiment_dir, plan_data, int(seed), arm)
        for seed in plan_data["seeds"] for arm in ARMS
    ]

    generation_rows = []
    run_rows = []
    for run in run_data:
        run_rows.append({
            "seed": run["seed"], "arm": run["arm"], "run_id": run["run_id"],
            "run_status": run["run_status"],
            "actual_sequence_evaluations": run["actual_sequence_evaluations"],
            "native_fitness_score_count": run.get("native_fitness_score_count", 0),
            "scored_candidate_count": sum(
                item.get("scored_candidate_count", 0) for item in run.get("generations", [])
            ),
            "audit_total_candidates": run.get("audit_total_candidates"),
            "generations_completed": run.get("generations_completed", 0),
            "initial_pool_sha256_plan": run["pool_sha256"],
            "initial_pool_sha256_audit": run.get("pool_sha256_audit"),
            "structure_queue_count_g1": len(run.get("structure_attempts", [])),
            "structure_attempt_count_g1": sum(bool(row.get("attempted")) for row in run.get("structure_attempts", [])),
            "structure_attempt_failures_g1": sum(
                bool(row.get("attempted")) and row.get("structure_available") is False
                for row in run.get("structure_attempts", [])
            ),
            "pool_hash_matches_plan": run.get("pool_hash_matches_plan"),
            "runtime_seconds": run.get("runtime_seconds"),
            "target_evidence_coverage_g1": (run.get("p4_generation1") or {}).get("target_evidence_coverage"),
            "target_eligibility_coverage_g1": (run.get("p4_generation1") or {}).get("selection_eligibility_coverage"),
            "target_abstentions_g1": (run.get("p4_generation1") or {}).get("abstention_count"),
            "pareto_reordered_g1": (run.get("p4_generation1") or {}).get("pareto_reordered_candidates"),
        })
        for item in run.get("generations", []):
            generation_rows.append({
                "seed": run["seed"], "arm": run["arm"], "run_id": run["run_id"],
                **{key: item.get(key) for key in (
                    "generation", "sequence_evaluations", "scored_candidate_count",
                    "native_fitness_score_coverage", "mean_native_fitness",
                    "best_native_fitness", "median_native_fitness", "unique_sequences",
                    "mean_pairwise_identity", "mutation_count", "proposal_operator_active",
                )},
                "pairwise_operator_diagnostics_json": json.dumps(
                    item.get("pairwise_operator_diagnostics"), sort_keys=True, default=str
                ),
                "arise_ale_guidance_json": json.dumps(
                    item.get("arise_ale_guidance"), sort_keys=True, default=str
                ),
                "common_sequence_endpoint_json": json.dumps(
                    item.get("common_sequence_endpoint"), sort_keys=True, default=str
                ),
            })

    write_csv(analysis_dir / "run_metrics.csv", run_rows)
    write_csv(analysis_dir / "generation_metrics.csv", generation_rows)

    structure_attempt_rows = [
        row for run in run_data for row in run.get("structure_attempts", [])
    ]
    target_candidate_rows = [
        row for run in run_data for row in run.get("target_candidates", [])
    ]
    write_csv(analysis_dir / "structural_attempts.csv", structure_attempt_rows)
    write_csv(analysis_dir / "target_candidate_evidence.csv", target_candidate_rows)

    comparisons = []
    contrasts = [
        ("P2-P1", "P1", "P2", "primary_incremental"),
        ("P3-P2", "P2", "P3", "primary_incremental"),
        ("P4-TA-P3", "P3", "P4-TA", "primary_incremental"),
        ("P3-P1", "P1", "P3", "secondary_cumulative"),
        ("P4-TA-P1", "P1", "P4-TA", "secondary_cumulative"),
        ("P4-TA-P2", "P2", "P4-TA", "secondary_cumulative"),
    ]
    by_key = {(run["seed"], run["arm"]): run for run in run_data}
    for name, control, treatment, comparison_role in contrasts:
        per_seed = []
        for seed in plan_data["seeds"]:
            left = by_key[(int(seed), control)]
            right = by_key[(int(seed), treatment)]
            left_gen = next((x for x in left.get("generations", []) if x["generation"] == plan_data["generations"]), None)
            right_gen = next((x for x in right.get("generations", []) if x["generation"] == plan_data["generations"]), None)
            valid = (
                left.get("run_status") == "COMPLETED"
                and right.get("run_status") == "COMPLETED"
                and left.get("pool_hash_matches_plan") is True
                and right.get("pool_hash_matches_plan") is True
                and left_gen is not None and right_gen is not None
                and left_gen.get("mean_native_fitness") is not None
                and right_gen.get("mean_native_fitness") is not None
                and left_gen.get("best_native_fitness") is not None
                and right_gen.get("best_native_fitness") is not None
            )
            if not valid:
                per_seed.append({"seed": int(seed), "status": "INCOMPLETE_PAIR", "mean_native_fitness_difference": None, "best_native_fitness_difference": None})
                continue
            per_seed.append({
                "seed": int(seed),
                "status": "PAIRED",
                "mean_native_fitness_difference": right_gen["mean_native_fitness"] - left_gen["mean_native_fitness"],
                "best_native_fitness_difference": right_gen["best_native_fitness"] - left_gen["best_native_fitness"],
            })
        valid_mean = [float(x["mean_native_fitness_difference"]) for x in per_seed if x["mean_native_fitness_difference"] is not None]
        valid_best = [float(x["best_native_fitness_difference"]) for x in per_seed if x["best_native_fitness_difference"] is not None]
        comparisons.append({
            "contrast": name,
            "comparison_role": comparison_role,
            "control": control,
            "treatment": treatment,
            "primary_endpoint": "final-generation mean native_fitness",
            "n_valid_seed_pairs": len(valid_mean),
            "per_seed": per_seed,
            "mean_paired_difference": mean(valid_mean) if valid_mean else None,
            "median_paired_difference": median(valid_mean) if valid_mean else None,
            "exact_two_sided_signflip_p_exploratory": _exact_two_sided_signflip(valid_mean),
            "mean_best_score_difference_secondary": mean(valid_best) if valid_best else None,
            "interpretation": "Descriptive paired pilot; two seeds provide limited inferential power.",
        })
    paired_report = {
        "schema_version": 1,
        "design": plan_data["design"],
        "seeds": plan_data["seeds"],
        "arms": list(ARMS),
        "primary_endpoint": "mean native C fitness in final generation",
        "planned_evaluations_per_run": plan_data["planned_sequence_evaluations_per_run"],
        "no_imputation": True,
        "comparisons": comparisons,
        "external_validation_scope": plan_data["external_family_validation_summary"],
    }
    (analysis_dir / "paired_comparisons.json").write_text(json.dumps(paired_report, indent=2, default=str) + "\n", encoding="utf-8")
    comparison_rows = []
    for contrast in comparisons:
        for row in contrast["per_seed"]:
            comparison_rows.append({"contrast": contrast["contrast"], **row})
    write_csv(analysis_dir / "paired_comparisons.csv", comparison_rows)

    abstention_rows = []
    for run in run_data:
        evidence = run.get("p4_generation1")
        if not evidence:
            continue
        for reason, count in evidence["abstention_reasons"].items():
            abstention_rows.append({"seed": run["seed"], "arm": run["arm"], "reason": reason, "count": count})
    write_csv(analysis_dir / "target_abstention_reasons.csv", abstention_rows)

    figures = _write_figures(analysis_dir, run_data, plan_data)
    summary = {
        "completed_runs": sum(row["run_status"] == "COMPLETED" for row in run_data),
        "incomplete_or_failed_runs": sum(row["run_status"] != "COMPLETED" for row in run_data),
        "planned_runs": len(run_data),
        "actual_sequence_evaluations": sum(row["actual_sequence_evaluations"] for row in run_data),
        "planned_sequence_evaluations": plan_data["planned_sequence_evaluations_total"],
        "structural_attempts": len(structure_attempt_rows),
        "structural_attempt_failures": sum(
            bool(row.get("attempted")) and row.get("structure_available") is False
            for row in structure_attempt_rows
        ),
        "external_validation_snapshot": plan_data["external_family_validation_summary"],
        "pool_hash_mismatches": [row["run_id"] for row in run_data if row.get("pool_hash_matches_plan") is False],
        "figures": figures,
        "paired_comparisons": str(analysis_dir / "paired_comparisons.json"),
        "interpretation": "Exploratory computational optimization comparison. Structural scores are not experimental affinity; PEPBI is an external benchmark, not direct HERV-K Env evidence.",
    }
    (analysis_dir / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0 if summary["incomplete_or_failed_runs"] == 0 and not summary["pool_hash_mismatches"] else 2


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)

    plan = sub.add_parser("plan", help="Create pools, manifests, and an unlaunched run plan")
    plan.add_argument("--experiment-dir", type=Path, required=True)
    plan.add_argument("--seeds", default="631,632")
    plan.add_argument("--population-size", type=int, default=100)
    plan.add_argument("--generations", type=int, default=2)
    plan.add_argument("--shortlist-size", type=int, default=2)
    plan.add_argument("--refill-cap", type=int, default=2)
    plan.add_argument("--dock-seeds", default="1917,2917")
    plan.add_argument("--haddock-ncores", type=int, default=4)
    plan.add_argument("--seed-sequence", default=DEFAULT_SEED_SEQUENCE)
    plan.add_argument("--historical-manifest", type=Path, default=DEFAULT_HISTORICAL_MANIFEST)
    plan.add_argument("--initial-pool-dir", type=Path, default=DEFAULT_POOL_DIR)
    plan.add_argument("--external-validation-summary", type=Path, default=DEFAULT_EXTERNAL_SUMMARY)
    plan.set_defaults(func=plan_experiment)

    run = sub.add_parser("run", help="Execute the registered 8-run matched experiment")
    run.add_argument("--experiment-dir", type=Path, required=True)
    run.set_defaults(func=run_experiment)

    analyze = sub.add_parser("analyze", help="Analyze completed/incomplete runs and generate four graphs")
    analyze.add_argument("--experiment-dir", type=Path, required=True)
    analyze.set_defaults(func=analyze_experiment)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return int(args.func(args))
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
