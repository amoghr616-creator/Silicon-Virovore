#!/usr/bin/env python3
"""Plan and analyze a matched P3 versus exploratory P4-TA experiment.

The generated launcher uses the same initial FASTA, seed, population size,
generation count, mutation mode, and structural-evaluation budget for both
arms. P3 receives target structures as a blinded evaluation only; P4-TA may
use the same type of target geometry for parent ordering after strict
sequence-match and repeatability checks. The experiment remains exploratory:
its geometry endpoint is not experimental binding affinity.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import shlex
import subprocess
import sys
from pathlib import Path
from statistics import mean
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.config import DEFAULT_SEED_SEQUENCE

AA = "ACDEFGHIKLMNPQRSTVWY"
DEFAULT_MODE = "hotspot_entropy_substitution_epistasis_exploratory_v2_matched"


def parse_seeds(value: str) -> list[int]:
    value = value.strip()
    if "-" in value and "," not in value:
        left, right = value.split("-", 1)
        seeds = list(range(int(left), int(right) + 1))
    else:
        seeds = [int(x.strip()) for x in value.split(",") if x.strip()]
    if not seeds or len(set(seeds)) != len(seeds):
        raise ValueError("Provide one or more unique integer seeds.")
    if any(seed < 0 or seed > 0xFFFFFFFF for seed in seeds):
        raise ValueError("Seeds must be between 0 and 2**32 - 1.")
    return seeds


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_pool(seed: int, pool_size: int, seed_sequence: str) -> list[str]:
    if pool_size < 1:
        raise ValueError("pool_size must be positive")
    sequence = seed_sequence.upper()
    if not sequence or set(sequence) - set(AA):
        raise ValueError("Seed sequence must contain canonical amino acids only")
    neighbors = []
    for i, old in enumerate(sequence):
        for replacement in AA:
            if replacement != old:
                neighbors.append(sequence[:i] + replacement + sequence[i + 1:])
    if pool_size > len(neighbors) + 1:
        raise ValueError(f"pool_size exceeds available unique single-mutant pool ({len(neighbors)+1})")
    rng = random.Random(seed)
    if pool_size == 1:
        return [sequence]
    selected = rng.sample(neighbors, pool_size - 1)
    return [sequence, *selected]


def write_fasta(path: Path, seed: int, sequences: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = "".join(
        f">matched_seed_{seed}_candidate_{index:03d}\n{sequence}\n"
        for index, sequence in enumerate(sequences, start=1)
    )
    path.write_text(text, encoding="utf-8")


def verify_fasta(path: Path, *, expected_count: int, expected_length: int) -> list[str]:
    sequences = []
    current = None
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current is not None:
                sequences.append("".join(current))
            current = []
        else:
            if current is None:
                raise ValueError(f"FASTA sequence before header: {path}")
            current.append(line.upper())
    if current is not None:
        sequences.append("".join(current))
    if len(sequences) != expected_count or len(set(sequences)) != expected_count:
        raise ValueError(f"Invalid FASTA count/uniqueness in {path}")
    if any(len(s) != expected_length or set(s) - set(AA) for s in sequences):
        raise ValueError(f"Invalid sequence length/residue in {path}")
    return sequences


def render_launcher(
    *, out_dir: Path, seeds: list[int], population_size: int,
    generations: int, shortlist_size: int, dock_seeds: tuple[int, int],
    mutation_mode: str, ncores: int,
) -> str:
    repo = str(ROOT)
    base = str(out_dir.resolve())
    pool_root = str((out_dir / "initial_pools").resolve())
    lines = [
        "#!/usr/bin/env bash",
        "set -Eeuo pipefail",
        f"cd {shlex.quote(repo)}",
        f"BASE={shlex.quote(base)}",
        f"POOL_ROOT={shlex.quote(pool_root)}",
        f"SHORTLIST={int(shortlist_size)}",
        f"DOCK_SEEDS={shlex.quote(','.join(map(str, dock_seeds)))}",
        'export SILICON_VIROVORE_P4_DOCK_SEEDS="$DOCK_SEEDS"',
        "mkdir -p \"$BASE/logs\" \"$BASE/runs\"",
        "test ! -e \"$BASE/LAUNCH_STARTED\" || { echo 'Refusing to relaunch an existing experiment directory.' >&2; exit 2; }",
        "touch \"$BASE/LAUNCH_STARTED\"",
        "unset SILICON_VIROVORE_P4_TARGET_MANIFEST",
        "export SILICON_VIROVORE_INITIAL_POOL_FASTA=''",
        "export SILICON_VIROVORE_P4_GENERATE_COMPLEXES=1",
        "export SILICON_VIROVORE_P4_STRUCTURE_LIMIT=\"$SHORTLIST\"",
        "export SILICON_VIROVORE_P4_MAX_DOCK_CANDIDATES=\"$SHORTLIST\"",
        f"export SILICON_VIROVORE_P4_DOCK_NCORES={int(ncores)}",
        "export SILICON_VIROVORE_ARISE_MEMORY_ENABLED=0",
        "export SILICON_VIROVORE_P4_TARGET_PARETO=0",
        "",
    ]
    for seed in seeds:
        pool = str((out_dir / "initial_pools" / f"seed_{seed}.fasta").resolve())
        for arm in ("P3", "P4-TA"):
            run_id = f"{arm.replace('-', '_')}_seed{seed}"
            log = str((out_dir / "logs" / f"{run_id}.log").resolve())
            run_dir = str((out_dir / "runs" / run_id).resolve())
            lines += [
                f"echo 'START {run_id} seed={seed} arm={arm}'",
                f"test ! -e {shlex.quote(run_dir)} || {{ echo 'Run directory exists; refusing overwrite: {run_dir}' >&2; exit 3; }}",
                f"export SILICON_VIROVORE_INITIAL_POOL_FASTA={shlex.quote(pool)}",
                "export SILICON_VIROVORE_P4_TARGET_PARETO=1" if arm == "P4-TA" else "export SILICON_VIROVORE_P4_TARGET_PARETO=0",
            ]
            command = [
                sys.executable, "run_pipeline.py",
                "--condition", "adaptive",
                "--population-size", str(population_size),
                "--generations", str(generations),
                "--seed", str(seed),
                "--run-id", run_id,
                "--output-dir", str((out_dir / "runs").resolve()),
                "--adaptive-enabled", "--arise-enabled", "--ale-enabled", "--hotspot-enabled",
                "--mutation-guidance-mode", mutation_mode,
                "--no-arise-memory-enabled", "--no-arise-memory-use-history",
                "--no-enable-vina-tier2", "--no-cache-enabled", "--no-generate-plots",
            ]
            if arm == "P4-TA":
                command.append("--p4-exploratory-target-optimization")
            lines += [
                "set +e",
                " ".join(shlex.quote(arg) for arg in command) + f" 2>&1 | tee {shlex.quote(log)}",
                "rc=${PIPESTATUS[0]}",
                "set -e",
                f"if [[ $rc -ne 0 ]]; then printf '%s\\t%s\\t%s\\tINCOMPLETE(exit=%s)\\n' '{seed}' '{arm}' '{run_id}' \"$rc\" >> \"$BASE/run_ledger.tsv\"; exit $rc; fi",
                f"python scripts/p4_target_coverage_audit.py --run-dir {shlex.quote(run_dir)} || {{ echo 'coverage audit failed for {run_id}' >&2; exit 4; }}",
                f'echo -e "{seed}\t{arm}\t{run_id}\tCOMPLETED" >> "$BASE/run_ledger.tsv"',
                f"echo 'COMPLETED {run_id}'",
                "",
            ]
    lines += [
        "echo 'Both arms finished for all registered seeds.'",
        "echo 'Run: python scripts/p4_matched_experiment.py analyze --experiment-dir \"$BASE\"'",
    ]
    return "\n".join(lines) + "\n"


def plan(args: argparse.Namespace) -> int:
    seeds = parse_seeds(args.seeds)
    out_dir = args.experiment_dir.expanduser().resolve()
    if out_dir.exists() and any(out_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty experiment directory: {out_dir}")
    pool_root = out_dir / "initial_pools"
    pool_root.mkdir(parents=True, exist_ok=True)
    pool_records = []
    for seed in seeds:
        pool = make_pool(seed, args.population_size, args.seed_sequence)
        path = pool_root / f"seed_{seed}.fasta"
        write_fasta(path, seed, pool)
        verified = verify_fasta(path, expected_count=args.population_size, expected_length=len(args.seed_sequence))
        pool_records.append({
            "seed": seed,
            "path": str(path),
            "sha256": sha256_file(path),
            "record_count": len(verified),
            "unique_sequence_count": len(set(verified)),
            "seed_sequence_included": args.seed_sequence.upper() in verified,
        })
    manifest = {
        "schema_version": 1,
        "design": "matched P3 versus exploratory P4-TA; identical per-seed initial pools",
        "arms": ["P3", "P4-TA"],
        "seeds": seeds,
        "population_size": args.population_size,
        "generations": args.generations,
        "sequence_evaluations_per_run_expected": args.population_size * args.generations,
        "shortlist_size_per_generation": args.shortlist_size,
        "haddock_independent_seeds": list(args.dock_seeds),
        "mutation_guidance_mode": args.mutation_mode,
        "seed_sequence": args.seed_sequence.upper(),
        "pool_generation": "deterministic-seeded sample of single-substitution neighbours of configured project seed",
        "pool_records": pool_records,
        "target_context_validation": "exploratory only; PEPBI external-transfer gate not passed",
        "independent_experimental_binding_validation": False,
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "experiment_plan.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    ledger = out_dir / "run_ledger.csv"
    with ledger.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["seed", "arm", "run_id", "initial_pool_sha256", "status"])
        for record in pool_records:
            for arm in ("P3", "P4-TA"):
                writer.writerow([record["seed"], arm, f"{arm.replace('-', '_')}_seed{record['seed']}", record["sha256"], "PLANNED"])
    ledger_tsv = out_dir / "run_ledger.tsv"
    with ledger_tsv.open("w", encoding="utf-8") as handle:
        handle.write("seed\tarm\trun_id\tstatus\n")
        for seed in seeds:
            for arm in ("P3", "P4-TA"):
                handle.write(f"{seed}\t{arm}\t{arm.replace('-', '_')}_seed{seed}\tPLANNED\n")
    launcher = out_dir / "launch_matched_p3_p4ta.sh"
    launcher.write_text(render_launcher(
        out_dir=out_dir,
        seeds=seeds,
        population_size=args.population_size,
        generations=args.generations,
        shortlist_size=args.shortlist_size,
        dock_seeds=args.dock_seeds,
        mutation_mode=args.mutation_mode,
        ncores=args.ncores,
    ), encoding="utf-8")
    launcher.chmod(0o755)
    print(f"Plan: {out_dir / 'experiment_plan.json'}")
    print(f"Pools: {pool_root}")
    print(f"Launcher: {launcher}")
    print("Planned runs:", len(seeds) * 2)
    print("Target structural executions per P4 run:", args.generations * args.shortlist_size * len(args.dock_seeds))
    print("P3 receives the same target-evaluation budget as blinded evaluation; no parent ordering change.")
    return 0


def read_generation_metrics(run_dir: Path) -> list[dict[str, Any]]:
    reports = sorted((run_dir / "p4_target_evidence").glob("generation_*.json"))
    metrics = []
    for report_path in reports:
        data = json.loads(report_path.read_text(encoding="utf-8"))
        rows = data.get("candidates", [])
        shortlist = [r for r in rows if (r.get("structure_shortlist") or {}).get("selected")]
        valid = [r for r in shortlist if r.get("status") == "success"]
        repeatable = [r for r in shortlist if (r.get("reliability") or {}).get("status") == "diagnostic_only"]
        geometry = [float(r["target_geometry_score_mean"]) for r in valid if r.get("target_geometry_score_mean") is not None]
        metrics.append({
            "generation": int(data.get("generation", 0)),
            "shortlist_size": len(shortlist),
            "valid_target_evidence_count": len(valid),
            "repeatability_valid_count": len(repeatable),
            "mean_target_geometry": mean(geometry) if geometry else None,
            "best_target_geometry": max(geometry) if geometry else None,
            "pareto_reordered_candidates": bool(data.get("pareto_reordered_candidates")),
            "exploratory_override_used": bool(data.get("exploratory_override_used")),
            "selection_validation_status": data.get("selection_validation_status"),
        })
    return metrics


def _read_run_status(run_dir: Path, expected_generations: int, expected_evaluations: int) -> tuple[str, str]:
    audit_path = run_dir / "experiment_audit.json"
    if not audit_path.is_file():
        return "INCOMPLETE", "missing experiment_audit.json"
    try:
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return "INCOMPLETE", f"audit unreadable: {exc}"
    summaries = audit.get("optimization", {}).get("generation_summaries", [])
    if len(summaries) != expected_generations:
        return "INCOMPLETE", f"expected {expected_generations} generation summaries, got {len(summaries)}"
    total = int(audit.get("execution", {}).get("total_candidates", 0) or 0)
    if total != expected_evaluations:
        return "INCOMPLETE", f"expected {expected_evaluations} evaluations, got {total}"
    history_path = run_dir / "candidate_history.json"
    if not history_path.is_file():
        return "INCOMPLETE", "missing candidate_history.json"
    return "PASS", "run-level execution records complete"


def exact_signflip_p(values: list[float]) -> float | None:
    n = len(values)
    if n == 0 or n > 20:
        return None
    observed = sum(values) / n
    total = 1 << n
    extreme = 0
    for mask in range(total):
        permuted = sum((1 if mask & (1 << i) else -1) * value for i, value in enumerate(values)) / n
        if abs(permuted) >= abs(observed) - 1e-12:
            extreme += 1
    return extreme / total


def bootstrap_ci(values: list[float], *, repetitions: int = 10000, seed: int = 1701) -> list[float] | None:
    if len(values) < 2:
        return None
    rng = random.Random(seed)
    means = []
    for _ in range(repetitions):
        sample = [values[rng.randrange(len(values))] for _ in values]
        means.append(mean(sample))
    means.sort()
    return [means[int(0.025 * (len(means) - 1))], means[int(0.975 * (len(means) - 1))]]


def analyze(args: argparse.Namespace) -> int:
    exp_dir = args.experiment_dir.expanduser().resolve()
    plan_data = json.loads((exp_dir / "experiment_plan.json").read_text(encoding="utf-8"))
    seeds = [int(x) for x in plan_data["seeds"]]
    generations = int(plan_data["generations"])
    expected = int(plan_data["sequence_evaluations_per_run_expected"])
    rows = []
    generation_rows = []
    paired_effects = []
    for seed in seeds:
        arm_data = {}
        for arm in ("P3", "P4-TA"):
            run_id = f"{arm.replace('-', '_')}_seed{seed}"
            run_dir = exp_dir / "runs" / run_id
            status, reason = _read_run_status(run_dir, generations, expected)
            gm = read_generation_metrics(run_dir) if (run_dir / "p4_target_evidence").exists() else []
            for generation_metric in gm:
                generation_rows.append({"seed": seed, "arm": arm, "run_id": run_id, **generation_metric})
            final = next((m for m in reversed(gm) if m["generation"] == generations), None)
            if status == "PASS" and len(gm) != generations:
                status, reason = "INCOMPLETE", f"expected {generations} target-evidence generation reports, got {len(gm)}"
            if status == "PASS":
                for generation_metric in gm:
                    if generation_metric["valid_target_evidence_count"] != generation_metric["shortlist_size"]:
                        status, reason = "INCOMPLETE", f"generation {generation_metric['generation']} lacks full shortlist target-evidence coverage"
                        break
                    if generation_metric["repeatability_valid_count"] != generation_metric["shortlist_size"]:
                        status, reason = "INCOMPLETE", f"generation {generation_metric['generation']} lacks full shortlist repeatability coverage"
                        break
                    if arm == "P4-TA" and not generation_metric["pareto_reordered_candidates"]:
                        status, reason = "INCOMPLETE", f"P4-TA did not apply Pareto ordering in generation {generation_metric['generation']}"
                        break
                    if arm == "P3" and generation_metric["pareto_reordered_candidates"]:
                        status, reason = "INCOMPLETE", f"control P3 unexpectedly applied target-aware ordering in generation {generation_metric['generation']}"
                        break
            row = {
                "seed": seed,
                "arm": arm,
                "run_id": run_id,
                "run_status": status,
                "status_reason": reason,
                "final_generation": generations if final else None,
                "final_shortlist_size": final["shortlist_size"] if final else None,
                "final_valid_target_evidence_count": final["valid_target_evidence_count"] if final else None,
                "final_repeatability_valid_count": final["repeatability_valid_count"] if final else None,
                "final_mean_target_geometry": final["mean_target_geometry"] if final else None,
                "final_best_target_geometry": final["best_target_geometry"] if final else None,
                "p4_pareto_reordered": final["pareto_reordered_candidates"] if final else False,
            }
            rows.append(row)
            arm_data[arm] = row
        if all(arm_data[a]["run_status"] == "PASS" for a in ("P3", "P4-TA")):
            p3_value = arm_data["P3"]["final_best_target_geometry"]
            p4_value = arm_data["P4-TA"]["final_best_target_geometry"]
            if p3_value is not None and p4_value is not None:
                p3_mean = arm_data["P3"]["final_mean_target_geometry"]
                p4_mean = arm_data["P4-TA"]["final_mean_target_geometry"]
                paired_effects.append({
                    "seed": seed,
                    "P3_final_best_geometry": p3_value,
                    "P4_TA_final_best_geometry": p4_value,
                    "P4_minus_P3_best_geometry_difference": p4_value - p3_value,
                    "P3_final_mean_geometry": p3_mean,
                    "P4_TA_final_mean_geometry": p4_mean,
                    "P4_minus_P3_mean_geometry_difference": (
                        p4_mean - p3_mean if p3_mean is not None and p4_mean is not None else None
                    ),
                })

    out_csv = exp_dir / "matched_run_metrics.csv"
    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else [])
        if rows:
            writer.writeheader()
            writer.writerows(rows)
    generation_csv = exp_dir / "matched_generation_metrics.csv"
    with generation_csv.open("w", newline="", encoding="utf-8") as handle:
        fields = list(generation_rows[0]) if generation_rows else ["seed", "arm", "run_id", "generation"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(generation_rows)
    effects_csv = exp_dir / "paired_seed_effects.csv"
    with effects_csv.open("w", newline="", encoding="utf-8") as handle:
        fields = list(paired_effects[0]) if paired_effects else ["seed", "P3_final_best_geometry", "P4_TA_final_best_geometry", "P4_minus_P3_best_geometry_difference", "P3_final_mean_geometry", "P4_TA_final_mean_geometry", "P4_minus_P3_mean_geometry_difference"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(paired_effects)

    differences = [float(r["P4_minus_P3_best_geometry_difference"]) for r in paired_effects]
    summary = {
        "experiment_dir": str(exp_dir),
        "planned_seed_count": len(seeds),
        "complete_paired_seed_count": len(paired_effects),
        "planned_paired_seed_count": len(seeds),
        "p4_minus_p3_mean_final_best_geometry_difference": mean(differences) if differences else None,
        "paired_seed_bootstrap_95ci": bootstrap_ci(differences),
        "paired_seed_exact_two_sided_signflip_p": exact_signflip_p(differences),
        "paired_seed_effects": paired_effects,
        "all_runs_complete": len(paired_effects) == len(seeds),
        "endpoint": "best target-geometry score in final-generation structural shortlist",
        "mean_final_geometry_difference": mean([float(r["P4_minus_P3_mean_geometry_difference"]) for r in paired_effects if r.get("P4_minus_P3_mean_geometry_difference") is not None]) if paired_effects else None,
        "mean_geometry_bootstrap_95ci": bootstrap_ci([float(r["P4_minus_P3_mean_geometry_difference"]) for r in paired_effects if r.get("P4_minus_P3_mean_geometry_difference") is not None]),
        "endpoint_is_independent_of_selection": False,
        "experimental_binding_validation": False,
        "status": "EXPLORATORY_COMPUTATIONAL_COMPARISON_NOT_BINDING_VALIDATION",
        "limitations": [
            "Target geometry is a computational structural endpoint, not experimental affinity.",
            "P4-TA uses target geometry to alter parent ordering; the final best geometry endpoint is not independent of that selection objective.",
            "Incomplete run pairs are not imputed and are excluded from paired effects.",
            "A separate untouched experimental HERV-K Env interaction measurement would be required for biological validation.",
        ],
        "matched_run_metrics_csv": str(out_csv),
        "paired_seed_effects_csv": str(effects_csv),
        "matched_generation_metrics_csv": str(generation_csv),
    }
    summary_path = exp_dir / "matched_experiment_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    if not summary["all_runs_complete"]:
        return 2
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p_plan = sub.add_parser("plan", help="Create identical per-seed FASTA pools and a launcher")
    p_plan.add_argument("--experiment-dir", type=Path, required=True)
    p_plan.add_argument("--seeds", default="631,632")
    p_plan.add_argument("--population-size", type=int, default=100)
    p_plan.add_argument("--generations", type=int, default=3)
    p_plan.add_argument("--shortlist-size", type=int, default=4)
    p_plan.add_argument("--dock-seeds", default="1917,2917")
    p_plan.add_argument("--mutation-mode", default=DEFAULT_MODE)
    p_plan.add_argument("--ncores", type=int, default=4)
    p_plan.add_argument("--seed-sequence", default=DEFAULT_SEED_SEQUENCE)
    p_plan.set_defaults(func=plan)
    p_analyze = sub.add_parser("analyze", help="Analyze only complete matched run pairs")
    p_analyze.add_argument("--experiment-dir", type=Path, required=True)
    p_analyze.set_defaults(func=analyze)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    if hasattr(args, "dock_seeds"):
        dock_seeds = tuple(int(x.strip()) for x in args.dock_seeds.split(",") if x.strip())
        if len(dock_seeds) != 2 or len(set(dock_seeds)) != 2:
            parser.error("--dock-seeds must contain exactly two distinct seeds")
        args.dock_seeds = dock_seeds
    if hasattr(args, "population_size") and args.population_size < 2:
        parser.error("--population-size must be >= 2")
    if hasattr(args, "generations") and args.generations < 1:
        parser.error("--generations must be >= 1")
    if hasattr(args, "shortlist_size") and not 1 <= args.shortlist_size <= args.population_size:
        parser.error("--shortlist-size must be between 1 and population-size")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
