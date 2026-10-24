from __future__ import annotations

import csv
import hashlib
import json
import math
import random
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from statistics import mean, stdev

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "prototype2/p2_baseline_vs_substitution_replication_20261022"
RUNS = OUT / "runs"
LOGS = OUT / "logs"
STATE = OUT / "state"
SUMMARY = OUT / "run_summary.csv"
CONTRASTS = OUT / "paired_contrasts.csv"
MANIFEST_PATH = OUT / "experiment_manifest.json"

POPULATION = 16
GENERATIONS = 10
MUTATION_RATE = 0.05
ELITE_COUNT = 2
SEED_SEQUENCE = "MKLAVFALLVFFAGSSDLIRR"
SEEDS = [20900021 + 1000 * i for i in range(200)]
BOOTSTRAP_REPLICATES = 20_000

METHOD_FLAGS = [
    "--no-adaptive-enabled",
    "--arise-enabled",
    "--no-ale-enabled",
    "--hotspot-enabled",
]

METHODS = {
    "adaptive_substitution": {
        "mode": "hotspot_entropy_substitution_epistasis_exploratory",
        "flags": list(METHOD_FLAGS),
    },
    "uniform_baseline": {
        "mode": "uniform",
        "flags": list(METHOD_FLAGS),
    },
}

SOURCE_FILES = [
    "run_pipeline.py",
    "src/config.py",
    "src/population_runner.py",
    "src/epistasis_entropy.py",
    "src/hotspot.py",
    "src/ranking.py",
    "src/ale_policy.py",
    "src/mutation_credit.py",
    "src/models.py",
    "c/bridge.py",
    "c/engine.c",
    "c/ga_loop.c",
    "c/hydropathy.c",
    "c/safety.c",
    "c/virovore.h",
]
CSV_FIELDS = [
    "seed", "method", "cli_mode", "final_best", "mean_mutation_count",
    "mean_pairwise_identity", "run_id", "log_file",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def score_from_summary(line: str, key: str) -> float | None:
    match = re.search(
        rf"['\"]{re.escape(key)}['\"]\s*:\s*"
        rf"(-?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)",
        line,
    )
    return float(match.group(1)) if match else None


def percentile(sorted_values: list[float], p: float) -> float:
    if not sorted_values:
        raise ValueError("percentile requires nonempty data")
    position = (len(sorted_values) - 1) * p
    low = math.floor(position)
    high = math.ceil(position)
    if low == high:
        return sorted_values[low]
    frac = position - low
    return sorted_values[low] * (1 - frac) + sorted_values[high] * frac


def save_rows(rows: list[dict[str, object]]) -> None:
    tmp = SUMMARY.with_suffix(".csv.tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(SUMMARY)


def prepare_memory_snapshot() -> tuple[Path, bool]:
    memory = ROOT / "results/arise_memory.json"
    snapshot = STATE / "arise_memory_original.json"
    existed_marker = STATE / "arise_memory_original_existed.txt"
    if not existed_marker.exists():
        STATE.mkdir(parents=True, exist_ok=True)
        if memory.exists():
            shutil.copy2(memory, snapshot)
            existed_marker.write_text("true\n")
        else:
            existed_marker.write_text("false\n")
    existed = existed_marker.read_text().strip() == "true"
    if existed and not snapshot.exists():
        raise SystemExit("ABORT: original ARISE-memory snapshot is missing.")
    return snapshot, existed


def restore_memory(snapshot: Path, existed: bool) -> None:
    memory = ROOT / "results/arise_memory.json"
    if existed:
        memory.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(snapshot, memory)
    elif memory.exists():
        memory.unlink()


def main() -> None:
    ROOT_RESOLVED = ROOT.resolve()
    import os
    os.chdir(ROOT_RESOLVED)

    if not (ROOT / ".venv/bin/python").exists() and Path(sys.executable).name != "python":
        print(f"Warning: runner Python is {sys.executable}; activate .venv before launch.", flush=True)
    if len(SEEDS) != 200 or len(set(SEEDS)) != 200:
        raise SystemExit("ABORT: expected 200 unique matched seeds.")
    if METHODS["adaptive_substitution"]["flags"] != METHODS["uniform_baseline"]["flags"]:
        raise SystemExit("ABORT: arm-level flags differ; only mutation-guidance mode may vary.")
    if set(SEEDS) & {20700021 + 1000 * i for i in range(200)}:
        raise SystemExit("ABORT: new seed list overlaps the prior 200-seed seed range.")

    OUT.mkdir(parents=True, exist_ok=True)
    RUNS.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    STATE.mkdir(parents=True, exist_ok=True)

    missing = [name for name in SOURCE_FILES if not (ROOT / name).is_file()]
    if missing:
        raise SystemExit("ABORT: required source files missing: " + ", ".join(missing))
    source_hashes = {name: sha256_file(ROOT / name) for name in SOURCE_FILES}
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=False,
    ).stdout.strip()
    expected_manifest = {
        "experiment": "Fresh 200-pair matched-seed replication: substitution-specific adaptive vs ordinary uniform baseline",
        "created_at": datetime.now().isoformat(),
        "seed_count": 200,
        "seeds": SEEDS,
        "population_size": POPULATION,
        "generations": GENERATIONS,
        "mutation_rate": MUTATION_RATE,
        "elite_count": ELITE_COUNT,
        "seed_sequence": SEED_SEQUENCE,
        "methods": METHODS,
        "design_status": "Fresh-seed replication specified after prior exploratory results were inspected; not a preregistered confirmatory experiment.",
        "seed_rule": "20900021 + 1000*i for i=0..199; fixed before this run and not replaced based on outcomes.",
        "primary_endpoint": "Generation-10 best sequence-level score; paired adaptive-minus-uniform difference.",
        "comparability": "Both modes use the same evaluation budget, run parameters, and selection pipeline. The adaptive proposal deliberately includes 10 single mutants and 5 double mutants per 16-candidate generation; ordinary uniform proposes single mutants. Thus this is an end-to-end search-algorithm comparison at equal candidate-evaluation budget, not an equal-mutation-event ablation.",
        "fitness_objective": "sequence-level objective; no target-binding inference",
        "docking_enabled": False,
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "git_commit": commit,
        "source_sha256": source_hashes,
        "runner_sha256": sha256_file(Path(__file__)),
    }

    if MANIFEST_PATH.exists():
        old_manifest = json.loads(MANIFEST_PATH.read_text())
        for key in expected_manifest:
            if key == "created_at":
                continue
            if old_manifest.get(key) != expected_manifest.get(key):
                raise SystemExit(
                    f"ABORT: existing manifest differs at {key}. "
                    "Do not mix results from different settings/code."
                )
    else:
        MANIFEST_PATH.write_text(json.dumps(expected_manifest, indent=2) + "\n")

    if not SUMMARY.exists():
        save_rows([])
    with SUMMARY.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        row["seed"] = int(row["seed"])
        row["final_best"] = float(row["final_best"])
        for field in ("mean_mutation_count", "mean_pairwise_identity"):
            row[field] = float(row[field]) if row.get(field) not in (None, "") else ""
    completed = {(int(row["seed"]), str(row["method"])) for row in rows}
    legal = {(seed, method) for seed in SEEDS for method in METHODS}
    if not completed.issubset(legal):
        raise SystemExit("ABORT: existing CSV has unexpected seeds or methods.")
    if len(completed) != len(rows):
        raise SystemExit("ABORT: duplicate (seed, method) rows in existing CSV.")

    snapshot, memory_existed = prepare_memory_snapshot()
    # Restore before every run to prevent any run's memory writes from changing
    # the starting memory state of a later condition.
    try:
        help_result = subprocess.run(
            [sys.executable, "run_pipeline.py", "--help"], cwd=ROOT,
            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            check=False,
        )
        if help_result.returncode != 0:
            raise SystemExit("ABORT: run_pipeline.py --help failed.")
        help_text = help_result.stdout
        if "hotspot_entropy_substitution_epistasis_exploratory" not in help_text:
            raise SystemExit("ABORT: current CLI does not expose the requested adaptive mode.")
        if "--arise-memory-use-history" not in help_text or "--arise-memory-enabled" not in help_text:
            raise SystemExit("ABORT: memory isolation CLI options are missing.")

        total = len(SEEDS) * len(METHODS)
        finished = len(completed)
        for seed in SEEDS:
            # Keep each pair together for easier auditability.
            for method in ("adaptive_substitution", "uniform_baseline"):
                key = (seed, method)
                if key in completed:
                    continue
                spec = METHODS[method]
                run_id = f"p2-rep-{method}-seed{seed}"
                log_path = LOGS / f"{method}_seed{seed}.log"
                output_dir = RUNS / f"{method}_seed{seed}"
                restore_memory(snapshot, memory_existed)

                command = [
                    sys.executable, "-u", "run_pipeline.py",
                    "--population-size", str(POPULATION),
                    "--generations", str(GENERATIONS),
                    "--mutation-rate", str(MUTATION_RATE),
                    "--elite-count", str(ELITE_COUNT),
                    "--seed", str(seed),
                    "--seed-sequence", SEED_SEQUENCE,
                    "--run-id", run_id,
                    "--output-dir", str(output_dir.relative_to(ROOT)),
                    "--mutation-guidance-mode", str(spec["mode"]),
                    "--no-arise-memory-enabled",
                    "--no-arise-memory-use-history",
                    "--no-enable-vina-tier2",
                    "--no-generate-plots",
                    *list(spec["flags"]),
                ]
                print(
                    f"[{finished + 1}/{total}] START {method} seed={seed}",
                    flush=True,
                )
                with log_path.open("w") as log_handle:
                    result = subprocess.run(
                        command, cwd=ROOT, stdout=log_handle,
                        stderr=subprocess.STDOUT, text=True, timeout=600,
                        check=False,
                    )
                if result.returncode != 0:
                    tail = "\n".join(log_path.read_text(errors="replace").splitlines()[-60:])
                    print("RUN FAILED; last log lines:\n" + tail, flush=True)
                    raise SystemExit(
                        f"Stopped at {method}, seed={seed}; resume by rerunning this same script."
                    )

                lines = log_path.read_text(errors="replace").splitlines()
                matches = [
                    line for line in lines
                    if f"Generation {GENERATIONS} summary:" in line
                ]
                if not matches:
                    raise SystemExit(f"No generation-{GENERATIONS} summary in {log_path}")
                last_summary = matches[-1]
                score = score_from_summary(last_summary, "best_overall")
                if score is None or not math.isfinite(score):
                    raise SystemExit(f"Could not parse final score in {log_path}")
                row = {
                    "seed": seed,
                    "method": method,
                    "cli_mode": spec["mode"],
                    "final_best": f"{score:.12f}",
                    "mean_mutation_count": score_from_summary(last_summary, "mean_mutation_count"),
                    "mean_pairwise_identity": score_from_summary(last_summary, "mean_pairwise_identity"),
                    "run_id": run_id,
                    "log_file": str(log_path.relative_to(ROOT)),
                }
                rows.append(row)
                save_rows(rows)
                completed.add(key)
                finished += 1
                print(f"DONE {method} seed={seed}: final_best={score:.6f}", flush=True)

        restore_memory(snapshot, memory_existed)
    finally:
        # This executes on normal exit and KeyboardInterrupt. If the terminal
        # is force-closed, the persisted snapshot is restored on the next run.
        restore_memory(snapshot, memory_existed)

    by_key = {(int(r["seed"]), str(r["method"])): float(r["final_best"]) for r in rows}
    differences = [
        by_key[(seed, "adaptive_substitution")] - by_key[(seed, "uniform_baseline")]
        for seed in SEEDS
    ]
    if len(differences) != 200:
        raise SystemExit(f"ABORT: only {len(differences)}/200 complete pairs are available.")

    rng = random.Random(271022)
    boot = []
    n = len(differences)
    for _ in range(BOOTSTRAP_REPLICATES):
        boot.append(mean(differences[rng.randrange(n)] for _ in range(n)))
    boot.sort()
    wins = sum(d > 1e-9 for d in differences)
    ties = sum(abs(d) <= 1e-9 for d in differences)
    losses = sum(d < -1e-9 for d in differences)
    contrast = {
        "adaptive": "adaptive_substitution",
        "baseline": "uniform_baseline",
        "n_matched_seeds": n,
        "mean_adaptive_final_best": mean(by_key[(s, "adaptive_substitution")] for s in SEEDS),
        "sd_adaptive_final_best": stdev(by_key[(s, "adaptive_substitution")] for s in SEEDS),
        "mean_baseline_final_best": mean(by_key[(s, "uniform_baseline")] for s in SEEDS),
        "sd_baseline_final_best": stdev(by_key[(s, "uniform_baseline")] for s in SEEDS),
        "mean_paired_delta_adaptive_minus_baseline": mean(differences),
        "sd_paired_delta": stdev(differences),
        "paired_bootstrap_ci95_lower": percentile(boot, 0.025),
        "paired_bootstrap_ci95_upper": percentile(boot, 0.975),
        "wins": wins,
        "ties": ties,
        "losses": losses,
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "interpretation_limit": "sequence-level score only; does not establish target binding",
    }
    with CONTRASTS.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(contrast))
        writer.writeheader()
        writer.writerow(contrast)
    (OUT / "final_summary.json").write_text(json.dumps(contrast, indent=2) + "\n")

    print("\n===== FINAL PAIRED RESULTS (200 SEEDS) =====")
    print(f"Adaptive mean (SD): {contrast['mean_adaptive_final_best']:.6f} ({contrast['sd_adaptive_final_best']:.6f})")
    print(f"Baseline mean (SD): {contrast['mean_baseline_final_best']:.6f} ({contrast['sd_baseline_final_best']:.6f})")
    print(f"Mean paired delta (adaptive - baseline): {contrast['mean_paired_delta_adaptive_minus_baseline']:+.6f}")
    print(f"Paired bootstrap 95% CI: [{contrast['paired_bootstrap_ci95_lower']:+.6f}, {contrast['paired_bootstrap_ci95_upper']:+.6f}]")
    print(f"Wins-ties-losses: {wins}-{ties}-{losses}")
    print(f"Saved: {SUMMARY}")
    print(f"Saved: {CONTRASTS}")
    print("Interpretation limit: sequence-level fitness only; no target-binding claim.")


if __name__ == "__main__":
    main()
