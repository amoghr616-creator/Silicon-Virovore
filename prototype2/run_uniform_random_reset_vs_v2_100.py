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
OUT = ROOT / "prototype2/p2_uniform_random_reset_vs_v2_100_20261023"
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
SEEDS = [30100023 + 1000 * i for i in range(100)]
BOOTSTRAP_REPLICATES = 20_000

METHOD_FLAGS = [
    "--no-adaptive-enabled",
    "--arise-enabled",
    "--no-ale-enabled",
    "--hotspot-enabled",
]

METHODS = {
    "adaptive_substitution": {
        "mode": "hotspot_entropy_substitution_epistasis_exploratory_v2",
        "condition": "adaptive",
        "flags": list(METHOD_FLAGS),
    },
    "uniform_baseline": {
        "mode": "uniform_random_reset",
        "condition": "baseline",
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
    "seed", "method", "cli_mode", "best_across_generations",
    "generation10_best", "mean_mutation_count",
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
    if len(SEEDS) != 100 or len(set(SEEDS)) != 100:
        raise SystemExit("ABORT: expected 100 unique matched seeds.")
    if METHODS["adaptive_substitution"]["flags"] != METHODS["uniform_baseline"]["flags"]:
        raise SystemExit("ABORT: arm-level flags differ; only mutation-guidance mode may vary.")
    historical_seeds = (
        {20700021 + 1000 * i for i in range(200)}
        | {21100021 + 1000 * i for i in range(200)}
    )
    if set(SEEDS) & historical_seeds:
        raise SystemExit("ABORT: new seed list overlaps a prior 200-seed range.")

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
        "experiment": "Fresh 100-pair matched-seed comparison: adaptive V2 versus fixed-parent uniform random-mutation baseline",
        "created_at": datetime.now().isoformat(),
        "seed_count": 100,
        "seeds": SEEDS,
        "population_size": POPULATION,
        "generations": GENERATIONS,
        "mutation_rate": MUTATION_RATE,
        "elite_count": ELITE_COUNT,
        "seed_sequence": SEED_SEQUENCE,
        "methods": METHODS,
        "design_status": "Fresh-seed replication specified after prior exploratory results were inspected; not a preregistered confirmatory experiment.",
        "seed_rule": "30100023 + 1000*i for i=0..99; fresh fixed seed range disjoint from the previous 200-seed ranges.",
        "primary_endpoint": "Maximum best_overall score across all 10 generations; paired adaptive-minus-uniform_random_reset difference. Generation-10 score is secondary.",
        "comparability": "Both arms use population 16, 10 generations, mutation_rate 0.05, the same paired seed, and the same sequence-level objective. Adaptive V2 uses its existing guided proposals and best-parent carry-forward. The baseline proposes uniform one-position mutants from the original seed each generation and injects no preserved parent. Both evaluate 160 candidates per run. This is an end-to-end algorithm comparison, not a mutation-guidance-only ablation.",
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
        row["best_across_generations"] = float(row["best_across_generations"])
        row["generation10_best"] = float(row["generation10_best"])
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
        if "hotspot_entropy_substitution_epistasis_exploratory_v2" not in help_text:
            raise SystemExit("ABORT: current CLI does not expose the requested adaptive mode.")
        if "uniform_random_reset" not in help_text:
            raise SystemExit("ABORT: current CLI does not expose uniform_random_reset.")
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
                run_id = f"p2-rreset100-{method}-seed{seed}"
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
                    "--condition", str(spec["condition"]),
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
                summaries = [
                    line for line in lines
                    if re.search(r"Generation \d+ summary:", line)
                ]
                if len(summaries) != GENERATIONS:
                    raise SystemExit(
                        f"Expected {GENERATIONS} generation summaries, "
                        f"found {len(summaries)} in {log_path}"
                    )

                parsed_scores = [
                    score_from_summary(line, "best_overall")
                    for line in summaries
                ]
                if any(
                    score is None or not math.isfinite(score)
                    for score in parsed_scores
                ):
                    raise SystemExit(
                        f"Could not parse every generation's best score in {log_path}"
                    )

                generation_scores = [float(score) for score in parsed_scores]
                score = max(generation_scores)
                generation10_best = generation_scores[-1]
                last_summary = summaries[-1]

                # Audit the baseline's defining behavior in every generation:
                # original seed stays fixed and every candidate is a single mutant.
                if method == "uniform_baseline":
                    seed_lines = [
                        line for line in lines
                        if re.search(r"INFO:__main__:Seed:\s*", line)
                    ]
                    if (
                        len(seed_lines) != GENERATIONS
                        or any(
                            not line.rstrip().endswith(SEED_SEQUENCE)
                            for line in seed_lines
                        )
                    ):
                        raise SystemExit(
                            f"Baseline did not retain the original reference seed "
                            f"for all generations; inspect {log_path}"
                        )

                    for generation_number, summary_line in enumerate(
                        summaries, start=1
                    ):
                        mean_mutations = score_from_summary(
                            summary_line, "mean_mutation_count"
                        )
                        if (
                            mean_mutations is None
                            or not math.isclose(
                                mean_mutations, 1.0, rel_tol=0, abs_tol=1e-6
                            )
                        ):
                            raise SystemExit(
                                f"Baseline generation {generation_number} does "
                                f"not show the expected one-mutation population; "
                                f"inspect {log_path}"
                            )
                row = {
                    "seed": seed,
                    "method": method,
                    "cli_mode": spec["mode"],
                    "best_across_generations": f"{score:.12f}",
                    "generation10_best": f"{generation10_best:.12f}",
                    "mean_mutation_count": score_from_summary(last_summary, "mean_mutation_count"),
                    "mean_pairwise_identity": score_from_summary(last_summary, "mean_pairwise_identity"),
                    "run_id": run_id,
                    "log_file": str(log_path.relative_to(ROOT)),
                }
                rows.append(row)
                save_rows(rows)
                completed.add(key)
                finished += 1
                print(f"DONE {method} seed={seed}: best_across_generations={score:.6f}", flush=True)

        restore_memory(snapshot, memory_existed)
    finally:
        # This executes on normal exit and KeyboardInterrupt. If the terminal
        # is force-closed, the persisted snapshot is restored on the next run.
        restore_memory(snapshot, memory_existed)

    by_key = {(int(r["seed"]), str(r["method"])): float(r["best_across_generations"]) for r in rows}
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
        "mean_adaptive_best_across_generations": mean(by_key[(s, "adaptive_substitution")] for s in SEEDS),
        "sd_adaptive_best_across_generations": stdev(by_key[(s, "adaptive_substitution")] for s in SEEDS),
        "mean_baseline_best_across_generations": mean(by_key[(s, "uniform_baseline")] for s in SEEDS),
        "sd_baseline_best_across_generations": stdev(by_key[(s, "uniform_baseline")] for s in SEEDS),
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
    print(f"Adaptive mean (SD): {contrast['mean_adaptive_best_across_generations']:.6f} ({contrast['sd_adaptive_best_across_generations']:.6f})")
    print(f"Baseline mean (SD): {contrast['mean_baseline_best_across_generations']:.6f} ({contrast['sd_baseline_best_across_generations']:.6f})")
    print(f"Mean paired delta (adaptive - baseline): {contrast['mean_paired_delta_adaptive_minus_baseline']:+.6f}")
    print(f"Paired bootstrap 95% CI: [{contrast['paired_bootstrap_ci95_lower']:+.6f}, {contrast['paired_bootstrap_ci95_upper']:+.6f}]")
    print(f"Wins-ties-losses: {wins}-{ties}-{losses}")
    print(f"Saved: {SUMMARY}")
    print(f"Saved: {CONTRASTS}")
    print("Interpretation limit: sequence-level fitness only; no target-binding claim.")


if __name__ == "__main__":
    main()
