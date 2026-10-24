from pathlib import Path
from statistics import mean, stdev
from datetime import datetime
import csv
import hashlib
import json
import os
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
os.chdir(ROOT)

OUT = ROOT / "prototype2/substitution_corrected_pilot_20261022"
if OUT.exists():
    raise SystemExit(f"Output directory already exists: {OUT}")

OUT.mkdir(parents=True)

POPULATION = 16
GENERATIONS = 10
MUTATION_RATE = 0.05
ELITE_COUNT = 2
SEEDS = [20600021 + 1000 * i for i in range(10)]

# Experiment labels map to actual CLI modes.
METHODS = [
    ("substitution_specific",
     "hotspot_entropy_substitution_epistasis_exploratory"),
    ("previous_entropy",
     "epistasis_entropy"),
    ("pairwise_uniform_control",
     "uniform_pair_matched"),
    ("uniform_baseline",
     "uniform"),
]

CSV_PATH = OUT / "run_summary.csv"
LOG_DIR = OUT / "logs"
LOG_DIR.mkdir()

def save_csv(rows):
    columns = [
        "seed", "method", "cli_mode", "final_best",
        "mean_mutation_count", "mean_pairwise_identity", "log_file"
    ]
    with CSV_PATH.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)

def summary_metric(line, name):
    match = re.search(
        rf"""['"]{re.escape(name)}['"]\s*:\s*
            (-?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)""",
        line,
        re.X,
    )
    return float(match.group(1)) if match else None

# Verify that the CLI option exists and check the listed mode choices
# whenever argparse exposes them.
help_result = subprocess.run(
    [sys.executable, "run_pipeline.py", "--help"],
    text=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
)
help_text = help_result.stdout
if "--mutation-guidance-mode" not in help_text:
    raise SystemExit("ABORT: mutation-guidance-mode CLI option not found.")

choice_match = re.search(
    r"--mutation-guidance-mode[^\n]*?\{([^}]+)\}", help_text
)
if choice_match:
    allowed = {x.strip() for x in choice_match.group(1).split(",")}
    invalid = [mode for _, mode in METHODS if mode not in allowed]
    if invalid:
        raise SystemExit(
            "ABORT: modes not listed by the CLI: " + ", ".join(invalid)
        )

manifest = {
    "created_at": datetime.now().isoformat(),
    "seed_count": len(SEEDS),
    "seeds": SEEDS,
    "population_size": POPULATION,
    "generations": GENERATIONS,
    "mutation_rate": MUTATION_RATE,
    "elite_count": ELITE_COUNT,
    "seed_sequence": "MKLAVFALLVFFAGSSDLIRR",
    "methods": [{"label": a, "cli_mode": b} for a, b in METHODS],
    "docking": "same defaults as the prior CLI replication",
    "git_commit": subprocess.run(
        ["git", "rev-parse", "HEAD"],
        text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL
    ).stdout.strip(),
}
(OUT / "experiment_manifest.json").write_text(
    json.dumps(manifest, indent=2) + "\n"
)
(OUT / "seed_list.txt").write_text(
    "\n".join(map(str, SEEDS)) + "\n"
)

# Preserve memory because the pipeline may write its memory file even when
# prior-memory use is disabled. All pilot conditions should start without
# learned history influencing the next run.
memory_path = ROOT / "results/arise_memory.json"
memory_existed = memory_path.exists()
memory_original = memory_path.read_bytes() if memory_existed else None

rows = []
try:
    total = len(SEEDS) * len(METHODS)
    completed = 0

    for seed in SEEDS:
        for label, cli_mode in METHODS:
            completed += 1
            logfile = LOG_DIR / f"{label}_seed{seed}.log"
            print(
                f"[{completed}/{total}] START {label} seed={seed}",
                flush=True,
            )

            command = [
                sys.executable, "run_pipeline.py",
                "--population-size", str(POPULATION),
                "--generations", str(GENERATIONS),
                "--mutation-rate", str(MUTATION_RATE),
                "--elite-count", str(ELITE_COUNT),
                "--seed", str(seed),
                "--mutation-guidance-mode", cli_mode,
            ]

            # Use unique run IDs if supported by this version of the CLI.
            if "--run-id" in help_text:
                command += ["--run-id", f"pilot-{label}-seed{seed}"]

            result = subprocess.run(
                command,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=600,
            )
            logfile.write_text(result.stdout)

            if result.returncode != 0:
                print("\nRUN FAILED. Last 50 log lines:", flush=True)
                print("\n".join(result.stdout.splitlines()[-50:]))
                raise SystemExit(
                    f"Pilot stopped at {label}, seed={seed}; "
                    f"log: {logfile}"
                )

            # Read the exact final-generation score from the generation
            # summary, not the rounded "Best Overall Score" log line.
            summary_lines = [
                line for line in result.stdout.splitlines()
                if f"Generation {GENERATIONS} summary:" in line
            ]
            if not summary_lines:
                raise SystemExit(
                    f"Could not find final-generation summary for "
                    f"{label}, seed={seed}. Inspect {logfile}"
                )

            line = summary_lines[-1]
            best = summary_metric(line, "best_overall")
            if best is None:
                raise SystemExit(
                    f"Could not parse final best score for "
                    f"{label}, seed={seed}. Inspect {logfile}"
                )

            row = {
                "seed": seed,
                "method": label,
                "cli_mode": cli_mode,
                "final_best": best,
                "mean_mutation_count": summary_metric(
                    line, "mean_mutation_count"
                ),
                "mean_pairwise_identity": summary_metric(
                    line, "mean_pairwise_identity"
                ),
                "log_file": str(logfile.relative_to(ROOT)),
            }
            rows.append(row)
            save_csv(rows)
            print(
                f"DONE {label} seed={seed}: final_best={best:.6f}",
                flush=True,
            )

finally:
    if memory_existed:
        memory_path.parent.mkdir(parents=True, exist_ok=True)
        memory_path.write_bytes(memory_original)
    elif memory_path.exists():
        memory_path.unlink()

if len(rows) != len(SEEDS) * len(METHODS):
    raise SystemExit(f"Incomplete pilot: {len(rows)}/{total} runs recorded.")

by_method = {
    name: [r for r in rows if r["method"] == name]
    for name, _ in METHODS
}

print("\n===== PILOT METHOD SUMMARY =====")
means = {}
for name, _ in METHODS:
    scores = [r["final_best"] for r in by_method[name]]
    means[name] = mean(scores)
    print(
        f"{name}: n={len(scores)}, mean={mean(scores):.6f}, "
        f"SD={stdev(scores):.6f}"
    )

primary = "substitution_specific"
comparators = [
    "uniform_baseline",
    "pairwise_uniform_control",
    "previous_entropy",
]

print("\n===== PAIRED COMPARISONS (NEW METHOD MINUS COMPARATOR) =====")
paired = {}
for comparator in comparators:
    a = {r["seed"]: r["final_best"] for r in by_method[primary]}
    b = {r["seed"]: r["final_best"] for r in by_method[comparator]}
    differences = [a[s] - b[s] for s in SEEDS]
    wins = sum(d > 1e-9 for d in differences)
    ties = sum(abs(d) <= 1e-9 for d in differences)
    losses = sum(d < -1e-9 for d in differences)
    paired[comparator] = {
        "mean_delta": mean(differences),
        "wins": wins,
        "ties": ties,
        "losses": losses,
    }
    print(
        f"{primary} vs {comparator}: "
        f"mean_delta={mean(differences):+.6f}; "
        f"W-T-L={wins}-{ties}-{losses}"
    )

# Screening flag only; this is not a significance test.
screen_pass = (
    paired["uniform_baseline"]["mean_delta"] > 0
    and paired["uniform_baseline"]["wins"] >= 7
    and paired["pairwise_uniform_control"]["mean_delta"] > 0
)
print("\n===== PILOT SCREEN =====")
print(
    "PASS: candidate merits a fresh 126-seed confirmation"
    if screen_pass
    else "NO-PASS: current pilot does not justify the 126-seed expansion"
)
print("This screen is exploratory, not a statistical significance claim.")
print("Output directory:", OUT)
print("Run summary:", CSV_PATH)
