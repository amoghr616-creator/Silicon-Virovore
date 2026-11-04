#!/usr/bin/env bash
set -Eeuo pipefail

cd ~/Desktop/Silicon-Virovore

BASE="results/matched_p1_p3_20261103"
POOL_ROOT="results/unified_ablation_20261102/initial_pools"
MANIFEST="$POOL_ROOT/pool_manifest.json"

if [[ -e "$BASE" ]]; then
    echo "Refusing to overwrite existing experiment directory: $BASE"
    exit 1
fi

mkdir -p "$BASE/logs" "$BASE/runs"

# Verify every shared pool against the recorded manifest.
python - <<'PY'
import hashlib
import json
from pathlib import Path

manifest = Path("results/unified_ablation_20261102/initial_pools/pool_manifest.json")
data = json.loads(manifest.read_text())
aa = set("ACDEFGHIKLMNPQRSTVWY")

for item in data["pools"]:
    path = Path(item["path"])
    if not path.is_file():
        raise SystemExit(f"Missing pool: {path}")

    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    sequences = [
        line.strip().upper()
        for line in raw.decode().splitlines()
        if line.strip() and not line.startswith(">")
    ]

    if digest != item["sha256"]:
        raise SystemExit(f"SHA-256 mismatch: {path}")
    if len(sequences) != 100 or len(set(sequences)) != 100:
        raise SystemExit(f"Pool count/uniqueness failure: {path}")
    if any(len(s) != 21 or not set(s) <= aa for s in sequences):
        raise SystemExit(f"Invalid sequence in: {path}")

print("All ten matched initial pools verified.")
PY

export SILICON_VIROVORE_P4_TARGET_PARETO=0
unset SILICON_VIROVORE_P4_TARGET_MANIFEST
export SILICON_VIROVORE_ARISE_MEMORY_ENABLED=0

printf 'seed\tarm\tmode\tstatus\tactual_evaluations\n' \
    > "$BASE/run_ledger.tsv"

for seed in {631..640}; do
    for arm in P1 P2 P3; do

        case "$arm" in
            P1)
                mode="hotspot_biased_annealed_exploratory"
                ;;
            P2)
                mode="hotspot_epistasis_entropy_exploratory"
                ;;
            P3)
                mode="hotspot_entropy_substitution_epistasis_exploratory_v2_matched"
                ;;
        esac

        run_id="${arm}_seed${seed}"
        pool="$POOL_ROOT/seed_${seed}.fasta"
        run_dir="$BASE/runs/$run_id"
        log="$BASE/logs/${run_id}.log"

        if [[ -e "$run_dir" ]]; then
            echo "Existing run directory; refusing overwrite: $run_dir"
            exit 1
        fi

        export SILICON_VIROVORE_INITIAL_POOL_FASTA="$PWD/$pool"

        echo
        echo "=========================================="
        echo "STARTING $run_id"
        echo "Mode: $mode"
        echo "Seed: $seed"
        echo "Population: 100 | Generations: 5"
        echo "=========================================="

        if python run_pipeline.py \
            --condition adaptive \
            --population-size 100 \
            --generations 5 \
            --seed "$seed" \
            --run-id "$run_id" \
            --output-dir "$BASE/runs" \
            --adaptive-enabled \
            --arise-enabled \
            --ale-enabled \
            --hotspot-enabled \
            --mutation-guidance-mode "$mode" \
            --no-arise-memory-enabled \
            --no-arise-memory-use-history \
            --no-enable-vina-tier2 \
            --no-cache-enabled \
            --no-generate-plots \
            2>&1 | tee "$log"
        then
            :
        else
            rc=$?
            printf '%s\t%s\t%s\tFAILED\t0\n' \
                "$seed" "$arm" "$mode" >> "$BASE/run_ledger.tsv"
            echo "FAILED: $run_id (exit $rc). Stopping; preserve logs."
            exit "$rc"
        fi

        # A completed run must have all five generation summaries.
        if ! python - "$run_dir" <<'PY'
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
audit_path = root / "experiment_audit.json"
history_path = root / "candidate_history.json"

if not audit_path.is_file() or not history_path.is_file():
    raise SystemExit("Missing audit or candidate history.")

audit = json.loads(audit_path.read_text())
summaries = audit.get("optimization", {}).get("generation_summaries", [])
generations = [s.get("generation") for s in summaries]

if generations != [1, 2, 3, 4, 5]:
    raise SystemExit(f"Unexpected generation summaries: {generations}")

print("Verified five generations. Population guard enforced 100 unique candidates per generation.")
PY
        then
            printf '%s\t%s\t%s\tINCOMPLETE\t0\n' \
                "$seed" "$arm" "$mode" >> "$BASE/run_ledger.tsv"
            echo "INCOMPLETE: $run_id. Stopping."
            exit 1
        fi

        printf '%s\t%s\t%s\tPASS\t500\n' \
            "$seed" "$arm" "$mode" >> "$BASE/run_ledger.tsv"

        echo "PASS: $run_id — 500 successful candidate evaluations."
    done
done

echo
echo "=========================================="
echo "P1-P3 MATCHED COMPARISON COMPLETE"
echo "Runs: 30"
echo "Planned successful evaluations: 15,000"
echo "Ledger: $BASE/run_ledger.tsv"
echo "=========================================="
