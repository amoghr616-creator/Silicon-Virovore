#!/usr/bin/env bash
set -Eeuo pipefail

cd ~/Desktop/Silicon-Virovore

BASE="results/matched_p1_p3_abstention_20261103"
POOL_ROOT="results/unified_ablation_20261102/initial_pools"
MANIFEST="$PWD/prototype4/results/real_structure_integration_check_20261101/historical_pose_manifest.json"

if [[ -e "$BASE" ]]; then
    echo "Refusing to overwrite: $BASE"
    exit 1
fi

if [[ ! -f "$MANIFEST" ]]; then
    echo "Target manifest missing."
    exit 1
fi

mkdir -p "$BASE/logs" "$BASE/runs"

# Verify matched pools.
python - <<'PY'
import hashlib
import json
from pathlib import Path

manifest = json.loads(Path(
    "results/unified_ablation_20261102/initial_pools/pool_manifest.json"
).read_text())

for item in manifest["pools"]:
    path = Path(item["path"])
    if not path.is_file():
        raise SystemExit(f"Missing pool: {path}")

    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    sequences = [
        x.strip().upper()
        for x in path.read_text().splitlines()
        if x.strip() and not x.startswith(">")
    ]

    assert digest == item["sha256"], f"Hash mismatch: {path}"
    assert len(sequences) == 100
    assert len(set(sequences)) == 100
    assert all(len(s) == 21 for s in sequences)

print("All 10 matched pools verified.")
PY

export SILICON_VIROVORE_P4_TARGET_MANIFEST="$MANIFEST"
export SILICON_VIROVORE_P4_TARGET_PARETO=0
export SILICON_VIROVORE_ARISE_MEMORY_ENABLED=0
export SILICON_VIROVORE_ARISE_MEMORY_USE_HISTORY=0

printf 'seed\tarm\tgeneration\tcandidates\ttarget_matches\tcoverage\tabstentions\tindependent_run_statuses\n' \
    > "$BASE/abstention_by_generation.tsv"

printf 'seed\tarm\tstatus\tevaluations\n' \
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
            echo "Refusing to overwrite $run_dir"
            exit 1
        fi

        export SILICON_VIROVORE_INITIAL_POOL_FASTA="$PWD/$pool"
        export SILICON_VIROVORE_RUN_ID="$run_id"

        echo
        echo "STARTING $run_id"

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
            printf '%s\t%s\tFAILED\t0\n' "$seed" "$arm" \
                >> "$BASE/run_ledger.tsv"
            echo "Run failed: $run_id"
            exit "$rc"
        fi

        # Validate all five generations and extract P4 abstention data.
        python - "$run_dir" "$pool" "$BASE/abstention_by_generation.tsv" "$seed" "$arm" <<'PY'
import csv
import json
import sys
from pathlib import Path
from collections import Counter

root = Path(sys.argv[1])
pool_path = Path(sys.argv[2])
ledger_path = Path(sys.argv[3])
seed, arm = sys.argv[4], sys.argv[5]

def read_json(path):
    if not path.is_file():
        raise RuntimeError(f"Missing required report: {path}")
    return json.loads(path.read_text())

pool_sequences = {
    line.strip().upper()
    for line in pool_path.read_text().splitlines()
    if line.strip() and not line.startswith(">")
}

rows = []

for generation in range(1, 6):
    target_path = root / "p4_target_evidence" / (
        f"generation_{generation:04d}.json"
    )
    assessment_path = root / "p4_integrated_assessment" / (
        f"generation_{generation:04d}.json"
    )

    target = read_json(target_path)
    assessment = read_json(assessment_path)

    if target["candidate_count"] != 100:
        raise RuntimeError(
            f"{arm} seed {seed} generation {generation}: "
            "candidate count is not 100."
        )

    if assessment["candidate_count"] != 100:
        raise RuntimeError("P4 assessment candidate count mismatch.")

    if target["pareto_requested"] is not False:
        raise RuntimeError("Pareto should be disabled in this comparison.")

    if target["pareto_reordered_candidates"] is not False:
        raise RuntimeError("Unexpected Pareto reordering.")

    if target["manifest"].get("status") != "loaded":
        raise RuntimeError("Target manifest was not successfully loaded.")

    if generation == 1:
        observed = {
            row["sequence"].upper()
            for row in target["candidates"]
        }
        if observed != pool_sequences:
            raise RuntimeError("Generation-1 pool mismatch.")

    gate = assessment["integrated_gate"]
    if gate["status"] != "candidate_level_adjudication_executed":
        raise RuntimeError("P4 candidate gate was not executed.")

    records = target.get("candidates", [])
    reliability_statuses = Counter(
        (row.get("reliability") or {}).get("status", "missing")
        for row in records
    )

    rows.append({
        "seed": seed,
        "arm": arm,
        "generation": generation,
        "candidates": target["candidate_count"],
        "target_matches": target["target_evidence_count"],
        "coverage": target["target_evidence_coverage"],
        "abstentions": gate["abstain_count"],
        "independent_run_statuses": json.dumps(
            reliability_statuses, sort_keys=True
        ),
    })

with open(ledger_path, "a", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(rows[0]),
        delimiter="\t",
    )
    writer.writerows(rows)

print(f"Verified {arm} seed {seed}: 5 generations, 500 candidates.")
PY

        printf '%s\t%s\tPASS\t500\n' "$seed" "$arm" \
            >> "$BASE/run_ledger.tsv"

        echo "PASS: $run_id"
    done
done

echo
echo "P1-P3 COMPARISON FINISHED"
echo "Runs: 30"
echo "Planned successful evaluations: 15,000"
echo "Run ledger: $BASE/run_ledger.tsv"
echo "P4 abstention data: $BASE/abstention_by_generation.tsv"
