#!/bin/bash

echo "10-seed panel started: $(date)"
pids=""
count=0

run_seed() {
    seed="$1"
    out="$PANEL/seed_$seed"
    log="$PANEL/seed_$seed.log"
    dlog="$PANEL/seed_$seed.dispatch.log"

    {
        echo "START seed=$seed time=$(date)"
        "$CLI" dock \
            --peptide "LKAAHEALLKFFWEWSDLLKR" \
            --receptor "$CROP" \
            --site -8.405 -5.823 1.523 \
            --box 32 \
            --long-checkpoint-threshold 22 \
            --n-samples 1 \
            --seed "$seed" \
            --output-dir "$out" > "$log" 2>&1
        rc=$?
        echo "$rc" > "$PANEL/seed_$seed.exitcode"
        echo "END seed=$seed exit=$rc time=$(date)"
    } >> "$dlog" 2>&1
}

for seed in 81101 81102 81103 81104 81105 81106 81107 81108 81109 81110; do
    if [ -f "$PANEL/seed_$seed.exitcode" ]; then
        echo "SKIP seed=$seed: exitcode record already exists"
        continue
    fi

    run_seed "$seed" &
    pid=$!
    pids="$pids $pid"
    count=$((count + 1))
    echo "QUEUED seed=$seed pid=$pid time=$(date)"

    if [ "$count" -ge 3 ]; then
        for p in $pids; do wait "$p" || true; done
        pids=""
        count=0
    fi
done

for p in $pids; do wait "$p" || true; done
echo "PANEL DISPATCH COMPLETE: $(date)"
