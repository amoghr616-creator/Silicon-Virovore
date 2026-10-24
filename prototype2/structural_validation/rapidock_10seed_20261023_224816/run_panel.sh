#!/bin/bash
set -u

echo "Panel dispatcher started: $(date)"
echo "Waiting for pilot PID $PILOT_PID to finish..."

while ps -p "$PILOT_PID" -o stat= 2>/dev/null \
    | grep -v 'Z' | grep -q '[A-Za-z]'; do
    sleep 5
done

echo "Pilot finished. Launching ten independent seeds: $(date)"

run_one() {
    seed="$1"
    out="$PANEL/seed_$seed"
    log="$PANEL/seed_$seed.log"

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
    echo "END seed=$seed exit=$rc time=$(date)"
    return 0
}

for seed in 81101 81102 81103 81104 81105 \
            81106 81107 81108 81109 81110; do
    run_one "$seed" &
    pids="$pids $!"
    count=$(( ${count:-0} + 1 ))

    if [ "$count" -eq 2 ]; then
        for pid in $pids; do wait "$pid"; done
        count=0
        pids=""
    fi
done

for pid in ${pids:-}; do wait "$pid"; done

echo "TEN-SEED PANEL DISPATCH COMPLETE: $(date)"
