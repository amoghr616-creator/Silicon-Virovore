#!/bin/bash

set -e

RECEPTOR="$HOME/Desktop/Silicon-Virovore/prototype2/9SYA_chainA_meeko_clean.pdb"
CALIB="$HOME/Desktop/Silicon-Virovore/hybridock-pep/data/calibration_v1_2_production_entropy.json"
OUT="$HOME/Desktop/Silicon-Virovore/prototype2/p2_heldout"

mkdir -p "$OUT"

run_pool () {
    METHOD="$1"
    FILE="$2"

    echo
    echo "============================================================"
    echo "$METHOD HELD-OUT POOL"
    echo "============================================================"

    tail -n +2 "$FILE" | while IFS=, read -r rank sequence; do

        echo
        echo "------------------------------------------------------------"
        echo "$METHOD  $sequence"
        echo "------------------------------------------------------------"

        hybridock-pep dock \
          --peptide "$sequence" \
          --receptor "$RECEPTOR" \
          --site -8.405 -5.823 1.523 \
          --box 52 \
          --n-samples 20 \
          --calibration "$CALIB" \
          --output-dir "$OUT/${METHOD}-${rank}"

    done
}

run_pool "Baseline" \
  "$HOME/Desktop/Silicon-Virovore/prototype2/baseline_gen10_pool.csv"

run_pool "Adaptive" \
  "$HOME/Desktop/Silicon-Virovore/prototype2/adaptive_gen10_pool.csv"

echo
echo "============================================================"
echo "HELD-OUT P2 EVALUATION COMPLETE"
echo "============================================================"
