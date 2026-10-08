#!/bin/bash

set -e

ROOT="$HOME/Desktop/Silicon-Virovore"
RECEPTOR="$ROOT/prototype2/9SYA_chainA_meeko_clean.pdb"
CAL="$ROOT/hybridock-pep/data/calibration_v1_2_production_entropy.json"
OUT="$ROOT/prototype2/p2_rapidock"

run_candidate() {
    NAME="$1"
    SEQ="$2"

    echo
    echo "============================================================"
    echo "$NAME  $SEQ"
    echo "============================================================"

    hybridock-pep dock \
      --peptide "$SEQ" \
      --receptor "$RECEPTOR" \
      --site -8.405 -5.823 1.523 \
      --box 52 \
      --n-samples 20 \
      --calibration "$CAL" \
      --output-dir "$OUT/$NAME"
}

run_candidate "P2-02" "DKLAVTALLVVFAESSDKLRR"
run_candidate "P2-03" "MKLAPFALLVVFAGASDWIRR"
run_candidate "P2-04" "MKLAVFALLVCFAESSDLVRR"
run_candidate "P2-05" "MKQAVFALLQVFAGSSDWIRR"
run_candidate "P2-06" "MKLAVFALRVFFAESSDAIRR"

echo
echo "ALL P2 CANDIDATES COMPLETE"
