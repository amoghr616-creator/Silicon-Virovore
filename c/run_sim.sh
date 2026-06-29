#!/bin/bash

# Clear previous data points to guarantee clean live run streams
rm -f metrics.csv evolution_consensus.csv virovore_engine

echo ""
echo "============================================================="
echo "[LAUNCH] Starting Silicon Virovore Simulation Pipeline..."
echo "============================================================="
echo ""

# 1. Compile backend assets natively via high optimization flags
echo "[STEP 1/3] Compiling optimized biophysical backend structures..."
clang -O3 engine.c hydropathy.c ga_loop.c -o virovore_engine

if [ $? -ne 0 ]; then
    echo "CRITICAL ERROR: Native compilation pipeline failed to construct executable."
    exit 1
fi

# 2. Run simulation engine to produce data metrics logs
echo "[STEP 2/3] Executing 10,000 epoch evolutionary sweep..."
./virovore_engine

if [ $? -ne 0 ]; then
    echo "CRITICAL ERROR: Simulation crashed or memory fault detected."
    exit 1
fi

# 3. Fire up visualization pipeline using virtual environment wrapper
echo "[STEP 3/3] Passing raw logs to Python 6-Panel Visualization Suite..."
/Users/centurion616/Desktop/Silicon-Virovore/.venv/bin/python dashboard.py

echo ""
echo "============================================================="
echo "[COMPLETE] 6-Panel Visual Suite generated as virovore_dashboard.png"
echo "============================================================="
echo ""