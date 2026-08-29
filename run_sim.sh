#!/bin/bash

# ============================================================
# Silicon Virovore
# Computational Peptide Engineering Platform
# Complete Pipeline Launcher
# ============================================================

set -e

PROJECT_ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$PROJECT_ROOT"

echo "PROJECT_ROOT = $PROJECT_ROOT"
pwd
echo ""
echo "Contents of c/:"
ls -l c
echo ""

BUILD_DIR="build"
LOG_DIR="logs"

mkdir -p "$BUILD_DIR"
mkdir -p "$LOG_DIR"

echo ""
echo "============================================================"
echo " Silicon Virovore"
echo " Computational Peptide Engineering Platform"
echo "============================================================"
echo ""

# ------------------------------------------------------------
# Clean Previous Outputs
# ------------------------------------------------------------

echo "[0/4] Cleaning previous outputs..."

rm -f metrics.csv
rm -f evolution_consensus.csv
rm -f virovore_dashboard.png
rm -f "$BUILD_DIR/virovore_engine"

echo "Done."
echo ""

# ------------------------------------------------------------
# Verify Source Files
# ------------------------------------------------------------

echo "[1/4] Verifying C source files..."

FILES=(
    "c/engine.c"
    "c/hydropathy.c"
    "c/ga_loop.c"
)

for file in "${FILES[@]}"
do
    if [ ! -f "$file" ]; then
        echo ""
        echo "ERROR: Missing source file:"
        echo "  $file"
        exit 1
    fi
done

echo "All source files located."
echo ""

# ------------------------------------------------------------
# Compile Native Backend
# ------------------------------------------------------------

echo "[2/4] Compiling native backend..."

clang \
    -O3 \
    c/engine.c \
    c/hydropathy.c \
    c/ga_loop.c \
    -o "$BUILD_DIR/virovore_engine"

echo "Compilation successful."
echo ""

# ------------------------------------------------------------
# Execute Simulation
# ------------------------------------------------------------

echo "[3/4] Running evolutionary simulation..."

"$BUILD_DIR/virovore_engine"

echo "Simulation finished."
echo ""

# ------------------------------------------------------------
# Launch Dashboard
# ------------------------------------------------------------

echo "[4/4] Launching visualization..."

if [ -f ".venv/bin/python" ]; then
    PYTHON=".venv/bin/python"
else
    PYTHON="python3"
fi

if [ -f "dashboard/dashboard.py" ]; then
    "$PYTHON" dashboard/dashboard.py
elif [ -f "dashboard.py" ]; then
    "$PYTHON" dashboard.py
else
    echo "Dashboard script not found."
fi

echo ""
echo "============================================================"
echo " Silicon Virovore completed successfully."
echo "============================================================"
echo ""