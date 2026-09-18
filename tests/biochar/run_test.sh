#!/usr/bin/env bash
# Integration test for the DSSAT biochar model.
#
# Usage:
#   ./run_test.sh [DSSAT_BINARY] [RUN_DIR]
#
# Arguments:
#   DSSAT_BINARY  path to dscsm048 (default: ../../build/dscsm048)
#   RUN_DIR       directory containing experiment input files for a
#                 one-season maize or wheat run with biochar applied
#                 (default: ./testdata if it exists)
#
# The script:
#   1. Runs DSSAT for one season using the provided input data.
#   2. Calls check_biochar_out.py on the resulting BIOCHAR.OUT.
#   3. Exits 0 on success, non-zero on failure.
#
# NOTE: This test requires a complete DSSAT input dataset (FILEX, .SOL,
# .WTH) with a *BIOCHAR section.  A minimal synthetic dataset can be
# generated following the instructions in Data/BIOCHAR_template.X.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

DSSAT_BIN="${1:-$REPO_ROOT/build/dscsm048}"
RUN_DIR="${2:-$SCRIPT_DIR/testdata}"
CHECKER="$SCRIPT_DIR/check_biochar_out.py"

# --------------------------------------------------------------------------
echo "=== DSSAT Biochar Integration Test ==="

# 1. Check prerequisites
if [ ! -f "$DSSAT_BIN" ]; then
    echo "ERROR: binary not found: $DSSAT_BIN"
    echo "  Build the model first: mkdir build && cd build && cmake .. && make"
    exit 2
fi

if [ ! -d "$RUN_DIR" ]; then
    echo "SKIP: test data directory not found: $RUN_DIR"
    echo "  Provide a directory with FILEX + .SOL + .WTH + *BIOCHAR section."
    exit 0
fi

# Find the batch file or FILEX in RUN_DIR
BATCH_FILE="$(find "$RUN_DIR" -maxdepth 1 -name 'DSSBatch*' | head -1)"
if [ -z "$BATCH_FILE" ]; then
    echo "ERROR: No DSSBatch file found in $RUN_DIR"
    exit 2
fi

# 2. Run DSSAT
echo "Running: $DSSAT_BIN MZCER048 B $BATCH_FILE"
cd "$RUN_DIR"
"$DSSAT_BIN" MZCER048 B "$(basename "$BATCH_FILE")" > dssat_run.log 2>&1
STATUS=$?
if [ $STATUS -ne 0 ]; then
    echo "ERROR: DSSAT run failed (exit $STATUS) — see $RUN_DIR/dssat_run.log"
    exit 1
fi

# 3. Validate BIOCHAR.OUT
BIOCHAR_OUT="$RUN_DIR/BIOCHAR.OUT"
if [ ! -f "$BIOCHAR_OUT" ]; then
    echo "ERROR: BIOCHAR.OUT was not produced"
    exit 1
fi

echo "Validating $BIOCHAR_OUT ..."
python3 "$CHECKER" "$BIOCHAR_OUT"

echo "=== TEST PASSED ==="
