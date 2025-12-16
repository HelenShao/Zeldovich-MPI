#!/bin/bash
# Quick test script for RNG consistency verification
# Tests that overlapping modes use the same random numbers across different N values

set -e  # Exit on error

# Get the project root directory (parent of v15_3_tests)
PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT_ROOT"

echo "===================================================================================="
echo "RNG Consistency Test"
echo "===================================================================================="
echo ""

# Check if executable exists
if [ ! -f "./hermitian_3d_matrix" ]; then
    echo "Error: Executable not found. Building..."
    make clean
    make
fi

# Check if parameter file exists
PARAM_FILE="examples/param_template.par"
if [ ! -f "$PARAM_FILE" ]; then
    echo "Warning: Parameter file $PARAM_FILE not found. Using default (no param file)."
    PARAM_FILE=""
fi

# Test configuration
N1=${1:-4}   # First N value (default: 4)
N2=${2:-6}   # Second N value (default: 6)
NP=${3:-2}   # Number of MPI ranks (default: 2)

echo "Test configuration:"
echo "  N1 = $N1"
echo "  N2 = $N2"
echo "  MPI ranks = $NP"
echo "  Parameter file = ${PARAM_FILE:-none}"
echo ""

# Create output directory
OUTPUT_DIR="v15_3_tests"
mkdir -p "$OUTPUT_DIR"

# Run first test
echo "Running test with N=$N1..."
if [ -n "$PARAM_FILE" ]; then
    mpirun -np $NP ./hermitian_3d_matrix $N1 "$PARAM_FILE" > "$OUTPUT_DIR/output_N${N1}.out" 2>&1
else
    mpirun -np $NP ./hermitian_3d_matrix $N1 > "$OUTPUT_DIR/output_N${N1}.out" 2>&1
fi

if [ $? -ne 0 ]; then
    echo "Error: Test with N=$N1 failed!"
    exit 1
fi

# Run second test
echo "Running test with N=$N2..."
if [ -n "$PARAM_FILE" ]; then
    mpirun -np $NP ./hermitian_3d_matrix $N2 "$PARAM_FILE" > "$OUTPUT_DIR/output_N${N2}.out" 2>&1
else
    mpirun -np $NP ./hermitian_3d_matrix $N2 > "$OUTPUT_DIR/output_N${N2}.out" 2>&1
fi

if [ $? -ne 0 ]; then
    echo "Error: Test with N=$N2 failed!"
    exit 1
fi

# Check if debug output was generated
if ! grep -q "\[RNG-DEBUG\]" "$OUTPUT_DIR/output_N${N1}.out"; then
    echo "Warning: No [RNG-DEBUG] lines found in output_N${N1}.out"
    echo "  Make sure DEBUG_RNG_CONSISTENCY=1 is set in src/config.h"
fi

if ! grep -q "\[RNG-DEBUG\]" "$OUTPUT_DIR/output_N${N2}.out"; then
    echo "Warning: No [RNG-DEBUG] lines found in output_N${N2}.out"
    echo "  Make sure DEBUG_RNG_CONSISTENCY=1 is set in src/config.h"
fi

# Compare results
echo ""
echo "Comparing results..."
echo "===================================================================================="
cd "$OUTPUT_DIR"
python compare_rng_consistency.py "output_N${N1}.out" "output_N${N2}.out"

echo ""
echo "===================================================================================="
echo "Test complete!"
echo "Output files:"
echo "  - $OUTPUT_DIR/output_N${N1}.out"
echo "  - $OUTPUT_DIR/output_N${N2}.out"
echo "===================================================================================="

