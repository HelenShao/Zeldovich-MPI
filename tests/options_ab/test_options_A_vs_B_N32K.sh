#!/bin/bash
# ====================================================================================
# COMPREHENSIVE TEST: Options A vs B for N=2048
# ====================================================================================
# This script tests both Option A and Option B with N=2048 and compares:
#   1. Performance (memory usage, execution time)
#   2. Correctness (output file comparison)
#   3. Scalability metrics
#
# Usage:
#   ./test_options_A_vs_B_N32K.sh [param_file] [num_ranks] [num_nodes]
#
# Example:
#   ./test_options_A_vs_B_N32K.sh examples/param_N2048.par 256 16
# ====================================================================================

set -e  # Exit on error

# Configuration
N=2048
PARAM_FILE="${1:-examples/param_N2048.par}"
NUM_RANKS="${2:-256}"   # 16×16 = 256 ranks (default for N=2048)
NUM_NODES="${3:-16}"    # 16 ranks per node (adjust based on your cluster)
MAIN_EXEC="./main"

# Derived configuration
GRID_X=16
GRID_Z=16
X_COUNT=$((N / GRID_X))  # 128 per rank

echo "========================================================================"
echo "OPTIONS A vs B COMPARISON TEST: N=$N"
echo "========================================================================"
echo "Configuration:"
echo "  N: $N"
echo "  Ranks: $NUM_RANKS (${GRID_X}×${GRID_Z})"
echo "  Nodes: $NUM_NODES"
echo "  X-count per rank: $X_COUNT"
echo "  Parameter file: $PARAM_FILE"
echo "========================================================================"
echo ""

# Check prerequisites
if [ ! -f "$MAIN_EXEC" ]; then
    echo "ERROR: $MAIN_EXEC not found. Please compile first."
    exit 1
fi

if [ ! -f "$PARAM_FILE" ]; then
    echo "ERROR: Parameter file $PARAM_FILE not found"
    echo "Creating minimal parameter file..."
    # Could create one here if needed
    exit 1
fi

# Create output directories
OUTPUT_DIR_A="output_option_a"
OUTPUT_DIR_B="output_option_b"
TIMING_DIR="timing_results"

mkdir -p "$OUTPUT_DIR_A" "$OUTPUT_DIR_B" "$TIMING_DIR"

echo "========================================================================"
echo "TEST 1: OPTION A (Transpose + WriteParticlesSlab_range)"
echo "========================================================================"

# Compile with Option A
echo "Compiling with Option A (USE_PARTICLE_OUTPUT_OPTION_B=0)..."
make clean
make CFLAGS="-DUSE_PARTICLE_OUTPUT_OPTION_B=0" -j$(nproc)

# Clean previous outputs
rm -rf rank_* "$OUTPUT_DIR_A"/*

# Run Option A with timing
echo ""
echo "Running Option A..."
echo "Command: mpirun -np $NUM_RANKS $MAIN_EXEC $N $PARAM_FILE"
echo ""

START_TIME_A=$(date +%s.%N)
/usr/bin/time -v mpirun -np $NUM_RANKS $MAIN_EXEC $N $PARAM_FILE 2>&1 | tee "${TIMING_DIR}/option_a_timing.log"
END_TIME_A=$(date +%s.%N)

# Move outputs
if [ -d "output" ]; then
    mv output/* "$OUTPUT_DIR_A/" 2>/dev/null || true
fi

# Extract timing information
ELAPSED_A=$(echo "$END_TIME_A - $START_TIME_A" | bc)
echo "Option A elapsed time: ${ELAPSED_A} seconds" | tee "${TIMING_DIR}/option_a_summary.txt"

# Count files
NUM_FILES_A=$(find "$OUTPUT_DIR_A" -name "ic_rank*" 2>/dev/null | wc -l)
echo "Option A: $NUM_FILES_A particle IC files created" | tee -a "${TIMING_DIR}/option_a_summary.txt"

# Get memory info from timing log if available
if grep -q "Maximum resident set size" "${TIMING_DIR}/option_a_timing.log"; then
    MAX_RSS_A=$(grep "Maximum resident set size" "${TIMING_DIR}/option_a_timing.log" | awk '{print $6}')
    echo "Option A: Max RSS per rank: ${MAX_RSS_A} KB" | tee -a "${TIMING_DIR}/option_a_summary.txt"
fi

echo ""
echo "Option A complete. Outputs in: $OUTPUT_DIR_A"
echo ""

# Wait a bit between tests
sleep 5

echo "========================================================================"
echo "TEST 2: OPTION B (Direct WriteParticlesSlab_range_from_zslab)"
echo "========================================================================"

# Compile with Option B
echo "Compiling with Option B (USE_PARTICLE_OUTPUT_OPTION_B=1)..."
make clean
make CFLAGS="-DUSE_PARTICLE_OUTPUT_OPTION_B=1" -j$(nproc)

# Clean previous outputs
rm -rf rank_* "$OUTPUT_DIR_B"/*

# Run Option B with timing
echo ""
echo "Running Option B..."
echo "Command: mpirun -np $NUM_RANKS $MAIN_EXEC $N $PARAM_FILE"
echo ""

START_TIME_B=$(date +%s.%N)
/usr/bin/time -v mpirun -np $NUM_RANKS $MAIN_EXEC $N $PARAM_FILE 2>&1 | tee "${TIMING_DIR}/option_b_timing.log"
END_TIME_B=$(date +%s.%N)

# Move outputs
if [ -d "output" ]; then
    mv output/* "$OUTPUT_DIR_B/" 2>/dev/null || true
fi

# Extract timing information
ELAPSED_B=$(echo "$END_TIME_B - $START_TIME_B" | bc)
echo "Option B elapsed time: ${ELAPSED_B} seconds" | tee "${TIMING_DIR}/option_b_summary.txt"

# Count files
NUM_FILES_B=$(find "$OUTPUT_DIR_B" -name "ic_rank*" 2>/dev/null | wc -l)
echo "Option B: $NUM_FILES_B particle IC files created" | tee -a "${TIMING_DIR}/option_b_summary.txt"

# Get memory info from timing log if available
if grep -q "Maximum resident set size" "${TIMING_DIR}/option_b_timing.log"; then
    MAX_RSS_B=$(grep "Maximum resident set size" "${TIMING_DIR}/option_b_timing.log" | awk '{print $6}')
    echo "Option B: Max RSS per rank: ${MAX_RSS_B} KB" | tee -a "${TIMING_DIR}/option_b_summary.txt"
fi

echo ""
echo "Option B complete. Outputs in: $OUTPUT_DIR_B"
echo ""

echo "========================================================================"
echo "COMPARISON SUMMARY"
echo "========================================================================"

# Calculate speedup
if [ -n "$ELAPSED_A" ] && [ -n "$ELAPSED_B" ] && [ "$ELAPSED_A" != "0" ]; then
    SPEEDUP=$(echo "scale=2; $ELAPSED_A / $ELAPSED_B" | bc)
    TIME_SAVED=$(echo "scale=2; $ELAPSED_A - $ELAPSED_B" | bc)
    echo "Time comparison:"
    echo "  Option A: ${ELAPSED_A} seconds"
    echo "  Option B: ${ELAPSED_B} seconds"
    echo "  Speedup: ${SPEEDUP}× (Option B is ${SPEEDUP}× faster)"
    echo "  Time saved: ${TIME_SAVED} seconds"
    echo ""
fi

# Memory comparison
if [ -n "$MAX_RSS_A" ] && [ -n "$MAX_RSS_B" ]; then
    MEM_SAVED_KB=$(echo "$MAX_RSS_A - $MAX_RSS_B" | bc)
    MEM_SAVED_GB=$(echo "scale=2; $MEM_SAVED_KB / 1048576" | bc)
    MEM_RATIO=$(echo "scale=2; $MAX_RSS_A / $MAX_RSS_B" | bc)
    echo "Memory comparison (per rank):"
    echo "  Option A: ${MAX_RSS_A} KB (~$(echo "scale=2; $MAX_RSS_A / 1048576" | bc) GB)"
    echo "  Option B: ${MAX_RSS_B} KB (~$(echo "scale=2; $MAX_RSS_B / 1048576" | bc) GB)"
    echo "  Memory saved: ${MEM_SAVED_KB} KB (~${MEM_SAVED_GB} GB per rank)"
    echo "  Memory ratio: ${MEM_RATIO}× (Option A uses ${MEM_RATIO}× more memory)"
    echo ""
fi

# File count comparison
echo "Output file comparison:"
echo "  Option A: $NUM_FILES_A files"
echo "  Option B: $NUM_FILES_B files"
if [ "$NUM_FILES_A" -eq "$NUM_FILES_B" ]; then
    echo "  ✓ File counts match"
else
    echo "  ⚠ WARNING: File counts differ!"
fi
echo ""

echo "========================================================================"
echo "NEXT STEPS: Correctness Verification"
echo "========================================================================"
echo "1. Compare output files:"
echo "   ./compare_option_outputs.sh $OUTPUT_DIR_A $OUTPUT_DIR_B"
echo ""
echo "2. Verify a sample of particles match between options"
echo ""
echo "3. Check timing logs for detailed metrics:"
echo "   cat ${TIMING_DIR}/option_a_summary.txt"
echo "   cat ${TIMING_DIR}/option_b_summary.txt"
echo ""
echo "========================================================================"

