#!/bin/bash
# Script to run zeldovich-PLT with PLT enabled (rescaling ENABLED) for PLT test (ZD_qdensity=1 or 2)
# Serial version: single-threaded (no MPI, no OpenMP)
# Testing k_cutoff implementation (k_cutoff=1.0, CornerModes=0)
# This should be run separately (not in PBS) to generate reference output

set -e

# load modules
module load frameworks
module load fftw/3.3.10

ZELDOVICH_DIR="/home/helenshao/InitialConditions/zeldovich-PLT"
PARAM_FILE="/home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_N1024/CPD-225/param_N1024_CPD_225.par"
OUTPUT_DIR="${ZELDOVICH_DIR}/output_N1024_CPD_225"
RNG_LOG_DIR="/home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_N1024/CPD-225/rng_logs"

echo "=========================================="
echo "Running zeldovich-PLT with PLT ENABLED (rescaling ENABLED) for CPD-225 test (N=1024)"
echo "Using ZD_qdensity=1 or 2 (density + displacement or density only)"
echo "Testing k_cutoff implementation (k_cutoff=1.0, CornerModes=0)"
echo "CPD-225 test: comparing with 32 nodes, 1 rank per node, CPD=225"
echo "=========================================="
echo "Parameter file: $PARAM_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "RNG debug output: $RNG_LOG_DIR/zeldovich_rng_debug.log"
echo ""

# Check if zeldovich executable exists (must be built first)
ZELDOVICH_EXEC="${ZELDOVICH_DIR}/build/zeldovich"
if [ ! -x "$ZELDOVICH_EXEC" ]; then
    echo "ERROR: zeldovich executable not found: $ZELDOVICH_EXEC"
    echo "Please build zeldovich-PLT first (see hermitian_3d_matrix_production/docs/zeldovich-PLT-build-notes.md):"
    echo "  module load fftw/3.3.10"
    echo "  cd $ZELDOVICH_DIR"
    echo "  meson setup build    # if not yet configured"
    echo "  meson compile -C build"
    exit 1
fi

# Check if parameter file exists
if [ ! -f "$PARAM_FILE" ]; then
    echo "ERROR: Parameter file not found: $PARAM_FILE"
    exit 1
fi

# Verify parameter file has PLT enabled
if grep -q "ZD_qPLT = 0" "$PARAM_FILE"; then
    echo "WARNING: Parameter file has ZD_qPLT = 0 (PLT disabled)"
    echo "This script expects PLT to be enabled. Please check the parameter file."
    exit 1
fi

# Verify parameter file has PLT rescaling enabled
if grep -q "ZD_qPLTrescale = 0" "$PARAM_FILE"; then
    echo "WARNING: Parameter file has ZD_qPLTrescale = 0 (PLT rescaling disabled)"
    echo "This script expects PLT rescaling to be enabled. Please check the parameter file."
    exit 1
fi

# Note: This script accepts ZD_qdensity=1 or ZD_qdensity=2
# ZD_qdensity=1: outputs density + displacement (narray=2)
# ZD_qdensity=2: outputs density only (narray=1)

# Create output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$RNG_LOG_DIR"
echo "Output will be written to: $OUTPUT_DIR"
echo "RNG debug will be written to: $RNG_LOG_DIR"
echo ""

# Update parameter file to use the output directory
TEMP_PARAM=$(mktemp)
cp "$PARAM_FILE" "$TEMP_PARAM"
sed -i "s|InitialConditionsDirectory = \".*\"|InitialConditionsDirectory = \"$OUTPUT_DIR\"|" "$TEMP_PARAM"

echo "Running zeldovich-PLT (serial, single-threaded)..."

# Run zeldovich, capture stderr (where RNG-DEBUG goes) separately
"$ZELDOVICH_EXEC" "$TEMP_PARAM" \
    > "$OUTPUT_DIR/zeldovich_dens.log" \
    2> "$RNG_LOG_DIR/zeldovich_rng_debug.log"

# Clean up temp file
rm -f "$TEMP_PARAM"

# Extract RNG debug lines
echo "Extracting RNG-DEBUG lines from log..."
grep "\[RNG-DEBUG\]" "$RNG_LOG_DIR/zeldovich_rng_debug.log" > "$RNG_LOG_DIR/zeldovich_rng_debug_filtered.txt" || true
RNG_LINES=$(wc -l < "$RNG_LOG_DIR/zeldovich_rng_debug_filtered.txt")
echo "Found $RNG_LINES RNG-DEBUG lines"

# Extract REAL-FFT-DEBUG lines
echo "Extracting REAL-FFT-DEBUG lines from log..."
grep "\[REAL-FFT-DEBUG\]" "$RNG_LOG_DIR/zeldovich_rng_debug.log" > "$RNG_LOG_DIR/zeldovich_real_fft_debug.txt" || true
REAL_FFT_LINES=$(wc -l < "$RNG_LOG_DIR/zeldovich_real_fft_debug.txt")
echo "Found $REAL_FFT_LINES REAL-FFT-DEBUG lines"

# Check for output files
IC_COUNT=$(find "$OUTPUT_DIR" -name "ic_*" 2>/dev/null | wc -l)
echo ""
echo "=========================================="
echo "zeldovich-PLT (PLT enabled, rescaling ENABLED) completed"
echo "Using ZD_qdensity=1 or 2 (density + displacement or density only)"
echo "Testing k_cutoff implementation"
echo "CPD-225 test (N=1024): comparing with 32 nodes, 1 rank per node, CPD=225"
echo "=========================================="
echo "Output directory: $OUTPUT_DIR"
echo "Particle IC files: $IC_COUNT"
echo "RNG debug lines: $RNG_LINES"
echo "REAL-FFT-DEBUG lines: $REAL_FFT_LINES"
echo ""
echo "RNG debug output files:"
echo "  - $RNG_LOG_DIR/zeldovich_rng_debug_filtered.txt"
echo "  - $RNG_LOG_DIR/hermitian_rng_debug_filtered.txt (from hermitian run)"
echo ""
echo "REAL-FFT-DEBUG output files:"
echo "  - $RNG_LOG_DIR/hermitian_real_fft_debug.txt (from hermitian run)"
echo "  - $RNG_LOG_DIR/zeldovich_real_fft_debug.txt (from zeldovich run)"
echo ""
echo "You can now compare RNG outputs:"
echo "  - diff $RNG_LOG_DIR/hermitian_rng_debug_filtered.txt $RNG_LOG_DIR/zeldovich_rng_debug_filtered.txt"
echo "  - Or use the comparison script in workflow.txt"
echo ""
echo "You can now compare REAL-FFT-DEBUG outputs:"
echo "  - See workflow.txt for comparison script"
echo ""


