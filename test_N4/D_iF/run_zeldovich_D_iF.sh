#!/bin/bash
# Script to run zeldovich-PLT with D+i*F test (ZD_qdensity=2, but code explicitly stores D+i*F)
# Serial version: single-threaded (no MPI, no OpenMP)
# This should be run separately (not in PBS) to generate reference output

set -e

# load modules
module load frameworks
module load fftw/3.3.10

ZELDOVICH_DIR="/home/helenshao/InitialConditions/zeldovich-PLT"
PARAM_FILE="/home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_N4/D_iF/param_N4_D_iF.par"
OUTPUT_DIR="${ZELDOVICH_DIR}/output_N4_D_iF"
RNG_LOG_DIR="/home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_N4/D_iF/rng_logs"

echo "=========================================="
echo "Running zeldovich-PLT with D+i*F test"
echo "Using ZD_qdensity=2, but code explicitly stores D+i*F"
echo "=========================================="
echo "Parameter file: $PARAM_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "RNG debug output: $RNG_LOG_DIR/zeldovich_rng_debug.log"
echo ""

# Check if zeldovich executable exists
ZELDOVICH_EXEC="${ZELDOVICH_DIR}/build/zeldovich"
if [ ! -f "$ZELDOVICH_EXEC" ]; then
    echo "ERROR: zeldovich executable not found: $ZELDOVICH_EXEC"
    echo "Please build zeldovich-PLT first:"
    echo "  cd $ZELDOVICH_DIR && meson compile -C build"
    exit 1
fi

# Check if parameter file exists
if [ ! -f "$PARAM_FILE" ]; then
    echo "ERROR: Parameter file not found: $PARAM_FILE"
    exit 1
fi

# Verify parameter file has PLT disabled
if grep -q "ZD_qPLT = 1" "$PARAM_FILE"; then
    echo "WARNING: Parameter file has ZD_qPLT = 1 (PLT enabled)"
    echo "This script expects PLT to be disabled. Please check the parameter file."
    exit 1
fi

# Verify parameter file has ZD_qdensity=2
if ! grep -q "ZD_qdensity = 2" "$PARAM_FILE"; then
    echo "WARNING: Parameter file does not have ZD_qdensity = 2"
    echo "This script expects ZD_qdensity=2. Please check the parameter file."
    exit 1
fi

# Create output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$RNG_LOG_DIR"
echo "Output will be written to: $OUTPUT_DIR"
echo "RNG debug will be written to: $RNG_LOG_DIR"
echo ""

# Update parameter file to use the new output directory
# Create a temporary parameter file with updated output directory
TEMP_PARAM=$(mktemp)
cp "$PARAM_FILE" "$TEMP_PARAM"
sed -i "s|InitialConditionsDirectory = \".*\"|InitialConditionsDirectory = \"$OUTPUT_DIR\"|" "$TEMP_PARAM"

# Rebuild zeldovich to ensure debug code is included
echo "Rebuilding zeldovich-PLT to include REAL-FFT-DEBUG code..."
cd "$ZELDOVICH_DIR"

# Use manual rebuild method (same as run_zeldovich_real_fft_comparison.sh)
echo "Using manual rebuild method..."
module load fftw/3.3.10 2>/dev/null || true

# Touch source to force rebuild
touch src/zeldovich.cpp

# Rebuild object file and library
echo "Rebuilding object file and library with DEBUG_RNG_CONSISTENCY flag..."
# Workaround for meson Python environment issues: touch build files to prevent regeneration
touch build/build.ninja meson.build 2>/dev/null || true

# Modify build.ninja to add DEBUG_RNG_CONSISTENCY flag to cpp_args
if [ -f "build/build.ninja" ]; then
    echo "Modifying build.ninja to add DEBUG_RNG_CONSISTENCY flag..."
    # Check if flag is already present
    if ! grep -q "DEBUG_RNG_CONSISTENCY" build/build.ninja; then
        # Add -DDEBUG_RNG_CONSISTENCY to cpp_args (most common variable name in meson builds)
        # Try different patterns that meson might use
        sed -i 's/^cpp_args = \(.*\)$/cpp_args = \1 -DDEBUG_RNG_CONSISTENCY/' build/build.ninja 2>/dev/null || true
        sed -i 's/^cpp_args =$/cpp_args = -DDEBUG_RNG_CONSISTENCY/' build/build.ninja 2>/dev/null || true
        # Also try adding to the compile command directly for the zeldovich.cpp object
        sed -i '/zeldovich\.cpp\.o:/s|\( \$in\)|\1 -DDEBUG_RNG_CONSISTENCY|' build/build.ninja 2>/dev/null || true
        # Alternative: add to the rule that compiles .cpp files
        sed -i 's|\(command = .*icpx.*\)|\1 -DDEBUG_RNG_CONSISTENCY|' build/build.ninja 2>/dev/null || true
        sed -i 's|\(command = .*g++.*\)|\1 -DDEBUG_RNG_CONSISTENCY|' build/build.ninja 2>/dev/null || true
        sed -i 's|\(command = .*clang++.*\)|\1 -DDEBUG_RNG_CONSISTENCY|' build/build.ninja 2>/dev/null || true
    else
        echo "DEBUG_RNG_CONSISTENCY flag already present in build.ninja"
    fi
fi

# Rebuild with modified flags
echo "Rebuilding with DEBUG_RNG_CONSISTENCY flag..."
if ninja -C build zeldovich.p/src_zeldovich.cpp.o libzeldovich.so 2>&1 | grep -v "symbolextractor\|ModuleNotFoundError\|mesonbuild\|FAILED.*build.ninja"; then
    OBJ_LIB_STATUS=0
else
    OBJ_LIB_STATUS=${PIPESTATUS[0]}
fi

if [ $OBJ_LIB_STATUS -ne 0 ]; then
    echo "WARNING: Build had errors, but continuing..."
fi

# Verify object file and library exist
if [ ! -f "build/zeldovich.p/src_zeldovich.cpp.o" ] || [ ! -f "build/libzeldovich.so" ]; then
    echo "ERROR: Object file or library not found!"
    exit 1
fi

# Manually link executable
echo "Linking executable manually..."
icpx -fopenmp -DFMT_SHARED build/zeldovich.p/src_zeldovich.cpp.o \
     -Lbuild -Lbuild/subprojects/fmt-11.2.0 \
     -lzeldovich -lfmt -lfftw3 -lfftw3_omp \
     -Wl,-rpath,build:build/subprojects/fmt-11.2.0 \
     -o build/zeldovich 2>&1 || {
    echo "WARNING: Manual link had errors, but continuing..."
}

if [ ! -f "./build/zeldovich" ]; then
    echo "ERROR: zeldovich executable not found after compilation!"
    exit 1
fi
echo "Build successful"

echo ""
echo "Running zeldovich-PLT (serial, single-threaded)..."

# Set LD_LIBRARY_PATH for shared libraries
export LD_LIBRARY_PATH="${ZELDOVICH_DIR}/build:${ZELDOVICH_DIR}/build/subprojects/fmt-11.2.0:${LD_LIBRARY_PATH}"

# Run zeldovich, capture stderr (where RNG-DEBUG and REAL-FFT-DEBUG go) separately
"$ZELDOVICH_EXEC" "$TEMP_PARAM" \
    > "$OUTPUT_DIR/zeldovich_D_iF.log" \
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
echo "zeldovich-PLT (D+i*F test) completed"
echo "Using ZD_qdensity=2, but code explicitly stores D+i*F"
echo "=========================================="
echo "Output directory: $OUTPUT_DIR"
echo "Particle IC files: $IC_COUNT"
echo "RNG debug lines: $RNG_LINES"
echo "REAL-FFT-DEBUG lines: $REAL_FFT_LINES"
echo ""
echo "REAL-FFT-DEBUG output files:"
echo "  - $RNG_LOG_DIR/hermitian_real_fft_debug.txt (from hermitian run)"
echo "  - $RNG_LOG_DIR/zeldovich_real_fft_debug.txt (from zeldovich run)"
echo ""
echo "You can now compare REAL-FFT-DEBUG outputs:"
echo "  - See workflow.txt for comparison script"
echo ""

