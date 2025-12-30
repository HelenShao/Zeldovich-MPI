#!/bin/bash
# Script to run zeldovich-PLT with PLT disabled for RNG comparison
# Serial version: single-threaded (no MPI, no OpenMP)
# Testing k_cutoff implementation (k_cutoff=1.0, CornerModes=0)
# This should be run separately (not in PBS) to generate reference output

set -e

# load modules
module load frameworks
module load fftw/3.3.10

ZELDOVICH_DIR="/home/helenshao/InitialConditions/zeldovich-PLT"
PARAM_FILE="/home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_rng/k2_cutoff/N64/param_N64_rng.par"
OUTPUT_DIR="${ZELDOVICH_DIR}/output_rng_k2cutoff_N64"
RNG_LOG_DIR="/home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_rng/k2_cutoff/N64/rng_logs"

echo "=========================================="
echo "Running zeldovich-PLT with PLT DISABLED (RNG comparison)"
echo "Testing k_cutoff implementation (k_cutoff=1.0, CornerModes=0)"
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
    echo "  cd $ZELDOVICH_DIR && make"
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
echo "Rebuilding zeldovich-PLT to include RNG debug code..."
cd "$ZELDOVICH_DIR"

# Load required modules for build
module load fftw/3.3.10 2>/dev/null || true

# Use ninja directly (recommended method per BUILD_AND_RUN.md)
# This works even if meson has Python interpreter issues
if [ -f "build/build.ninja" ] && command -v ninja >/dev/null 2>&1; then
    echo "Rebuilding with ninja (full rebuild to ensure debug code is included)..."
    # Force full rebuild to ensure executable is updated
    # Touch the source file to force rebuild, then rebuild everything
    touch src/zeldovich.cpp
    
    # Rebuild object file, library, and executable
    # Note: Symbol extraction may fail (non-critical per BUILD_AND_RUN.md), but executable is still built
    echo "Rebuilding object file and library (per BUILD_AND_RUN.md troubleshooting)..."
    # Per BUILD_AND_RUN.md, build object file and library separately to avoid symbol extraction error
    if ninja -C build zeldovich.p/src_zeldovich.cpp.o libzeldovich.so 2>&1 | grep -v "symbolextractor\|ModuleNotFoundError"; then
        OBJ_LIB_STATUS=0
    else
        OBJ_LIB_STATUS=${PIPESTATUS[0]}
    fi
    
    if [ $OBJ_LIB_STATUS -ne 0 ]; then
        echo "WARNING: Object file/library build had errors, but continuing..."
    fi
    
    echo "Linking executable manually (per BUILD_AND_RUN.md troubleshooting)..."
    # Per BUILD_AND_RUN.md, manually link the executable to avoid symbol extraction step
    BUILD_STATUS=0
    
    # Manual link per BUILD_AND_RUN.md troubleshooting section
    if [ -f "build/zeldovich.p/src_zeldovich.cpp.o" ] && [ -f "build/libzeldovich.so" ]; then
        echo "Linking executable with icpx..."
        if icpx -fopenmp -DFMT_SHARED build/zeldovich.p/src_zeldovich.cpp.o \
             -Lbuild -Lbuild/subprojects/fmt-11.2.0 \
             -lzeldovich -lfmt -lfftw3 -lfftw3_omp \
             -Wl,-rpath,build:build/subprojects/fmt-11.2.0 \
             -o build/zeldovich 2>&1; then
            BUILD_STATUS=0
        else
            BUILD_STATUS=$?
            echo "WARNING: Manual link had errors (exit status: $BUILD_STATUS)"
        fi
    else
        echo "ERROR: Object file or library not found. Cannot link."
        echo "  Object file: build/zeldovich.p/src_zeldovich.cpp.o - $([ -f build/zeldovich.p/src_zeldovich.cpp.o ] && echo 'EXISTS' || echo 'NOT FOUND')"
        echo "  Library: build/libzeldovich.so - $([ -f build/libzeldovich.so ] && echo 'EXISTS' || echo 'NOT FOUND')"
        echo "Please rebuild manually:"
        echo "  cd $ZELDOVICH_DIR && module load fftw/3.3.10"
        echo "  ninja -C build zeldovich.p/src_zeldovich.cpp.o libzeldovich.so"
        exit 1
    fi
    
    # Check if executable exists
    if [ ! -f "$ZELDOVICH_EXEC" ]; then
        echo "ERROR: Executable not found after manual link. Build failed."
        echo "Please rebuild manually:"
        echo "  cd $ZELDOVICH_DIR && module load fftw/3.3.10"
        echo "  ninja -C build zeldovich.p/src_zeldovich.cpp.o libzeldovich.so"
        echo "  icpx -fopenmp -DFMT_SHARED build/zeldovich.p/src_zeldovich.cpp.o \\"
        echo "       -Lbuild -Lbuild/subprojects/fmt-11.2.0 \\"
        echo "       -lzeldovich -lfmt -lfftw3 -lfftw3_omp \\"
        echo "       -Wl,-rpath,build:build/subprojects/fmt-11.2.0 \\"
        echo "       -o build/zeldovich"
        exit 1
    elif [ $BUILD_STATUS -eq 0 ]; then
        echo "Build completed successfully"
    else
        echo "WARNING: Build had errors, but executable exists"
    fi
    
    # Verify executable was updated
    if [ -f "src/zeldovich.cpp" ] && [ -f "$ZELDOVICH_EXEC" ]; then
        if [ "src/zeldovich.cpp" -nt "$ZELDOVICH_EXEC" ]; then
            echo "WARNING: Executable is still older than source after rebuild"
            echo "Trying one more rebuild..."
            # Force rebuild by removing object file
            rm -f build/zeldovich.p/src_zeldovich.cpp.o
            ninja -C build zeldovich 2>&1 | grep -v "symbolextractor\|ModuleNotFoundError" || true
            
            # Check again
            if [ "src/zeldovich.cpp" -nt "$ZELDOVICH_EXEC" ]; then
                echo "ERROR: Executable is still older than source after rebuild!"
                echo "Please rebuild manually:"
                echo "  cd $ZELDOVICH_DIR && module load fftw/3.3.10 && ninja -C build"
                exit 1
            fi
        fi
        echo "Rebuild successful - executable is up to date"
    fi
else
    echo "WARNING: ninja not available or build.ninja not found"
    echo "Attempting to use meson..."
    if command -v meson >/dev/null 2>&1; then
        meson compile -C build 2>&1 || {
            echo "WARNING: meson rebuild failed, but executable may still be valid"
        }
    else
        echo "WARNING: Neither ninja nor meson available. Using existing executable."
    fi
fi

# Check if executable exists (even if build had warnings)
if [ ! -f "$ZELDOVICH_EXEC" ]; then
    echo "ERROR: zeldovich executable not found: $ZELDOVICH_EXEC"
    echo "Please build zeldovich-PLT first:"
    echo "  cd $ZELDOVICH_DIR"
    echo "  module load fftw/3.3.10"
    echo "  ninja -C build"
    exit 1
fi

# Verify executable is recent (check if source is newer)
if [ -f "src/zeldovich.cpp" ] && [ "src/zeldovich.cpp" -nt "$ZELDOVICH_EXEC" ]; then
    echo "WARNING: zeldovich.cpp is newer than executable. Rebuild may have failed."
    echo "You may need to rebuild manually:"
    echo "  cd $ZELDOVICH_DIR && module load fftw/3.3.10 && ninja -C build"
else
    echo "Build check complete - executable is up to date"
fi

echo "Running zeldovich-PLT (serial, single-threaded)..."

# Run zeldovich, capture stderr (where RNG-DEBUG goes) separately
"$ZELDOVICH_EXEC" "$TEMP_PARAM" \
    > "$OUTPUT_DIR/zeldovich_rng.log" \
    2> "$RNG_LOG_DIR/zeldovich_rng_debug.log"

# Clean up temp file
rm -f "$TEMP_PARAM"

# Extract RNG debug lines
echo "Extracting RNG-DEBUG lines from log..."
grep "\[RNG-DEBUG\]" "$RNG_LOG_DIR/zeldovich_rng_debug.log" > "$RNG_LOG_DIR/zeldovich_rng_debug_filtered.txt" || true
RNG_LINES=$(wc -l < "$RNG_LOG_DIR/zeldovich_rng_debug_filtered.txt")
echo "Found $RNG_LINES RNG-DEBUG lines"

# Check for output files
IC_COUNT=$(find "$OUTPUT_DIR" -name "ic_*" 2>/dev/null | wc -l)
echo ""
echo "=========================================="
echo "zeldovich-PLT (RNG comparison) completed"
echo "Testing k_cutoff implementation"
echo "=========================================="
echo "Output directory: $OUTPUT_DIR"
echo "Particle IC files: $IC_COUNT"
echo "RNG debug lines: $RNG_LINES"
echo ""
echo "RNG debug output files:"
echo "  - $RNG_LOG_DIR/zeldovich_rng_debug_filtered.txt"
echo "  - $RNG_LOG_DIR/hermitian_rng_debug_filtered.txt (from hermitian run)"
echo ""
echo "You can now compare RNG outputs:"
echo "  - diff $RNG_LOG_DIR/hermitian_rng_debug_filtered.txt $RNG_LOG_DIR/zeldovich_rng_debug_filtered.txt"
echo "  - Or use the comparison script in workflow.txt"
echo ""
echo "Expected: All modes with k² >= 64 (k2_cutoff) should have D=0 in both codes"

