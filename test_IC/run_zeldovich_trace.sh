#!/bin/bash
# Script to run zeldovich-PLT with matrix dump enabled for trace comparison
# Serial version: single-threaded (no MPI, no OpenMP)
# Testing matrix values before and after FFT (N=8)

# Don't use set -e, we'll check exit codes manually

# Load modules
module load frameworks
module load fftw/3.3.10

ZELDOVICH_DIR="/home/helenshao/InitialConditions/zeldovich-PLT"
PARAM_FILE="/home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_IC/param_N8_trace.par"
OUTPUT_DIR="${ZELDOVICH_DIR}/output_trace"
TRACE_LOG_DIR="/home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_IC/zeldovich_output"

echo "=========================================="
echo "Running zeldovich-PLT with matrix trace (N=8)"
echo "=========================================="
echo "Parameter file: $PARAM_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Trace log directory: $TRACE_LOG_DIR"
echo ""

# Check if zeldovich executable exists
ZELDOVICH_EXEC="${ZELDOVICH_DIR}/build/zeldovich"

# Rebuild zeldovich to ensure dump code is included
echo "Rebuilding zeldovich-PLT to include matrix dump code..."
cd "$ZELDOVICH_DIR"

# Load required modules for build
module load fftw/3.3.10 2>/dev/null || true

# Check if meson was configured with DUMP_MATRIX flags
NEEDS_RECONFIGURE=0
if [ -f "build/build.ninja" ]; then
    if ! grep -q "DUMP_MATRIX" build/build.ninja 2>/dev/null; then
        echo "WARNING: build.ninja does not contain DUMP_MATRIX flags"
        echo "Reconfiguring meson with DUMP_MATRIX flags..."
        
        # Try meson command first
        if meson setup build --buildtype=release -Dcpp_args='-DDUMP_MATRIX_BEFORE_FFT=1 -DDUMP_MATRIX_AFTER_FFT=1' 2>&1; then
            NEEDS_RECONFIGURE=1
        else
            # If meson command fails (e.g., bad Python interpreter), try using python3 directly
            echo "meson command failed, trying with python3..."
            if command -v python3 >/dev/null 2>&1; then
                if python3 -m mesonbuild.mesonmain setup build --buildtype=release -Dcpp_args='-DDUMP_MATRIX_BEFORE_FFT=1 -DDUMP_MATRIX_AFTER_FFT=1' 2>&1; then
                    NEEDS_RECONFIGURE=1
                else
                    echo "ERROR: meson setup failed even with python3"
                    echo "Please reconfigure manually using one of these methods:"
                    echo ""
                    echo "Method 1: Use python3 directly:"
                    echo "  cd $ZELDOVICH_DIR"
                    echo "  module load fftw/3.3.10"
                    echo "  python3 -m mesonbuild.mesonmain setup build --buildtype=release -Dcpp_args='-DDUMP_MATRIX_BEFORE_FFT=1 -DDUMP_MATRIX_AFTER_FFT=1'"
                    echo ""
                    echo "Method 2: Manually edit build/build.ninja to add -DDUMP_MATRIX_BEFORE_FFT=1 -DDUMP_MATRIX_AFTER_FFT=1 to cpp_args"
                    echo ""
                    echo "Continuing with existing build (matrix dumps may not work)..."
                fi
            else
                echo "ERROR: meson setup failed and python3 not available"
                echo "Please reconfigure manually:"
                echo "  cd $ZELDOVICH_DIR"
                echo "  module load fftw/3.3.10"
                echo "  python3 -m mesonbuild.mesonmain setup build --buildtype=release -Dcpp_args='-DDUMP_MATRIX_BEFORE_FFT=1 -DDUMP_MATRIX_AFTER_FFT=1'"
                echo ""
                echo "Continuing with existing build (matrix dumps may not work)..."
            fi
        fi
    fi
else
    echo "build.ninja not found. Setting up meson with DUMP_MATRIX flags..."
    if meson setup build --buildtype=release -Dcpp_args='-DDUMP_MATRIX_BEFORE_FFT=1 -DDUMP_MATRIX_AFTER_FFT=1' 2>&1; then
        NEEDS_RECONFIGURE=1
    else
        # Try with python3
        if command -v python3 >/dev/null 2>&1; then
            echo "meson command failed, trying with python3..."
            if python3 -m mesonbuild.mesonmain setup build --buildtype=release -Dcpp_args='-DDUMP_MATRIX_BEFORE_FFT=1 -DDUMP_MATRIX_AFTER_FFT=1' 2>&1; then
                NEEDS_RECONFIGURE=1
            else
                echo "ERROR: meson setup failed"
                exit 1
            fi
        else
            echo "ERROR: meson setup failed and python3 not available"
            exit 1
        fi
    fi
fi

# Use ninja directly (recommended method per BUILD_AND_RUN.md)
# This works even if meson has Python interpreter issues
if [ -f "build/build.ninja" ] && command -v ninja >/dev/null 2>&1; then
    echo "Rebuilding with ninja (full rebuild to ensure dump code is included)..."
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
        echo "  meson setup build --buildtype=release -Dcpp_args='-DDUMP_MATRIX_BEFORE_FFT=1 -DDUMP_MATRIX_AFTER_FFT=1'"
        echo "  ninja -C build zeldovich.p/src_zeldovich.cpp.o libzeldovich.so"
        exit 1
    fi
    
    # Check if executable exists
    if [ ! -f "$ZELDOVICH_EXEC" ]; then
        echo "ERROR: Executable not found after manual link. Build failed."
        echo "Please rebuild manually:"
        echo "  cd $ZELDOVICH_DIR && module load fftw/3.3.10"
        echo "  meson setup build --buildtype=release -Dcpp_args='-DDUMP_MATRIX_BEFORE_FFT=1 -DDUMP_MATRIX_AFTER_FFT=1'"
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
                echo "  cd $ZELDOVICH_DIR && module load fftw/3.3.10"
                echo "  meson setup build --buildtype=release -Dcpp_args='-DDUMP_MATRIX_BEFORE_FFT=1 -DDUMP_MATRIX_AFTER_FFT=1'"
                echo "  ninja -C build"
                exit 1
            fi
        fi
        echo "Rebuild successful - executable is up to date"
    fi
else
    echo "WARNING: ninja not available or build.ninja not found"
    echo "Attempting to use meson..."
    if command -v meson >/dev/null 2>&1; then
        echo "Reconfiguring meson with DUMP_MATRIX flags..."
        meson setup build --buildtype=release -Dcpp_args='-DDUMP_MATRIX_BEFORE_FFT=1 -DDUMP_MATRIX_AFTER_FFT=1' 2>&1 || {
            echo "WARNING: meson setup failed, trying compile anyway..."
        }
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
    echo "  meson setup build --buildtype=release -Dcpp_args='-DDUMP_MATRIX_BEFORE_FFT=1 -DDUMP_MATRIX_AFTER_FFT=1'"
    echo "  ninja -C build"
    exit 1
fi

# Verify executable is recent (check if source is newer)
if [ -f "src/zeldovich.cpp" ] && [ "src/zeldovich.cpp" -nt "$ZELDOVICH_EXEC" ]; then
    echo "WARNING: zeldovich.cpp is newer than executable. Rebuild may have failed."
    echo "You may need to rebuild manually:"
    echo "  cd $ZELDOVICH_DIR && module load fftw/3.3.10"
    echo "  meson setup build --buildtype=release -Dcpp_args='-DDUMP_MATRIX_BEFORE_FFT=1 -DDUMP_MATRIX_AFTER_FFT=1'"
    echo "  ninja -C build"
else
    echo "Build check complete - executable is up to date"
fi

# Check if executable was compiled with dump flags (basic check)
if ! strings "$ZELDOVICH_EXEC" | grep -q "DUMP_MATRIX"; then
    echo ""
    echo "WARNING: zeldovich executable may not have been compiled with DUMP_MATRIX flags"
    echo "To ensure matrix dumps work, reconfigure and rebuild:"
    echo "  cd $ZELDOVICH_DIR"
    echo "  module load fftw/3.3.10"
    echo "  meson setup build --buildtype=release -Dcpp_args='-DDUMP_MATRIX_BEFORE_FFT=1 -DDUMP_MATRIX_AFTER_FFT=1'"
    echo "  ninja -C build"
    echo ""
    echo "Continuing anyway, but matrix dumps may not be generated..."
fi

# Check if parameter file exists
if [ ! -f "$PARAM_FILE" ]; then
    echo "ERROR: Parameter file not found: $PARAM_FILE"
    exit 1
fi

# Create output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$TRACE_LOG_DIR"

# Run zeldovich-PLT
# Note: zeldovich-PLT needs to be compiled with -DDUMP_MATRIX_BEFORE_FFT=1 -DDUMP_MATRIX_AFTER_FFT=1
# Matrix dump files will be written to the current working directory
echo "Running zeldovich-PLT..."
cd "$ZELDOVICH_DIR"
"$ZELDOVICH_EXEC" "$PARAM_FILE" 2>&1 | tee "$TRACE_LOG_DIR/zeldovich_N8_trace.log"

# Move matrix dump files to trace log directory if they exist
echo ""
echo "Checking for matrix dump files..."
if [ -f "matrix_before_fft.txt" ]; then
    mv "matrix_before_fft.txt" "$TRACE_LOG_DIR/" 2>/dev/null || true
    echo "  Found matrix_before_fft.txt"
fi
if [ -f "matrix_after_fft.txt" ]; then
    mv "matrix_after_fft.txt" "$TRACE_LOG_DIR/" 2>/dev/null || true
    echo "  Found matrix_after_fft.txt"
fi

echo ""
echo "=========================================="
echo "zeldovich-PLT trace completed"
echo "Log saved to: $TRACE_LOG_DIR/zeldovich_N8_trace.log"
if [ -f "$TRACE_LOG_DIR/matrix_before_fft.txt" ] || [ -f "$TRACE_LOG_DIR/matrix_after_fft.txt" ]; then
    echo "Matrix dumps saved to: $TRACE_LOG_DIR/"
else
    echo "WARNING: Matrix dump files not found. Make sure zeldovich-PLT was compiled with:"
    echo "  -DDUMP_MATRIX_BEFORE_FFT=1 -DDUMP_MATRIX_AFTER_FFT=1"
fi
echo "=========================================="

