#!/bin/bash
# Script to run hermitian_3d_matrix with matrix dump enabled for trace comparison
# Serial version: 1 rank, 1 thread (N=8)

# Don't use set -e, we'll check exit codes manually

# Load modules
echo "Loading frameworks module..."
module load frameworks
echo "Loading FFTW module..."
module load fftw/3.3.10

# Test parameters
N=8
PARAM_FILE="param_N8_trace.par"
OUTPUT_DIR="hermitian_output"

echo "=========================================="
echo "Running hermitian_3d_matrix with matrix trace (N=8)"
echo "=========================================="
echo "Parameter file: $PARAM_FILE"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Check if executable exists (use absolute path)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
EXEC="$PROJECT_ROOT/hermitian_3d_matrix"

# Rebuild executables with correct flags
cd "$PROJECT_ROOT"
make clean
make CFLAGS="-DUSE_DOUBLE_PRECISION -DPRINT_DETAILED_SLICES=1 -DPRINT_Z_SLABS=0 -DSKIP_FILE_WRITE=0 -DDUMP_MATRIX_BEFORE_FFT=1 -DDUMP_MATRIX_AFTER_FFT=1"
make reassembly CFLAGS="-DUSE_DOUBLE_PRECISION -DPRINT_DETAILED_SLICES=1 -DPRINT_Z_SLABS=0 -DSKIP_FILE_WRITE=0 -DDUMP_MATRIX_BEFORE_FFT=1 -DDUMP_MATRIX_AFTER_FFT=1"
BUILD_EXIT_CODE=$?

if [ $BUILD_EXIT_CODE -ne 0 ]; then
    echo "ERROR: Build failed with exit code $BUILD_EXIT_CODE"
    exit 1
fi

# Return to script directory
cd "$SCRIPT_DIR"

# Verify executable exists
if [ ! -f "$EXEC" ]; then
    echo "ERROR: Executable not found after build: $EXEC"
    exit 1
fi
echo "Using executable: $EXEC"

# Check if parameter file exists
if [ ! -f "$PARAM_FILE" ]; then
    echo "ERROR: Parameter file not found: $PARAM_FILE"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run hermitian_3d_matrix
echo "Running hermitian_3d_matrix..."
echo "Command: mpiexec -n 1 $EXEC $N $PARAM_FILE"
echo "Working directory: $(pwd)"
echo "Executable exists: $([ -f "$EXEC" ] && echo "YES" || echo "NO")"
echo "Parameter file exists: $([ -f "$PARAM_FILE" ] && echo "YES" || echo "NO")"
echo ""

# Use absolute path for executable (like PBS script does)
EXECUTABLE_PATH="$PROJECT_ROOT/hermitian_3d_matrix"
PARAM_FILE_ABS="$(cd "$SCRIPT_DIR" && pwd)/$PARAM_FILE"

# Change to output directory so matrix dump files are written there
# (The program writes to current working directory)
cd "$SCRIPT_DIR/$OUTPUT_DIR"

# Run with explicit output redirection
# Use absolute paths for both executable and parameter file
echo "Changed to directory: $(pwd)"
echo "Running: mpiexec -n 1 $EXECUTABLE_PATH $N $PARAM_FILE_ABS"
mpiexec -n 1 "$EXECUTABLE_PATH" $N "$PARAM_FILE_ABS" > "hermitian_N8_trace.log" 2>&1
EXIT_CODE=$?

# Change back to script directory
cd "$SCRIPT_DIR"

# Also display to terminal
cat "$OUTPUT_DIR/hermitian_N8_trace.log"

if [ $EXIT_CODE -ne 0 ]; then
    echo ""
    echo "ERROR: hermitian_3d_matrix exited with code $EXIT_CODE"
    if [ $EXIT_CODE -eq 127 ]; then
        echo ""
        echo "Exit code 127 means 'command not found'. This usually indicates:"
        echo "  1. mpiexec cannot find or execute the program"
        echo "  2. You may need to run this in a PBS job or interactive compute node"
        echo "  3. The MPI environment may not be properly configured for interactive use"
        echo ""
        echo "Try running this script inside a PBS job:"
        echo "  Create a PBS script (e.g., test_IC_trace.pbs) and submit it with: qsub test_IC_trace.pbs"
        echo "  OR request an interactive session:"
        echo "    qsub -I -l select=1:ncpus=1:mpiprocs=1 -l walltime=01:00:00"
    fi
    echo "Check log file: $OUTPUT_DIR/hermitian_N8_trace.log"
    if [ -f "$OUTPUT_DIR/hermitian_N8_trace.log" ]; then
        echo "Last 20 lines of log:"
        tail -20 "$OUTPUT_DIR/hermitian_N8_trace.log"
    fi
    exit $EXIT_CODE
fi

# Matrix dump files should already be in OUTPUT_DIR (since we cd'd there)
# Just verify they exist
echo ""
echo "Checking for matrix dump files in $OUTPUT_DIR/..."
ls -lh "$OUTPUT_DIR/matrix_*.txt" 2>/dev/null || echo "Warning: Matrix dump files not found in $OUTPUT_DIR/"

# Check for .bin files (if SKIP_FILE_WRITE=0 was used)
echo ""
echo "Checking for .bin files in $OUTPUT_DIR/rank_0/..."
BIN_COUNT=$(find "$OUTPUT_DIR/rank_0" -name "*.bin" -type f 2>/dev/null | wc -l)
if [ "$BIN_COUNT" -gt 0 ]; then
    echo "Found $BIN_COUNT .bin files in $OUTPUT_DIR/rank_0/"
    ls -lh "$OUTPUT_DIR/rank_0"/*.bin 2>/dev/null | head -5
    echo ""
    echo "NOTE: .bin files are saved to: $OUTPUT_DIR/rank_0/i*_slab_N8.bin"
    echo "      (relative to the working directory when the executable runs)"
else
    echo "No .bin files found (executable was compiled with SKIP_FILE_WRITE=1)"
    echo ""
    echo "To generate .bin files, recompile with:"
    echo "  make CFLAGS=\"-DUSE_DOUBLE_PRECISION -DSKIP_FILE_WRITE=0\""
    echo ""
    echo "Then .bin files will be written to: $OUTPUT_DIR/rank_0/i*_slab_N8.bin"
fi

echo ""
echo "=========================================="
echo "hermitian_3d_matrix trace completed"
echo "Log saved to: $OUTPUT_DIR/hermitian_N8_trace.log"
echo "Matrix dumps saved to: $OUTPUT_DIR/"
if [ "$BIN_COUNT" -gt 0 ]; then
    echo ".bin files saved to: $OUTPUT_DIR/rank_0/"
fi
echo "=========================================="

