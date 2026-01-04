#!/bin/bash

# Load modules
echo "Loading frameworks module..."
module load frameworks
echo "Loading FFTW module..."
module load fftw/3.3.10

# Job parameters - SERIAL: 1 rank, 1 thread
NNODES=1
NRANKS_PER_NODE=1
NTHREADS_PER_RANK=1
NTOTRANKS=1

# OpenMP settings
export OMP_NUM_THREADS=$NTHREADS_PER_RANK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

# Test parameters
N=8
PARAM_FILE="/home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_IC/param_N8_trace.par"
BIN_OUTPUT_DIR="test_IC/hermitian_output"  # .bin files are in rank_0/ subdirectory
PARTICLE_OUTPUT_DIR="test_IC/particle_ics"

# Create output directories
mkdir -p test_IC
mkdir -p "$BIN_OUTPUT_DIR"
mkdir -p "$PARTICLE_OUTPUT_DIR"

# ====================================================================================
# STEP: Reassemble and write particle ICs
# ====================================================================================
echo ""
echo "=========================================="
echo "Reassembling and writing particle ICs (N=8)"
echo "=========================================="

cd /home/helenshao/InitialConditions/hermitian_3d_matrix_production

I_START=0
I_END=$N

# Check if .bin files exist
BIN_COUNT=$(find "$BIN_OUTPUT_DIR/rank_0" -name "*.bin" 2>/dev/null | wc -l)
if [ "$BIN_COUNT" -eq 0 ]; then
    echo "ERROR: No .bin files found in $BIN_OUTPUT_DIR/rank_0/"
    echo "You need to run hermitian_3d_matrix with SKIP_FILE_WRITE=0 to generate .bin files"
    echo "Recompile with: make CFLAGS=\"-DUSE_DOUBLE_PRECISION -DSKIP_FILE_WRITE=0\""
    exit 1
fi

echo "Found $BIN_COUNT .bin files"
echo "Running write_particles_from_reassembled_mpi..."
./write_particles_from_reassembled_mpi \
    "$BIN_OUTPUT_DIR" \
    $N \
    $NTOTRANKS \
    "$PARAM_FILE" \
    $I_START \
    $I_END \
    2>&1 | tee test_IC/N8_trace_reassembly.log

# Move output files
# Note: The parameter file specifies InitialConditionsDirectory = "./output_trace"
# so files are written to output_trace/, not output/
if [ -d "output_trace" ]; then
    mkdir -p "$PARTICLE_OUTPUT_DIR"
    mv output_trace/ic_* "$PARTICLE_OUTPUT_DIR/" 2>/dev/null || true
    mv output_trace/dens* "$PARTICLE_OUTPUT_DIR/" 2>/dev/null || true
    echo "Moved IC files from output_trace/ to $PARTICLE_OUTPUT_DIR/"
elif [ -d "output" ]; then
    # Fallback: also check output/ directory in case parameter file uses different setting
    mkdir -p "$PARTICLE_OUTPUT_DIR"
    mv output/ic_* "$PARTICLE_OUTPUT_DIR/" 2>/dev/null || true
    mv output/dens* "$PARTICLE_OUTPUT_DIR/" 2>/dev/null || true
    echo "Moved IC files from output/ to $PARTICLE_OUTPUT_DIR/"
fi

IC_COUNT=$(find "$PARTICLE_OUTPUT_DIR" -name "ic_*" 2>/dev/null | wc -l)
echo "Generated $IC_COUNT particle IC files"

