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

# Test parameters - N=1024 PLT test
N=1024
PARAM_FILE="/home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_rng/k2_cutoff/N1024_PLT/param_N1024_PLT.par"
BIN_OUTPUT_DIR="test_rng/k2_cutoff/N1024_PLT/bin_files"
PARTICLE_OUTPUT_DIR="test_rng/k2_cutoff/N1024_PLT/particle_ics"
RNG_LOG_DIR="test_rng/k2_cutoff/N1024_PLT/rng_logs"

# Create output directories
mkdir -p test_rng/k2_cutoff/N1024_PLT
mkdir -p "$BIN_OUTPUT_DIR"
mkdir -p "$PARTICLE_OUTPUT_DIR"
mkdir -p "$RNG_LOG_DIR"

# ====================================================================================
# STEP 3: Reassemble and write particle ICs
# ====================================================================================
echo ""
echo "=========================================="
echo "STEP 3: Reassembling and writing particle ICs"
echo "=========================================="

cd /home/helenshao/InitialConditions/hermitian_3d_matrix_production

I_START=0
I_END=$N

echo "Running write_particles_from_reassembled_mpi..."
./write_particles_from_reassembled_mpi \
    "$BIN_OUTPUT_DIR" \
    $N \
    $NTOTRANKS \
    "$PARAM_FILE" \
    $I_START \
    $I_END \
    2>&1 | tee test_rng/k2_cutoff/N1024_PLT/N1024_PLT_reassembly.log

# Move output files
if [ -d "output" ]; then
    mkdir -p "$PARTICLE_OUTPUT_DIR"
    mv output/ic_* "$PARTICLE_OUTPUT_DIR/" 2>/dev/null || true
    mv output/dens* "$PARTICLE_OUTPUT_DIR/" 2>/dev/null || true
fi

IC_COUNT=$(find "$PARTICLE_OUTPUT_DIR" -name "ic_*" 2>/dev/null | wc -l)
echo "Generated $IC_COUNT particle IC files"
