#!/bin/bash
#SBATCH --job-name=dd_test_${TEST_NAME}
#SBATCH --nodes=${NODES}
#SBATCH --ntasks=${RANKS}
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-${TEST_NAME}.out
#SBATCH --error=slurm-${TEST_NAME}.err

echo "========================================"
echo "Domain Decomposition Test: ${TEST_NAME}"
echo "N=${N}, Ranks=${RANKS}"
echo "========================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Nodes: $SLURM_NODELIST"
echo "Start time: $(date)"
echo ""

# Navigate to test directory
cd /scratch/gpfs/hshao/C_Bible/FFT_Hermitian/fftw_openmp/mpi_mock/hermitian_3d_matrix_production

# Load modules
module purge
module load openmpi/gcc/4.1.6
module load fftw/gcc/openmpi-4.1.6/3.3.10

# Set environment
export OMP_NUM_THREADS=8
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

# Build if needed (with ENABLE_ABACUS_VALIDATION=0)
BINARY="./hermitian_3d_matrix"
if [ ! -f "$BINARY" ]; then
    echo "Building binary with ENABLE_ABACUS_VALIDATION=0..."
    make clean
    make CFLAGS="-DENABLE_ABACUS_VALIDATION=0"
    echo ""
fi

# Navigate to test directory
cd tests/domain_decomposition

# Create results directory
RESULTS_DIR="results"
mkdir -p "$RESULTS_DIR"

# Run test
echo "Running test: N=${N}, Ranks=${RANKS}"
echo ""

# Run with srun (SLURM's MPI launcher)
# Note: Code distributes N/2+1 Y-slice pairs across available ranks
# Each rank processes multiple pairs if num_ranks < N/2+1
srun -n ${RANKS} ../../hermitian_3d_matrix ${N} > "${RESULTS_DIR}/${TEST_NAME}.log" 2>&1
EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ Test ${TEST_NAME} completed successfully"
else
    echo "✗ Test ${TEST_NAME} failed with exit code $EXIT_CODE"
fi

# Validate Hermitian symmetry
VALIDATE_SCRIPT="validate_hermitian.sh"
if [ -f "$VALIDATE_SCRIPT" ] && [ -f "${RESULTS_DIR}/${TEST_NAME}.log" ]; then
    echo ""
    echo "Validating Hermitian symmetry..."
    bash "$VALIDATE_SCRIPT" "${RESULTS_DIR}/${TEST_NAME}.log" | tail -10
fi

echo ""
echo "End time: $(date)"
echo "========================================"
echo "Test complete. Output: ${RESULTS_DIR}/${TEST_NAME}.log"
echo ""

