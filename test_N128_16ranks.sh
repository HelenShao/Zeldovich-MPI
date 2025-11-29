#!/bin/bash
# Simple test script for N=128, 16 ranks on login node

echo "========================================"
echo "Testing: N=128, 16 ranks"
echo "========================================"
echo "Date: $(date)"
echo ""

# Navigate to directory
cd /scratch/gpfs/hshao/C_Bible/FFT_Hermitian/fftw_openmp/mpi_mock/hermitian_3d_matrix_production

# Load modules
module purge
module load openmpi/gcc/4.1.6
module load fftw/gcc/openmpi-4.1.6/3.3.10

# Set environment (use fewer threads for login node)
export OMP_NUM_THREADS=4
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

# Build binary (always rebuild to ensure latest changes)
BINARY="./hermitian_3d_matrix"
echo "Building binary with ENABLE_ABACUS_VALIDATION=0..."
make clean
make CFLAGS="-DENABLE_ABACUS_VALIDATION=0"
echo ""

# Expected grid decomposition
echo "Expected grid decomposition:"
echo "  grid_x = 4, grid_z = 4"
echo "  X-dimension: 4 chunks (128/4 = 32 elements per chunk)"
echo "  Z-dimension: 4 chunks (128/4 = 32 elements per chunk)"
echo "  Each rank: 32 X-values × 128 Y-values × 32 Z-values"
echo ""

# Run test
echo "Running: mpirun -np 16 ./hermitian_3d_matrix 128"
echo ""

mpirun -np 16 ./hermitian_3d_matrix 128

echo ""
echo "Test complete!"
echo "========================================"
