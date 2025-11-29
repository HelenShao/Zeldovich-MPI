#!/bin/bash
# Test script for N=256, 14 ranks (Version 15 with PCG RNG and nskip tracking)

echo "========================================"
echo "Testing Version 15: N=256, 14 ranks"
echo "Features: PCG RNG with nskip tracking"
echo "========================================"
echo "Date: $(date)"
echo ""

# Navigate to directory
cd /scratch/gpfs/hshao/C_Bible/FFT_Hermitian/fftw_openmp/mpi_mock/hermitian_3d_matrix_production

# Load modules
module purge
module load openmpi/gcc/4.1.6
module load fftw/gcc/openmpi-4.1.6/3.3.10

# Set environment
export OMP_NUM_THREADS=4
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

# Build binary (always rebuild to ensure latest changes)
BINARY="./hermitian_3d_matrix"
echo "Building binary with latest changes..."
make clean
make CFLAGS="-DENABLE_ABACUS_VALIDATION=0"
echo ""

# Expected grid decomposition (14 = 2 × 7, so likely 2×7 or 7×2)
echo "Expected grid decomposition:"
echo "  For 14 ranks (14 = 2 × 7):"
echo "  - grid_x = 2, grid_z = 7 (or 7×2)"
echo "  - X-dimension: 2 chunks (256/2 = 128 elements per chunk)"
echo "  - Z-dimension: 7 chunks (256/7 = 36 remainder 4)"
echo "  - First 4 ranks in Z: 37 Z-values each"
echo "  - Last 3 ranks in Z: 36 Z-values each"
echo "  - Each rank: 128 X-values × 256 Y-values × (36-37) Z-values"
echo ""

# Run test
echo "Running: mpirun -np 14 ./hermitian_3d_matrix 256"
echo ""

mpirun -np 14 ./hermitian_3d_matrix 256

echo ""
echo "Test complete!"
echo "========================================"

