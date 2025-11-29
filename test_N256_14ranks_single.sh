#!/bin/bash
# Test script for N=256, 14 ranks - SINGLE PRECISION
# Version 15.1.2: Tests that matrix output is purely real

echo "========================================"
echo "Testing: N=256, 14 ranks (SINGLE PRECISION)"
echo "Version: 15.1.2"
echo "Date: $(date)"
echo "========================================"
echo ""

# Navigate to directory
cd /scratch/gpfs/hshao/C_Bible/InitialConditions/hermitian_3d_matrix_production

# Load modules
module purge
module load openmpi/gcc/4.1.6
module load fftw/gcc/openmpi-4.1.6/3.3.10

# Set environment
export OMP_NUM_THREADS=4
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

# Build binary with SINGLE precision (default)
echo "Building binary (SINGLE PRECISION)..."
make clean
make CFLAGS="-DENABLE_ABACUS_VALIDATION=0"
echo ""

# Run test
echo "Running: mpirun -np 14 ./hermitian_3d_matrix 256"
echo "Expected: Matrix output should be purely real (imaginary parts ≈ 0)"
echo ""

mpirun -np 14 ./hermitian_3d_matrix 256 2>&1 | tee v15_1_2_test/test_N256_14ranks_single.log

echo ""
echo "Test complete!"
echo "Check test_N256_14ranks_single.log for output"
echo "========================================"

