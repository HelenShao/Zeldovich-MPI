#!/bin/bash
#SBATCH --job-name=hermitian_test_16ranks
#SBATCH --nodes=2
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=4
#SBATCH --time=00:10:00
#SBATCH --output=v15_1_2_test/test_N256_16ranks_2nodes_double_%j.out
#SBATCH --error=v15_1_2_test/test_N256_16ranks_2nodes_double_%j.err

# Test script for N=256, 16 ranks across 2 nodes - DOUBLE PRECISION
# Version 15.1.2: Tests that matrix output is purely real

echo "========================================"
echo "Testing: N=256, 16 ranks across 2 nodes (DOUBLE PRECISION)"
echo "Version: 15.1.2"
echo "Date: $(date)"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Nodes: $SLURM_JOB_NUM_NODES"
echo "Tasks per node: $SLURM_NTASKS_PER_NODE"
echo "Total tasks: $SLURM_NTASKS"
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

# Build binary with DOUBLE precision
echo "Building binary (DOUBLE PRECISION)..."
make clean
make CFLAGS="-DUSE_DOUBLE_PRECISION -DENABLE_ABACUS_VALIDATION=0"
echo ""

# Run test with SLURM's srun
echo "Running: srun -n 16 ./hermitian_3d_matrix 256"
echo "Expected: Matrix output should be purely real (imaginary parts ≈ 0)"
echo ""

srun -n 16 ./hermitian_3d_matrix 256 2>&1 | tee v15_1_2_test/test_N256_16ranks_2nodes_double.log

echo ""
echo "Test complete!"
echo "Check v15_1_2_test/test_N256_16ranks_2nodes_double.log for output"
echo "========================================"

