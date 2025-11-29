#!/bin/bash
# Interactive script to run test with 16 ranks across 2 nodes using salloc
# This script uses salloc to allocate resources, then runs the test

echo "========================================"
echo "Requesting 2 nodes with 16 tasks (8 per node)"
echo "This will allocate resources interactively"
echo "========================================"
echo ""

# Request 2 nodes, 16 tasks total (8 per node), 4 CPUs per task
# Time limit: 10 minutes
salloc --nodes=2 --ntasks=16 --ntasks-per-node=8 --cpus-per-task=4 --time=00:10:00 bash << 'EOF'
# Script runs inside the allocated resources

echo "========================================"
echo "Resources allocated!"
echo "Nodes: $SLURM_JOB_NUM_NODES"
echo "Tasks: $SLURM_NTASKS"
echo "Tasks per node: $SLURM_NTASKS_PER_NODE"
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

# Run test
echo "========================================"
echo "Testing: N=256, 16 ranks across 2 nodes (DOUBLE PRECISION)"
echo "Version: 15.1.2"
echo "Date: $(date)"
echo "========================================"
echo ""
echo "Running: srun -n 16 ./hermitian_3d_matrix 256"
echo "Expected: Matrix output should be purely real (imaginary parts ≈ 0)"
echo ""

srun -n 16 ./hermitian_3d_matrix 256 2>&1 | tee v15_1_2_test/test_N256_16ranks_2nodes_double.log

echo ""
echo "Test complete!"
echo "Check v15_1_2_test/test_N256_16ranks_2nodes_double.log for output"
echo "========================================"

EOF

echo ""
echo "salloc session ended"
echo "========================================"

