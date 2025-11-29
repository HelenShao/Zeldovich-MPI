#!/bin/bash
# Test v15.2: Fallback mode (no parameter file)
# Should use local RNG implementation like v15.1.2

#SBATCH --job-name=test_v15_2_fallback
#SBATCH --output=/scratch/gpfs/hshao/C_Bible/InitialConditions/hermitian_3d_matrix_production/v15_2_tests/test_v15_2_fallback.log
#SBATCH --error=/scratch/gpfs/hshao/C_Bible/InitialConditions/hermitian_3d_matrix_production/v15_2_tests/test_v15_2_fallback.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00

module purge
module load openmpi/gcc/4.1.6
module load fftw/gcc/openmpi-4.1.6/3.3.10

cd /scratch/gpfs/hshao/C_Bible/InitialConditions/hermitian_3d_matrix_production

N=256
N_RANKS=14

echo "=========================================="
echo "Test: v15.2 Fallback Mode (no param file)"
echo "N=$N, ranks=$N_RANKS"
echo "Should use local RNG (like v15.1.2)"
echo "=========================================="
echo ""

srun ./hermitian_3d_matrix $N 2>&1 | tee v15_2_tests/test_v15_2_fallback_output.log

echo ""
echo "Checking if output is purely real..."
grep -i "imaginary\|complex\|real" v15_2_tests/test_v15_2_fallback_output.log | tail -5

