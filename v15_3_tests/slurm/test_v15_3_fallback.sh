#!/bin/bash
# Test v15.3: Fallback mode (no parameter file)
# Should use local RNG implementation like v15.1.2

#SBATCH --job-name=v15_3_fallback
#SBATCH --output=/scratch/gpfs/hshao/C_Bible/InitialConditions/hermitian_3d_matrix_production/v15_3_tests/test_v15_3_fallback.log
#SBATCH --error=/scratch/gpfs/hshao/C_Bible/InitialConditions/hermitian_3d_matrix_production/v15_3_tests/test_v15_3_fallback.err
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
echo "Test: v15.3 Fallback Mode (no param file)"
echo "N=$N, ranks=$N_RANKS"
echo "Should use local RNG (like v15.1.2)"
echo "=========================================="
echo ""

# Build with single precision
echo "Building with single precision..."
make clean
if ! make 2>&1 | tee /tmp/build_single.log | tail -10; then
    echo "ERROR: Build failed!"
    cat /tmp/build_single.log | tail -20
    exit 1
fi

# Verify binary exists
if [ ! -f "./hermitian_3d_matrix" ]; then
    echo "ERROR: Binary not found after build!"
    exit 1
fi
echo "Build successful!"

srun ./hermitian_3d_matrix $N 2>&1 | tee v15_3_tests/test_v15_3_fallback_output.log

echo ""
echo "Checking if output is purely real..."
grep -i "imaginary\|complex\|real" v15_3_tests/test_v15_3_fallback_output.log | tail -5

