#!/bin/bash
# Test v15.3: Power spectrum mode (with parameter file)
# Should use zeldovich-PLT PowerSpectrum and RNG

#SBATCH --job-name=v15_3_ps
#SBATCH --output=/scratch/gpfs/hshao/C_Bible/InitialConditions/hermitian_3d_matrix_production/v15_3_tests/test_v15_3_ps.log
#SBATCH --error=/scratch/gpfs/hshao/C_Bible/InitialConditions/hermitian_3d_matrix_production/v15_3_tests/test_v15_3_ps.err
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
PARAM_FILE="examples/param_template.par"

echo "=========================================="
echo "Test: v15.3 Power Spectrum Mode"
echo "N=$N, ranks=$N_RANKS"
echo "Parameter file: $PARAM_FILE"
echo "Should use zeldovich-PLT PowerSpectrum"
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

# Check if parameter file exists
if [ ! -f "$PARAM_FILE" ]; then
    echo "ERROR: Parameter file not found: $PARAM_FILE"
    exit 1
fi

srun ./hermitian_3d_matrix $N $PARAM_FILE 2>&1 | tee v15_3_tests/test_v15_3_ps_output.log

echo ""
echo "Checking if output is purely real..."
grep -i "imaginary\|complex\|real" v15_3_tests/test_v15_3_ps_output.log | tail -5

echo ""
echo "Checking power spectrum initialization..."
grep -i "power.*spectrum\|PowerSpectrum\|powerlaw" v15_3_tests/test_v15_3_ps_output.log | head -10

