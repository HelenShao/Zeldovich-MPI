#!/bin/bash
# Test v15.2: Power spectrum mode (with parameter file) - DOUBLE PRECISION
# Should use zeldovich-PLT PowerSpectrum and RNG

#SBATCH --job-name=test_v15_2_ps_dbl
#SBATCH --output=/scratch/gpfs/hshao/C_Bible/InitialConditions/hermitian_3d_matrix_production/v15_2_tests/test_v15_2_ps_double.log
#SBATCH --error=/scratch/gpfs/hshao/C_Bible/InitialConditions/hermitian_3d_matrix_production/v15_2_tests/test_v15_2_ps_double.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00

module purge
module load openmpi/gcc/4.1.6
module load fftw/gcc/openmpi-4.1.6/3.3.10

# Verify modules are loaded
echo "Module environment:"
module list
echo "FFTW paths:"
echo "  FFTW_INCLUDE: ${FFTW_INCLUDE:-not set}"
echo "  PKG_CONFIG_PATH: ${PKG_CONFIG_PATH:-not set}"

cd /scratch/gpfs/hshao/C_Bible/InitialConditions/hermitian_3d_matrix_production

N=256
N_RANKS=14
PARAM_FILE="examples/param_template.par"

echo "=========================================="
echo "Test: v15.2 Power Spectrum Mode - DOUBLE PRECISION"
echo "N=$N, ranks=$N_RANKS"
echo "Parameter file: $PARAM_FILE"
echo "Should use zeldovich-PLT PowerSpectrum"
echo "=========================================="
echo ""

# Rebuild with double precision
echo "Rebuilding with double precision..."
make clean
if ! make CFLAGS="-DUSE_DOUBLE_PRECISION" 2>&1 | tee /tmp/build_double.log | tail -10; then
    echo "ERROR: Build failed!"
    cat /tmp/build_double.log | tail -20
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

srun ./hermitian_3d_matrix $N $PARAM_FILE 2>&1 | tee v15_2_tests/test_v15_2_ps_double_output.log

echo ""
echo "Checking if output is purely real..."
grep -i "Array.*Max\|RESULT\|purely real" v15_2_tests/test_v15_2_ps_double_output.log | tail -10

echo ""
echo "Checking power spectrum initialization..."
grep -i "power.*spectrum\|PowerSpectrum\|powerlaw\|Precision:" v15_2_tests/test_v15_2_ps_double_output.log | head -10

