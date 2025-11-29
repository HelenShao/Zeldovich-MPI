#!/bin/bash
# Test v15.2: Power spectrum mode (with parameter file)
# Should use zeldovich-PLT PowerSpectrum and RNG

#SBATCH --job-name=test_v15_2_ps
#SBATCH --output=/scratch/gpfs/hshao/C_Bible/InitialConditions/hermitian_3d_matrix_production/v15_2_tests/test_v15_2_ps.log
#SBATCH --error=/scratch/gpfs/hshao/C_Bible/InitialConditions/hermitian_3d_matrix_production/v15_2_tests/test_v15_2_ps.err
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
echo "Test: v15.2 Power Spectrum Mode"
echo "N=$N, ranks=$N_RANKS"
echo "Parameter file: $PARAM_FILE"
echo "Should use zeldovich-PLT PowerSpectrum"
echo "=========================================="
echo ""

# Check if parameter file exists
if [ ! -f "$PARAM_FILE" ]; then
    echo "ERROR: Parameter file not found: $PARAM_FILE"
    exit 1
fi

srun ./hermitian_3d_matrix $N $PARAM_FILE 2>&1 | tee v15_2_tests/test_v15_2_ps_output.log

echo ""
echo "Checking if output is purely real..."
grep -i "imaginary\|complex\|real" v15_2_tests/test_v15_2_ps_output.log | tail -5

echo ""
echo "Checking power spectrum initialization..."
grep -i "power.*spectrum\|PowerSpectrum\|powerlaw" v15_2_tests/test_v15_2_ps_output.log | head -10

