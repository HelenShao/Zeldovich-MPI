#!/bin/bash
# Test v15.2: Power spectrum mode (with parameter file) - DOUBLE PRECISION
# Uses pre-built double precision binary

#SBATCH --job-name=test_v15_2_ps_dbl2
#SBATCH --output=/scratch/gpfs/hshao/C_Bible/InitialConditions/hermitian_3d_matrix_production/v15_2_tests/test_v15_2_ps_double2.log
#SBATCH --error=/scratch/gpfs/hshao/C_Bible/InitialConditions/hermitian_3d_matrix_production/v15_2_tests/test_v15_2_ps_double2.err
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
echo "Test: v15.2 Power Spectrum Mode - DOUBLE PRECISION"
echo "N=$N, ranks=$N_RANKS"
echo "Parameter file: $PARAM_FILE"
echo "=========================================="
echo ""

# Check if binary exists
if [ ! -f "./hermitian_3d_matrix" ]; then
    echo "ERROR: Binary not found! Please build with: make CFLAGS=\"-DUSE_DOUBLE_PRECISION\""
    exit 1
fi

# Check if parameter file exists
if [ ! -f "$PARAM_FILE" ]; then
    echo "ERROR: Parameter file not found: $PARAM_FILE"
    exit 1
fi

srun ./hermitian_3d_matrix $N $PARAM_FILE 2>&1 | tee v15_2_tests/test_v15_2_ps_double2_output.log

echo ""
echo "Checking if output is purely real..."
grep -i "Array.*Max\|RESULT\|purely real" v15_2_tests/test_v15_2_ps_double2_output.log | tail -10

echo ""
echo "Checking precision info..."
grep -i "Precision:" v15_2_tests/test_v15_2_ps_double2_output.log | head -3

