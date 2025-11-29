#!/bin/bash
# Test v15.2: Power spectrum mode (with parameter file) - DOUBLE PRECISION
# Medium array test: N=256, 64 ranks, 8 nodes

#SBATCH --job-name=test_v15_2_ps_dbl_N256_64
#SBATCH --output=/scratch/gpfs/hshao/C_Bible/InitialConditions/hermitian_3d_matrix_production/v15_2_tests/test_v15_2_ps_double_N256_64ranks_8nodes.log
#SBATCH --error=/scratch/gpfs/hshao/C_Bible/InitialConditions/hermitian_3d_matrix_production/v15_2_tests/test_v15_2_ps_double_N256_64ranks_8nodes.err
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00

module purge
module load openmpi/gcc/4.1.6
module load fftw/gcc/openmpi-4.1.6/3.3.10

cd /scratch/gpfs/hshao/C_Bible/InitialConditions/hermitian_3d_matrix_production

N=256
N_RANKS=64
PARAM_FILE="examples/param_template.par"

echo "=========================================="
echo "Test: v15.2 Power Spectrum Mode - DOUBLE PRECISION"
echo "Medium Array Test: N=$N, ranks=$N_RANKS, nodes=8"
echo "Parameter file: $PARAM_FILE"
echo "=========================================="
echo ""

# MPI/UCX WORKAROUND:
# Disable hcoll and libnbc collectives (avoids UCX/hcoll teardown issues)
export OMPI_MCA_coll_hcoll_enable=0
export OMPI_MCA_coll=^hcoll,libnbc

# Check if binary exists (built without ASan, double precision)
if [ ! -f "./hermitian_3d_matrix" ]; then
    echo "ERROR: Binary not found! Please build with: make clean && make CFLAGS=\"-DUSE_DOUBLE_PRECISION -O3 -march=native\""
    exit 1
fi

# Check if parameter file exists
if [ ! -f "$PARAM_FILE" ]; then
    echo "ERROR: Parameter file not found: $PARAM_FILE"
    exit 1
fi

echo "Starting computation..."
echo "Memory estimate for N=$N (double precision):"
echo "  Full matrix: ~$(echo "scale=2; $N * $N * $N * 4 * 16 / 1024 / 1024 / 1024" | bc) GB"
echo ""

srun -N 8 -n $N_RANKS ./hermitian_3d_matrix $N $PARAM_FILE 2>&1 | tee v15_2_tests/test_v15_2_ps_double_N256_64ranks_8nodes_output.log

echo ""
echo "Checking if output is purely real..."
grep -i "Array.*Max\|RESULT\|purely real" v15_2_tests/test_v15_2_ps_double_N256_64ranks_8nodes_output.log | tail -10

echo ""
echo "Checking precision info..."
grep -i "Precision:" v15_2_tests/test_v15_2_ps_double_N256_64ranks_8nodes_output.log | head -3

echo ""
echo "Checking timing..."
grep -i "TIMING\|Stage.*time\|Total.*time" v15_2_tests/test_v15_2_ps_double_N256_64ranks_8nodes_output.log | tail -10

echo ""
echo "Done."


