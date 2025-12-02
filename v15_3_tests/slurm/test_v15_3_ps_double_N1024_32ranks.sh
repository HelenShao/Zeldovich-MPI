#!/bin/bash
# Test v15.3: Power spectrum mode (with parameter file) - DOUBLE PRECISION
# Large array test: N=1024, 32 ranks, 2 nodes

#SBATCH --job-name=v15_3_ps_dbl_N1024_32
#SBATCH --output=/scratch/gpfs/hshao/C_Bible/InitialConditions/hermitian_3d_matrix_production/v15_3_tests/test_v15_3_ps_double_N1024_32ranks.log
#SBATCH --error=/scratch/gpfs/hshao/C_Bible/InitialConditions/hermitian_3d_matrix_production/v15_3_tests/test_v15_3_ps_double_N1024_32ranks.err
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00

module purge
module load openmpi/gcc/4.1.6
module load fftw/gcc/openmpi-4.1.6/3.3.10

cd /scratch/gpfs/hshao/C_Bible/InitialConditions/hermitian_3d_matrix_production

N=1024
N_RANKS=32
PARAM_FILE="examples/param_template.par"

echo "=========================================="
echo "Test: v15.3 Power Spectrum Mode - DOUBLE PRECISION"
echo "Large Array Test: N=$N, ranks=$N_RANKS"
echo "Parameter file: $PARAM_FILE"
echo "=========================================="
echo ""

# MPI/UCX WORKAROUND:
# The application completes all computation and verification successfully,
# but OpenMPI + UCX + hcoll occasionally crash during MPI_Finalize with
# heap corruption inside UCX/hcoll. Disable hcoll and libnbc collectives
# to route around the buggy path.
export OMPI_MCA_coll_hcoll_enable=0
export OMPI_MCA_coll=^hcoll,libnbc

# Build with double precision
echo "Building with double precision..."
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

# Update parameter file for N=1024 if needed
# Note: The parameter file has CPD=256, but we're using N=1024 from command line
# This is fine - the code will use N from command line and warn if mismatch

echo "Starting computation..."
echo "Memory estimate for N=$N (double precision):"
echo "  Full matrix: ~$(echo "scale=2; $N * $N * $N * 4 * 16 / 1024 / 1024 / 1024" | bc) GB"
echo ""

srun ./hermitian_3d_matrix $N $PARAM_FILE 2>&1 | tee v15_3_tests/test_v15_3_ps_double_N1024_32ranks_output.log

echo ""
echo "Checking if output is purely real..."
grep -i "Array.*Max\|RESULT\|purely real" v15_3_tests/test_v15_3_ps_double_N1024_32ranks_output.log | tail -10

echo ""
echo "Checking precision info..."
grep -i "Precision:" v15_3_tests/test_v15_3_ps_double_N1024_32ranks_output.log | head -3

echo ""
echo "Checking timing..."
grep -i "TIMING\|Stage.*time\|Total.*time" v15_3_tests/test_v15_3_ps_double_N1024_32ranks_output.log | tail -10

