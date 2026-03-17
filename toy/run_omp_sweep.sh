#!/bin/bash
# OMP thread sweep for toy_zloop_omp_scaling (or any executable that respects OMP_NUM_THREADS).
#
# Usage:
#   ./run_omp_sweep.sh <exe> <N> <narray> <param_file> [num_y_repeats] [thread_list] [mpi_ranks]
#
# Examples:
#   ./run_omp_sweep.sh ./toy_zloop_omp_scaling 1024 2 param.par 10
#   ./run_omp_sweep.sh ./toy_zloop_omp_scaling 1024 2 param.par 10 "1 2 4 8 16"
#   ./run_omp_sweep.sh ./toy_zloop_omp_scaling 1024 2 param.par 10 "1 2 4 8 16" 9
#
# mpi_ranks: if 2 or more, run with mpiexec -n <mpi_ranks> to mimic production multi-node.
#            if 0 or 1 or omitted, run exe directly (single process).
#
# Pure OMP control (no MPI_Init / no mpiexec), using the toy:
#   cd /home/helenshao/InitialConditions/hermitian_3d_matrix_production/toy
#   make toy-zloop
#   export TOY_DISABLE_MPI=1
#   ./run_omp_sweep.sh ./toy_zloop_omp_scaling 1024 2 \
#       /home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_N1024/CPD-225/MPI-9/mode3/1D/param_N1024_CPD_225.par \
#       10 "1 2 4 8 16"
# Output: one block per thread count; Z-loop time in "[Rank r] [PTimerWall]" lines.

set -e
EXE="${1:?Usage: $0 <exe> <N> <narray> <param_file> [num_y_repeats] [thread_list] [mpi_ranks]}"
N="${2:?}"
NARRAY="${3:?}"
PARAM="${4:?}"
REPEATS="${5:-10}"
THREADS="${6:-1 2 4 8 16}"
MPI_RANKS="${7:-1}"

if [ "$MPI_RANKS" -ge 2 ] 2>/dev/null; then
  USE_MPI=1
else
  USE_MPI=0
  MPI_RANKS=1
fi

echo "OMP sweep: exe=$EXE N=$N narray=$NARRAY param=$PARAM num_y_repeats=$REPEATS MPI_ranks=$MPI_RANKS"
echo "Thread list: $THREADS"
echo "---"

for t in $THREADS; do
  echo "========== OMP_NUM_THREADS=$t =========="
  export OMP_NUM_THREADS=$t
  if [ "$USE_MPI" -eq 1 ]; then
    mpiexec -n "$MPI_RANKS" "$EXE" "$N" "$NARRAY" "$PARAM" "$REPEATS"
  else
    "$EXE" "$N" "$NARRAY" "$PARAM" "$REPEATS"
  fi
  echo ""
done

echo "--- sweep done ---"
