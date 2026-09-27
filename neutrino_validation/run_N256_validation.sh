#!/bin/bash
# N = 256 external-white-noise validation, run on an interactive compute node.
# Steps (each can be skipped by setting the matching SKIP_* variable to 1):
#   1. monofonIC at GridRes = 256  -> white-noise slabs + first-order delta HDF5
#   2. Zeldovich-MPI with LOAD_EXTERNAL_NOISE + DUMP_D_SLABS -> D_slabs/
#   3. validate_external_noise.py  -> validation_results.txt
# Every step writes its output to a log file in $RUN_DIR.
#
# Usage:
#   MONO_EXE=$HOME/monofonic/build/monofonIC \
#   ZMPI_EXE=<repo>/build_extnoise/src/Zeldovich_MPI \
#   PK_FILE=/path/to/pk_cb_z31.txt \
#   bash run_N256_validation.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_DIR="${RUN_DIR:-$HERE/run256}"
MONO_EXE="${MONO_EXE:?set MONO_EXE to the monofonIC executable}"
ZMPI_EXE="${ZMPI_EXE:?set ZMPI_EXE to the Zeldovich_MPI executable}"
PK_FILE="${PK_FILE:?set PK_FILE to the P(k) file at InitialRedshift}"

N=256
BOX=673.2
MONO_RANKS="${MONO_RANKS:-4}"
ZMPI_RANKS="${ZMPI_RANKS:-4}"
NTHREADS="${NTHREADS:-8}"
MPIEXEC="${MPIEXEC:-mpiexec}"
# Aurora-style binding; override MPI_OPTS for other clusters (e.g. MPI_OPTS="" for plain mpirun)
MPI_OPTS="${MPI_OPTS:---ppn ${ZMPI_RANKS} --depth=${NTHREADS} --cpu-bind depth}"

export OMP_NUM_THREADS=$NTHREADS OMP_PLACES=cores OMP_PROC_BIND=close

mkdir -p "$RUN_DIR"/{mono,external_noise}
cd "$RUN_DIR"

# ---------------------------------------------------------------- 1. monofonIC
if [[ "${SKIP_MONO:-0}" != 1 ]]; then
    echo "[1] monofonIC GridRes=$N (log: $RUN_DIR/mono/monofonic.log)"
    cp "$HERE/monofonic_N256_test.conf" mono/
    (cd mono && $MPIEXEC -n "$MONO_RANKS" $MPI_OPTS "$MONO_EXE" monofonic_N256_test.conf > monofonic.log 2>&1)
    mv mono/output_k_space_slab_y.* external_noise/
fi
nslab=$(ls external_noise/output_k_space_slab_y.* 2>/dev/null | wc -l)
echo "    external_noise/: $nslab slab files (expect $((N / 2 + 1)))"
echo "    slab size: $(stat -c %s external_noise/output_k_space_slab_y.0 2>/dev/null || stat -f %z external_noise/output_k_space_slab_y.0) bytes (expect $((N * N * 16)) for double)"

# ---------------------------------------------------------------- 2. Zeldovich-MPI
if [[ "${SKIP_ZMPI:-0}" != 1 ]]; then
    echo "[2] Zeldovich-MPI N=$N (log: $RUN_DIR/zmpi_run.log)"
    sed "s|REPLACE_WITH/pk_cb_z31.txt|${PK_FILE}|" "$HERE/param_N256_extnoise.par" > param_N256_extnoise.par
    rm -rf D_slabs
    $MPIEXEC -n "$ZMPI_RANKS" $MPI_OPTS "$ZMPI_EXE" param_N256_extnoise.par > zmpi_run.log 2>&1
fi
echo "    D_slabs/: $(ls D_slabs/D_slab_y.* 2>/dev/null | wc -l) files (expect $((N / 2 + 1)))"

# ---------------------------------------------------------------- 3. validation
echo "[3] validate_external_noise.py (results: $RUN_DIR/validation_results.txt)"
python3 "$HERE/validate_external_noise.py" \
    --slab-dir external_noise --N0 "$N" \
    --ours-slabs D_slabs --N "$N" \
    --box "$BOX" --pk "$PK_FILE" \
    --mono-hdf5 mono/mono_N256_delta.hdf5 \
    --out validation_results.txt > validate.log 2>&1
cat validation_results.txt
