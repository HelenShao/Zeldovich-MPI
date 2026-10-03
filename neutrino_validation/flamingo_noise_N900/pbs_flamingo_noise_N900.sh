#!/bin/bash
#PBS -N flam_noise900
#PBS -A Abacus
#PBS -q debug
#PBS -l select=1
#PBS -l place=scatter
#PBS -l walltime=01:00:00
#PBS -l filesystems=home:flare
#PBS -j oe
#PBS -o /home/helenshao/InitialConditions/abacus/subprojects/zeldovich_mpi/neutrino_validation/flamingo_noise_N900/pbs_flamingo_noise_N900.log
#
# One monofonIC run with the FLAMINGO Panphasia descriptor at GridRes 900 (N0 = 1024), box 673.2 Mpc/h.
# Output: $OUT/noise/output_k_space_slab_y.{0..512} (~8.6 GB), logs, provenance and checksums in $OUT.
# If the existing monofonIC binary does not resolve its libraries on this compute image, it is rebuilt
# (same source, separate build dir) as in pbs_N256_validation.sh.
# Submit: qsub pbs_flamingo_noise_N900.sh
set -eo pipefail

VAL=/home/helenshao/InitialConditions/abacus/subprojects/zeldovich_mpi/neutrino_validation
HERE=$VAL/flamingo_noise_N900
CONF=FLAMINGO_noise_N900.conf
OUT=/flare/Abacus/helenshao/ICs/flamingo_panphasia_N1024
N0=1024
export MONO_DIR=$HOME/monofonic
RANKS_PER_NODE=12
NTHREADS=8

module load frameworks fftw/3.3.10 hdf5/1.14.6
export GSL_ROOT_DIR=/home/helenshao/env
export CC=mpicc CXX=mpicxx FC=mpifort
export OMP_NUM_THREADS=$NTHREADS OMP_PLACES=cores OMP_PROC_BIND=close
NNODES=$(wc -l < "$PBS_NODEFILE")
image=$(basename "$(readlink -f /opt/aurora/default)")

MONO_EXE=$MONO_DIR/build/monofonIC
if ldd "$MONO_EXE" | grep -q "not found"; then
    export BUILD_DIR=build_$image
    MONO_EXE=$MONO_DIR/$BUILD_DIR/monofonIC
    if [[ ! -x $MONO_EXE ]] || ldd "$MONO_EXE" | grep -q "not found"; then
        echo "[build] monofonIC for image $image -> $HERE/build_monofonic_$image.log"
        bash "$VAL/build_monofonic.sh" > "$HERE/build_monofonic_$image.log" 2>&1
    fi
fi

[[ -e $OUT/noise ]] && { echo "$OUT/noise exists; move it aside first" >&2; exit 2; }
mkdir -p "$OUT/noise"
cp "$HERE/$CONF" "$OUT/"
{
    date
    echo "host $(hostname)  job ${PBS_JOBID:-none}  nodes $NNODES  ranks $((NNODES * RANKS_PER_NODE))  image $image"
    echo "monofonIC: $MONO_EXE, $(git -C "$MONO_DIR" rev-parse HEAD) ($(git -C "$MONO_DIR" branch --show-current))"
} > "$OUT/provenance.txt"

cd "$OUT/noise"
echo "=== monofonIC start $(date)"
t0=$SECONDS
set +e
mpiexec -n $((NNODES * RANKS_PER_NODE)) --ppn $RANKS_PER_NODE --depth=$NTHREADS --cpu-bind depth \
    "$MONO_EXE" ../$CONF > ../monofonic.log 2>&1
rc=$?
set -e
echo "=== monofonIC rc=$rc after $((SECONDS - t0)) s $(date)" | tee -a ../provenance.txt

nexp=$((N0 / 2 + 1))
size=$((N0 * N0 * 16))
nslab=$(ls output_k_space_slab_y.* 2>/dev/null | wc -l)
nbad=$(find . -name "output_k_space_slab_y.*" ! -size ${size}c | wc -l)
echo "slabs: $nslab (expect $nexp), wrong size: $nbad (expect 0, size $size B)" | tee -a ../provenance.txt
grep -iE "level|N0|panphasia|descriptor" ../monofonic.log | head -20 >> ../provenance.txt || true
if [[ $nslab -eq $nexp && $nbad -eq 0 ]]; then
    sha256sum output_k_space_slab_y.* > ../noise_sha256.txt
    chmod a-w output_k_space_slab_y.* && chmod a-w .
    echo "slabs complete; checksums in $OUT/noise_sha256.txt (rc=$rc is expected from the testing-mode power spectrum)"
else
    echo "slabs incomplete" >&2
    exit 1
fi
echo "done $(date)"
