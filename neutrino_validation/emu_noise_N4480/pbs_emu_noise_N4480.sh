#!/bin/bash
#PBS -N emu_noise4480
#PBS -A Abacus
#PBS -q capacity
#PBS -l select=16
#PBS -l place=scatter
#PBS -l walltime=08:00:00
#PBS -l filesystems=home:flare
#PBS -j oe
#PBS -o /home/helenshao/InitialConditions/abacus/subprojects/zeldovich_mpi/neutrino_validation/emu_noise_N4480/pbs_emu_noise_N4480.log
#
# One monofonIC run for the emulator box: GridRes 4480 (N0 = 8192), box 1215 Mpc/h, own descriptor.
# Output: $OUT/noise/output_k_space_slab_y.{0..4096} (~4.4 TB), $OUT/mono_N4480_white_noise.hdf5 (~720 GB),
# logs, provenance and checksums in $OUT.
# 256 ranks (16 per node): the slab writer's int index needs < 64 y-planes of the 8192 grid per rank.
# Memory: the 8192^3 Panphasia field alone is ~4.4 TB; 16 nodes give ~17 TB.
# If the existing monofonIC binary does not resolve its libraries on this compute image, it is rebuilt
# (same source, separate build dir) as in pbs_N256_validation.sh.
# Submit: qsub pbs_emu_noise_N4480.sh
set -eo pipefail

VAL=/home/helenshao/InitialConditions/abacus/subprojects/zeldovich_mpi/neutrino_validation
HERE=$VAL/emu_noise_N4480
CONF=EMU1215_noise_N4480.conf
OUT=/flare/Abacus/helenshao/ICs/emu1215_panphasia_N8192
N0=8192
export MONO_DIR=$HOME/monofonic
RANKS_PER_NODE=16
NTHREADS=6

module load frameworks fftw/3.3.10 hdf5/1.14.6
export GSL_ROOT_DIR=/home/helenshao/env
export CC=mpicc CXX=mpicxx FC=mpifort
export OMP_NUM_THREADS=$NTHREADS OMP_PLACES=cores OMP_PROC_BIND=close
NNODES=$(wc -l < "$PBS_NODEFILE")
NRANKS=$((NNODES * RANKS_PER_NODE))
image=$(basename "$(readlink -f /opt/aurora/default)")

if (( (N0 + NRANKS - 1) / NRANKS >= 64 )); then
    echo "$NRANKS ranks give >= 64 y-planes per rank of N0 = $N0: the slab writer's int index overflows" >&2
    exit 2
fi

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
    echo "host $(hostname)  job ${PBS_JOBID:-none}  nodes $NNODES  ranks $NRANKS  image $image"
    echo "monofonIC: $MONO_EXE, $(git -C "$MONO_DIR" rev-parse HEAD) ($(git -C "$MONO_DIR" branch --show-current))"
} > "$OUT/provenance.txt"

cd "$OUT/noise"
echo "=== monofonIC start $(date)"
t0=$SECONDS
set +e
mpiexec -n $NRANKS --ppn $RANKS_PER_NODE --depth=$NTHREADS --cpu-bind depth \
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
    ls output_k_space_slab_y.* | xargs -P 32 -n 16 sha256sum | sort -k2 -V > ../noise_sha256.txt
    chmod a-w output_k_space_slab_y.* && chmod a-w .
    echo "slabs complete; checksums in $OUT/noise_sha256.txt"
else
    echo "slabs incomplete" >&2
    exit 1
fi
exit $rc
