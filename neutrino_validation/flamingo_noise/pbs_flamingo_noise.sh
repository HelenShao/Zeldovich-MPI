#!/bin/bash
#PBS -N flamingo_noise
#PBS -A Abacus
#PBS -q debug-scaling
#PBS -l select=2
#PBS -l place=scatter
#PBS -l walltime=01:00:00
#PBS -l filesystems=home:flare
#PBS -j oe
#PBS -o /home/helenshao/InitialConditions/abacus/zeldovich_mpi/neutrino_validation/flamingo_noise/pbs_flamingo_noise.log
#
# One monofonIC run with the FLAMINGO L1000N1800 Panphasia descriptor (GridRes 1800, N0 = 2048).
# Output: $OUT/noise/output_k_space_slab_y.{0..1024} (~69 GB), $OUT/mono_N1800_white_noise.hdf5,
# logs and checksums in $OUT. monofonIC was built on the default (26.26.0) image, hence debug-scaling.
# Submit: qsub pbs_flamingo_noise.sh
set -eo pipefail

HERE=/home/helenshao/InitialConditions/abacus/subprojects/zeldovich_mpi/neutrino_validation/flamingo_noise
OUT=/flare/Abacus/helenshao/ICs/flamingo_panphasia_N2048
MONO_EXE=$HOME/monofonic/build/monofonIC
RANKS_PER_NODE=12
NTHREADS=8

module load frameworks fftw/3.3.10 hdf5/1.14.6
export OMP_NUM_THREADS=$NTHREADS OMP_PLACES=cores OMP_PROC_BIND=close
NNODES=$(wc -l < "$PBS_NODEFILE")

[[ -e $OUT/noise ]] && { echo "$OUT/noise exists; move it aside first" >&2; exit 2; }
mkdir -p "$OUT/noise"
cp "$HERE/FLAMINGO_noise_N1800.conf" "$OUT/"
{
    date
    echo "host $(hostname)  job ${PBS_JOBID:-none}  nodes $NNODES  ranks $((NNODES * RANKS_PER_NODE))"
    echo "monofonIC: $(git -C "$HOME/monofonic" rev-parse HEAD) ($(git -C "$HOME/monofonic" branch --show-current))"
} > "$OUT/provenance.txt"

cd "$OUT/noise"
echo "=== monofonIC start $(date)"
t0=$SECONDS
mpiexec -n $((NNODES * RANKS_PER_NODE)) --ppn $RANKS_PER_NODE --depth=$NTHREADS --cpu-bind depth \
    "$MONO_EXE" ../FLAMINGO_noise_N1800.conf > ../monofonic.log 2>&1
echo "=== monofonIC done after $((SECONDS - t0)) s $(date)"

nslab=$(ls output_k_space_slab_y.* | wc -l)
nbad=$(find . -name "output_k_space_slab_y.*" ! -size $((2048 * 2048 * 16))c | wc -l)
echo "slabs: $nslab (expect 1025), wrong size: $nbad (expect 0, size $((2048 * 2048 * 16)) B)" | tee -a ../provenance.txt
grep -iE "level|N0|panphasia|descriptor" ../monofonic.log | head -20 >> ../provenance.txt || true
sha256sum output_k_space_slab_y.* > ../noise_sha256.txt
chmod a-w output_k_space_slab_y.* && chmod a-w .
echo "done $(date)"
