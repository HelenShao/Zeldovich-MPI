#!/bin/bash
#PBS -N nu_crop1800
#PBS -A Abacus
#PBS -q debug-scaling
#PBS -l select=2
#PBS -l place=scatter
#PBS -l walltime=01:00:00
#PBS -l filesystems=home:flare
#PBS -j oe
#PBS -o /home/helenshao/InitialConditions/abacus/zeldovich_mpi/neutrino_validation/crop1800/pbs_crop1800.log
#
# Crop test (b): standalone Zeldovich-MPI at N = 1800 on the FLAMINGO noise (N0 = 2048) with a
# flat P(k), then crop_check_1800.py (per-plane k-space identity + real-space block vs monofonIC's
# white_noise). Outputs (~260 GB: ICs + D slabs) go to $RUN on Flare.
# Needs build_rt_load (roundtrip256) and the FLAMINGO noise run. Submit: qsub crop1800/pbs_crop1800.sh
set -eo pipefail

VAL=/home/helenshao/InitialConditions/abacus/subprojects/zeldovich_mpi/neutrino_validation
HERE=$VAL/crop1800
ZMPI_EXE=$VAL/../build_rt_load/src/Zeldovich_MPI
NOISE=/flare/Abacus/helenshao/ICs/flamingo_panphasia_N2048
RUN=/flare/Abacus/helenshao/ICs/crop1800_flat
RANKS_PER_NODE=12
NTHREADS=8

module load frameworks fftw/3.3.10 hdf5/1.14.6
export NUMEXPR_MAX_THREADS=256
export OMP_NUM_THREADS=$NTHREADS OMP_PLACES=cores OMP_PROC_BIND=close
export ZD_EXTERNAL_NOISE_DIR=$NOISE/noise
NNODES=$(wc -l < "$PBS_NODEFILE")

[[ -x $ZMPI_EXE ]] || { echo "missing $ZMPI_EXE" >&2; exit 1; }
[[ $(ls "$NOISE"/noise/output_k_space_slab_y.* | wc -l) == 1025 ]] || { echo "FLAMINGO noise incomplete in $NOISE/noise" >&2; exit 1; }
[[ -e $RUN ]] && { echo "$RUN exists; move it aside first" >&2; exit 2; }
mkdir -p "$RUN"
cp "$HERE/param_N1800_flat.par" "$RUN/"
cd "$RUN"
echo "host $(hostname)  job ${PBS_JOBID:-none}  nodes $NNODES  $(date)"

t0=$SECONDS
mpiexec -n $((NNODES * RANKS_PER_NODE)) --ppn $RANKS_PER_NODE --depth=$NTHREADS --cpu-bind depth \
    "$ZMPI_EXE" param_N1800_flat.par > zmpi_run.log 2>&1
echo "Zeldovich-MPI done after $((SECONDS - t0)) s"
grep -m1 "external noise" zmpi_run.log || true

python3 "$HERE/crop_check_1800.py" --noise-dir "$NOISE/noise" --mono-h5 "$NOISE/mono_N1800_white_noise.hdf5" \
    --ref-noise-dir "$VAL/crop300/run/external_noise" --ref-N0 512 \
    --out "$HERE/crop1800_results.txt" > crop_check.log 2>&1 || true
cat "$HERE/crop1800_results.txt"
echo "done $(date)"
