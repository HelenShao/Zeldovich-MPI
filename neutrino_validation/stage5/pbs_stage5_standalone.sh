#!/bin/bash
#PBS -N nu_st5_sa
#PBS -A Abacus
#PBS -q debug
#PBS -l select=1
#PBS -l place=scatter
#PBS -l walltime=00:30:00
#PBS -l filesystems=home:flare
#PBS -j oe
#PBS -o /home/helenshao/InitialConditions/abacus/subprojects/zeldovich_mpi/neutrino_validation/stage5/pbs_stage5_standalone.log

# Stage 5, standalone half: build_rt_load (LOAD_EXTERNAL_NOISE) reads the FLAMINGO noise and writes
# cosm202 ICs at 300^3 from the hand-written par. Needs roundtrip256's builds and the flamingo_noise run.

set -o pipefail
HERE=/home/helenshao/InitialConditions/abacus/subprojects/zeldovich_mpi/neutrino_validation/stage5
REPO=/home/helenshao/InitialConditions/abacus/subprojects/zeldovich_mpi
EXE=$REPO/build_rt_load/src/Zeldovich_MPI
NOISE=/flare/Abacus/helenshao/ICs/flamingo_panphasia_N2048/noise
RUN=/flare/Abacus/helenshao/stage5_c202_N300/standalone

module load frameworks fftw/3.3.10
[[ -x $EXE ]] || { echo "missing $EXE (run roundtrip256 build first)" >&2; exit 1; }
[[ $(ls "$NOISE"/output_k_space_slab_y.* | wc -l) -eq 1025 ]] || { echo "FLAMINGO noise incomplete in $NOISE" >&2; exit 1; }
if [[ -e $RUN ]]; then
    echo "$RUN exists; move it aside first" >&2
    exit 2
fi
mkdir -p "$RUN"
cp "$HERE/param_stage5_c202_N300.par" "$RUN/"
{
    date; hostname
    echo "exe: $EXE"
    echo "zeldovich_mpi: $(git -C "$REPO" rev-parse HEAD) (+ uncommitted, see git diff)"
    echo "noise: $NOISE"
} > "$RUN/provenance.txt"

cd "$RUN"
export ZD_EXTERNAL_NOISE_DIR=$NOISE
export OMP_NUM_THREADS=16 OMP_PLACES=cores OMP_PROC_BIND=close
t0=$SECONDS
mpiexec -n 4 --ppn 4 --depth=16 --cpu-bind depth "$EXE" param_stage5_c202_N300.par > zmpi.log 2>&1
rc=$?
echo "Zeldovich_MPI rc=$rc after $((SECONDS - t0)) s"
grep -m1 "\[external noise\] reading" zmpi.log
echo "ic files: $(find ic -name 'ic_*' | wc -l)"
exit $rc
