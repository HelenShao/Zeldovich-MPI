#!/bin/bash
#PBS -N nu_rv2_seed
#PBS -A Abacus
#PBS -q debug
#PBS -l select=1
#PBS -l place=scatter
#PBS -l walltime=00:45:00
#PBS -l filesystems=home:flare
#PBS -j oe
#PBS -o /home/helenshao/InitialConditions/abacus/subprojects/zeldovich_mpi/neutrino_validation/roundtrip_rockv2/pbs_seed_dump.log

# rockV2 round trip, part 1 (debug image, same as the rt256 standalone builds): standalone
# Zeldovich-MPI at the rockV2_n800 settings with DUMP_NOISE_SLABS -> seed noise w + reference ICs.
# Needs build_rt_seed from roundtrip256/build_rt256.sh.

set -o pipefail
HERE=/home/helenshao/InitialConditions/abacus/subprojects/zeldovich_mpi/neutrino_validation/roundtrip_rockv2
REPO=/home/helenshao/InitialConditions/abacus/subprojects/zeldovich_mpi
EXE=$REPO/build_rt_seed/src/Zeldovich_MPI
RUN=/flare/Abacus/helenshao/roundtrip_rockv2/seed_standalone

module load frameworks fftw/3.3.10
[[ -x $EXE ]] || { echo "missing $EXE (run roundtrip256 build first)" >&2; exit 1; }
if [[ -e $RUN ]]; then
    echo "$RUN exists; move it aside first" >&2
    exit 2
fi
mkdir -p "$RUN"
cp "$HERE/param_rockv2_seed.par" "$RUN/"
{
    date; hostname
    echo "exe: $EXE"
    echo "zeldovich_mpi: $(git -C "$REPO" rev-parse HEAD) (+ uncommitted, see git diff)"
} > "$RUN/provenance.txt"

cd "$RUN"
export OMP_NUM_THREADS=16 OMP_PLACES=cores OMP_PROC_BIND=close
t0=$SECONDS
mpiexec -n 4 --ppn 4 --depth=16 --cpu-bind depth "$EXE" param_rockv2_seed.par > zmpi.log 2>&1
rc=$?
echo "Zeldovich_MPI rc=$rc after $((SECONDS - t0)) s"
echo "noise slabs: $(ls noise_slabs/output_k_space_slab_y.* 2>/dev/null | wc -l), ic files: $(find ic -name 'ic_*' | wc -l)"
exit $rc
