#!/bin/bash
# Run the N = 256 round-trip test after build_rt256.sh (compute node).
#   run_plain/: seed ICs, no dumps (reference for "dumping does not perturb the RNG")
#   run_seed/:  seed ICs + noise_slabs/ (w = D/sqrt(P)) + D_slabs/
#   run_load/:  LOAD_EXTERNAL_NOISE from run_seed/noise_slabs (symlinked as external_noise) + D_slabs/
# then compare_rt256.py -> results.txt. Refuses to overwrite existing run directories.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
RANKS="${RANKS:-4}"          # NumZRanks = 2 -> 2 x 2 rank grid
NTHREADS="${NTHREADS:-16}"
export OMP_NUM_THREADS=$NTHREADS OMP_PLACES=cores OMP_PROC_BIND=close

cd "$HERE"
for v in plain seed load; do
    [[ -e run_$v ]] && { echo "run_$v exists; move it aside first" >&2; exit 2; }
done

for v in plain seed load; do
    mkdir run_$v
    cp param_rt256.par run_$v/
    [[ $v == load ]] && ln -s ../run_seed/noise_slabs run_load/external_noise
    echo "=== $v start $(date)"
    t0=$SECONDS
    (cd run_$v && mpiexec -n "$RANKS" --ppn "$RANKS" --depth="$NTHREADS" --cpu-bind depth \
        "$REPO/build_rt_$v/src/Zeldovich_MPI" param_rt256.par > zmpi.log 2>&1)
    echo "=== $v rc=0 after $((SECONDS - t0)) s"
done

python3 compare_rt256.py --N 256 --out results.txt > compare.log 2>&1
cat results.txt
