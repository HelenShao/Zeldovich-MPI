#!/bin/bash
# Build the three Zeldovich-MPI variants of the N = 256 round-trip test (compute node only:
# meson/ninja on a login node fails on libxpmem). All use FFTW_ESTIMATE so the FFT plans are identical.
#   plain: seed path, no dumps                  -> $REPO/build_rt_plain
#   seed:  seed path + DUMP_NOISE_SLABS + D dump -> $REPO/build_rt_seed
#   load:  LOAD_EXTERNAL_NOISE + D dump          -> $REPO/build_rt_load
# Otherwise the options of the production embed (abacus 9575f04 src/meson.build): single precision,
# particle_output_mode 3, parallelize_z_loop 1 (standalone: no wisdom, since ESTIMATE plans are deterministic).
#
# Usage: [MESON=/path/to/meson] bash build_rt256.sh > build_rt256.log 2>&1
set -euo pipefail

MESON="${MESON:-meson}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO"
git branch --show-current
git log --oneline -1
git status --short src

COMMON=(
  -Dproduction_mode=true
  -Dparticle_output_mode=3
  -Dskip_file_write=0
  -Dparallelize_z_loop=1
  -Ddebug_prints=false
  -Ddebug_rng_consistency=false
  -Ddebug_rng_skip=false
  -Duse_fftw_wisdom=false
)
declare -A ARGS=(
  [plain]="-DFFTW_PLANNER_FLAGS=FFTW_ESTIMATE"
  [seed]="-DFFTW_PLANNER_FLAGS=FFTW_ESTIMATE,-DDUMP_NOISE_SLABS=1,-DDUMP_D_SLABS=1"
  [load]="-DFFTW_PLANNER_FLAGS=FFTW_ESTIMATE,-DLOAD_EXTERNAL_NOISE=1,-DDUMP_D_SLABS=1"
)

for v in plain seed load; do
  B="$REPO/build_rt_$v"
  echo "=== $v: ${ARGS[$v]}"
  if [[ -f "$B/meson-private/coredata.dat" ]]; then
    "$MESON" setup "$B" "${COMMON[@]}" "-Dextra_cpp_args=${ARGS[$v]}" --reconfigure
  else
    "$MESON" setup "$B" "${COMMON[@]}" "-Dextra_cpp_args=${ARGS[$v]}"
  fi
  "$MESON" compile -C "$B" Zeldovich_MPI
  ls -l "$B/src/Zeldovich_MPI"
done
