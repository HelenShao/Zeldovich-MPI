#!/bin/bash
# Build Zeldovich-MPI (branch neutrinos) with the external-noise reader and the D dump.
# Run from anywhere after loading the usual modules (e.g. module load frameworks fftw/3.3.10).
#
# Usage: bash build_zmpi_extnoise.sh > build_zmpi_extnoise.log 2>&1
# Result: $BUILD_DIR/src/Zeldovich_MPI
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$REPO_ROOT/build_extnoise}"
DUMP_D="${DUMP_D:-1}"   # set DUMP_D=0 for production runs (the dump is ~47 GB at N = 1800)

cd "$REPO_ROOT"
git branch --show-current
git log --oneline -1

SETUP=(
  -Dproduction_mode=true
  -Dparticle_output_mode=3
  -Dskip_file_write=0
  -Dparallelize_z_loop=1
  -Ddebug_prints=false
  -Ddebug_rng_consistency=false
  -Ddebug_rng_skip=false
  "-Dextra_cpp_args=-DLOAD_EXTERNAL_NOISE=1,-DDUMP_D_SLABS=${DUMP_D}"
)

if [[ -f "$BUILD_DIR/meson-private/coredata.dat" ]]; then
  meson setup "$BUILD_DIR" "${SETUP[@]}" --reconfigure
else
  meson setup "$BUILD_DIR" "${SETUP[@]}"
fi
meson compile -C "$BUILD_DIR"

ls -l "$BUILD_DIR/src/Zeldovich_MPI"
