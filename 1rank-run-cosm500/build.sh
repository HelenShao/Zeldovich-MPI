#!/bin/bash
# Meson build for Cosm500 1-rank (mode 3, FFTW wisdom, float precision).
#
# Requires: meson, ninja, MPI C/C++ wrappers (mpicc/mpicxx or CC/CXX), FFTW with OpenMP.
# Override: REPO_ROOT, MESON_BUILD, FFTW_ROOT (-Dfftw_root=...)
set -euo pipefail

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${REPO_ROOT:-$(cd "${BUNDLE_DIR}/../.." && pwd)}"
MESON_BUILD="${MESON_BUILD:-${REPO_ROOT}/build/cosm500_1d_1rank}"

# shellcheck source=/dev/null
source "${BUNDLE_DIR}/env.sh"

if ! command -v meson >/dev/null 2>&1 || ! command -v ninja >/dev/null 2>&1; then
  echo "ERROR: meson and ninja must be on PATH (or set ENV_SCRIPT)." >&2
  exit 1
fi

cd "${REPO_ROOT}"

MESON_SETUP=(
  -Dproduction_mode=true
  -Dparticle_output_mode=3
  -Dskip_file_write=0
  -Dparallelize_z_loop=1
  -Ddebug_prints=false
  -Ddebug_rng_consistency=true
  -Ddebug_rng_skip=false
  -Duse_fftw_wisdom=true
  -Dextra_cpp_args=-g
)

if [[ -n "${FFTW_ROOT:-}" ]]; then
  MESON_SETUP+=(-Dfftw_root="${FFTW_ROOT}")
fi

if [[ -f "${MESON_BUILD}/meson-private/coredata.dat" ]]; then
  meson setup "${MESON_BUILD}" "${MESON_SETUP[@]}" --reconfigure
else
  meson setup "${MESON_BUILD}" "${MESON_SETUP[@]}"
fi

meson compile -C "${MESON_BUILD}"

export ZMPI_EXE="${MESON_BUILD}/src/Zeldovich_MPI"
export WISDOM_EXE="${MESON_BUILD}/src/wisdom_rank0"

for exe in "${ZMPI_EXE}" "${WISDOM_EXE}"; do
  if [[ ! -x "${exe}" ]]; then
    echo "Build failed: missing ${exe}" >&2
    exit 1
  fi
done

echo "Build OK:"
echo "  ZMPI_EXE=${ZMPI_EXE}"
echo "  WISDOM_EXE=${WISDOM_EXE}"
