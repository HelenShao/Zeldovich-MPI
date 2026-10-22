#!/bin/bash
# Build Kernel_MPI (double precision) in its own Meson build directory.
#
# Kernel_MPI is only built when the build is configured with -Ddouble_precision=true,
# so it lives in a separate build dir from the (float) Zeldovich_MPI production build.
#
# Usage:  tests/kernel_regression/build_kernel.sh            # -> build/kernel_double/src/Kernel_MPI
# Override: KERNEL_BUILD=<dir>, ENV_SCRIPT=<site env to source (modules, PATH, PKG_CONFIG_PATH)>
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${HERE}/../.." && pwd)"
KERNEL_BUILD="${KERNEL_BUILD:-${REPO_ROOT}/build/kernel_double}"

if [[ -n "${ENV_SCRIPT:-}" && -f "${ENV_SCRIPT}" ]]; then
  # shellcheck source=/dev/null
  source "${ENV_SCRIPT}"
fi

MESON_SETUP=(
  -Ddouble_precision=true
  -Dproduction_mode=true
  -Ddebug_prints=false
  -Ddebug_rng_consistency=false
  -Ddebug_rng_skip=false
  -Duse_fftw_wisdom=false
)

cd "${REPO_ROOT}"
if [[ -f "${KERNEL_BUILD}/meson-private/coredata.dat" ]]; then
  meson setup "${KERNEL_BUILD}" "${MESON_SETUP[@]}" --reconfigure
else
  meson setup "${KERNEL_BUILD}" "${MESON_SETUP[@]}"
fi
meson compile -C "${KERNEL_BUILD}" Kernel_MPI

echo "Build OK: ${KERNEL_BUILD}/src/Kernel_MPI"
echo "Quick regression: python3 ${HERE}/run_regression.py --build-dir ${KERNEL_BUILD}"
echo "  (or: meson test -C ${KERNEL_BUILD} kernel_regression_quick)"
