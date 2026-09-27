#!/bin/bash
# Optional site-specific setup. Source a cluster env script if provided, otherwise
# assume CC, CXX, meson, ninja, mpiexec, and FFTW are already on PATH.
#
# Override before sourcing:
#   ENV_SCRIPT=/path/to/your/site-env.sh
#   DATA_DIR=/path/to/pk-and-plt-files
#   IC_OUTPUT_DIR=/path/to/write/ics
#   NTHREADS_PER_RANK=8
set -euo pipefail

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -n "${ENV_SCRIPT:-}" && -f "${ENV_SCRIPT}" ]]; then
  # shellcheck source=/dev/null
  source "${ENV_SCRIPT}"
fi

export DATA_DIR="${DATA_DIR:-${BUNDLE_DIR}/data}"
export IC_OUTPUT_DIR="${IC_OUTPUT_DIR:-${BUNDLE_DIR}/output}"
export PK_FILE="${PK_FILE:-${DATA_DIR}/camb_planck15.dat}"
export PLT_FILE="${PLT_FILE:-${DATA_DIR}/eigmodes128}"

export NTHREADS_PER_RANK="${NTHREADS_PER_RANK:-${OMP_NUM_THREADS:-1}}"
export MPIEXEC="${MPIEXEC:-mpiexec}"

# Optional: set MPIEXEC_EXTRA="--bind-to core" (or site-specific flags) before sourcing run.sh
export MPIEXEC_EXTRA="${MPIEXEC_EXTRA:-}"
