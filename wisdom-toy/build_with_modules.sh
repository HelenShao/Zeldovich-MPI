#!/usr/bin/env bash
# Load Princeton-style (frameworks + fftw) modules, then run make in this directory.
# Usage: ./build_with_modules.sh          # same as make all
#        ./build_with_modules.sh clean
#        ./build_with_modules.sh clean all

set -eo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Lmod init may reference optional shell vars; do not use nounset here.
if ! type module &>/dev/null; then
    if [[ -f /usr/share/lmod/lmod/init/bash ]]; then
        # shellcheck source=/dev/null
        source /usr/share/lmod/lmod/init/bash
    elif [[ -f /etc/profile.d/modules.sh ]]; then
        # shellcheck source=/dev/null
        source /etc/profile.d/modules.sh
    fi
fi

module load frameworks
module load fftw/3.3.10

make -C "$ROOT" "${@:-all}"
