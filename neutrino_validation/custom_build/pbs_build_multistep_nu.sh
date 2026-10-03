#!/bin/bash
#PBS -N build_ms_nu
#PBS -A Abacus
#PBS -q capacity
#PBS -l select=1
#PBS -l place=scatter
#PBS -l walltime=01:30:00
#PBS -l filesystems=home:flare
#PBS -j oe
#PBS -o /home/helenshao/InitialConditions/abacus/zeldovich_mpi/neutrino_validation/custom_build/pbs_build_multistep_nu.log
#
# Compute-node half of the custom multistep build (after prepare_checkout.sh). Same recipe as
# AbacusAurora job/hashrun.sh (uv sync --no-editable; meson setup build; meson compile -C build)
# on the 26.181.0 compute image (default since 2026-09-30), the stack the production binary links, plus
#   zeldovich_mpi:extra_cpp_args = LOAD_EXTERNAL_NOISE=1, EXTERNAL_NOISE_DIR=<FLAMINGO noise>
# The noise directory can be overridden at run time with ZD_EXTERNAL_NOISE_DIR.
set -eo pipefail

DEST=/home/helenshao/abacus-store/9575f04_zmpi-neutrinos
NOISE_DIR=/flare/Abacus/helenshao/ICs/flamingo_panphasia_N2048/noise

[[ -e /opt/aurora/26.181.0/oneapi/compiler/latest/lib/libsycl.so.9 ]] \
    || { echo "wrong system image on $(hostname): /opt/aurora/26.181.0 missing" >&2; exit 3; }

cd "$DEST"
. ./env/aurora.sh
echo "host $(hostname)  job ${PBS_JOBID:-none}  $(date)"
module list 2>&1

# c-blosc2's CMake fetches lz4 from GitHub (FetchContent); compute nodes reach it only via the ALCF proxy
export http_proxy=http://proxy.alcf.anl.gov:3128 https_proxy=http://proxy.alcf.anl.gov:3128
export HTTP_PROXY=$http_proxy HTTPS_PROXY=$https_proxy
uv sync --no-editable --offline
meson setup $([[ -f build/build.ninja ]] && echo --reconfigure) build "-Dzeldovich_mpi:extra_cpp_args=['-DLOAD_EXTERNAL_NOISE=1','-DEXTERNAL_NOISE_DIR=\"$NOISE_DIR\"']"
meson configure build | grep -E "extra_cpp_args|particle_output_mode|use_fftw_wisdom|double_precision" || true
meson compile -C build
ls -l build/src/Multistep/multistep* 2>/dev/null || find build -maxdepth 3 -name "multistep*" -type f
echo "done $(date)"
