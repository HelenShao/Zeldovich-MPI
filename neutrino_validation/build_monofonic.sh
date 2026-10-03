#!/bin/bash
# Clone and build monofonIC with Willem's output_k_space_slabs branch.
# Needs: CMake, MPI, FFTW3 (double, with MPI + OpenMP), GSL, HDF5, a Fortran compiler.
# Panphasia licence: register at http://icc.dur.ac.uk/Panphasia.php before using the ICs.
#
# Usage: MONO_DIR=$HOME/monofonic [BUILD_DIR=build] bash build_monofonic.sh > build_monofonic.log 2>&1
set -euo pipefail

MONO_DIR="${MONO_DIR:-$HOME/monofonic}"
BUILD_DIR="${BUILD_DIR:-build}"

if [[ ! -d "$MONO_DIR/.git" ]]; then
    git clone --branch output_k_space_slabs https://github.com/wullm/monofonic.git "$MONO_DIR"
fi
cd "$MONO_DIR"
git checkout output_k_space_slabs
git log --oneline -1

# CLASS / zwindstroom / FastDF are only needed for transfer functions and neutrino particles,
# not for the white noise; turning them off avoids the FetchContent downloads.
cmake -S . -B "$BUILD_DIR" \
      -DCMAKE_BUILD_TYPE=Release \
      -DENABLE_MPI=ON \
      -DENABLE_PANPHASIA=ON \
      -DENABLE_CLASS=OFF \
      -DENABLE_ZWINDSTROOM=OFF \
      -DENABLE_FASTDF=OFF \
      -DCODE_PRECISION=DOUBLE
cmake --build "$BUILD_DIR" -j 16

ls -l "$BUILD_DIR/monofonIC"
