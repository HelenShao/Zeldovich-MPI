#!/bin/bash
#PBS -N nu_rt256
#PBS -A Abacus
#PBS -q debug
#PBS -l select=1
#PBS -l place=scatter
#PBS -l walltime=01:00:00
#PBS -l filesystems=home:flare
#PBS -j oe
#PBS -o /home/helenshao/InitialConditions/abacus/zeldovich_mpi/neutrino_validation/roundtrip256/pbs_rt256.log
#
# N = 256 seed-noise round trip: build the three variants, run them, compare (results.txt).
# Submit: qsub neutrino_validation/roundtrip256/pbs_rt256.sh
set -euo pipefail

RT=/home/helenshao/InitialConditions/abacus/subprojects/zeldovich_mpi/neutrino_validation/roundtrip256
set +u
module load frameworks fftw/3.3.10
set -u
export CC=mpicc CXX=mpicxx FC=mpifort
export NUMEXPR_MAX_THREADS=256
# ~/.local/bin/meson is a broken uv-tool shim; use the meson of the custom abacus checkout venv
export MESON=/home/helenshao/abacus-store/9575f04_zmpi-neutrinos/.venv/bin/meson

cd "$RT"
echo "host $(hostname)  job ${PBS_JOBID:-none}  $(date)"
echo "[build] -> $RT/build_rt256.log"
bash build_rt256.sh > build_rt256.log 2>&1
echo "[run]"
bash run_rt256.sh
echo "done $(date)"
