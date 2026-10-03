#!/bin/bash
#PBS -N nu_val256
#PBS -A Abacus
#PBS -q debug
#PBS -l select=1
#PBS -l walltime=01:00:00
#PBS -l filesystems=home
#PBS -j oe
#PBS -o /home/helenshao/InitialConditions/abacus/zeldovich_mpi/neutrino_validation/pbs_N256_validation.log
#
# Steps 1, 2 and 4 of README.md on one compute node (builds must not run on login nodes).
# P(k): CLASS P_cb at z = 31 for abacus_cosm202 (FLAMINGO Planck, 60 meV), made with pk/cosm202_z31.ini.
# Submit: qsub neutrino_validation/pbs_N256_validation.sh
set -euo pipefail

ZMPI_ROOT=/home/helenshao/InitialConditions/abacus/subprojects/zeldovich_mpi
VAL=$ZMPI_ROOT/neutrino_validation
export MONO_DIR="${MONO_DIR:-$HOME/monofonic}"
export PK_FILE="${PK_FILE:-$VAL/pk/cosm202_z31.z1_pk_cb.dat}"

set +u
module load frameworks fftw/3.3.10 hdf5/1.14.6
set -u
export GSL_ROOT_DIR=/home/helenshao/env
export CC=mpicc CXX=mpicxx FC=mpifort
export NUMEXPR_MAX_THREADS=256

cd "$ZMPI_ROOT"
echo "host $(hostname)  job ${PBS_JOBID:-none}  $(date)"
[[ -s "$PK_FILE" ]] || { echo "missing PK_FILE $PK_FILE"; exit 1; }

echo "[build] monofonIC -> $VAL/build_monofonic.log"
bash "$VAL/build_monofonic.sh" > "$VAL/build_monofonic.log" 2>&1

echo "[build] Zeldovich-MPI extnoise -> $VAL/build_zmpi_extnoise.log"
bash "$VAL/build_zmpi_extnoise.sh" > "$VAL/build_zmpi_extnoise.log" 2>&1

MONO_EXE="$MONO_DIR/build/monofonIC" \
ZMPI_EXE="$ZMPI_ROOT/build_extnoise/src/Zeldovich_MPI" \
bash "$VAL/run_N256_validation.sh"
echo "done $(date)"
