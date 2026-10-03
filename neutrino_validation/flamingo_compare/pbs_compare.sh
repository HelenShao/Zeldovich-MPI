#!/bin/bash
#PBS -N flam_cmp
#PBS -A Abacus
#PBS -q debug
#PBS -l select=1
#PBS -l walltime=01:00:00
#PBS -l filesystems=home:flare
#PBS -j oe

# Compare a cosm202 subsample redshift with a FLAMINGO snapshot (compare_subsample_flamingo.py).
#   qsub -v Z=0.500,FLAM=/path/to/flamingo_0068.hdf5 -o .../flamingo_compare/cmp_z0.500.pbs.log pbs_compare.sh
#   qsub -v Z=4.000,ABACUS_ONLY=1 ...    # test the Abacus loader only
set -eo pipefail
: "${Z:?qsub -v Z=0.500,...}"

CODE=/home/helenshao/abacus-store/9575f04_zmpi-neutrinos
HERE=/home/helenshao/InitialConditions/abacus/subprojects/zeldovich_mpi/neutrino_validation/flamingo_compare
SUB=/flare/Abacus/helenshao/AbacusAurora_FLAMINGO_nu_c${COSM:-202}_N1800/subsamples/z$Z

cd "$CODE" && . ./env/aurora.sh
export NUMEXPR_MAX_THREADS=256 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
cd "$HERE"

tag=c${COSM:-202}_z$Z
if [[ -n ${ABACUS_ONLY:-} ]]; then
    python -u compare_subsample_flamingo.py --abacus "$SUB" --abacus-only --out "abacus_only_$tag.txt"
else
    : "${FLAM:?qsub -v FLAM=<FLAMINGO snapshot .hdf5 or dir>}"
    python -u compare_subsample_flamingo.py --abacus "$SUB" --flamingo "$FLAM" ${OVERLOAD:+--overload $OVERLOAD} \
        --out "compare_$tag.txt" --plot "compare_$tag.png"
fi
