#!/bin/bash
#PBS -N fl_noise2500_sha
#PBS -A Abacus
#PBS -q capacity
#PBS -l select=1
#PBS -l walltime=02:00:00
#PBS -l filesystems=home:flare
#PBS -j oe
#PBS -o /home/helenshao/InitialConditions/abacus/subprojects/zeldovich_mpi/neutrino_validation/flamingo_noise_N2500/pbs_flamingo_noise_N2500_sha256.log
#
# Checksums and read-only protection for the FLAMINGO N0 = 4096 noise slabs written by pbs_flamingo_noise_N2500.sh.
#   qsub -W depend=afterok:<noise job id> pbs_flamingo_noise_N2500_sha256.sh
set -eo pipefail

OUT=/flare/Abacus/helenshao/ICs/flamingo_panphasia_N4096
N0=4096
cd "$OUT/noise"

nexp=$((N0 / 2 + 1))
nslab=$(ls output_k_space_slab_y.* | wc -l)
[[ $nslab -eq $nexp ]] || { echo "slabs: $nslab, expect $nexp" >&2; exit 1; }

echo "=== sha256 start $(date)"
ls output_k_space_slab_y.* | xargs -P 64 -n 16 sha256sum | sort -k2 -V > ../noise_sha256.txt
echo "=== sha256 done $(date): $(wc -l < ../noise_sha256.txt) entries"
chmod a-w output_k_space_slab_y.* && chmod a-w .
echo "slabs read-only; checksums in $OUT/noise_sha256.txt"
