#!/bin/bash
#PBS -N flamingo_nu
#PBS -A Abacus
#PBS -q capacity
#PBS -l select=2
#PBS -l place=scatter
#PBS -l walltime=24:00:00
#PBS -l filesystems=home:flare
#PBS -j oe

# FLAMINGO-phase-matched neutrino production run (PRODUCTION_NEUTRINO_ICS.md): 1800^3 abacus_cosmNNN
# from z = 99 to 0 (standard CLASS_power, ZD_Pk_norm = 8) with embedded Zeldovich-MPI reading the
# FLAMINGO Panphasia noise. The z = 31 runs are in /flare/Abacus/helenshao/*_N1800_z31start.
#   qsub -v COSM=202 -o .../production/flamingo_c202_z99.pbs.log pbs_flamingo_nu.sh
# Resubmitting the same COSM resumes from the last checkpoint (abacus.run without --clean).
# Optional: -v COSM=202,EXTRA="ZD_CornerModes=1" forwards KEY=VAL overrides to abacus.run.

set -o pipefail
: "${COSM:?qsub -v COSM=202 (one of 202..206)}"

CODE=/home/helenshao/abacus-store/9575f04_zmpi-neutrinos   # 9575f04 + zeldovich_mpi neutrinos, LOAD_EXTERNAL_NOISE=1
AA=/home/helenshao/InitialConditions/AbacusAurora
PAR2=$AA/Simulations/FLAMINGO_nu_c${COSM}.par2
PROD=/home/helenshao/InitialConditions/abacus/subprojects/zeldovich_mpi/neutrino_validation/production
NOISE=/flare/Abacus/helenshao/ICs/flamingo_panphasia_N2048/noise
SIMNAME=AbacusAurora_FLAMINGO_nu_c${COSM}_N1800

cd "$CODE" && . ./env/aurora.sh
unset ZD_EXTERNAL_NOISE_DIR   # use the compiled-in FLAMINGO noise directory

if [[ ! -e /opt/aurora/26.181.0/oneapi/compiler/latest/lib/libsycl.so.9 ]]; then
    echo "wrong system image on $(hostname): /opt/aurora/26.181.0 missing" >&2
    exit 3
fi
[[ -r $PAR2 ]] || { echo "missing $PAR2" >&2; exit 2; }
[[ -r $NOISE/output_k_space_slab_y.0 && -r $NOISE/output_k_space_slab_y.1024 ]] \
    || { echo "FLAMINGO noise slabs missing under $NOISE" >&2; exit 2; }

RUN=$PROD/c${COSM}
mkdir -p "$RUN"
cd "$RUN"
attempt_tag=${PBS_JOBID%%.*}

{
    date
    echo "job: $PBS_JOBID on $(sort -u "$PBS_NODEFILE" | tr '\n' ' ')"
    echo "code: $CODE"
    echo "par2: $PAR2 (AbacusAurora $(git -C "$AA" rev-parse --short HEAD))"
    echo "zeldovich_mpi: $(git -C "$CODE/subprojects/zeldovich_mpi" rev-parse HEAD 2>/dev/null) (+ uncommitted neutrinos-branch edits)"
    echo "overrides: ${EXTRA:-none}"
    if [[ -d /flare/Abacus/helenshao/$SIMNAME ]]; then echo "existing working dir: resuming"; else echo "fresh start"; fi
    module list 2>&1
} > "provenance.$attempt_tag.txt"

# Halt cleanly (and restartably) 20 min before walltime, as multisim.pbs does.
wall=$(qstat -f "$PBS_JOBID" 2>/dev/null | awk -F'= ' '/Resource_List.walltime/{print $2}')
IFS=: read -r wh wm ws <<< "${wall:-24:00:00}"
export ABACUS_JOB_HALT_TIME=$(( $(date +%s) + (10#$wh * 60 + 10#$wm - 20) * 60 ))
export ABACUS_JOB_HALT_FILE=/flare/Abacus/helenshao/HALT.$attempt_tag
trap 'rm -f "$ABACUS_JOB_HALT_FILE"' EXIT
echo "halt at $(date -d @"$ABACUS_JOB_HALT_TIME"); touch $ABACUS_JOB_HALT_FILE to stop early"

sort -u "$PBS_NODEFILE" > "hostfile.$attempt_tag"
overrides=()
for o in ${EXTRA:-}; do overrides+=("$o"); done

echo "=== onesim start: $(date) ==="
t0=$SECONDS
"$AA/job/onesim.sh" "$PAR2" "$RUN/hostfile.$attempt_tag" ${overrides[@]+"${overrides[@]}"} \
    < /dev/null > "sim.$attempt_tag.out" 2>&1
rc=$?
echo "=== onesim rc=$rc after $((SECONDS - t0)) s: $(date) ==="

if grep -q "\[external noise\] reading $NOISE/" "sim.$attempt_tag.out"; then
    echo "external noise: OK ($(grep -m1 '\[external noise\] reading' "sim.$attempt_tag.out"))"
elif grep -q "Embedded Zeldovich\|Zeldovich_MPI" "sim.$attempt_tag.out"; then
    echo "WARNING: ICs were generated but no '[external noise] reading $NOISE/' line was found" >&2
else
    echo "external noise: no IC generation in this attempt (resumed from a checkpoint)"
fi
exit $rc
