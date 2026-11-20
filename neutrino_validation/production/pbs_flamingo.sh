#!/bin/bash
#PBS -N flamingo
#PBS -A Abacus
#PBS -q capacity
#PBS -l select=16
#PBS -l place=scatter
#PBS -l walltime=06:00:00
#PBS -l filesystems=home:flare
#PBS -j oe
#PBS -o /home/helenshao/InitialConditions/abacus/subprojects/zeldovich_mpi/neutrino_validation/production/flamingo.pbs.log

# FLAMINGO-matched boxes at AuroraFlagship mass resolution (2026-11-17): NP = 2500^3 in FLAMINGO's 1 Gpc box, any FLAMINGO
# cosmology (Simulations/flamingo_nu_N2500/AbacusAurora_FLAMINGO_<label>_N2500.par2, Cosmologies/FLAMINGO_cosm/flamingo_<label>/),
# z = 99 (2LPT + PLT) to 0, embedded Zeldovich-MPI reading FLAMINGO's Panphasia noise at N0 = 4096 (GridRes 2500).
# 16 nodes = the capacity-queue maximum (0.98e9 particles/node; CPD 567, layout 8 x 4 from FLAMINGO.par2 NumZRanks = 4). Estimated ~2.5-3 h
# (1800^3 took 7.0 h on 2 nodes = 8.7e-6 node-s per particle); 6 h walltime, halts cleanly 20 min before and resumes on resubmit.
# Adapted from pbs_flamingo_m9.sh (now in flamingo_1800_local/).
#   qsub -v SIMNAME=AbacusAurora_FLAMINGO_<label>_N2500 -o $PWD/flamingo_<label>_N2500.pbs.log pbs_flamingo.sh
# Resubmitting resumes from the last checkpoint (abacus.run without --clean).

set -o pipefail

CODE=/home/helenshao/abacus-store/9575f04_zmpi-neutrinos   # 9575f04 + zeldovich_mpi neutrinos, LOAD_EXTERNAL_NOISE=1
AA=${AA:-/home/helenshao/InitialConditions/AbacusAurora-flamingo}   # git worktree of AbacusAurora on branch flamingo-n2500 (the main checkout may be on another branch)
[[ -n ${SIMNAME:-} ]] || { echo "set SIMNAME with qsub -v SIMNAME=AbacusAurora_FLAMINGO_<label>_N2500" >&2; exit 2; }
PAR2=$AA/Simulations/flamingo_nu_N2500/$SIMNAME.par2
PROD=/home/helenshao/InitialConditions/abacus/subprojects/zeldovich_mpi/neutrino_validation/production
NOISE=/flare/Abacus/helenshao/ICs/flamingo_panphasia_N4096/noise

cd "$CODE" && . ./env/aurora.sh
# The build's compiled-in noise directory is the N0 = 2048 set: this override selects the N0 = 4096 set (required for NP = 2500^3).
export ZD_EXTERNAL_NOISE_DIR=$NOISE

if [[ ! -e /opt/aurora/26.181.0/oneapi/compiler/latest/lib/libsycl.so.9 ]]; then
    echo "wrong system image on $(hostname): /opt/aurora/26.181.0 missing" >&2
    exit 3
fi
[[ -r $PAR2 ]] || { echo "missing $PAR2" >&2; exit 2; }
COSM=$(sed -n 's|^#include "\.\./\.\./Cosmologies/\(.*\)/cosm.def"|\1|p' "$PAR2")
[[ -n $COSM && -s $AA/Cosmologies/$COSM/CLASS_power ]] || { echo "missing CLASS_power for cosmology '$COSM' in $PAR2" >&2; exit 2; }
# Checkpoint containers (DAOS) are not supported by this build; aurora.def needs the variable defined.
export ABACUS_CHECKPOINT_CONTAINERS=0
[[ -r $NOISE/output_k_space_slab_y.0 && -r $NOISE/output_k_space_slab_y.2048 && ! -e $NOISE/output_k_space_slab_y.2049 ]] \
    || { echo "N0 = 4096 noise slabs (0..2048) not found under $NOISE" >&2; exit 2; }

# Pre-generated ICs (ExternalICs = 1): check they are there and laid out for this run (4 z-ranks of CPD files each).
EXTERNAL_ICS=0
if grep -Eq '^ExternalICs *= *1' "$PAR2"; then
    EXTERNAL_ICS=1
    ICDIR=$(sed -n 's|^InitialConditionsDirectory *= *"\(.*\)".*|\1|p' "$PAR2")
    nz=$(find "$ICDIR" -maxdepth 1 -type d -name "z*" 2>/dev/null | wc -l)
    (( nz == 4 )) || { echo "ExternalICs = 1 but $ICDIR has $nz z-rank directories (expect 4 for NumZRanks = 4)" >&2; exit 2; }
fi

RUN=$PROD/${RUNDIR:-${SIMNAME#AbacusAurora_FLAMINGO_}}
mkdir -p "$RUN"
cd "$RUN"
attempt_tag=${PBS_JOBID%%.*}

{
    date
    echo "job: $PBS_JOBID on $(sort -u "$PBS_NODEFILE" | tr '\n' ' ')"
    echo "code: $CODE"
    echo "par2: $PAR2 (AbacusAurora $(git -C "$AA" rev-parse --short HEAD) + uncommitted FLAMINGO_N2500 / FLAMINGO_cosm files); cosmology $COSM"
    echo "zeldovich_mpi: $(git -C "$CODE/subprojects/zeldovich_mpi" rev-parse HEAD 2>/dev/null)"
    echo "ZD_EXTERNAL_NOISE_DIR: $ZD_EXTERNAL_NOISE_DIR"
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

echo "=== onesim start: $(date) ==="
t0=$SECONDS
"$AA/job/onesim.sh" "$PAR2" "$RUN/hostfile.$attempt_tag" < /dev/null > "sim.$attempt_tag.out" 2>&1
rc=$?
echo "=== onesim rc=$rc after $((SECONDS - t0)) s: $(date) ==="

if (( EXTERNAL_ICS )); then
    echo "ICs: pre-generated (ExternalICs = 1), read from $ICDIR; no IC generation in this run"
elif grep -q "\[external noise\] reading $NOISE/" "sim.$attempt_tag.out"; then
    echo "external noise: OK ($(grep -m1 '\[external noise\] reading' "sim.$attempt_tag.out"))"
elif grep -q "\[external noise\] reading" "sim.$attempt_tag.out"; then
    echo "ERROR: ICs read the WRONG noise: $(grep -m1 '\[external noise\] reading' "sim.$attempt_tag.out")" >&2
    rc=4
elif grep -q "Embedded Zeldovich\|Zeldovich_MPI" "sim.$attempt_tag.out"; then
    echo "WARNING: ICs were generated but no '[external noise] reading $NOISE/' line was found" >&2
else
    echo "external noise: no IC generation in this attempt (resumed from a checkpoint)"
fi
exit $rc
