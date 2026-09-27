#!/usr/bin/env bash
# HAL single-socket OMP scaling sweep for fftw_2d_thread_sweep.
# Assumes this script and fftw_2d_thread_sweep live in the same folder.

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXE="${SCRIPT_DIR}/fftw_2d_thread_sweep"

# You can override these at launch, e.g. N=4096 REPEATS=20 ./run_hal_single_socket_sweep.sh
N="${N:-2048}"
NARRAY="${NARRAY:-4}"
REPEATS="${REPEATS:-30}"
THREADS=(${THREADS:-1 2 4 7 10 12 14})

RUN_ROOT="${SCRIPT_DIR}/hal-single-socket-sweep"
AFFINITY="0-13"
NUMAPIN=(numactl --cpunodebind=0 --membind=0)

[[ -x "${EXE}" ]] || {
    echo "Missing executable: ${EXE}"
    echo "Build first (example): make -C \"${SCRIPT_DIR}\" fftw_2d_thread_sweep"
    exit 1
}

mkdir -p "${RUN_ROOT}"

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_DISPLAY_AFFINITY=false
export GOMP_CPU_AFFINITY="${AFFINITY}"

echo "=== HAL single-socket sweep ==="
echo "Host: $(hostname)"
echo "Run root: ${RUN_ROOT}"
echo "Threads: ${THREADS[*]}"
echo "Config: N=${N} NARRAY=${NARRAY} REPEATS=${REPEATS}"
echo "Affinity: GOMP_CPU_AFFINITY=${GOMP_CPU_AFFINITY}"
echo "NUMA pin: ${NUMAPIN[*]}"
echo

for T in "${THREADS[@]}"; do
    RUNDIR="${RUN_ROOT}/omp-${T}"
    LOG="${RUNDIR}/N${N}_fftw2d_OMP${T}_timings.log"
    mkdir -p "${RUNDIR}"

    export OMP_NUM_THREADS="${T}"
    export FFTW_WISDOM_FILE="${RUNDIR}/fftw_wisdom_float.wisdom"

    {
        echo "=== FFTW 2D single-socket sweep ==="
        echo "Date: $(date -Is)"
        echo "Host: $(hostname)"
        echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}"
        echo "OMP_PLACES=${OMP_PLACES} OMP_PROC_BIND=${OMP_PROC_BIND}"
        echo "GOMP_CPU_AFFINITY=${GOMP_CPU_AFFINITY}"
        echo "FFTW_WISDOM_FILE=${FFTW_WISDOM_FILE}"
        echo "N=${N} NARRAY=${NARRAY} REPEATS=${REPEATS}"
        echo "NUMAPIN=${NUMAPIN[*]}"
        echo
    } | tee "${LOG}"

    "${NUMAPIN[@]}" "${EXE}" "${N}" "${NARRAY}" "${T}" "${REPEATS}" 2>&1 | tee -a "${LOG}"

    {
        echo
        echo "=== Done OMP=${T} ==="
        echo "Date: $(date -Is)"
        echo "Log: ${LOG}"
    } | tee -a "${LOG}"
done

echo
echo "Sweep complete: ${RUN_ROOT}"
