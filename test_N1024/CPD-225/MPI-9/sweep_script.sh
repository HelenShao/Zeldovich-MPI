#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/home/helenshao/InitialConditions/hermitian_3d_matrix_production"
WORKFLOW_ROOT="$PROJECT_ROOT/test_N1024/CPD-225/MPI-9"
PARAM_SOURCE="$WORKFLOW_ROOT/omp-8/param_N1024_CPD_225.par"
PLOT_SOURCE="$WORKFLOW_ROOT/plot_omp_stage_scaling.py"

EXPERIMENT_DIR=""
CHANGE_SUMMARY=""
HYPOTHESIS=""
THREAD_LIST="1 2 4 8 16"
POLL_SECONDS=60
PREPARE_ONLY=0

STATUS_LOG=""
COMMANDS_FILE=""
LAST_QSTAT_OUTPUT=""

usage() {
    cat <<'EOF'
Usage:
  sweep_script.sh --experiment-dir DIR --change-summary "text" [options]

Required:
  --experiment-dir DIR          Destination directory for one experiment run
  --change-summary TEXT         Short description of the code change being tested

Optional:
  --hypothesis TEXT             Expected timing effect to record in README
  --threads "1 2 4 8 16"        Space-separated OMP thread counts
  --poll-seconds N              Queue polling interval in seconds (default: 60)
  --prepare-only                Set up experiment files only; do not submit jobs
  --help                        Show this message
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

abspath() {
    readlink -m -- "$1"
}

timestamp() {
    date "+%Y-%m-%d %H:%M:%S"
}

log() {
    local message="$*"
    local line
    line="[$(timestamp)] $message"
    echo "$line"
    if [[ -n "$STATUS_LOG" ]]; then
        printf "%s\n" "$line" >> "$STATUS_LOG"
    fi
}

append_command() {
    printf "%s\n" "$*" >> "$COMMANDS_FILE"
}

run_dir_for() {
    local threads="$1"
    printf "%s/runs/omp-%s" "$EXPERIMENT_DIR" "$threads"
}

job_name_for() {
    local threads="$1"
    printf "N1024_CPD225_MPI9_OMP%s" "$threads"
}

archived_out_for() {
    local threads="$1"
    printf "%s/%s.out" "$(run_dir_for "$threads")" "$(job_name_for "$threads")"
}

archived_bin_log_for() {
    local threads="$1"
    printf "%s/%s_bin_generation.log" "$(run_dir_for "$threads")" "$(job_name_for "$threads")"
}

archived_reassembly_log_for() {
    local threads="$1"
    printf "%s/%s_reassembly.log" "$(run_dir_for "$threads")" "$(job_name_for "$threads")"
}

archived_rng_dir_for() {
    local threads="$1"
    printf "%s/rng_logs" "$(run_dir_for "$threads")"
}

source_pbs_for() {
    local threads="$1"
    case "$threads" in
        1) printf "%s/N1024_CPD_225_MPI_9.pbs" "$WORKFLOW_ROOT" ;;
        2) printf "%s/omp-2/N1024_CPD225_MPI9_OMP2.pbs" "$WORKFLOW_ROOT" ;;
        4) printf "%s/omp-4/N1024_CPD225_MPI9_OMP4.pbs" "$WORKFLOW_ROOT" ;;
        8) printf "%s/omp-8/N1024_CPD225_MPI9_OMP8.pbs" "$WORKFLOW_ROOT" ;;
        16) printf "%s/omp-16/N1024_CPD225_MPI9_OMP16.pbs" "$WORKFLOW_ROOT" ;;
        *) die "No reference PBS script configured for OMP=$threads" ;;
    esac
}

local_pbs_for() {
    local threads="$1"
    printf "%s/%s" "$(run_dir_for "$threads")" "$(basename "$(source_pbs_for "$threads")")"
}

source_run_dir_for() {
    local threads="$1"
    case "$threads" in
        1) printf "%s" "$WORKFLOW_ROOT" ;;
        *) printf "%s/omp-%s" "$WORKFLOW_ROOT" "$threads" ;;
    esac
}

source_out_for() {
    local threads="$1"
    case "$threads" in
        1) printf "%s/N1024_CPD_225_MPI_9.out" "$(source_run_dir_for "$threads")" ;;
        *) printf "%s/N1024_CPD225_MPI9_OMP%s.out" "$(source_run_dir_for "$threads")" "$threads" ;;
    esac
}

source_bin_log_for() {
    local threads="$1"
    case "$threads" in
        1) printf "%s/N1024_CPD_225_MPI_9_bin_generation.log" "$(source_run_dir_for "$threads")" ;;
        *) printf "%s/N1024_CPD225_MPI9_OMP%s_bin_generation.log" "$(source_run_dir_for "$threads")" "$threads" ;;
    esac
}

source_reassembly_log_for() {
    local threads="$1"
    case "$threads" in
        1) printf "%s/N1024_CPD_225_MPI_9_reassembly.log" "$(source_run_dir_for "$threads")" ;;
        *) printf "%s/N1024_CPD225_MPI9_OMP%s_reassembly.log" "$(source_run_dir_for "$threads")" "$threads" ;;
    esac
}

source_rng_dir_for() {
    local threads="$1"
    printf "%s/rng_logs" "$(source_run_dir_for "$threads")"
}

local_out_for() {
    local threads="$1"
    printf "%s/%s" "$(run_dir_for "$threads")" "$(basename "$(source_out_for "$threads")")"
}

local_bin_log_for() {
    local threads="$1"
    printf "%s/%s" "$(run_dir_for "$threads")" "$(basename "$(source_bin_log_for "$threads")")"
}

local_reassembly_log_for() {
    local threads="$1"
    printf "%s/%s" "$(run_dir_for "$threads")" "$(basename "$(source_reassembly_log_for "$threads")")"
}

local_rng_dir_for() {
    local threads="$1"
    printf "%s/rng_logs" "$(run_dir_for "$threads")"
}

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --experiment-dir)
                [[ $# -ge 2 ]] || die "--experiment-dir requires a value"
                EXPERIMENT_DIR="$2"
                shift 2
                ;;
            --change-summary)
                [[ $# -ge 2 ]] || die "--change-summary requires a value"
                CHANGE_SUMMARY="$2"
                shift 2
                ;;
            --hypothesis)
                [[ $# -ge 2 ]] || die "--hypothesis requires a value"
                HYPOTHESIS="$2"
                shift 2
                ;;
            --threads)
                [[ $# -ge 2 ]] || die "--threads requires a value"
                THREAD_LIST="$2"
                shift 2
                ;;
            --poll-seconds)
                [[ $# -ge 2 ]] || die "--poll-seconds requires a value"
                POLL_SECONDS="$2"
                shift 2
                ;;
            --prepare-only)
                PREPARE_ONLY=1
                shift
                ;;
            --help|-h)
                usage
                exit 0
                ;;
            *)
                die "Unknown argument: $1"
                ;;
        esac
    done

    [[ -n "$EXPERIMENT_DIR" ]] || die "--experiment-dir is required"
    [[ -n "$CHANGE_SUMMARY" ]] || die "--change-summary is required"
    [[ "$POLL_SECONDS" =~ ^[0-9]+$ ]] || die "--poll-seconds must be an integer"
    [[ -f "$PARAM_SOURCE" ]] || die "Missing parameter file: $PARAM_SOURCE"
    [[ -f "$PLOT_SOURCE" ]] || die "Missing plotting script: $PLOT_SOURCE"
    EXPERIMENT_DIR="$(abspath "$EXPERIMENT_DIR")"
}

normalize_threads() {
    local token
    local raw_threads=()
    while read -r token; do
        [[ -n "$token" ]] || continue
        [[ "$token" =~ ^[0-9]+$ ]] || die "Invalid thread count: $token"
        source_pbs_for "$token" >/dev/null
        raw_threads+=("$token")
    done < <(printf "%s\n" $THREAD_LIST)

    [[ ${#raw_threads[@]} -gt 0 ]] || die "At least one OMP thread count is required"
    mapfile -t THREADS < <(printf "%s\n" "${raw_threads[@]}" | sort -n | awk '!seen[$0]++')
}

setup_layout() {
    mkdir -p "$EXPERIMENT_DIR/logs" "$EXPERIMENT_DIR/metadata" "$EXPERIMENT_DIR/scripts" "$EXPERIMENT_DIR/runs" "$EXPERIMENT_DIR/results"
    STATUS_LOG="$EXPERIMENT_DIR/logs/sweep_status.log"
    COMMANDS_FILE="$EXPERIMENT_DIR/metadata/commands.txt"
    touch "$STATUS_LOG" "$COMMANDS_FILE"

    if [[ "${SWEEP_STDOUT_REDIRECTED:-0}" != "1" ]]; then
        export SWEEP_STDOUT_REDIRECTED=1
        exec > >(tee -a "$EXPERIMENT_DIR/logs/sweep_stdout_stderr.txt") 2>&1
    fi

    log "Experiment directory: $EXPERIMENT_DIR"
}

copy_static_inputs() {
    cp "$PARAM_SOURCE" "$EXPERIMENT_DIR/metadata/param_N1024_CPD_225.par"
    cp "$PLOT_SOURCE" "$EXPERIMENT_DIR/scripts/plot_omp_stage_scaling.py"

    local self_path
    self_path="$(readlink -f "${BASH_SOURCE[0]}")"
    if [[ "$(abspath "$self_path")" != "$(abspath "$EXPERIMENT_DIR/scripts/sweep_script.sh")" ]]; then
        cp "$self_path" "$EXPERIMENT_DIR/scripts/sweep_script.sh"
    fi

    chmod +x "$EXPERIMENT_DIR/scripts/sweep_script.sh"
}

prepare_run_directories() {
    local threads source_pbs dest_dir local_pbs source_run_dir source_rel_dir
    for threads in "${THREADS[@]}"; do
        dest_dir="$(run_dir_for "$threads")"
        mkdir -p "$dest_dir"
        source_pbs="$(source_pbs_for "$threads")"
        local_pbs="$(local_pbs_for "$threads")"
        cp "$source_pbs" "$local_pbs"

        source_run_dir="$(source_run_dir_for "$threads")"
        case "$threads" in
            1) source_rel_dir="test_N1024/CPD-225/MPI-9" ;;
            *) source_rel_dir="test_N1024/CPD-225/MPI-9/omp-$threads" ;;
        esac

        python3 - "$local_pbs" "$PROJECT_ROOT" "$source_run_dir" "$source_rel_dir" "$dest_dir" "$PARAM_SOURCE" "$EXPERIMENT_DIR/metadata/param_N1024_CPD_225.par" <<'PY'
from pathlib import Path
import re
import sys

pbs_path = Path(sys.argv[1])
project_root = sys.argv[2]
source_run_dir = sys.argv[3]
source_rel_dir = sys.argv[4]
dest_dir = sys.argv[5]
param_source = sys.argv[6]
param_dest = sys.argv[7]

text = pbs_path.read_text()
text = text.replace(source_run_dir, dest_dir)
text = re.sub(
    rf'(?<!{re.escape(project_root + "/")}){re.escape(source_rel_dir)}',
    dest_dir,
    text,
)
text = text.replace(param_source, param_dest)
text = text.replace(f"{dest_dir}/param_N1024_CPD_225.par", param_dest)
pbs_path.write_text(text)
PY
    done
}

write_change_summary() {
    {
        echo "Experiment directory: $EXPERIMENT_DIR"
        echo "Date: $(date)"
        echo "Code change summary: $CHANGE_SUMMARY"
        if [[ -n "$HYPOTHESIS" ]]; then
            echo "Hypothesis: $HYPOTHESIS"
        fi
    } > "$EXPERIMENT_DIR/README_change_summary.txt"
}

write_manifest() {
    {
        echo "experiment_dir=$EXPERIMENT_DIR"
        echo "created_at=$(date --iso-8601=seconds)"
        echo "code_change_summary=$CHANGE_SUMMARY"
        echo "hypothesis=${HYPOTHESIS:-none}"
        echo "parameter_file=$EXPERIMENT_DIR/metadata/param_N1024_CPD_225.par"
        echo "omp_threads=${THREADS[*]}"
        echo "pbs_mode=experiment_local_copies"
        echo "plotting_script=$EXPERIMENT_DIR/scripts/plot_omp_stage_scaling.py"
        echo "project_root=$PROJECT_ROOT"
        echo "workflow_root=$WORKFLOW_ROOT"
    } > "$EXPERIMENT_DIR/metadata/experiment_manifest.txt"
}

capture_environment() {
    {
        echo "Date: $(date)"
        echo "Hostname: $(hostname)"
        echo "Shell: ${SHELL:-unknown}"
        echo "PWD: $(pwd)"
        echo
        echo "[relevant env vars]"
        echo "OMP_NUM_THREADS=${OMP_NUM_THREADS:-unset}"
        echo "OMP_PLACES=${OMP_PLACES:-unset}"
        echo "OMP_PROC_BIND=${OMP_PROC_BIND:-unset}"
        echo "OMP_MAX_ACTIVE_LEVELS=${OMP_MAX_ACTIVE_LEVELS:-unset}"
        echo
        echo "[module list]"
        if type module >/dev/null 2>&1; then
            module list 2>&1 || true
        else
            echo "module command not available in current shell"
        fi
        echo
        echo "[make check-omp]"
        if make -C "$PROJECT_ROOT" check-omp; then
            true
        else
            echo "make check-omp failed"
        fi
    } > "$EXPERIMENT_DIR/metadata/environment.txt"
}

write_commands_header() {
    {
        echo "# Workflow commands"
        echo "queue_poll_command=qstat -u \$USER"
        echo "plot_command=python3 \"$EXPERIMENT_DIR/scripts/plot_omp_stage_scaling.py\" --base-dir \"$EXPERIMENT_DIR\" --runs-dir \"$EXPERIMENT_DIR/runs\" --results \"$EXPERIMENT_DIR/results\""
        echo
        local threads
        for threads in "${THREADS[@]}"; do
            echo "qsub \"$(local_pbs_for "$threads")\""
        done
    } > "$COMMANDS_FILE"
}

queue_is_empty() {
    local attempt output
    for attempt in 1 2 3; do
        if output="$(qstat -u "$USER" 2>&1)"; then
            LAST_QSTAT_OUTPUT="$output"
            if [[ -z "${output//[[:space:]]/}" ]]; then
                return 0
            fi
            if printf "%s\n" "$output" | grep -qiE "no jobs|no active jobs|no pending jobs"; then
                return 0
            fi
            if printf "%s\n" "$output" | grep -Eq '^[0-9][^[:space:]]*[[:space:]].*[[:space:]]Q[[:space:]]'; then
                return 1
            fi
            return 0
        fi

        LAST_QSTAT_OUTPUT="$output"
        if (( attempt < 3 )); then
            log "qstat -u \$USER failed on attempt $attempt; retrying in 10 seconds"
            sleep 10
        fi
    done

    return 2
}

wait_for_queue_empty() {
    local reason="$1"
    local queue_status
    while true; do
        if queue_is_empty; then
            queue_status=0
        else
            queue_status=$?
        fi

        if [[ $queue_status -eq 0 ]]; then
            log "No queued jobs: $reason"
            return 0
        fi

        if [[ $queue_status -eq 2 ]]; then
            log "Queue polling failed repeatedly"
            printf "%s\n" "$LAST_QSTAT_OUTPUT" >> "$STATUS_LOG"
            return 1
        fi

        log "Queued job still present: $reason"
        if [[ -n "${LAST_QSTAT_OUTPUT//[[:space:]]/}" ]]; then
            printf "%s\n" "$LAST_QSTAT_OUTPUT" >> "$STATUS_LOG"
        fi
        sleep "$POLL_SECONDS"
    done
}

wait_for_file() {
    local path="$1"
    local attempts="$2"
    local delay="$3"
    local try
    for ((try = 1; try <= attempts; try++)); do
        if [[ -f "$path" ]]; then
            return 0
        fi
        sleep "$delay"
    done
    return 1
}

validate_local_run() {
    local threads="$1"
    local out_path bin_log reassembly_log rng_dir rng_log

    out_path="$(local_out_for "$threads")"
    bin_log="$(local_bin_log_for "$threads")"
    reassembly_log="$(local_reassembly_log_for "$threads")"
    rng_dir="$(local_rng_dir_for "$threads")"
    rng_log="$rng_dir/hermitian_rng_debug.log"

    wait_for_file "$out_path" 12 10 || {
        log "Validation failed for OMP=$threads: missing PBS output file $out_path"
        return 1
    }
    wait_for_file "$bin_log" 12 10 || {
        log "Validation failed for OMP=$threads: missing bin log $bin_log"
        return 1
    }
    wait_for_file "$reassembly_log" 12 10 || {
        log "Validation failed for OMP=$threads: missing reassembly log $reassembly_log"
        return 1
    }

    grep -q "TIMING SUMMARY" "$bin_log" || {
        log "Validation failed for OMP=$threads: TIMING SUMMARY missing"
        return 1
    }
    grep -q "Stage 1 (Y-slice generation + 2D FFT)" "$bin_log" || {
        log "Validation failed for OMP=$threads: Stage 1 timing missing"
        return 1
    }
    grep -q "Stage 3 (Communication: Alltoallv)" "$bin_log" || {
        log "Validation failed for OMP=$threads: communication timing missing"
        return 1
    }
    grep -q "Stage 4 (Streaming: Unpack+FFT+Write)" "$bin_log" || {
        log "Validation failed for OMP=$threads: streaming timing missing"
        return 1
    }
    grep -q "Total time (including all stages)" "$bin_log" || {
        log "Validation failed for OMP=$threads: total timing missing"
        return 1
    }

    if grep -qiE "ERROR|MPI_Abort|Traceback|Segmentation fault|core dumped" "$out_path" "$bin_log" "$reassembly_log" "$rng_log" 2>/dev/null; then
        log "Validation failed for OMP=$threads: fatal error string detected"
        return 1
    fi

    log "Validation passed for OMP=$threads"
}

run_plotting() {
    local plot_script="$EXPERIMENT_DIR/scripts/plot_omp_stage_scaling.py"
    log "Running plotting step"
    python3 "$plot_script" --base-dir "$EXPERIMENT_DIR" --runs-dir "$EXPERIMENT_DIR/runs" --results "$EXPERIMENT_DIR/results"
}

run_sweep() {
    local threads local_pbs qsub_output
    for threads in "${THREADS[@]}"; do
        wait_for_queue_empty "before submitting OMP=$threads" || return 1

        local_pbs="$(local_pbs_for "$threads")"
        log "Submitting OMP=$threads with $local_pbs"
        qsub_output="$(qsub "$local_pbs")" || {
            log "qsub failed for OMP=$threads"
            return 1
        }
        log "Submitted OMP=$threads: $qsub_output"

        wait_for_queue_empty "waiting for OMP=$threads to finish" || return 1
        validate_local_run "$threads" || return 1
    done

    run_plotting || {
        log "Plotting failed"
        return 1
    }
}

main() {
    parse_args "$@"
    normalize_threads
    setup_layout
    copy_static_inputs
    prepare_run_directories
    write_change_summary
    write_manifest
    capture_environment
    write_commands_header

    log "Prepared OMP sweep for thread counts: ${THREADS[*]}"
    if [[ $PREPARE_ONLY -eq 1 ]]; then
        log "Prepare-only mode enabled; not submitting jobs"
        exit 0
    fi

    run_sweep
    log "Sweep completed successfully"
}

main "$@"
