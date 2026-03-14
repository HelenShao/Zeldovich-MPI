#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="/home/helenshao/InitialConditions/hermitian_3d_matrix_production"
THREADS=(1 2 4 8 16 32)
PARAM_SOURCE="$SCRIPT_DIR/omp-8/param_N1024_CPD_225.par"

source_pbs() {
    case $1 in
        1) echo "$SCRIPT_DIR/N1024_CPD_225_MPI_9.pbs" ;;
        *) echo "$SCRIPT_DIR/omp-$1/N1024_CPD225_MPI9_OMP${1}.pbs" ;;
    esac
}

source_absdir() {
    case $1 in
        1) echo "$SCRIPT_DIR" ;;
        *) echo "$SCRIPT_DIR/omp-$1" ;;
    esac
}

source_reldir() {
    case $1 in
        1) echo "test_N1024/CPD-225/MPI-9" ;;
        *) echo "test_N1024/CPD-225/MPI-9/omp-$1" ;;
    esac
}

rewrite_pbs() {
    local local_pbs="$1" src_abs="$2" src_rel="$3" dest_dir="$4" param_dest="$5"

    if [[ "$dest_dir" == "${PROJECT_ROOT}/"* ]]; then
        local dest_rel="${dest_dir#${PROJECT_ROOT}/}"
        sed -i "s|${src_rel}|${dest_rel}|g" "$local_pbs"
    else
        sed -i "s|${src_abs}|${dest_dir}|g" "$local_pbs"
        sed -i "s|${src_rel}|${dest_dir}|g" "$local_pbs"
    fi
    sed -i 's|PARAM_FILE="[^"]*"|PARAM_FILE="'"${param_dest}"'"|' "$local_pbs"
}

cmd_setup() {
    local dir="" summary=""
    while [[ $# -gt 0 ]]; do
        case $1 in
            --dir)     dir="$2"; shift 2 ;;
            --summary) summary="$2"; shift 2 ;;
            *)         echo "Unknown option: $1" >&2; exit 1 ;;
        esac
    done
    [[ -z "$dir" ]] && { echo "Usage: $0 setup --dir DIR [--summary TEXT]" >&2; exit 1; }
    [[ "$dir" != /* ]] && dir="$SCRIPT_DIR/$dir"

    mkdir -p "$dir"/{results,metadata,scripts}
    for t in "${THREADS[@]}"; do mkdir -p "$dir/runs/omp-$t"; done

    cp "$PARAM_SOURCE" "$dir/metadata/param_N1024_CPD_225.par"
    cp "$SCRIPT_DIR/plot_omp_stage_scaling.py" "$dir/scripts/"
    [[ -n "$summary" ]] && printf '%s\n' "$summary" > "$dir/README_change_summary.txt"

    local param_dest="$dir/metadata/param_N1024_CPD_225.par"
    for t in "${THREADS[@]}"; do
        local src dest_dir local_pbs
        src="$(source_pbs "$t")"
        dest_dir="$dir/runs/omp-$t"
        local_pbs="$dest_dir/$(basename "$src")"

        cp "$src" "$local_pbs"
        rewrite_pbs "$local_pbs" "$(source_absdir "$t")" "$(source_reldir "$t")" "$dest_dir" "$param_dest"
    done

    echo ""
    echo "Experiment directory: $dir"
    echo ""
    echo "Submit jobs (wait for Q to clear between each):"
    for t in "${THREADS[@]}"; do
        echo "  qsub $dir/runs/omp-$t/$(basename "$(source_pbs "$t")")"
    done
    echo ""
    echo "After all jobs finish:"
    echo "  $0 analyze --dir $dir"
}

cmd_analyze() {
    local dir=""
    while [[ $# -gt 0 ]]; do
        case $1 in
            --dir) dir="$2"; shift 2 ;;
            *)     echo "Unknown option: $1" >&2; exit 1 ;;
        esac
    done
    [[ -z "$dir" ]] && { echo "Usage: $0 analyze --dir DIR" >&2; exit 1; }
    [[ "$dir" != /* ]] && dir="$SCRIPT_DIR/$dir"

    python3 "$dir/scripts/plot_omp_stage_scaling.py" \
        --runs-dir "$dir/runs" \
        --results "$dir/results"
}

case "${1:-}" in
    setup)   shift; cmd_setup "$@" ;;
    analyze) shift; cmd_analyze "$@" ;;
    *)       echo "Usage: $0 {setup|analyze} [options]" >&2; exit 1 ;;
esac
