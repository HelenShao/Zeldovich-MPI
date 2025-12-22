#!/bin/bash
# ====================================================================================
# SUBMIT ALL TESTS: Options A vs B
# ====================================================================================
# Submits all PBS jobs for testing Options A and B at different scales.
# After jobs complete, run compare_option_outputs.sh to verify results.
#
# Usage:
#   ./submit_all_tests.sh [--small-only] [--medium-only]
#
# Examples:
#   ./submit_all_tests.sh              # Submit all tests
#   ./submit_all_tests.sh --small-only # Only N=256 tests
#   ./submit_all_tests.sh --medium-only # Only N=2048 tests
# ====================================================================================

set -e

SMALL_ONLY=false
MEDIUM_ONLY=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --small-only)
            SMALL_ONLY=true
            shift
            ;;
        --medium-only)
            MEDIUM_ONLY=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--small-only] [--medium-only]"
            exit 1
            ;;
    esac
done

echo "========================================================================"
echo "SUBMITTING TESTS: Options A vs B"
echo "========================================================================"
echo ""

# Get absolute path of script directory
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# Change to project root (one level up from tests/options_ab/)
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

echo "Project root: $PROJECT_ROOT"
echo "Script directory: $SCRIPT_DIR"
echo ""

# PBS script directory (absolute path)
PBS_SCRIPT_DIR="$PROJECT_ROOT/tests/options_ab"

# Submit small-scale tests (N=256)
if [ "$MEDIUM_ONLY" = false ]; then
    echo "Submitting small-scale tests (N=256)..."
    
    JOB_A_256=$(qsub "$PBS_SCRIPT_DIR/test_option_a_N256.pbs")
    echo "  Option A: $JOB_A_256"
    
    JOB_B_256=$(qsub "$PBS_SCRIPT_DIR/test_option_b_N256.pbs")
    echo "  Option B: $JOB_B_256"
    
    echo "  Small-scale jobs submitted"
    echo ""
fi

# Submit medium-scale tests (N=2048)
if [ "$SMALL_ONLY" = false ]; then
    echo "Submitting medium-scale tests (N=2048)..."
    
    JOB_A_2048=$(qsub "$PBS_SCRIPT_DIR/test_option_a_N2048.pbs")
    echo "  Option A: $JOB_A_2048"
    
    JOB_B_2048=$(qsub "$PBS_SCRIPT_DIR/test_option_b_N2048.pbs")
    echo "  Option B: $JOB_B_2048"
    
    echo "  Medium-scale jobs submitted"
    echo ""
fi

echo "========================================================================"
echo "ALL JOBS SUBMITTED"
echo "========================================================================"
echo ""
echo "Job IDs:"
if [ "$MEDIUM_ONLY" = false ]; then
    echo "  N=256 Option A: $JOB_A_256"
    echo "  N=256 Option B: $JOB_B_256"
fi
if [ "$SMALL_ONLY" = false ]; then
    echo "  N=2048 Option A: $JOB_A_2048"
    echo "  N=2048 Option B: $JOB_B_2048"
fi
echo ""
echo "Monitor jobs with:"
echo "  qstat -u $USER"
echo ""
echo "After jobs complete, compare outputs:"
echo "  cd tests/options_ab"
echo "  ./compare_option_outputs.sh output_option_a output_option_b"
echo ""

