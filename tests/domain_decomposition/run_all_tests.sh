#!/bin/bash
# run_all_tests.sh - Master test runner for domain decomposition tests
#
# Usage: ./run_all_tests.sh [category]
#   category: exact|remainder|prime|power2|nonsquare|small|abacus|edge|all
#   If omitted, runs all tests

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="$SCRIPT_DIR/results"
BINARY="$SCRIPT_DIR/../../hermitian_3d_matrix"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Create results directory
mkdir -p "$RESULTS_DIR"

# Function to run a test script
run_test_script() {
    local script=$1
    local category=$2
    
    if [ -f "$script" ]; then
        echo -e "${GREEN}Running $category tests...${NC}"
        bash "$script" "$BINARY" "$RESULTS_DIR"
        echo -e "${GREEN}✓ $category tests complete${NC}\n"
    else
        echo -e "${YELLOW}Warning: $script not found, skipping${NC}\n"
    fi
}

# Parse category argument
CATEGORY=${1:-all}

echo "=========================================="
echo "Domain Decomposition Test Suite"
echo "=========================================="
echo "Binary: $BINARY"
echo "Results: $RESULTS_DIR"
echo "Category: $CATEGORY"
echo ""
echo "IMPORTANT: Tests assume ENABLE_ABACUS_VALIDATION=0"
echo "           All outputs must be purely real (Hermitian symmetry)"
echo "           SKIPPED: Abacus tests (N=6075), tests with N > 1024,"
echo "                    and Advisor's test (N=1024, 175 ranks - needs many nodes)"
echo "=========================================="
echo ""

# Check if binary exists
if [ ! -f "$BINARY" ]; then
    echo -e "${RED}Error: Binary not found at $BINARY${NC}"
    echo "Please build the code first: cd ../.. && make"
    exit 1
fi

# Run tests based on category
case "$CATEGORY" in
    exact|all)
        run_test_script "$SCRIPT_DIR/test_exact_division.sh" "Exact Division"
        ;;
    remainder|all)
        run_test_script "$SCRIPT_DIR/test_remainders.sh" "Remainders"
        ;;
    prime|all)
        run_test_script "$SCRIPT_DIR/test_prime_ranks.sh" "Prime Ranks"
        ;;
    power2|all)
        run_test_script "$SCRIPT_DIR/test_power_of_2.sh" "Power-of-2"
        ;;
    nonsquare|all)
        run_test_script "$SCRIPT_DIR/test_non_square.sh" "Non-Square"
        ;;
    small|all)
        run_test_script "$SCRIPT_DIR/test_small_grids.sh" "Small Grids"
        ;;
    abacus|all)
        echo -e "${YELLOW}Skipping Abacus tests (N=6075 too large)${NC}\n"
        # run_test_script "$SCRIPT_DIR/test_abacus.sh" "Abacus"
        ;;
    nonpower2|all)
        run_test_script "$SCRIPT_DIR/test_non_power_of_2.sh" "Non-Power-of-2 (N=175)"
        ;;
    advisor|all)
        echo -e "${YELLOW}Skipping Advisor's Recommendation test (N=1024, 175 ranks - requires many MPI nodes)${NC}\n"
        # run_test_script "$SCRIPT_DIR/test_advisor_recommendation.sh" "Advisor's Recommendation (N=1024, 175 ranks)"
        ;;
    edge|all)
        run_test_script "$SCRIPT_DIR/test_edge_cases.sh" "Edge Cases"
        ;;
    *)
        echo -e "${RED}Error: Unknown category '$CATEGORY'${NC}"
        echo "Valid categories: exact, remainder, prime, power2, nonsquare, small, abacus, nonpower2, advisor, edge, all"
        exit 1
        ;;
esac

echo "=========================================="
echo "All tests complete!"
echo "Results saved in: $RESULTS_DIR"
echo "=========================================="

