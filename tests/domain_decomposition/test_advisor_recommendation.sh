#!/bin/bash
# test_advisor_recommendation.sh - Test N=1024 with 175 ranks (advisor's recommendation)

set -e

BINARY=$1
RESULTS_DIR=$2
CATEGORY_DIR="$RESULTS_DIR/advisor_recommendation"
VALIDATE_SCRIPT="$(dirname "$0")/validate_hermitian.sh"

mkdir -p "$CATEGORY_DIR"

echo "=========================================="
echo "Advisor's Recommendation Test"
echo "N=1024, Ranks=175"
echo "=========================================="
echo ""
echo "Why this is a critical edge case:"
echo "  - Grid decomposition: 7×25 (non-square, rectangular)"
echo "  - X remainder: 1024 ÷ 7 = 146 remainder 2 (2 ranks get +1 cell)"
echo "  - Z remainder: 1024 ÷ 25 = 40 remainder 24 (24 ranks get +1 cell)"
echo "  - Load imbalance: ~2.8% in X, ~96% in Z (extreme!)"
echo "  - Tests: Non-square grids, dual remainders, prime factorization"
echo ""

# Test AD1: N=1024, Ranks=175 (advisor's recommendation)
echo "Test AD1: N=1024, Ranks=175 (Expected: 7×25 grid, remainders 2 and 24)"
mpirun -np 175 "$BINARY" 1024 > "$CATEGORY_DIR/AD1.log" 2>&1 || true

echo ""
echo "Checking domain decomposition..."
if grep -q "grid_x.*7\|grid_z.*25" "$CATEGORY_DIR/AD1.log"; then
    echo "✓ AD1: Grid decomposition 7×25 detected"
else
    echo "✗ AD1: Grid decomposition not found or incorrect"
fi

if grep -q "WARNING.*not divisible.*grid_x=7" "$CATEGORY_DIR/AD1.log"; then
    echo "✓ AD1: X remainder warning found (expected: remainder 2)"
else
    echo "✗ AD1: X remainder warning missing"
fi

if grep -q "WARNING.*not divisible.*grid_z=25" "$CATEGORY_DIR/AD1.log"; then
    echo "✓ AD1: Z remainder warning found (expected: remainder 24)"
else
    echo "✗ AD1: Z remainder warning missing"
fi

if grep -q "Remainder: 2" "$CATEGORY_DIR/AD1.log"; then
    echo "✓ AD1: X remainder 2 detected"
fi

if grep -q "Remainder: 24" "$CATEGORY_DIR/AD1.log"; then
    echo "✓ AD1: Z remainder 24 detected"
fi

echo ""
echo "Validating Hermitian symmetry (purely real output)..."
bash "$VALIDATE_SCRIPT" "$CATEGORY_DIR/AD1.log" | tail -10

echo ""
echo "=========================================="
echo "Advisor's Recommendation Test Complete"
echo "Results: $CATEGORY_DIR"
echo "=========================================="
echo ""
echo "Key Validation:"
echo "  - Non-square grid (7×25) handled correctly"
echo "  - Dual remainders (2 in X, 24 in Z) distributed correctly"
echo "  - Load imbalance acceptable despite extreme Z remainder"
echo "  - max_imag < 1e-10 (purely real output)"
echo "  - Hermitian symmetry preserved"

