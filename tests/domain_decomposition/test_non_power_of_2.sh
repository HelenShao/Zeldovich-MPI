#!/bin/bash
# test_non_power_of_2.sh - Test N=175 (advisor's recommendation)

set -e

BINARY=$1
RESULTS_DIR=$2
CATEGORY_DIR="$RESULTS_DIR/non_power_of_2"
VALIDATE_SCRIPT="$(dirname "$0")/validate_hermitian.sh"

mkdir -p "$CATEGORY_DIR"

echo "=========================================="
echo "Non-Power-of-2 Tests (N=175)"
echo "Advisor's Recommendation: N³=1024³, CPD=175³"
echo "=========================================="
echo ""
echo "Why N=175 is important:"
echo "  - Factorization: 175 = 5² × 7 (small primes, good for FFT)"
echo "  - Not power-of-2: Tests mixed-radix FFT"
echo "  - Similar to Abacus: 6075 = 3⁵ × 5² (small prime factors)"
echo "  - Tests remainder handling in domain decomposition"
echo ""

# Test NP1: N=175, Ranks=16 (remainder 3 in both)
echo "Test NP1: N=175, Ranks=16 (Expected: 4×4 grid, remainder 3)"
mpirun -np 16 "$BINARY" 175 > "$CATEGORY_DIR/NP1.log" 2>&1 || true
grep -q "WARNING.*not divisible.*grid_x=4" "$CATEGORY_DIR/NP1.log" && echo "✓ NP1: X remainder warning found" || echo "✗ NP1: X remainder warning missing"
grep -q "Remainder: 3" "$CATEGORY_DIR/NP1.log" && echo "✓ NP1: Remainder 3 detected" || echo "✗ NP1: Remainder 3 not detected"
bash "$VALIDATE_SCRIPT" "$CATEGORY_DIR/NP1.log" | tail -5

# Test NP2: N=175, Ranks=25 (perfect square, exact division)
echo ""
echo "Test NP2: N=175, Ranks=25 (Expected: 5×5 grid, exact division)"
mpirun -np 25 "$BINARY" 175 > "$CATEGORY_DIR/NP2.log" 2>&1 || true
grep -q "Exact division\|Abacus-compatible" "$CATEGORY_DIR/NP2.log" && echo "✓ NP2: Exact division detected" || echo "✗ NP2: Exact division not detected"
grep -q "WARNING.*not divisible" "$CATEGORY_DIR/NP2.log" && echo "✗ NP2: Unexpected warning" || echo "✓ NP2: No warnings"
bash "$VALIDATE_SCRIPT" "$CATEGORY_DIR/NP2.log" | tail -5

# Test NP3: N=175, Ranks=7 (prime factor, exact division)
echo ""
echo "Test NP3: N=175, Ranks=7 (Expected: 1×7 grid, exact division)"
mpirun -np 7 "$BINARY" 175 > "$CATEGORY_DIR/NP3.log" 2>&1 || true
grep -q "Exact division\|Abacus-compatible" "$CATEGORY_DIR/NP3.log" && echo "✓ NP3: Exact division detected" || echo "✗ NP3: Exact division not detected"
bash "$VALIDATE_SCRIPT" "$CATEGORY_DIR/NP3.log" | tail -5

# Test NP4: N=175, Ranks=35 (matches factorization: 5×7)
echo ""
echo "Test NP4: N=175, Ranks=35 (Expected: 5×7 grid, exact division)"
mpirun -np 35 "$BINARY" 175 > "$CATEGORY_DIR/NP4.log" 2>&1 || true
grep -q "Exact division\|Abacus-compatible" "$CATEGORY_DIR/NP4.log" && echo "✓ NP4: Exact division detected" || echo "✗ NP4: Exact division not detected"
bash "$VALIDATE_SCRIPT" "$CATEGORY_DIR/NP4.log" | tail -5

# Test NP5: N=175, Ranks=8 (remainder 3 in both)
echo ""
echo "Test NP5: N=175, Ranks=8 (Expected: 2×4 grid, remainder 3)"
mpirun -np 8 "$BINARY" 175 > "$CATEGORY_DIR/NP5.log" 2>&1 || true
grep -q "WARNING.*not divisible" "$CATEGORY_DIR/NP5.log" && echo "✓ NP5: Remainder warnings found" || echo "✗ NP5: Remainder warnings missing"
bash "$VALIDATE_SCRIPT" "$CATEGORY_DIR/NP5.log" | tail -5

# Test NP6: N=175, Ranks=13 (prime ranks, remainder 2 in Z)
echo ""
echo "Test NP6: N=175, Ranks=13 (Expected: 1×13 grid, remainder 2 in Z)"
mpirun -np 13 "$BINARY" 175 > "$CATEGORY_DIR/NP6.log" 2>&1 || true
grep -q "WARNING.*not divisible.*grid_z=13" "$CATEGORY_DIR/NP6.log" && echo "✓ NP6: Z remainder warning found" || echo "✗ NP6: Z remainder warning missing"
bash "$VALIDATE_SCRIPT" "$CATEGORY_DIR/NP6.log" | tail -5

echo ""
echo "=========================================="
echo "Non-Power-of-2 Tests Complete"
echo "Results: $CATEGORY_DIR"
echo "=========================================="
echo ""
echo "Key Validation:"
echo "  - All tests should show max_imag < 1e-10 (purely real output)"
echo "  - Hermitian symmetry must be preserved despite non-power-of-2 size"
echo "  - FFT performance may differ from power-of-2 cases"


