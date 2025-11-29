#!/bin/bash
# test_exact_division.sh - Test exact division cases (Category 1)
# Note: Run with ENABLE_ABACUS_VALIDATION=0

set -e

BINARY=$1
RESULTS_DIR=$2
CATEGORY_DIR="$RESULTS_DIR/exact_division"
VALIDATE_SCRIPT="$(dirname "$0")/validate_hermitian.sh"

mkdir -p "$CATEGORY_DIR"

echo "=========================================="
echo "Exact Division Tests"
echo "=========================================="

# Test E1: N=1024, Ranks=16
echo "Test E1: N=1024, Ranks=16 (Expected: 4×4 grid, 256×256 cells/rank)"
mpirun -np 16 "$BINARY" 1024 > "$CATEGORY_DIR/E1.log" 2>&1 || true
grep -q "Exact division\|Abacus-compatible" "$CATEGORY_DIR/E1.log" && echo "✓ E1: Info message found" || echo "✗ E1: Info message missing"
grep -q "WARNING.*not divisible" "$CATEGORY_DIR/E1.log" && echo "✗ E1: Unexpected warning" || echo "✓ E1: No warnings"

# Test E2: N=512, Ranks=8
echo "Test E2: N=512, Ranks=8 (Expected: 2×4 grid, 256×128 cells/rank)"
mpirun -np 8 "$BINARY" 512 > "$CATEGORY_DIR/E2.log" 2>&1 || true
grep -q "Exact division\|Abacus-compatible" "$CATEGORY_DIR/E2.log" && echo "✓ E2: Info message found" || echo "✗ E2: Info message missing"

# Test E3: N=256, Ranks=4
echo "Test E3: N=256, Ranks=4 (Expected: 2×2 grid, 128×128 cells/rank)"
mpirun -np 4 "$BINARY" 256 > "$CATEGORY_DIR/E3.log" 2>&1 || true
grep -q "Exact division\|Abacus-compatible" "$CATEGORY_DIR/E3.log" && echo "✓ E3: Info message found" || echo "✗ E3: Info message missing"

# Test E4: N=2048, Ranks=32 - SKIPPED (N > 1024)
# echo "Test E4: N=2048, Ranks=32 (Expected: 4×8 grid, 512×256 cells/rank)"
# mpirun -np 32 "$BINARY" 2048 > "$CATEGORY_DIR/E4.log" 2>&1 || true
# grep -q "Exact division\|Abacus-compatible" "$CATEGORY_DIR/E4.log" && echo "✓ E4: Info message found" || echo "✗ E4: Info message missing"
echo "Test E4: SKIPPED (N=2048 > 1024)"

# Test E5: N=4096, Ranks=64 - SKIPPED (N > 1024)
# echo "Test E5: N=4096, Ranks=64 (Expected: 8×8 grid, 512×512 cells/rank)"
# mpirun -np 64 "$BINARY" 4096 > "$CATEGORY_DIR/E5.log" 2>&1 || true
# grep -q "Exact division\|Abacus-compatible" "$CATEGORY_DIR/E5.log" && echo "✓ E5: Info message found" || echo "✗ E5: Info message missing"
echo "Test E5: SKIPPED (N=4096 > 1024)"

# Test E6: N=1000, Ranks=16
echo "Test E6: N=1000, Ranks=16 (Expected: 4×4 grid, 250×250 cells/rank)"
mpirun -np 16 "$BINARY" 1000 > "$CATEGORY_DIR/E6.log" 2>&1 || true
grep -q "Exact division\|Abacus-compatible" "$CATEGORY_DIR/E6.log" && echo "✓ E6: Info message found" || echo "✗ E6: Info message missing"

echo ""
echo "Validating Hermitian symmetry (purely real output)..."
for test in E1 E2 E3 E6; do
    if [ -f "$CATEGORY_DIR/${test}.log" ]; then
        echo "  Checking $test..."
        bash "$VALIDATE_SCRIPT" "$CATEGORY_DIR/${test}.log" | tail -5
    fi
done
echo "  E4, E5: SKIPPED (N > 1024)"

echo "=========================================="
echo "Exact Division Tests Complete"
echo "Results: $CATEGORY_DIR"
echo "=========================================="

