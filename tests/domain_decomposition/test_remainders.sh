#!/bin/bash
# test_remainders.sh - Test remainder cases (Category 2)

set -e

BINARY=$1
RESULTS_DIR=$2
CATEGORY_DIR="$RESULTS_DIR/remainders"

mkdir -p "$CATEGORY_DIR"

echo "=========================================="
echo "Remainder Tests"
echo "=========================================="

# Test R1: N=1001, Ranks=16 (remainder 1 in both)
echo "Test R1: N=1001, Ranks=16 (Expected: 4×4 grid, remainder 1)"
mpirun -np 16 "$BINARY" 1001 > "$CATEGORY_DIR/R1.log" 2>&1 || true
grep -q "WARNING.*not divisible.*grid_x=4" "$CATEGORY_DIR/R1.log" && echo "✓ R1: X remainder warning found" || echo "✗ R1: X remainder warning missing"
grep -q "WARNING.*not divisible.*grid_z=4" "$CATEGORY_DIR/R1.log" && echo "✓ R1: Z remainder warning found" || echo "✗ R1: Z remainder warning missing"
grep -q "Remainder: 1" "$CATEGORY_DIR/R1.log" && echo "✓ R1: Remainder value correct" || echo "✗ R1: Remainder value incorrect"

# Test R2: N=1025, Ranks=16 (remainder 1 in both)
echo "Test R2: N=1025, Ranks=16 (Expected: 4×4 grid, remainder 1)"
mpirun -np 16 "$BINARY" 1025 > "$CATEGORY_DIR/R2.log" 2>&1 || true
grep -q "WARNING.*not divisible" "$CATEGORY_DIR/R2.log" && echo "✓ R2: Warnings found" || echo "✗ R2: Warnings missing"

# Test R3: N=1003, Ranks=16 (remainder 3 in both)
echo "Test R3: N=1003, Ranks=16 (Expected: 4×4 grid, remainder 3)"
mpirun -np 16 "$BINARY" 1003 > "$CATEGORY_DIR/R3.log" 2>&1 || true
grep -q "Remainder: 3" "$CATEGORY_DIR/R3.log" && echo "✓ R3: Remainder 3 detected" || echo "✗ R3: Remainder 3 not detected"

# Test R4: N=1024, Ranks=13 (prime ranks, remainder 10 in Z)
echo "Test R4: N=1024, Ranks=13 (Expected: 1×13 grid, remainder 10 in Z)"
mpirun -np 13 "$BINARY" 1024 > "$CATEGORY_DIR/R4.log" 2>&1 || true
grep -q "WARNING.*not divisible.*grid_z=13" "$CATEGORY_DIR/R4.log" && echo "✓ R4: Z remainder warning found" || echo "✗ R4: Z remainder warning missing"
grep -q "Remainder: 10" "$CATEGORY_DIR/R4.log" && echo "✓ R4: Remainder 10 detected" || echo "✗ R4: Remainder 10 not detected"

# Test R5: N=1000, Ranks=13 (prime ranks, remainder 12 in Z)
echo "Test R5: N=1000, Ranks=13 (Expected: 1×13 grid, remainder 12 in Z)"
mpirun -np 13 "$BINARY" 1000 > "$CATEGORY_DIR/R5.log" 2>&1 || true
grep -q "Remainder: 12" "$CATEGORY_DIR/R5.log" && echo "✓ R5: Remainder 12 detected" || echo "✗ R5: Remainder 12 not detected"

# Test R6: N=1001, Ranks=13 (prime ranks, no remainder in Z)
echo "Test R6: N=1001, Ranks=13 (Expected: 1×13 grid, no remainder in Z)"
mpirun -np 13 "$BINARY" 1001 > "$CATEGORY_DIR/R6.log" 2>&1 || true
grep -q "WARNING.*not divisible.*grid_z=13" "$CATEGORY_DIR/R6.log" && echo "✗ R6: Unexpected Z warning" || echo "✓ R6: No Z remainder (1001/13=77 exact)"

echo "=========================================="
echo "Remainder Tests Complete"
echo "Results: $CATEGORY_DIR"
echo "=========================================="

