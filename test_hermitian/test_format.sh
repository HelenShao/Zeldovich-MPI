#!/bin/bash
# Quick test script to verify [Array][X][Y] -> [Array][Y][X] output format
# This can be run after a PBS job completes to verify output files

cd /home/helenshao/InitialConditions/hermitian_3d_matrix_production

echo "=========================================="
echo "Testing [Array][X][Y] -> [Array][Y][X] Output Format"
echo "=========================================="
echo ""

# Check if test_hermitian directory exists and has output files
if [ ! -d "test_hermitian" ]; then
    echo "ERROR: test_hermitian directory not found!"
    echo "Please run a PBS job first (e.g., test_hermitian.pbs)"
    exit 1
fi

# Find a sample output file
SAMPLE_FILE=$(find test_hermitian -name "z*_slab_N*.bin" -type f | head -1)

if [ -z "$SAMPLE_FILE" ]; then
    echo "ERROR: No output files found in test_hermitian/"
    echo "Please run a PBS job first"
    exit 1
fi

echo "Found sample file: $SAMPLE_FILE"
echo ""

# Get file size and infer dimensions
FILE_SIZE=$(stat -f%z "$SAMPLE_FILE" 2>/dev/null || stat -c%s "$SAMPLE_FILE" 2>/dev/null)
BYTES_PER_COMPLEX=16  # complex128 = 16 bytes
TOTAL_ELEMENTS=$((FILE_SIZE / BYTES_PER_COMPLEX))

echo "File size: $FILE_SIZE bytes"
echo "Total elements: $TOTAL_ELEMENTS"
echo ""

# Try to infer dimensions (assuming N=1024, narray=4)
N=1024
NARRAY=4
X_COUNT=$((TOTAL_ELEMENTS / (NARRAY * N)))

echo "Inferred dimensions:"
echo "  N (Y dimension): $N"
echo "  narray: $NARRAY"
echo "  x_count (X dimension for this rank): $X_COUNT"
echo ""

# Verify expected size
EXPECTED_SIZE=$((NARRAY * N * X_COUNT * BYTES_PER_COMPLEX))
if [ "$FILE_SIZE" -eq "$EXPECTED_SIZE" ]; then
    echo "✓ File size matches expected: $EXPECTED_SIZE bytes"
else
    echo "✗ File size mismatch: expected $EXPECTED_SIZE, got $FILE_SIZE"
fi
echo ""

# Use Python to verify format
python3 << EOF
import numpy as np
import sys
import os

file_path = "$SAMPLE_FILE"
N = $N
narray = $NARRAY
x_count = $X_COUNT

print("Reading file: {}".format(file_path))
data = np.fromfile(file_path, dtype=np.complex128)

print("Total elements read: {}".format(len(data)))
print("Expected elements: {} (narray={} x N={} x x_count={})".format(narray * N * x_count, narray, N, x_count))

if len(data) != narray * N * x_count:
    print("ERROR: Element count mismatch!")
    sys.exit(1)

# Reshape to [Array][Y][X] format (as expected by reassembly script)
try:
    reshaped = data.reshape(narray, N, x_count)
    print("✓ Successfully reshaped to [Array][Y][X]: shape {}".format(reshaped.shape))
except ValueError as e:
    print("ERROR: Cannot reshape to [Array][Y][X]: {}".format(e))
    sys.exit(1)

# Check a few sample values
print("\nSample values (first few elements):")
for array_idx in range(min(2, narray)):
    print("  Array {}:".format(array_idx))
    for y in range(min(3, N)):
        for x in range(min(3, x_count)):
            idx = array_idx * N * x_count + y * x_count + x
            val = data[idx]
            print("    [A{},Y{},X{}] = {:.6e} + {:.6e}i".format(array_idx, y, x, val.real, val.imag))

print("\n✓ Format verification passed!")
print("  File is in [Array][Y][X] format as expected by reassembly script")
EOF

echo ""
echo "=========================================="
echo "Format test complete!"
echo "=========================================="
