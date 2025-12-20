#!/bin/bash
# Compare RNG consistency between N=256 and N=512 outputs
# Usage: ./compare_rng.sh

cd /home/helenshao/InitialConditions/hermitian_3d_matrix_production/conj-test/rng-consistency

echo "=========================================="
echo "RNG Consistency Comparison"
echo "Branch: zeldovich-conj"
echo "Comparing N=256 vs N=512"
echo "=========================================="
echo ""

# Check if output files exist
if [ ! -f "test_rng_N256.out" ]; then
    echo "ERROR: test_rng_N256.out not found!"
    echo "Run test_rng_N256.pbs first"
    exit 1
fi

if [ ! -f "test_rng_N512.out" ]; then
    echo "ERROR: test_rng_N512.out not found!"
    echo "Run test_rng_N512.pbs first"
    exit 1
fi

# Run comparison
python3 compare_rng_consistency.py test_rng_N256.out test_rng_N512.out > compare_256_vs_512.txt 2>&1

echo ""
echo "Comparison complete!"
echo "Results saved to: compare_256_vs_512.txt"
echo ""
echo "Summary:"
tail -20 compare_256_vs_512.txt

