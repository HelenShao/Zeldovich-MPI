#!/bin/bash
# Compare final output matrices between N=256 and N=512
# This tests RNG consistency by comparing real-space output in overlapping region

cd /home/helenshao/InitialConditions/hermitian_3d_matrix_production/conj-test/rng-consistency

echo "=========================================="
echo "Final Output Matrix Comparison"
echo "Branch: zeldovich-conj"
echo "Comparing N=256 vs N=512 final outputs"
echo "=========================================="
echo ""

# Determine output directories (where rank_*/z*_slab_N*.bin files are)
# These should be in the working directory where the jobs ran
N256_DIR="/home/helenshao/InitialConditions/hermitian_3d_matrix_production"
N512_DIR="/home/helenshao/InitialConditions/hermitian_3d_matrix_production"

# Check if files exist (use path pattern)
N256_FILES=$(find "$N256_DIR" -path "*/rank_*/z*_slab_N256.bin" 2>/dev/null | wc -l)
N512_FILES=$(find "$N512_DIR" -path "*/rank_*/z*_slab_N512.bin" 2>/dev/null | wc -l)

if [ "$N256_FILES" -eq 0 ]; then
    echo "ERROR: No N=256 output files found!"
    echo "Looking in: $N256_DIR/rank_*/z*_slab_N256.bin"
    echo "Make sure test_rng_N256.pbs completed and file writes are enabled"
    exit 1
fi

if [ "$N512_FILES" -eq 0 ]; then
    echo "ERROR: No N=512 output files found!"
    echo "Looking in: $N512_DIR/rank_*/z*_slab_N512.bin"
    echo "Make sure test_rng_N512.pbs completed and file writes are enabled"
    exit 1
fi

echo "Found $N256_FILES N=256 files and $N512_FILES N=512 files"
echo ""

# Run comparison
python3 compare_final_outputs.py \
    --N256-dir "$N256_DIR" \
    --N512-dir "$N512_DIR" \
    --z-start 0 \
    --z-end 256 \
    --N-overlap 256 \
    > compare_final_matrices.txt 2>&1

echo ""
echo "Comparison complete!"
echo "Results saved to: compare_final_matrices.txt"
echo ""
echo "Summary:"
tail -30 compare_final_matrices.txt

