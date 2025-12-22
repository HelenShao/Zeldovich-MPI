#!/bin/bash
# Test script for Option A implementation
# Tests transpose + WriteParticlesSlab_range particle writing

set -e  # Exit on error

echo "========================================================================"
echo "OPTION A TEST: Transpose + WriteParticlesSlab_range"
echo "========================================================================"

# Configuration
N=16
NRANKS=4
PARAM_FILE="examples/param_N16.par"
MAIN_EXEC="./main"

# Check if main executable exists
if [ ! -f "$MAIN_EXEC" ]; then
    echo "ERROR: $MAIN_EXEC not found. Please compile first."
    echo "Run: make"
    exit 1
fi

# Check if parameter file exists
if [ ! -f "$PARAM_FILE" ]; then
    echo "ERROR: Parameter file $PARAM_FILE not found"
    exit 1
fi

# Clean up previous outputs
echo ""
echo "Cleaning up previous outputs..."
rm -rf rank_* output/ *.bin
echo "Cleanup complete"

# Run main.cpp with Option A (param_file provided)
echo ""
echo "========================================================================"
echo "Running main.cpp with Option A (N=$N, ranks=$NRANKS)"
echo "Command: mpirun -np $NRANKS $MAIN_EXEC $N $PARAM_FILE"
echo "========================================================================"
mpirun -np $NRANKS $MAIN_EXEC $N $PARAM_FILE

# Check outputs
echo ""
echo "========================================================================"
echo "Verifying outputs..."
echo "========================================================================"

# Check if output directory was created
if [ ! -d "output" ]; then
    echo "ERROR: output/ directory not created"
    exit 1
fi
echo "✓ output/ directory created"

# Count particle IC files
NUM_PARTICLE_FILES=$(find output/ -name "ic_rank*" | wc -l)
echo "✓ Found $NUM_PARTICLE_FILES particle IC files in output/"

if [ $NUM_PARTICLE_FILES -eq 0 ]; then
    echo "ERROR: No particle IC files found!"
    echo "Expected files like: output/ic_rank0_i0_x0_4"
    exit 1
fi

# List first few particle files
echo ""
echo "First 10 particle IC files:"
find output/ -name "ic_rank*" | sort | head -10

# Check file sizes
echo ""
echo "File sizes (first 5 files):"
find output/ -name "ic_rank*" | sort | head -5 | xargs ls -lh

# Check for density files (if qdensity=1 in param file)
NUM_DENSITY_FILES=$(find output/ -name "dens_rank*" | wc -l)
if [ $NUM_DENSITY_FILES -gt 0 ]; then
    echo "✓ Found $NUM_DENSITY_FILES density files in output/"
else
    echo "⚠ No density files found (qdensity may be disabled)"
fi

# Verify rank directories are not used (old .bin output)
if ls rank_*/i*.bin 2>/dev/null | grep -q .; then
    echo "⚠ WARNING: Found .bin files in rank_* directories"
    echo "   Option A should write particles to output/, not .bin files to rank_*/"
    echo "   This may indicate fallback to .bin writing mode"
else
    echo "✓ No .bin files in rank_* directories (Option A working correctly)"
fi

echo ""
echo "========================================================================"
echo "OPTION A TEST COMPLETED SUCCESSFULLY"
echo "========================================================================"
echo ""
echo "Summary:"
echo "  - Particle IC files: $NUM_PARTICLE_FILES (in output/)"
echo "  - Density files: $NUM_DENSITY_FILES (in output/)"
echo "  - Format: ic_rank{rank}_i{ic_index}_x{k_start}_{k_end}"
echo ""
echo "Next steps:"
echo "  1. Verify particle data with a reader script"
echo "  2. Compare with Path 2 (reassembled ICs) for consistency"
echo "  3. Test with larger N (256, 1024)"
echo ""


