#!/bin/bash
# Script to run particle IC comparison for N=8 test (steps 8-11 from workflow)
# This script automates the particle IC comparison process

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=========================================="
echo "Particle IC Comparison (N=8, No PLT)"
echo "=========================================="
echo ""

# ====================================================================================
# STEP 8: Reassemble .bin files and write particle ICs (hermitian)
# ====================================================================================

echo "STEP 8: Reassembling .bin files and writing particle ICs (hermitian)"
echo "----------------------------------------------------------------------"

# Check if .bin files exist
BIN_COUNT=$(find hermitian_output/rank_0 -name "*.bin" -type f 2>/dev/null | wc -l)

if [ "$BIN_COUNT" -eq 0 ]; then
    echo "WARNING: No .bin files found in hermitian_output/rank_0/"
    echo "The trace test was run with SKIP_FILE_WRITE=1, so .bin files were not generated."
    echo ""
    echo "To generate .bin files, you need to:"
    echo "  1. Recompile without SKIP_FILE_WRITE:"
    echo "     cd /home/helenshao/InitialConditions/hermitian_3d_matrix_production"
    echo "     make clean"
    echo "     make CFLAGS=\"-DUSE_DOUBLE_PRECISION -DSKIP_FILE_WRITE=0\""
    echo "  2. Run hermitian_3d_matrix again to generate .bin files"
    echo ""
    echo "Skipping step 8 (reassembly) - no .bin files available"
    SKIP_REASSEMBLY=1
else
    echo "Found $BIN_COUNT .bin files"
    echo "Running write_ic_trace.sh..."
    bash write_ic_trace.sh
    SKIP_REASSEMBLY=0
fi

# ====================================================================================
# STEP 9: Read particle ICs
# ====================================================================================

echo ""
echo "STEP 9: Reading particle ICs"
echo "----------------------------------------------------------------------"

# Create particle_ics directory if it doesn't exist
mkdir -p particle_ics

# hermitian_3d_matrix particle ICs
if [ "$SKIP_REASSEMBLY" -eq 0 ]; then
    HERMITIAN_IC_COUNT=$(find particle_ics -name "ic_*" -type f 2>/dev/null | wc -l)
    if [ "$HERMITIAN_IC_COUNT" -gt 0 ]; then
        echo "Reading hermitian particle ICs ($HERMITIAN_IC_COUNT files)..."
        for f in ./particle_ics/ic_*; do
            if [ -f "$f" ]; then
                python3 /home/helenshao/InitialConditions/hermitian_3d_matrix_production/tests/options_ab/read_particles.py "$f" 10 >> ./particle_ics/hermitian_N8_trace_particle_values.txt 2>&1
            fi
        done
        echo "Hermitian particle values saved to: particle_ics/hermitian_N8_trace_particle_values.txt"
    else
        echo "WARNING: No hermitian particle IC files found in particle_ics/"
    fi
else
    echo "Skipping hermitian particle IC reading (no .bin files available)"
fi

# zeldovich-PLT particle ICs
ZELDOVICH_OUTPUT_DIR="/home/helenshao/InitialConditions/zeldovich-PLT/output_trace"
ZELDOVICH_IC_COUNT=$(find "$ZELDOVICH_OUTPUT_DIR" -name "ic_*" -type f 2>/dev/null | wc -l)

if [ "$ZELDOVICH_IC_COUNT" -gt 0 ]; then
    echo "Reading zeldovich particle ICs ($ZELDOVICH_IC_COUNT files)..."
    for f in "$ZELDOVICH_OUTPUT_DIR"/ic_*; do
        if [ -f "$f" ]; then
            python3 /home/helenshao/InitialConditions/hermitian_3d_matrix_production/tests/options_ab/read_particles.py "$f" 10 >> ./particle_ics/zeldovich_N8_trace_particle_values.txt 2>&1
        fi
    done
    echo "Zeldovich particle values saved to: particle_ics/zeldovich_N8_trace_particle_values.txt"
else
    echo "WARNING: No zeldovich particle IC files found in $ZELDOVICH_OUTPUT_DIR"
    echo "You may need to run zeldovich-PLT first to generate particle ICs"
    echo "Run: ./run_zeldovich_trace.sh"
fi

# ====================================================================================
# STEP 10: Compare particle ICs
# ====================================================================================

echo ""
echo "STEP 10: Comparing particle ICs"
echo "----------------------------------------------------------------------"

if [ "$SKIP_REASSEMBLY" -eq 0 ] && [ "$ZELDOVICH_IC_COUNT" -gt 0 ]; then
    echo "Running compare_ic_files.py..."
    python3 /home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_assembly/compare_ic_files.py \
        ./particle_ics \
        "$ZELDOVICH_OUTPUT_DIR" \
        > ic_comparison.txt 2>&1
    
    echo "Comparison report saved to: ic_comparison.txt"
    echo ""
    echo "Summary:"
    head -20 ic_comparison.txt
else
    echo "Skipping comparison - missing particle IC files"
    if [ "$SKIP_REASSEMBLY" -eq 1 ]; then
        echo "  - Hermitian ICs: Not available (no .bin files)"
    fi
    if [ "$ZELDOVICH_IC_COUNT" -eq 0 ]; then
        echo "  - Zeldovich ICs: Not available (run zeldovich-PLT first)"
    fi
fi

# ====================================================================================
# STEP 11: Visualize particle fields
# ====================================================================================

echo ""
echo "STEP 11: Visualizing particle fields"
echo "----------------------------------------------------------------------"

if [ "$SKIP_REASSEMBLY" -eq 0 ] && [ "$ZELDOVICH_IC_COUNT" -gt 0 ]; then
    mkdir -p visualizations
    echo "Running visualize_particle_comparison.py..."
    python3 /home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_assembly/visualize_particle_comparison.py \
        ./particle_ics \
        "$ZELDOVICH_OUTPUT_DIR" \
        ./visualizations
    
    echo "Visualizations saved to: ./visualizations/"
else
    echo "Skipping visualization - missing particle IC files"
fi

echo ""
echo "=========================================="
echo "Particle IC comparison completed"
echo "=========================================="
echo ""
echo "Files created:"
echo "  - particle_ics/hermitian_N8_trace_particle_values.txt (if .bin files were available)"
echo "  - particle_ics/zeldovich_N8_trace_particle_values.txt (if zeldovich ICs exist)"
echo "  - ic_comparison.txt (if both IC sets are available)"
echo "  - visualizations/ (if both IC sets are available)"

