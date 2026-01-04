#!/bin/bash
# Script to analyze core dump with gdb
# Usage: ./analyze_core.sh

cd bin_files
CORE_FILE=$(ls -t core.* 2>/dev/null | head -1)

if [ -z "$CORE_FILE" ]; then
    echo "ERROR: No core dump found"
    exit 1
fi

EXECUTABLE="/home/helenshao/InitialConditions/hermitian_3d_matrix_production/hermitian_3d_matrix"

echo "Found core dump: $CORE_FILE"
echo "Executable: $EXECUTABLE"
echo ""

if [ ! -f "$EXECUTABLE" ]; then
    echo "ERROR: Executable not found at $EXECUTABLE"
    echo "Please compile first, then run this script"
    exit 1
fi

echo "To analyze with gdb, run:"
echo "  module load gdb/15.2"
echo "  cd /home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_N4/rng_match/bin_files"
echo "  gdb $EXECUTABLE $CORE_FILE"
echo ""
echo "Then in gdb, type:"
echo "  bt          # Full backtrace with line numbers"
echo "  bt full     # Full backtrace with local variables"
echo "  frame 0     # Jump to crashing frame"
echo "  list        # Show source code at crash"
echo "  info registers  # CPU state at crash"
