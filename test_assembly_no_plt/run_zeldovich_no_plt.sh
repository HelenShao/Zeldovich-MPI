#!/bin/bash
# Script to run zeldovich-PLT with PLT disabled for comparison
# This should be run separately (not in PBS) to generate reference output

set -e

ZELDOVICH_DIR="/home/helenshao/InitialConditions/zeldovich-PLT"
PARAM_FILE="/home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_assembly_no_plt/param_N16_no_plt.par"
OUTPUT_DIR="${ZELDOVICH_DIR}/output_no_plt"

echo "=========================================="
echo "Running zeldovich-PLT with PLT DISABLED"
echo "=========================================="
echo "Parameter file: $PARAM_FILE"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Check if zeldovich executable exists
ZELDOVICH_EXEC="${ZELDOVICH_DIR}/build/zeldovich"
if [ ! -f "$ZELDOVICH_EXEC" ]; then
    echo "ERROR: zeldovich executable not found: $ZELDOVICH_EXEC"
    echo "Please build zeldovich-PLT first:"
    echo "  cd $ZELDOVICH_DIR && make"
    exit 1
fi

# Check if parameter file exists
if [ ! -f "$PARAM_FILE" ]; then
    echo "ERROR: Parameter file not found: $PARAM_FILE"
    exit 1
fi

# Verify parameter file has PLT disabled
if grep -q "ZD_qPLT = 1" "$PARAM_FILE"; then
    echo "WARNING: Parameter file has ZD_qPLT = 1 (PLT enabled)"
    echo "This script expects PLT to be disabled. Please check the parameter file."
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"
echo "Output will be written to: $OUTPUT_DIR"
echo ""

# Update parameter file to use the new output directory
# Create a temporary parameter file with updated output directory
TEMP_PARAM=$(mktemp)
cp "$PARAM_FILE" "$TEMP_PARAM"
sed -i "s|InitialConditionsDirectory = \".*\"|InitialConditionsDirectory = \"$OUTPUT_DIR\"|" "$TEMP_PARAM"

echo "Running zeldovich-PLT..."
cd "$ZELDOVICH_DIR"
"$ZELDOVICH_EXEC" "$TEMP_PARAM" 2>&1 | tee "$OUTPUT_DIR/zeldovich_no_plt.log"

# Clean up temp file
rm -f "$TEMP_PARAM"

# Check for output files
IC_COUNT=$(find "$OUTPUT_DIR" -name "ic_*" 2>/dev/null | wc -l)
echo ""
echo "=========================================="
echo "zeldovich-PLT (PLT disabled) completed"
echo "=========================================="
echo "Output directory: $OUTPUT_DIR"
echo "Particle IC files: $IC_COUNT"
echo ""
echo "You can now compare this output with:"
echo "  - hermitian_3d_matrix_production/test_assembly_no_plt/particle_ics/ (PLT disabled)"
echo "  - hermitian_3d_matrix_production/test_assembly/particle_ics/ (PLT enabled)"
echo "  - zeldovich-PLT/output/ (PLT enabled, if exists)"

