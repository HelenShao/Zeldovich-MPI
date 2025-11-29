#!/bin/bash
# validate_hermitian.sh - Validate that output is purely real (Hermitian symmetry)

LOG_FILE=$1

if [ -z "$LOG_FILE" ]; then
    echo "Usage: validate_hermitian.sh <log_file>"
    exit 1
fi

echo "=========================================="
echo "Hermitian Symmetry Validation"
echo "=========================================="
echo "Checking log: $LOG_FILE"
echo ""

# Check for verify_real_space_symmetry output
if grep -q "Verifying final real-space imaginary parts" "$LOG_FILE"; then
    echo "✓ Found real-space symmetry verification"
    
    # Extract max_imag
    MAX_IMAG=$(grep "max_imag" "$LOG_FILE" | sed 's/.*max_imag = \([0-9.eE+-]*\).*/\1/')
    if [ -n "$MAX_IMAG" ]; then
        echo "  max_imag = $MAX_IMAG"
        # Check if max_imag < 1e-10 (using awk for floating point comparison)
        if awk "BEGIN {exit !($MAX_IMAG < 1e-10)}"; then
            echo "  ✓ max_imag < 1e-10 (PASS)"
        else
            echo "  ✗ max_imag >= 1e-10 (FAIL)"
        fi
    fi
    
    # Extract rms_imag
    RMS_IMAG=$(grep "rms_imag" "$LOG_FILE" | sed 's/.*rms_imag = \([0-9.eE+-]*\).*/\1/')
    if [ -n "$RMS_IMAG" ]; then
        echo "  rms_imag = $RMS_IMAG"
        if awk "BEGIN {exit !($RMS_IMAG < 1e-12)}"; then
            echo "  ✓ rms_imag < 1e-12 (PASS)"
        else
            echo "  ✗ rms_imag >= 1e-12 (FAIL)"
        fi
    fi
    
    # Extract error count
    ERRORS=$(grep "errors = [0-9]" "$LOG_FILE" | sed 's/.*errors = \([0-9]*\).*/\1/')
    if [ -n "$ERRORS" ]; then
        echo "  errors = $ERRORS"
        if [ "$ERRORS" -eq 0 ]; then
            echo "  ✓ No errors (PASS)"
        else
            echo "  ✗ $ERRORS errors found (FAIL)"
        fi
    fi
    
    # Check final status
    if grep -q "OK: Final matrix is (numerically) real" "$LOG_FILE"; then
        echo "  ✓ Final status: OK (PASS)"
    elif grep -q "NOT OK: Matrix has non-zero imaginary parts" "$LOG_FILE"; then
        echo "  ✗ Final status: NOT OK (FAIL)"
    fi
else
    echo "✗ Real-space symmetry verification not found in log"
    echo "  (May need to enable RECONSTRUCT_GLOBAL_FOR_VERIFICATION=1)"
fi

echo ""
echo "=========================================="

