#!/bin/bash
# ====================================================================================
# COMPARE OUTPUT FILES: Options A vs B
# ====================================================================================
# Compares particle IC files from Option A and Option B to verify correctness.
# Checks file counts, sizes, and optionally particle data.
#
# Usage:
#   ./compare_option_outputs.sh output_option_a output_option_b [--verify-data]
# ====================================================================================

set -e

OUTPUT_DIR_A="${1:-output_option_a}"
OUTPUT_DIR_B="${2:-output_option_b}"
VERIFY_DATA="${3:-}"

echo "========================================================================"
echo "COMPARING OUTPUTS: Option A vs Option B"
echo "========================================================================"
echo "Option A directory: $OUTPUT_DIR_A"
echo "Option B directory: $OUTPUT_DIR_B"
echo "========================================================================"
echo ""

# Check directories exist
if [ ! -d "$OUTPUT_DIR_A" ]; then
    echo "ERROR: $OUTPUT_DIR_A not found"
    exit 1
fi

if [ ! -d "$OUTPUT_DIR_B" ]; then
    echo "ERROR: $OUTPUT_DIR_B not found"
    exit 1
fi

# Count files (support per-rank format ic_rank*_i*_k* and full-range format ic_*)
NUM_FILES_A=$(find "$OUTPUT_DIR_A" -type f \( -name "ic_rank*" -o -name "ic_[0-9]*" \) | wc -l)
NUM_FILES_B=$(find "$OUTPUT_DIR_B" -type f \( -name "ic_rank*" -o -name "ic_[0-9]*" \) | wc -l)

echo "File counts:"
echo "  Option A: $NUM_FILES_A particle IC files"
echo "  Option B: $NUM_FILES_B particle IC files"

if [ "$NUM_FILES_A" -ne "$NUM_FILES_B" ]; then
    echo "  ⚠ WARNING: File counts differ!"
    echo ""
    echo "Missing files in Option A:"
    comm -23 <(find "$OUTPUT_DIR_B" -type f \( -name "ic_rank*" -o -name "ic_[0-9]*" \) 2>/dev/null | xargs -n1 basename 2>/dev/null | sort) \
             <(find "$OUTPUT_DIR_A" -type f \( -name "ic_rank*" -o -name "ic_[0-9]*" \) 2>/dev/null | xargs -n1 basename 2>/dev/null | sort) | head -10
    echo ""
    echo "Missing files in Option B:"
    comm -13 <(find "$OUTPUT_DIR_B" -type f \( -name "ic_rank*" -o -name "ic_[0-9]*" \) 2>/dev/null | xargs -n1 basename 2>/dev/null | sort) \
             <(find "$OUTPUT_DIR_A" -type f \( -name "ic_rank*" -o -name "ic_[0-9]*" \) 2>/dev/null | xargs -n1 basename 2>/dev/null | sort) | head -10
    exit 1
else
    echo "  ✓ File counts match"
fi
echo ""

# Compare file sizes
echo "Comparing file sizes..."
SIZE_MISMATCH=0
for file_a in "$OUTPUT_DIR_A"/ic_rank* "$OUTPUT_DIR_A"/ic_[0-9]*; do
    if [ ! -f "$file_a" ]; then
        continue  # Skip if no files match the pattern
    fi
    file_b="$OUTPUT_DIR_B/$(basename "$file_a")"
    if [ -f "$file_b" ]; then
        size_a=$(stat -f%z "$file_a" 2>/dev/null || stat -c%s "$file_a" 2>/dev/null)
        size_b=$(stat -f%z "$file_b" 2>/dev/null || stat -c%s "$file_b" 2>/dev/null)
        if [ "$size_a" -ne "$size_b" ]; then
            echo "  ⚠ Size mismatch: $(basename "$file_a") (A: $size_a, B: $size_b)"
            SIZE_MISMATCH=$((SIZE_MISMATCH + 1))
        fi
    fi
done

if [ $SIZE_MISMATCH -eq 0 ]; then
    echo "  ✓ All file sizes match"
else
    echo "  ⚠ WARNING: $SIZE_MISMATCH files have size mismatches"
fi
echo ""

# Print sample outputs from each option
echo "========================================================================"
echo "SAMPLE PARTICLE DATA COMPARISON"
echo "========================================================================"
echo ""

# Check if Python is available for reading particle files
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
READ_PARTICLES_PY="$SCRIPT_DIR/read_particles.py"

if [ ! -f "$READ_PARTICLES_PY" ] || ! command -v python3 &> /dev/null; then
    echo "  ⚠ Python3 or read_particles.py not found - skipping detailed particle display"
    echo "  (File counts and sizes were compared above)"
    echo ""
else
    # Function to show particles using Python script
    show_particles() {
        local file="$1"
        local num_particles="${2:-3}"
        python3 "$READ_PARTICLES_PY" "$file" "$num_particles" 2>/dev/null || echo "  Error reading particles from $file"
    }

    # Find first few matching files to compare
    echo "Finding matching files to compare..."
    MATCHING_FILES=$(comm -12 \
        <(find "$OUTPUT_DIR_A" -type f \( -name "ic_rank*" -o -name "ic_[0-9]*" \) 2>/dev/null | xargs -n1 basename 2>/dev/null | sort) \
        <(find "$OUTPUT_DIR_B" -type f \( -name "ic_rank*" -o -name "ic_[0-9]*" \) 2>/dev/null | xargs -n1 basename 2>/dev/null | sort) | head -3)
    
    if [ -z "$MATCHING_FILES" ]; then
        echo "  ⚠ No matching files found to compare"
        echo ""
    else
        echo "Comparing first 3 matching files:"
        echo ""
        
        FILE_NUM=0
        for filename in $MATCHING_FILES; do
            FILE_NUM=$((FILE_NUM + 1))
            file_a="$OUTPUT_DIR_A/$filename"
            file_b="$OUTPUT_DIR_B/$filename"
            
            echo "----------------------------------------------------------------------"
            echo "File $FILE_NUM: $filename"
            echo "----------------------------------------------------------------------"
            echo ""
            echo "Option A:"
            show_particles "$file_a" 5
            echo ""
            
            echo "Option B:"
            show_particles "$file_b" 5
            echo ""
            
            # Quick binary comparison
            if cmp -s "$file_a" "$file_b"; then
                echo "  ✓ Files are identical (binary match)"
            else
                echo "  ⚠ Files differ (binary mismatch)"
                # Show first difference location
                local diff_pos=$(cmp "$file_a" "$file_b" 2>&1 | grep -oP 'differ: byte \K[0-9]+' || echo "unknown")
                if [ "$diff_pos" != "unknown" ]; then
                    local particle_num=$((diff_pos / 36))
                    echo "    First difference at byte $diff_pos (particle ~$particle_num)"
                fi
            fi
            echo ""
        done
    fi
fi

# Verify data if requested
if [ "$VERIFY_DATA" = "--verify-data" ]; then
    echo "========================================================================"
    echo "FULL BINARY VERIFICATION"
    echo "========================================================================"
    echo "Verifying particle data (comparing all files)..."
    VERIFY_COUNT=0
    FILE_COUNT=0
    for file_a in "$OUTPUT_DIR_A"/ic_rank* "$OUTPUT_DIR_A"/ic_[0-9]*; do
        if [ ! -f "$file_a" ]; then
            continue
        fi
        file_b="$OUTPUT_DIR_B/$(basename "$file_a")"
        if [ -f "$file_b" ]; then
            if cmp -s "$file_a" "$file_b"; then
                VERIFY_COUNT=$((VERIFY_COUNT + 1))
            else
                echo "  ⚠ Data mismatch: $(basename "$file_a")"
            fi
            FILE_COUNT=$((FILE_COUNT + 1))
        fi
    done
    echo "  ✓ $VERIFY_COUNT/$FILE_COUNT files have identical data"
    echo ""
fi

# Summary
echo "========================================================================"
echo "COMPARISON COMPLETE"
echo "========================================================================"
if [ $SIZE_MISMATCH -eq 0 ]; then
    echo "✓ All files match in size"
    echo "✓ Options A and B produce consistent outputs"
else
    echo "⚠ Some files have size differences - investigate further"
fi
echo ""

