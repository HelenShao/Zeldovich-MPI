#!/bin/bash
# ====================================================================================
# COMPARE TIMING: Options A vs B
# ====================================================================================
# Extracts and compares timing information from PBS job output files.
# Shows elapsed time, memory usage, and calculates speedup.
#
# Usage:
#   ./compare_timing.sh [N]
#
# Examples:
#   ./compare_timing.sh 256    # Compare N=256 tests
#   ./compare_timing.sh 2048   # Compare N=2048 tests
#   ./compare_timing.sh        # Compare all available tests
# ====================================================================================

set -e

N="${1:-}"

echo "========================================================================"
echo "TIMING COMPARISON: Options A vs B"
echo "========================================================================"
echo ""

# Function to extract and display timing from a PBS output file
extract_timing() {
    local out_file="$1"
    local option_name="$2"
    
    if [ ! -f "$out_file" ]; then
        echo "  ⚠ $out_file not found" >&2
        echo ""
        return 1
    fi
    
    # Extract elapsed time (try multiple patterns)
    # Look for "Total time (including all stages): X.XXXXXX s" or "Elapsed time: X seconds"
    local elapsed=$(grep -E "Total time.*including|Elapsed time|Total 3D FFT time" "$out_file" 2>/dev/null | tail -1 | awk '{for(i=1;i<=NF;i++) if($i~/^[0-9]+\.[0-9]+$/) {print $i; break} else if($i~/^[0-9]+$/ && $(i+1)~/^s/) {print $i; break}}' | head -1)
    
    # Extract file count from "Total files written: 512 (across all ranks)" or "Output files: X"
    local files=$(grep -E "Total files written|Output files:" "$out_file" 2>/dev/null | tail -1 | awk '{for(i=1;i<=NF;i++) if($i=="written:") {print $(i+1); break} else if($i=="files:" && $(i+1)~/^[0-9]+$/) {print $(i+1); break}}' | head -1)
    
    # Extract compilation time (if available)
    local compile_time=$(grep -i "compilation" "$out_file" | grep -i "time" | tail -1 || echo "")
    
    echo "  $option_name:" >&2
    if [ -n "$elapsed" ]; then
        echo "    Elapsed time: ${elapsed} seconds" >&2
    else
        echo "    Elapsed time: Not found" >&2
    fi
    
    if [ -n "$files" ]; then
        echo "    Output files: $files" >&2
    fi
    
    # Try to extract memory from PBS output (if available)
    local max_rss=$(grep -i "max.*rss\|maximum.*resident" "$out_file" | tail -1 | awk '{print $(NF-1), $NF}' || echo "N/A")
    if [ "$max_rss" != "N/A" ] && [ -n "$max_rss" ]; then
        echo "    Memory: $max_rss" >&2
    fi
    
    echo "" >&2
    
    # Return elapsed time for calculation (to stdout)
    echo "$elapsed"
}

# Compare N=256 if requested or if no N specified
if [ -z "$N" ] || [ "$N" = "256" ]; then
    echo "========================================================================"
    echo "N=256 Comparison"
    echo "========================================================================"
    echo ""
    
    OUT_A_256="test_option_a_N256.out"
    OUT_B_256="test_option_b_N256.out"
    
    if [ ! -f "$OUT_A_256" ] && [ ! -f "$OUT_B_256" ]; then
        echo "  ⚠ No N=256 output files found"
        echo ""
    else
        # Extract timing (function prints to stderr, returns elapsed time to stdout)
        ELAPSED_A_256=$(extract_timing "$OUT_A_256" "Option A")
        ELAPSED_B_256=$(extract_timing "$OUT_B_256" "Option B")
        
        # Calculate speedup if both times are available
        if [ -n "$ELAPSED_A_256" ] && [ -n "$ELAPSED_B_256" ] && [ "$ELAPSED_A_256" != "Not found" ] && [ "$ELAPSED_B_256" != "Not found" ]; then
            # Check if both are valid numbers (including decimals)
            if echo "$ELAPSED_A_256" | grep -qE '^[0-9]+\.?[0-9]*$' && echo "$ELAPSED_B_256" | grep -qE '^[0-9]+\.?[0-9]*$'; then
                if [ "$(echo "$ELAPSED_A_256 > 0" | bc 2>/dev/null || echo 0)" = "1" ] && [ "$(echo "$ELAPSED_B_256 > 0" | bc 2>/dev/null || echo 0)" = "1" ]; then
                SPEEDUP_256=$(echo "scale=2; $ELAPSED_A_256 / $ELAPSED_B_256" | bc 2>/dev/null || echo "N/A")
                TIME_SAVED_256=$(echo "scale=2; $ELAPSED_A_256 - $ELAPSED_B_256" | bc 2>/dev/null || echo "N/A")
                
                    if [ "$SPEEDUP_256" != "N/A" ]; then
                        echo "  Summary:"
                        echo "    Speedup: ${SPEEDUP_256}× (Option B is ${SPEEDUP_256}× faster)"
                        if [ "$TIME_SAVED_256" != "N/A" ]; then
                            echo "    Time saved: ${TIME_SAVED_256} seconds"
                        fi
                        echo ""
                    fi
                fi
            fi
        fi
    fi
fi

# Compare N=2048 if requested or if no N specified
if [ -z "$N" ] || [ "$N" = "2048" ]; then
    echo "========================================================================"
    echo "N=2048 Comparison"
    echo "========================================================================"
    echo ""
    
    OUT_A_2048="test_option_a_N2048.out"
    OUT_B_2048="test_option_b_N2048.out"
    
    if [ ! -f "$OUT_A_2048" ] && [ ! -f "$OUT_B_2048" ]; then
        echo "  ⚠ No N=2048 output files found"
        echo ""
    else
        # Extract timing (function prints to stderr, returns elapsed time to stdout)
        ELAPSED_A_2048=$(extract_timing "$OUT_A_2048" "Option A")
        ELAPSED_B_2048=$(extract_timing "$OUT_B_2048" "Option B")
        
        # Calculate speedup if both times are available
        if [ -n "$ELAPSED_A_2048" ] && [ -n "$ELAPSED_B_2048" ] && [ "$ELAPSED_A_2048" != "Not found" ] && [ "$ELAPSED_B_2048" != "Not found" ]; then
            # Check if both are valid numbers (including decimals)
            if echo "$ELAPSED_A_2048" | grep -qE '^[0-9]+\.?[0-9]*$' && echo "$ELAPSED_B_2048" | grep -qE '^[0-9]+\.?[0-9]*$'; then
                if [ "$(echo "$ELAPSED_A_2048 > 0" | bc 2>/dev/null || echo 0)" = "1" ] && [ "$(echo "$ELAPSED_B_2048 > 0" | bc 2>/dev/null || echo 0)" = "1" ]; then
                SPEEDUP_2048=$(echo "scale=2; $ELAPSED_A_2048 / $ELAPSED_B_2048" | bc 2>/dev/null || echo "N/A")
                TIME_SAVED_2048=$(echo "scale=2; $ELAPSED_A_2048 - $ELAPSED_B_2048" | bc 2>/dev/null || echo "N/A")
                
                    if [ "$SPEEDUP_2048" != "N/A" ]; then
                        echo "  Summary:"
                        echo "    Speedup: ${SPEEDUP_2048}× (Option B is ${SPEEDUP_2048}× faster)"
                        if [ "$TIME_SAVED_2048" != "N/A" ]; then
                            echo "    Time saved: ${TIME_SAVED_2048} seconds"
                        fi
                        echo ""
                    fi
                fi
            fi
        fi
    fi
fi

echo "========================================================================"
echo "TIMING COMPARISON COMPLETE"
echo "========================================================================"
echo ""
echo "Note: For detailed memory and resource usage, check PBS output files:"
echo "  cat test_option_*_N*.out"
echo ""

