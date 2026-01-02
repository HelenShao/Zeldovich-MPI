#!/usr/bin/env python3
"""
Generate summary report from trace comparisons.

Usage:
    python3 generate_trace_report.py <comparisons_dir> <output_file>
"""

import sys
import os
from pathlib import Path

def read_comparison_file(filename):
    """Read comparison file and extract key statistics."""
    stats = {}
    
    if not os.path.exists(filename):
        return None
    
    with open(filename, 'r') as f:
        content = f.read()
        
        # Extract statistics
        for line in content.split('\n'):
            if 'Total elements in hermitian:' in line:
                stats['hermitian_total'] = int(line.split(':')[1].strip())
            elif 'Total elements in zeldovich:' in line:
                stats['zeldovich_total'] = int(line.split(':')[1].strip())
            elif 'Common elements:' in line:
                stats['common'] = int(line.split(':')[1].strip())
            elif 'Exact matches' in line:
                stats['exact_matches'] = int(line.split(':')[1].strip())
            elif 'Close matches' in line:
                stats['close_matches'] = int(line.split(':')[1].strip())
            elif 'Different values:' in line:
                stats['different'] = int(line.split(':')[1].strip())
            elif 'Maximum differences:' in line:
                # Next lines contain max diff info
                pass
            elif 'Real part:' in line and 'max_diff_re' not in stats:
                stats['max_diff_re'] = float(line.split(':')[1].strip())
            elif 'Imaginary part:' in line and 'max_diff_im' not in stats:
                stats['max_diff_im'] = float(line.split(':')[1].strip())
    
    return stats

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 generate_trace_report.py <comparisons_dir> <output_file>")
        sys.exit(1)
    
    comparisons_dir = sys.argv[1]
    output_file = sys.argv[2]
    
    before_fft_file = os.path.join(comparisons_dir, "before_fft_comparison.txt")
    after_fft_file = os.path.join(comparisons_dir, "after_fft_comparison.txt")
    
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("Trace Comparison Summary Report\n")
        f.write("=" * 80 + "\n\n")
        
        # Before FFT comparison
        f.write("BEFORE FFT (Fourier Space, after generation):\n")
        f.write("-" * 80 + "\n")
        before_stats = read_comparison_file(before_fft_file)
        if before_stats:
            f.write(f"  Total elements (hermitian): {before_stats.get('hermitian_total', 'N/A')}\n")
            f.write(f"  Total elements (zeldovich): {before_stats.get('zeldovich_total', 'N/A')}\n")
            f.write(f"  Common elements: {before_stats.get('common', 'N/A')}\n")
            f.write(f"  Exact matches: {before_stats.get('exact_matches', 'N/A')}\n")
            f.write(f"  Close matches: {before_stats.get('close_matches', 'N/A')}\n")
            f.write(f"  Different values: {before_stats.get('different', 'N/A')}\n")
            if 'max_diff_re' in before_stats:
                f.write(f"  Max difference (real): {before_stats['max_diff_re']:.6e}\n")
            if 'max_diff_im' in before_stats:
                f.write(f"  Max difference (imag): {before_stats['max_diff_im']:.6e}\n")
        else:
            f.write("  Comparison file not found or could not be parsed.\n")
        f.write("\n")
        
        # After FFT comparison
        f.write("AFTER FFT (Real Space, after 3D FFT):\n")
        f.write("-" * 80 + "\n")
        after_stats = read_comparison_file(after_fft_file)
        if after_stats:
            f.write(f"  Total elements (hermitian): {after_stats.get('hermitian_total', 'N/A')}\n")
            f.write(f"  Total elements (zeldovich): {after_stats.get('zeldovich_total', 'N/A')}\n")
            f.write(f"  Common elements: {after_stats.get('common', 'N/A')}\n")
            f.write(f"  Exact matches: {after_stats.get('exact_matches', 'N/A')}\n")
            f.write(f"  Close matches: {after_stats.get('close_matches', 'N/A')}\n")
            f.write(f"  Different values: {after_stats.get('different', 'N/A')}\n")
            if 'max_diff_re' in after_stats:
                f.write(f"  Max difference (real): {after_stats['max_diff_re']:.6e}\n")
            if 'max_diff_im' in after_stats:
                f.write(f"  Max difference (imag): {after_stats['max_diff_im']:.6e}\n")
        else:
            f.write("  Comparison file not found or could not be parsed.\n")
        f.write("\n")
        
        # Summary
        f.write("SUMMARY:\n")
        f.write("-" * 80 + "\n")
        if before_stats and after_stats:
            if before_stats.get('exact_matches', 0) == before_stats.get('common', 0):
                f.write("  ✓ Matrices match exactly before FFT\n")
            else:
                f.write("  ✗ Matrices differ before FFT\n")
            
            if after_stats.get('exact_matches', 0) == after_stats.get('common', 0):
                f.write("  ✓ Matrices match exactly after FFT\n")
            else:
                f.write("  ✗ Matrices differ after FFT\n")
        else:
            f.write("  Could not generate summary (comparison files missing)\n")
    
    print(f"Summary report written to: {output_file}")

if __name__ == "__main__":
    main()

