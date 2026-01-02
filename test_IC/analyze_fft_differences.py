#!/usr/bin/env python3
"""
Analyze differences in matrix values after FFT between hermitian and zeldovich codes.

This script provides detailed analysis of where and why values differ:
- Patterns by array index
- Patterns by coordinate
- Scaling factors
- Statistical analysis
- Specific problematic locations

Usage:
    python3 analyze_fft_differences.py <hermitian_after_fft> <zeldovich_after_fft> <hermitian_before_fft> <zeldovich_before_fft> <output_file>
"""

import sys
import re
import numpy as np
from collections import defaultdict

def parse_matrix_file(filename):
    """Parse matrix dump file and return dictionary keyed by coordinates."""
    matrix = {}
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            
            # Parse format: Y=<y> X=<x> Z=<z> Array=<array> Re=<real> Im=<imag>
            # or: Z=<z> Y=<y> X=<x> Array=<array> Re=<real> Im=<imag>
            match = re.search(r'Y[=\s]+(\d+).*?X[=\s]+(\d+).*?Z[=\s]+(\d+).*?Array[=\s]+(\d+).*?Re[=\s]+([+-]?\d+\.?\d*[eE]?[+-]?\d*).*?Im[=\s]+([+-]?\d+\.?\d*[eE]?[+-]?\d*)', line)
            if not match:
                # Try Z-first format
                match = re.search(r'Z[=\s]+(\d+).*?Y[=\s]+(\d+).*?X[=\s]+(\d+).*?Array[=\s]+(\d+).*?Re[=\s]+([+-]?\d+\.?\d*[eE]?[+-]?\d*).*?Im[=\s]+([+-]?\d+\.?\d*[eE]?[+-]?\d*)', line)
                if match:
                    z, y, x, array_idx, re_val, im_val = match.groups()
                    key = (int(y), int(x), int(z), int(array_idx))
                    matrix[key] = (float(re_val), float(im_val))
            else:
                y, x, z, array_idx, re_val, im_val = match.groups()
                key = (int(y), int(x), int(z), int(array_idx))
                matrix[key] = (float(re_val), float(im_val))
    
    return matrix

def analyze_differences(hermitian_after, zeldovich_after, hermitian_before, zeldovich_before, output_file):
    """Analyze differences in detail."""
    
    common_keys = set(hermitian_after.keys()) & set(zeldovich_after.keys())
    
    # Collect differences by array
    diff_by_array = defaultdict(list)
    diff_by_coord = {}
    ratios = []
    
    # Statistical analysis
    all_diff_re = []
    all_diff_im = []
    all_ratio_re = []
    all_ratio_im = []
    
    for key in sorted(common_keys):
        y, x, z, array_idx = key
        
        h_re_after, h_im_after = hermitian_after[key]
        z_re_after, z_im_after = zeldovich_after[key]
        
        diff_re = abs(h_re_after - z_re_after)
        diff_im = abs(h_im_after - z_im_after)
        
        diff_by_array[array_idx].append((diff_re, diff_im, key))
        diff_by_coord[key] = {
            'hermitian': (h_re_after, h_im_after),
            'zeldovich': (z_re_after, z_im_after),
            'diff': (diff_re, diff_im)
        }
        
        all_diff_re.append(diff_re)
        all_diff_im.append(diff_im)
        
        # Calculate ratios (avoid division by zero)
        if abs(z_re_after) > 1e-15:
            ratio_re = h_re_after / z_re_after
            all_ratio_re.append(ratio_re)
        if abs(z_im_after) > 1e-15:
            ratio_im = h_im_after / z_im_after
            all_ratio_im.append(ratio_im)
        
        # Compare with before FFT if available
        if hermitian_before and zeldovich_before:
            if key in hermitian_before and key in zeldovich_before:
                h_re_before, h_im_before = hermitian_before[key]
                z_re_before, z_im_before = zeldovich_before[key]
                
                # Check if differences increased after FFT
                diff_re_before = abs(h_re_before - z_re_before)
                diff_im_before = abs(h_im_before - z_im_before)
                
                diff_by_coord[key]['before'] = {
                    'hermitian': (h_re_before, h_im_before),
                    'zeldovich': (z_re_before, z_im_before),
                    'diff': (diff_re_before, diff_im_before)
                }
    
    # Write analysis report
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("Detailed FFT Difference Analysis\n")
        f.write("=" * 80 + "\n\n")
        
        # Overall statistics
        f.write("OVERALL STATISTICS:\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total elements compared: {len(common_keys)}\n")
        f.write(f"Mean difference (real): {np.mean(all_diff_re):.6e}\n")
        f.write(f"Mean difference (imag): {np.mean(all_diff_im):.6e}\n")
        f.write(f"Max difference (real): {np.max(all_diff_re):.6e}\n")
        f.write(f"Max difference (imag): {np.max(all_diff_im):.6e}\n")
        f.write(f"Std deviation (real): {np.std(all_diff_re):.6e}\n")
        f.write(f"Std deviation (imag): {np.std(all_diff_im):.6e}\n")
        
        if all_ratio_re:
            f.write(f"\nRatio statistics (hermitian/zeldovich):\n")
            f.write(f"  Mean ratio (real): {np.mean(all_ratio_re):.6e}\n")
            f.write(f"  Median ratio (real): {np.median(all_ratio_re):.6e}\n")
            f.write(f"  Std ratio (real): {np.std(all_ratio_re):.6e}\n")
        if all_ratio_im:
            f.write(f"  Mean ratio (imag): {np.mean(all_ratio_im):.6e}\n")
            f.write(f"  Median ratio (imag): {np.median(all_ratio_im):.6e}\n")
            f.write(f"  Std ratio (imag): {np.std(all_ratio_im):.6e}\n")
        f.write("\n")
        
        # Analysis by array
        f.write("DIFFERENCES BY ARRAY:\n")
        f.write("-" * 80 + "\n")
        for array_idx in sorted(diff_by_array.keys()):
            diffs = diff_by_array[array_idx]
            mean_diff_re = np.mean([d[0] for d in diffs])
            mean_diff_im = np.mean([d[1] for d in diffs])
            max_diff_re = max([d[0] for d in diffs])
            max_diff_im = max([d[1] for d in diffs])
            
            # Find location of max difference
            max_re_key = max(diffs, key=lambda x: x[0])[2]
            max_im_key = max(diffs, key=lambda x: x[1])[2]
            
            f.write(f"Array {array_idx}:\n")
            f.write(f"  Elements: {len(diffs)}\n")
            f.write(f"  Mean diff (real): {mean_diff_re:.6e}, (imag): {mean_diff_im:.6e}\n")
            f.write(f"  Max diff (real): {max_diff_re:.6e} at Y={max_re_key[0]}, X={max_re_key[1]}, Z={max_re_key[2]}\n")
            f.write(f"  Max diff (imag): {max_diff_im:.6e} at Y={max_im_key[0]}, X={max_im_key[1]}, Z={max_im_key[2]}\n\n")
        
        # Top 20 largest differences
        f.write("TOP 20 LARGEST DIFFERENCES:\n")
        f.write("-" * 80 + "\n")
        sorted_diffs = sorted(diff_by_coord.items(), 
                             key=lambda x: max(x[1]['diff'][0], x[1]['diff'][1]), 
                             reverse=True)
        
        for i, (key, data) in enumerate(sorted_diffs[:20]):
            y, x, z, array_idx = key
            h_re, h_im = data['hermitian']
            z_re, z_im = data['zeldovich']
            diff_re, diff_im = data['diff']
            
            f.write(f"{i+1}. Y={y}, X={x}, Z={z}, Array={array_idx}:\n")
            f.write(f"   Hermitian: Re={h_re:.6e}, Im={h_im:.6e}\n")
            f.write(f"   Zeldovich: Re={z_re:.6e}, Im={z_im:.6e}\n")
            f.write(f"   Difference: Re={diff_re:.6e}, Im={diff_im:.6e}\n")
            
            if 'before' in data:
                h_re_b, h_im_b = data['before']['hermitian']
                z_re_b, z_im_b = data['before']['zeldovich']
                diff_re_b, diff_im_b = data['before']['diff']
                f.write(f"   Before FFT: H(Re={h_re_b:.6e}, Im={h_im_b:.6e}), Z(Re={z_re_b:.6e}, Im={z_im_b:.6e}), Diff(Re={diff_re_b:.6e}, Im={diff_im_b:.6e})\n")
                f.write(f"   Difference increase: Re={diff_re - diff_re_b:.6e}, Im={diff_im - diff_im_b:.6e}\n")
            f.write("\n")
        
        # Check for systematic patterns
        f.write("SYSTEMATIC PATTERNS:\n")
        f.write("-" * 80 + "\n")
        
        # Check if differences are larger for certain coordinates
        coord_diffs = defaultdict(list)
        for key, data in diff_by_coord.items():
            y, x, z, array_idx = key
            coord_diffs[(y, x, z)].append(max(data['diff'][0], data['diff'][1]))
        
        if coord_diffs:
            max_coord = max(coord_diffs.items(), key=lambda x: np.mean(x[1]))
            f.write(f"Coordinates with largest average differences: Y={max_coord[0][0]}, X={max_coord[0][1]}, Z={max_coord[0][2]}\n")
            f.write(f"  Average difference: {np.mean(max_coord[1]):.6e}\n")
            f.write(f"  Arrays at this location: {len(max_coord[1])}\n\n")
        
        # Check for zero/near-zero values
        f.write("ZERO/NEAR-ZERO ANALYSIS:\n")
        f.write("-" * 80 + "\n")
        hermitian_zeros = sum(1 for k, v in hermitian_after.items() if abs(v[0]) < 1e-10 and abs(v[1]) < 1e-10)
        zeldovich_zeros = sum(1 for k, v in zeldovich_after.items() if abs(v[0]) < 1e-10 and abs(v[1]) < 1e-10)
        f.write(f"Near-zero values (|val| < 1e-10): Hermitian={hermitian_zeros}, Zeldovich={zeldovich_zeros}\n")
        
        # Check where one is zero and the other isn't
        zero_mismatches = []
        for key in common_keys:
            h_re, h_im = hermitian_after[key]
            z_re, z_im = zeldovich_after[key]
            h_zero = abs(h_re) < 1e-10 and abs(h_im) < 1e-10
            z_zero = abs(z_re) < 1e-10 and abs(z_im) < 1e-10
            if h_zero != z_zero:
                zero_mismatches.append((key, h_zero, z_zero))
        
        f.write(f"Zero mismatches (one zero, other not): {len(zero_mismatches)}\n")
        if zero_mismatches:
            f.write("First 10 zero mismatches:\n")
            for key, h_zero, z_zero in zero_mismatches[:10]:
                y, x, z, array_idx = key
                h_re, h_im = hermitian_after[key]
                z_re, z_im = zeldovich_after[key]
                f.write(f"  Y={y}, X={x}, Z={z}, Array={array_idx}: H(zero={h_zero}, Re={h_re:.6e}, Im={h_im:.6e}), Z(zero={z_zero}, Re={z_re:.6e}, Im={z_im:.6e})\n")
    
    print(f"Analysis report written to: {output_file}")

def main():
    if len(sys.argv) != 6:
        print("Usage: python3 analyze_fft_differences.py <hermitian_after_fft> <zeldovich_after_fft> <hermitian_before_fft> <zeldovich_before_fft> <output_file>")
        print("  (hermitian_before_fft and zeldovich_before_fft can be empty strings if not available)")
        sys.exit(1)
    
    hermitian_after_file = sys.argv[1]
    zeldovich_after_file = sys.argv[2]
    hermitian_before_file = sys.argv[3] if sys.argv[3] else None
    zeldovich_before_file = sys.argv[4] if sys.argv[4] else None
    output_file = sys.argv[5]
    
    print(f"Loading hermitian after FFT from: {hermitian_after_file}")
    hermitian_after = parse_matrix_file(hermitian_after_file)
    
    print(f"Loading zeldovich after FFT from: {zeldovich_after_file}")
    zeldovich_after = parse_matrix_file(zeldovich_after_file)
    
    hermitian_before = None
    zeldovich_before = None
    
    if hermitian_before_file:
        print(f"Loading hermitian before FFT from: {hermitian_before_file}")
        hermitian_before = parse_matrix_file(hermitian_before_file)
    
    if zeldovich_before_file:
        print(f"Loading zeldovich before FFT from: {zeldovich_before_file}")
        zeldovich_before = parse_matrix_file(zeldovich_before_file)
    
    print("Analyzing differences...")
    analyze_differences(hermitian_after, zeldovich_after, hermitian_before, zeldovich_before, output_file)
    
    print("Analysis completed!")

if __name__ == "__main__":
    main()

