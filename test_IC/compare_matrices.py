#!/usr/bin/env python3
"""
Compare matrix dumps between hermitian_3d_matrix_production and zeldovich-PLT.

This script compares matrix values at the same stage (before or after FFT)
and reports differences.

Usage:
    python3 compare_matrices.py <hermitian_matrix_file> <zeldovich_matrix_file> <output_file>
"""

import sys
import re
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

def compare_matrices(hermitian_matrix, zeldovich_matrix, output_file):
    """Compare two matrices and write comparison report."""
    hermitian_keys = set(hermitian_matrix.keys())
    zeldovich_keys = set(zeldovich_matrix.keys())
    
    common_keys = hermitian_keys & zeldovich_keys
    only_hermitian = hermitian_keys - zeldovich_keys
    only_zeldovich = zeldovich_keys - hermitian_keys
    
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("Matrix Comparison Report\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"Total elements in hermitian: {len(hermitian_keys)}\n")
        f.write(f"Total elements in zeldovich: {len(zeldovich_keys)}\n")
        f.write(f"Common elements: {len(common_keys)}\n")
        f.write(f"Only in hermitian: {len(only_hermitian)}\n")
        f.write(f"Only in zeldovich: {len(only_zeldovich)}\n\n")
        
        # Compare common elements
        exact_matches = 0
        close_matches = 0
        different = 0
        max_diff_re = 0.0
        max_diff_im = 0.0
        max_diff_key = None
        
        differences = []
        
        for key in sorted(common_keys):
            h_re, h_im = hermitian_matrix[key]
            z_re, z_im = zeldovich_matrix[key]
            
            diff_re = abs(h_re - z_re)
            diff_im = abs(h_im - z_im)
            
            if diff_re < 1e-15 and diff_im < 1e-15:
                exact_matches += 1
            elif diff_re < 1e-5 and diff_im < 1e-5:
                close_matches += 1
            else:
                different += 1
                differences.append((key, h_re, h_im, z_re, z_im, diff_re, diff_im))
            
            if diff_re > max_diff_re or diff_im > max_diff_im:
                max_diff_re = max(max_diff_re, diff_re)
                max_diff_im = max(max_diff_im, diff_im)
                max_diff_key = key
        
        f.write("Value Comparison (common elements):\n")
        f.write(f"  Exact matches (diff < 1e-15): {exact_matches}\n")
        f.write(f"  Close matches (diff < 1e-5): {close_matches}\n")
        f.write(f"  Different values: {different}\n\n")
        
        if max_diff_key:
            f.write(f"Maximum differences:\n")
            f.write(f"  Real part: {max_diff_re:.6e}\n")
            f.write(f"  Imaginary part: {max_diff_im:.6e}\n")
            f.write(f"  Location: Y={max_diff_key[0]}, X={max_diff_key[1]}, Z={max_diff_key[2]}, Array={max_diff_key[3]}\n\n")
        
        # Show first 20 differences
        if differences:
            f.write("First 20 elements with differences:\n")
            for i, (key, h_re, h_im, z_re, z_im, diff_re, diff_im) in enumerate(differences[:20]):
                f.write(f"  Y={key[0]}, X={key[1]}, Z={key[2]}, Array={key[3]}:\n")
                f.write(f"    hermitian: Re={h_re:.6e}, Im={h_im:.6e}\n")
                f.write(f"    zeldovich: Re={z_re:.6e}, Im={z_im:.6e}\n")
                f.write(f"    Diff: Re={diff_re:.6e}, Im={diff_im:.6e}\n\n")
        
        # Show elements only in one code
        if only_hermitian:
            f.write(f"\nElements only in hermitian (first 10):\n")
            for key in sorted(only_hermitian)[:10]:
                re_val, im_val = hermitian_matrix[key]
                f.write(f"  Y={key[0]}, X={key[1]}, Z={key[2]}, Array={key[3]}: Re={re_val:.6e}, Im={im_val:.6e}\n")
        
        if only_zeldovich:
            f.write(f"\nElements only in zeldovich (first 10):\n")
            for key in sorted(only_zeldovich)[:10]:
                re_val, im_val = zeldovich_matrix[key]
                f.write(f"  Y={key[0]}, X={key[1]}, Z={key[2]}, Array={key[3]}: Re={re_val:.6e}, Im={im_val:.6e}\n")
    
    print(f"Comparison report written to: {output_file}")

def main():
    if len(sys.argv) != 4:
        print("Usage: python3 compare_matrices.py <hermitian_matrix_file> <zeldovich_matrix_file> <output_file>")
        sys.exit(1)
    
    hermitian_file = sys.argv[1]
    zeldovich_file = sys.argv[2]
    output_file = sys.argv[3]
    
    print(f"Loading hermitian matrix from: {hermitian_file}")
    hermitian_matrix = parse_matrix_file(hermitian_file)
    
    print(f"Loading zeldovich matrix from: {zeldovich_file}")
    zeldovich_matrix = parse_matrix_file(zeldovich_file)
    
    print("Comparing matrices...")
    compare_matrices(hermitian_matrix, zeldovich_matrix, output_file)
    
    print("Comparison completed!")

if __name__ == "__main__":
    main()

