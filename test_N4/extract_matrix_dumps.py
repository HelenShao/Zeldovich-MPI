#!/usr/bin/env python3
"""
Extract matrix dumps from hermitian_3d_matrix_production or zeldovich-PLT.

This script looks for matrix dump files directly (matrix_before_fft.txt, matrix_after_fft.txt)
or extracts from log files if dump files don't exist.

Usage:
    python3 extract_matrix_dumps.py <log_file_or_dir> <output_dir>
    
    If <log_file_or_dir> is a directory, it will look for matrix_*.txt files there.
    If it's a file, it will try to extract from the log or look for dump files in the same directory.
"""

import sys
import re
import os
from pathlib import Path

def copy_or_extract_matrix_dump(input_path, output_file, dump_type="before"):
    """Copy existing dump file or extract from log."""
    input_dir = os.path.dirname(input_path) if os.path.isfile(input_path) else input_path
    
    # First, try to find existing dump file
    if dump_type == "before":
        dump_filename = "matrix_before_fft.txt"
        pattern = r'Y[=\s]+(\d+).*?X[=\s]+(\d+).*?Z[=\s]+(\d+).*?Array[=\s]+(\d+).*?Re[=\s]+([+-]?\d+\.?\d*[eE]?[+-]?\d*).*?Im[=\s]+([+-]?\d+\.?\d*[eE]?[+-]?\d*)'
    else:
        dump_filename = "matrix_after_fft.txt"
        pattern = r'Z[=\s]+(\d+).*?Y[=\s]+(\d+).*?X[=\s]+(\d+).*?Array[=\s]+(\d+).*?Re[=\s]+([+-]?\d+\.?\d*[eE]?[+-]?\d*).*?Im[=\s]+([+-]?\d+\.?\d*[eE]?[+-]?\d*)'
    
    # Look for dump file in input directory
    dump_file_path = os.path.join(input_dir, dump_filename)
    if os.path.exists(dump_file_path):
        # Copy the dump file
        with open(dump_file_path, 'r') as f_in:
            content = f_in.read()
        with open(output_file, 'w') as f_out:
            f_out.write(content)
        
        # Count non-header lines
        lines = [l for l in content.split('\n') if l and not l.startswith('#')]
        print(f"Copied {len(lines)} matrix elements {dump_type} FFT from {dump_file_path} to {output_file}")
        return len(lines)
    
    # If dump file doesn't exist, try to extract from log file
    if os.path.isfile(input_path) and input_path.endswith('.log'):
        matrix_data = []
        with open(input_path, 'r') as f:
            for line in f:
                match = re.search(pattern, line)
                if match:
                    if dump_type == "before":
                        y, x, z, array_idx, re_val, im_val = match.groups()
                        matrix_data.append(f"Y={y} X={x} Z={z} Array={array_idx} Re={re_val} Im={im_val}\n")
                    else:
                        z, y, x, array_idx, re_val, im_val = match.groups()
                        matrix_data.append(f"Z={z} Y={y} X={x} Array={array_idx} Re={re_val} Im={im_val}\n")
        
        if matrix_data:
            with open(output_file, 'w') as f:
                if dump_type == "before":
                    f.write("# Matrix values before FFT (Fourier space)\n")
                    f.write("# Format: Y=<y> X=<x> Z=<z> Array=<array_idx> Re=<real> Im=<imag>\n")
                else:
                    f.write("# Matrix values after FFT (Real space)\n")
                    f.write("# Format: Z=<z> Y=<y> X=<x> Array=<array_idx> Re=<real> Im=<imag>\n")
                for line in sorted(matrix_data):
                    f.write(line)
            print(f"Extracted {len(matrix_data)} matrix elements {dump_type} FFT from log to {output_file}")
            return len(matrix_data)
    
    print(f"Warning: No matrix dump found {dump_type} FFT. Created empty file: {output_file}")
    # Create empty file with header
    with open(output_file, 'w') as f:
        if dump_type == "before":
            f.write("# Matrix values before FFT (Fourier space)\n")
            f.write("# Format: Y=<y> X=<x> Z=<z> Array=<array_idx> Re=<real> Im=<imag>\n")
        else:
            f.write("# Matrix values after FFT (Real space)\n")
            f.write("# Format: Z=<z> Y=<y> X=<x> Array=<array_idx> Re=<real> Im=<imag>\n")
    return 0

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 extract_matrix_dumps.py <log_file_or_dir> <output_dir>")
        print("  If <log_file_or_dir> is a directory, looks for matrix_*.txt files there")
        print("  If it's a file, extracts from log or looks for dump files in same directory")
        sys.exit(1)
    
    input_path = sys.argv[1]
    output_dir = sys.argv[2]
    
    if not os.path.exists(input_path):
        print(f"Error: Input path not found: {input_path}")
        sys.exit(1)
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Extract/copy matrices
    copy_or_extract_matrix_dump(input_path, os.path.join(output_dir, "matrix_before_fft.txt"), "before")
    copy_or_extract_matrix_dump(input_path, os.path.join(output_dir, "matrix_after_fft.txt"), "after")
    
    print("Matrix extraction completed!")

if __name__ == "__main__":
    main()

