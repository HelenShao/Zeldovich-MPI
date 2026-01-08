#!/usr/bin/env python3
"""
Compare real(FFT(D + i*F)) and imag(FFT(D + i*F)) values from debug output for D+i*F test.
In this test, D+i*F is explicitly stored in just_density mode, so the result should be purely real.
"""

import re
import sys
import os

def parse_real_fft_debug(line):
    """Parse [REAL-FFT-DEBUG] line."""
    pattern = r'\[REAL-FFT-DEBUG\] N=(\d+) just_density=(\d+) Z=(\d+) \(x,y\)=\((\d+),(\d+)\): real\(Array0\)=([\d.e+-]+) imag\(Array0\)=([\d.e+-]+)'
    match = re.search(pattern, line)
    if match:
        return {
            'N': int(match.group(1)),
            'just_density': int(match.group(2)),
            'Z': int(match.group(3)),
            'x': int(match.group(4)),
            'y': int(match.group(5)),
            'real': float(match.group(6)),
            'imag': float(match.group(7))
        }
    return None

def load_debug_file(filename):
    """Load and parse debug file."""
    values = {}
    if not os.path.exists(filename):
        print(f"Warning: File not found: {filename}", file=sys.stderr)
        return values
    with open(filename, 'r') as f:
        for line in f:
            parsed = parse_real_fft_debug(line)
            if parsed:
                key = (parsed['Z'], parsed['x'], parsed['y'])
                values[key] = parsed
    return values

def compare_files(file1, file2, output_file=None):
    """Compare two debug files."""
    values1 = load_debug_file(file1)
    values2 = load_debug_file(file2)
    
    if output_file:
        out = open(output_file, 'w')
    else:
        out = sys.stdout
    
    # Find common coordinates
    common_keys = set(values1.keys()) & set(values2.keys())
    only_file1 = set(values1.keys()) - set(values2.keys())
    only_file2 = set(values2.keys()) - set(values1.keys())
    
    print("=" * 70, file=out)
    print("D+i*F Test: REAL-FFT-DEBUG Comparison", file=out)
    print("=" * 70, file=out)
    print(f"File 1 (Hermitian): {file1}", file=out)
    print(f"  Total entries: {len(values1)}", file=out)
    print(f"File 2 (Zeldovich): {file2}", file=out)
    print(f"  Total entries: {len(values2)}", file=out)
    print(f"Common coordinates: {len(common_keys)}", file=out)
    print(f"Only in file 1: {len(only_file1)}", file=out)
    print(f"Only in file 2: {len(only_file2)}", file=out)
    print("", file=out)
    
    if len(common_keys) == 0:
        print("ERROR: No common coordinates found!", file=out)
        if output_file:
            out.close()
        return
    
    # Compare values
    max_real_diff = 0.0
    max_imag_diff = 0.0
    max_real_rel_diff = 0.0
    max_imag_rel_diff = 0.0
    num_errors = 0
    tolerance = 1e-10
    
    print("=" * 70, file=out)
    print("Comparison Results (D+i*F test: result should be purely REAL)", file=out)
    print("=" * 70, file=out)
    print(f"{'Z':>3} {'x':>3} {'y':>3} | {'Hermitian Real':>15} {'Hermitian Imag':>15} | {'Zeldovich Real':>15} {'Zeldovich Imag':>15} | {'Real Diff':>12} {'Imag Diff':>12}", file=out)
    print("-" * 70, file=out)
    
    for key in sorted(common_keys):
        v1 = values1[key]
        v2 = values2[key]
        
        real_diff = abs(v1['real'] - v2['real'])
        imag_diff = abs(v1['imag'] - v2['imag'])
        
        # Relative differences
        real_rel = real_diff / (abs(v1['real']) + 1e-20) if abs(v1['real']) > 1e-20 else real_diff
        imag_rel = imag_diff / (abs(v1['imag']) + 1e-20) if abs(v1['imag']) > 1e-20 else imag_diff
        
        max_real_diff = max(max_real_diff, real_diff)
        max_imag_diff = max(max_imag_diff, imag_diff)
        max_real_rel_diff = max(max_real_rel_diff, real_rel)
        max_imag_rel_diff = max(max_imag_rel_diff, imag_rel)
        
        if real_diff > tolerance or imag_diff > tolerance:
            num_errors += 1
            marker = " ***"
        else:
            marker = ""
        
        print(f"{v1['Z']:3d} {v1['x']:3d} {v1['y']:3d} | {v1['real']:15.10e} {v1['imag']:15.10e} | {v2['real']:15.10e} {v2['imag']:15.10e} | {real_diff:12.6e} {imag_diff:12.6e}{marker}", file=out)
    
    print("", file=out)
    print("=" * 70, file=out)
    print("Summary", file=out)
    print("=" * 70, file=out)
    print(f"Max absolute real difference: {max_real_diff:.6e}", file=out)
    print(f"Max absolute imag difference: {max_imag_diff:.6e}", file=out)
    print(f"Max relative real difference: {max_real_rel_diff:.6e}", file=out)
    print(f"Max relative imag difference: {max_imag_rel_diff:.6e}", file=out)
    print(f"Number of coordinates with differences > {tolerance}: {num_errors}", file=out)
    print("", file=out)
    
    # Check if result is purely real (imaginary should be ~0)
    print("=" * 70, file=out)
    print("Purely Real Check (D+i*F test)", file=out)
    print("=" * 70, file=out)
    imag_tolerance = 1e-10
    hermitian_imag_max = max(abs(v['imag']) for v in values1.values())
    zeldovich_imag_max = max(abs(v['imag']) for v in values2.values())
    print(f"Max |imag(Array0)| in Hermitian: {hermitian_imag_max:.6e}", file=out)
    print(f"Max |imag(Array0)| in Zeldovich: {zeldovich_imag_max:.6e}", file=out)
    if hermitian_imag_max < imag_tolerance and zeldovich_imag_max < imag_tolerance:
        print("✓ Both codes produce purely real results (imag ≈ 0)", file=out)
    else:
        print("✗ Results are NOT purely real (imaginary part is non-zero)", file=out)
    print("", file=out)
    
    if output_file:
        out.close()
        print(f"Comparison saved to: {output_file}")

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    rng_log_dir = os.path.join(script_dir, "rng_logs")
    
    hermitian_file = os.path.join(rng_log_dir, "hermitian_real_fft_debug.txt")
    zeldovich_file = os.path.join(rng_log_dir, "zeldovich_real_fft_debug.txt")
    output_file = os.path.join(script_dir, "real_fft_comparison_D_iF.txt")
    
    if len(sys.argv) > 1:
        hermitian_file = sys.argv[1]
    if len(sys.argv) > 2:
        zeldovich_file = sys.argv[2]
    if len(sys.argv) > 3:
        output_file = sys.argv[3]
    
    compare_files(hermitian_file, zeldovich_file, output_file)

