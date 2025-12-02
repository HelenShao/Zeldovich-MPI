#!/usr/bin/env python3
"""
Compare RNG consistency across different N values.

This script extracts [RNG-DEBUG] lines from output files and compares
the D, F, G, H values for overlapping grid points.

Usage:
    python compare_rng_consistency.py output_N4.out output_N6.out
"""

import sys
import re
from collections import defaultdict

def parse_rng_debug_line(line):
    """Parse a [RNG-DEBUG] line and extract values."""
    # Format: [RNG-DEBUG] N=4 Y=0 (x,z)=(0,0): k=(0,0,0) k2=0.000000 | D=(...,...) F=(...,...) G=(...,...) H=(...,...)
    pattern = r'\[RNG-DEBUG\] N=(\d+) Y=(\d+) \(x,z\)=\((\d+),(\d+)\): k=\((-?\d+),(-?\d+),(-?\d+)\) k2=([\d.]+) \| D=\(([\d.e+-]+),([\d.e+-]+)\) F=\(([\d.e+-]+),([\d.e+-]+)\) G=\(([\d.e+-]+),([\d.e+-]+)\) H=\(([\d.e+-]+),([\d.e+-]+)\)'
    
    match = re.match(pattern, line)
    if not match:
        return None
    
    return {
        'N': int(match.group(1)),
        'Y': int(match.group(2)),
        'x': int(match.group(3)),
        'z': int(match.group(4)),
        'kx': int(match.group(5)),
        'ky': int(match.group(6)),
        'kz': int(match.group(7)),
        'k2': float(match.group(8)),
        'D': (float(match.group(9)), float(match.group(10))),
        'F': (float(match.group(11)), float(match.group(12))),
        'G': (float(match.group(13)), float(match.group(14))),
        'H': (float(match.group(15)), float(match.group(16))),
    }

def extract_rng_debug(filepath):
    """Extract all [RNG-DEBUG] lines from a file."""
    results = {}
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if '[RNG-DEBUG]' in line:
                    parsed = parse_rng_debug_line(line.strip())
                    if parsed:
                        key = (parsed['x'], parsed['Y'], parsed['z'])
                        results[key] = parsed
    except FileNotFoundError:
        print(f"Error: File not found: {filepath}", file=sys.stderr)
        sys.exit(1)
    
    return results

def compare_values(val1, val2, name, tol=1e-10):
    """Compare two values and return match status."""
    if isinstance(val1, tuple):
        # Compare complex numbers
        diff_re = abs(val1[0] - val2[0])
        diff_im = abs(val1[1] - val2[1])
        match = diff_re < tol and diff_im < tol
        return match, f"re_diff={diff_re:.2e} im_diff={diff_im:.2e}"
    else:
        diff = abs(val1 - val2)
        match = diff < tol
        return match, f"diff={diff:.2e}"

def main():
    if len(sys.argv) < 3:
        print("Usage: python compare_rng_consistency.py output_N4.out output_N6.out", file=sys.stderr)
        sys.exit(1)
    
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    
    print("=" * 80)
    print("RNG Consistency Comparison")
    print("=" * 80)
    print()
    
    # Extract debug lines from both files
    data1 = extract_rng_debug(file1)
    data2 = extract_rng_debug(file2)
    
    if not data1:
        print(f"Warning: No [RNG-DEBUG] lines found in {file1}")
        print("Make sure to compile with -DDEBUG_RNG_CONSISTENCY=1")
        sys.exit(1)
    
    if not data2:
        print(f"Warning: No [RNG-DEBUG] lines found in {file2}")
        print("Make sure to compile with -DDEBUG_RNG_CONSISTENCY=1")
        sys.exit(1)
    
    # Find overlapping coordinates
    coords1 = set(data1.keys())
    coords2 = set(data2.keys())
    overlapping = coords1 & coords2
    
    print(f"File 1 ({file1}): {len(data1)} debug entries")
    print(f"File 2 ({file2}): {len(data2)} debug entries")
    print(f"Overlapping coordinates: {len(overlapping)}")
    print()
    
    if not overlapping:
        print("No overlapping coordinates found!")
        print(f"File 1 coordinates: {sorted(coords1)}")
        print(f"File 2 coordinates: {sorted(coords2)}")
        sys.exit(1)
    
    # Compare overlapping coordinates
    print("Comparing overlapping coordinates:")
    print("-" * 80)
    
    all_match = True
    for coord in sorted(overlapping):
        d1 = data1[coord]
        d2 = data2[coord]
        
        print(f"\nCoordinate (x,y,z) = {coord}:")
        print(f"  N1={d1['N']}, N2={d2['N']}")
        print(f"  k1=({d1['kx']},{d1['ky']},{d1['kz']}), k2=({d2['kx']},{d2['ky']},{d2['kz']})")
        print(f"  k2_1={d1['k2']:.6f}, k2_2={d2['k2']:.6f}")
        
        # Note: k-vectors will differ for different N, so we expect different k2
        # But D, F, G, H should match if RNG is consistent (same random numbers used)
        
        # Compare D
        match_d, msg_d = compare_values(d1['D'], d2['D'], 'D')
        print(f"  D: {'MATCH' if match_d else 'DIFFER'}: {msg_d}")
        print(f"    D1 (N={d1['N']}) = ({d1['D'][0]:.3g}, {d1['D'][1]:.3g})")
        print(f"    D2 (N={d2['N']}) = ({d2['D'][0]:.3g}, {d2['D'][1]:.3g})")
        if not match_d:
            all_match = False
        
        # Compare F
        match_f, msg_f = compare_values(d1['F'], d2['F'], 'F')
        print(f"  F: {'MATCH' if match_f else 'DIFFER'}: {msg_f}")
        print(f"    F1 (N={d1['N']}) = ({d1['F'][0]:.3g}, {d1['F'][1]:.3g})")
        print(f"    F2 (N={d2['N']}) = ({d2['F'][0]:.3g}, {d2['F'][1]:.3g})")
        if not match_f:
            all_match = False
        
        # Compare G
        match_g, msg_g = compare_values(d1['G'], d2['G'], 'G')
        print(f"  G: {'MATCH' if match_g else 'DIFFER'}: {msg_g}")
        print(f"    G1 (N={d1['N']}) = ({d1['G'][0]:.3g}, {d1['G'][1]:.3g})")
        print(f"    G2 (N={d2['N']}) = ({d2['G'][0]:.3g}, {d2['G'][1]:.3g})")
        if not match_g:
            all_match = False
        
        # Compare H
        match_h, msg_h = compare_values(d1['H'], d2['H'], 'H')
        print(f"  H: {'MATCH' if match_h else 'DIFFER'}: {msg_h}")
        print(f"    H1 (N={d1['N']}) = ({d1['H'][0]:.3g}, {d1['H'][1]:.3g})")
        print(f"    H2 (N={d2['N']}) = ({d2['H'][0]:.3g}, {d2['H'][1]:.3g})")
        if not match_h:
            all_match = False
    
    print()
    print("=" * 80)
    if all_match:
        print("RESULT: All overlapping coordinates match! RNG consistency verified.")
    else:
        print("RESULT: Some values differ. This may be expected if:")
        print("  - k-vectors differ (different N means different k-space sampling)")
        print("  - Power spectrum P(k) is evaluated at different k values")
        print("  - RNG consistency ensures same random numbers, but final values differ due to k-dependence")
    print("=" * 80)

if __name__ == '__main__':
    main()

