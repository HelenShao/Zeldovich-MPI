#!/usr/bin/env python3
"""
In-depth analysis of D value differences between hermitian and zeldovich codes.
"""

import re
from collections import defaultdict

def parse_comparison_file(filename):
    """Parse the comparison file and extract all differences."""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    differences = []
    missing_hermitian = []
    missing_zeldovich = []
    in_table = False
    
    for i, line in enumerate(lines):
        if 'coordinates missing in hermitian' in line:
            # Parse missing coordinates
            j = i + 1
            while j < len(lines) and lines[j].strip() and not lines[j].startswith('WARNING'):
                match = re.search(r'Y=(\d+)\s+\(x,z\)=\((\d+),(\d+)\)', lines[j])
                if match:
                    missing_hermitian.append((int(match.group(1)), int(match.group(2)), int(match.group(3))))
                j += 1
        
        if 'coordinates missing in zeldovich' in line:
            # Parse missing coordinates
            j = i + 1
            while j < len(lines) and lines[j].strip() and not lines[j].startswith('Found'):
                match = re.search(r'Y=(\d+)\s+\(x,z\)=\((\d+),(\d+)\)', lines[j])
                if match:
                    missing_zeldovich.append((int(match.group(1)), int(match.group(2)), int(match.group(3))))
                j += 1
        
        if 'Y    x    z' in line:
            in_table = True
            continue
        
        if in_table and line.strip() and not line.startswith('-') and not line.startswith('Summary'):
            parts = line.split()
            if len(parts) >= 10:
                try:
                    Y = int(parts[0])
                    x = int(parts[1])
                    z = int(parts[2])
                    kx = int(parts[3])
                    ky = int(parts[4])
                    kz = int(parts[5])
                    
                    # Parse D values
                    hermitian_match = re.search(r'\(([^,]+),([^)]+)\)', parts[6] + ' ' + parts[7])
                    zeldovich_match = re.search(r'\(([^,]+),([^)]+)\)', parts[8] + ' ' + parts[9])
                    
                    if hermitian_match and zeldovich_match:
                        D_h_real = float(hermitian_match.group(1))
                        D_h_imag = float(hermitian_match.group(2))
                        D_z_real = float(zeldovich_match.group(1))
                        D_z_imag = float(zeldovich_match.group(2))
                        
                        is_zero_mismatch = '[ZERO-MISMATCH]' in line
                        
                        differences.append({
                            'Y': Y, 'x': x, 'z': z,
                            'kx': kx, 'ky': ky, 'kz': kz,
                            'D_h': (D_h_real, D_h_imag),
                            'D_z': (D_z_real, D_z_imag),
                            'zero_mismatch': is_zero_mismatch
                        })
                except Exception as e:
                    pass
    
    return differences, missing_hermitian, missing_zeldovich

def analyze_patterns(differences, missing_hermitian, missing_zeldovich):
    """Perform comprehensive pattern analysis."""
    
    print("=" * 80)
    print("IN-DEPTH DIFFERENCE ANALYSIS")
    print("=" * 80)
    
    print(f"\n1. OVERALL STATISTICS")
    print(f"   Total differences: {len(differences)}")
    print(f"   Zero/non-zero mismatches: {sum(1 for d in differences if d['zero_mismatch'])}")
    print(f"   Missing in hermitian: {len(missing_hermitian)}")
    print(f"   Missing in zeldovich: {len(missing_zeldovich)}")
    
    # Group by Y
    print(f"\n2. DIFFERENCES BY Y SLICE")
    by_y = defaultdict(list)
    for d in differences:
        by_y[d['Y']].append(d)
    for y in sorted(by_y.keys()):
        zero_mismatches = sum(1 for d in by_y[y] if d['zero_mismatch'])
        print(f"   Y={y}: {len(by_y[y])} differences ({zero_mismatches} zero-mismatches)")
    
    # Y=0 specific analysis
    print(f"\n3. Y=0 SLICE ANALYSIS")
    y0_diffs = [d for d in differences if d['Y'] == 0]
    print(f"   Total Y=0 differences: {len(y0_diffs)}")
    
    if y0_diffs:
        # Check which ones are zero mismatches
        y0_zero_mismatch = [d for d in y0_diffs if d['zero_mismatch']]
        y0_nonzero = [d for d in y0_diffs if not d['zero_mismatch']]
        print(f"   Zero-mismatches: {len(y0_zero_mismatch)}")
        print(f"   Non-zero differences: {len(y0_nonzero)}")
        
        # Analyze k-space patterns for Y=0
        print(f"\n   Y=0 k-space patterns:")
        by_kx_y0 = defaultdict(int)
        by_kz_y0 = defaultdict(int)
        for d in y0_diffs:
            by_kx_y0[d['kx']] += 1
            by_kz_y0[d['kz']] += 1
        
        print(f"   Top kx values: {sorted(by_kx_y0.items(), key=lambda x: x[1], reverse=True)[:5]}")
        print(f"   Top kz values: {sorted(by_kz_y0.items(), key=lambda x: x[1], reverse=True)[:5]}")
        
        # Check if Y=0 differences are mostly at high k
        high_k_y0 = [d for d in y0_diffs if abs(d['kx']) >= 6 or abs(d['kz']) >= 6]
        print(f"   High k (|kx|>=6 or |kz|>=6): {len(high_k_y0)} out of {len(y0_diffs)}")
    
    # Missing coordinates analysis
    print(f"\n4. MISSING COORDINATES ANALYSIS")
    print(f"   Missing in hermitian (Y=0): {sum(1 for m in missing_hermitian if m[0] == 0)}")
    print(f"   Missing in zeldovich (Y=8): {sum(1 for m in missing_zeldovich if m[0] == 8)}")
    
    # Group missing by Y
    missing_h_by_y = defaultdict(list)
    missing_z_by_y = defaultdict(list)
    for m in missing_hermitian:
        missing_h_by_y[m[0]].append(m)
    for m in missing_zeldovich:
        missing_z_by_y[m[0]].append(m)
    
    print(f"   Missing coordinates by Y (hermitian):")
    for y in sorted(missing_h_by_y.keys()):
        print(f"     Y={y}: {len(missing_h_by_y[y])} coordinates")
    
    print(f"   Missing coordinates by Y (zeldovich):")
    for y in sorted(missing_z_by_y.keys()):
        print(f"     Y={y}: {len(missing_z_by_y[y])} coordinates")
    
    # k-space patterns
    print(f"\n5. K-SPACE PATTERNS")
    by_kx = defaultdict(int)
    by_ky = defaultdict(int)
    by_kz = defaultdict(int)
    for d in differences:
        by_kx[d['kx']] += 1
        by_ky[d['ky']] += 1
        by_kz[d['kz']] += 1
    
    print(f"   Top 10 kx values: {sorted(by_kx.items(), key=lambda x: x[1], reverse=True)[:10]}")
    print(f"   Top 10 ky values: {sorted(by_ky.items(), key=lambda x: x[1], reverse=True)[:10]}")
    print(f"   Top 10 kz values: {sorted(by_kz.items(), key=lambda x: x[1], reverse=True)[:10]}")
    
    # Check for Nyquist frequency patterns
    print(f"\n6. NYQUIST FREQUENCY ANALYSIS")
    # For N=16, Nyquist is at k=8
    nyquist_k = 8
    nyquist_diffs = [d for d in differences if abs(d['kx']) == nyquist_k or 
                     abs(d['ky']) == nyquist_k or abs(d['kz']) == nyquist_k]
    print(f"   Differences at Nyquist frequencies (|k|=8): {len(nyquist_diffs)}")
    
    # Check zero mismatches at Nyquist
    nyquist_zero_mismatch = [d for d in nyquist_diffs if d['zero_mismatch']]
    print(f"   Zero-mismatches at Nyquist: {len(nyquist_zero_mismatch)}")
    
    # High k analysis
    print(f"\n7. HIGH K ANALYSIS")
    high_k_diffs = [d for d in differences if abs(d['kx']) >= 6 or 
                    abs(d['ky']) >= 6 or abs(d['kz']) >= 6]
    print(f"   Differences at high k (|k|>=6): {len(high_k_diffs)}")
    high_k_zero = [d for d in high_k_diffs if d['zero_mismatch']]
    print(f"   Zero-mismatches at high k: {len(high_k_zero)}")
    
    # Sample actual differences (non-zero-mismatch)
    print(f"\n8. SAMPLE NON-ZERO-MISMATCH DIFFERENCES (first 10)")
    nonzero_diffs = [d for d in differences if not d['zero_mismatch']]
    for i, d in enumerate(nonzero_diffs[:10]):
        print(f"   {i+1}. Y={d['Y']} (x,z)=({d['x']},{d['z']}) k=({d['kx']},{d['ky']},{d['kz']})")
        print(f"      Hermitian: ({d['D_h'][0]:.6e}, {d['D_h'][1]:.6e})")
        print(f"      Zeldovich:  ({d['D_z'][0]:.6e}, {d['D_z'][1]:.6e})")
    
    # Check for coordinate-specific patterns
    print(f"\n9. COORDINATE PATTERNS")
    # Check if differences cluster at certain x,z values
    by_xz = defaultdict(int)
    for d in differences:
        by_xz[(d['x'], d['z'])] += 1
    
    print(f"   Coordinates with most differences:")
    for (x, z), count in sorted(by_xz.items(), key=lambda x: x[1], reverse=True)[:10]:
        print(f"     (x,z)=({x},{z}): {count} differences")
    
    print("\n" + "=" * 80)

if __name__ == '__main__':
    import sys
    
    filename = 'rng_logs/d_value_comparison.txt'
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    
    differences, missing_hermitian, missing_zeldovich = parse_comparison_file(filename)
    analyze_patterns(differences, missing_hermitian, missing_zeldovich)

