#!/usr/bin/env python3
"""
Quantitative analysis of RNG consistency test results
Extracts statistics from compare_256_vs_512_res512.txt
"""

import re
import sys

def analyze_rng_consistency(filename):
    """Analyze RNG consistency comparison results."""
    
    with open(filename, 'r') as f:
        content = f.read()
    
    lines = content.split('\n')
    
    # Statistics
    stats = {
        'D': {'matches': 0, 'differs': 0, 'diffs': []},
        'F': {'matches': 0, 'differs': 0, 'diffs': []},
        'G': {'matches': 0, 'differs': 0, 'diffs': []},
        'H': {'matches': 0, 'differs': 0, 'diffs': []},
    }
    
    coordinates = 0
    
    for i, line in enumerate(lines):
        # Count coordinates
        if 'Coordinate (x,y,z)' in line:
            coordinates += 1
        
        # Parse match/difference for each field
        for field in ['D', 'F', 'G', 'H']:
            pattern = rf'{field}:\s+(MATCH|DIFFER)'
            match = re.search(pattern, line)
            if match:
                if match.group(1) == 'MATCH':
                    stats[field]['matches'] += 1
                else:
                    stats[field]['differs'] += 1
                    
                    # Extract difference values
                    diff_match = re.search(r're_diff=([\d.e+-]+)\s+im_diff=([\d.e+-]+)', line)
                    if diff_match:
                        re_diff = float(diff_match.group(1))
                        im_diff = float(diff_match.group(2))
                        total_diff = (re_diff**2 + im_diff**2)**0.5
                        stats[field]['diffs'].append(total_diff)
    
    # Print results
    print("=" * 80)
    print("RNG Consistency Test Results")
    print("=" * 80)
    print(f"\nTest configuration:")
    print(f"  Comparing: N=256 vs N=512")
    print(f"  Total coordinates compared: {coordinates}")
    print(f"  Tolerance: 1e-10 (double precision)")
    print(f"  Spline resolution: 512")
    print()
    
    print("Match Statistics (per coordinate):")
    print("-" * 80)
    
    total_matches = 0
    total_checks = 0
    
    for field in ['D', 'F', 'G', 'H']:
        matches = stats[field]['matches']
        differs = stats[field]['differs']
        total = matches + differs
        match_rate = 100.0 * matches / total if total > 0 else 0
        
        print(f"\n{field} ({'Density' if field=='D' else f'{field}-displacement'}):")
        print(f"  Matches: {matches}/{total} ({match_rate:.2f}%)")
        print(f"  Differences: {differs}/{total} ({100-match_rate:.2f}%)")
        
        if stats[field]['diffs']:
            diffs = stats[field]['diffs']
            print(f"  Difference statistics (when non-zero):")
            print(f"    Count: {len(diffs)}")
            print(f"    Mean: {sum(diffs)/len(diffs):.6e}")
            print(f"    Median: {sorted(diffs)[len(diffs)//2]:.6e}")
            print(f"    Max: {max(diffs):.6e}")
            print(f"    Min: {min(diffs):.6e}")
        
        total_matches += matches
        total_checks += total
    
    print("\n" + "=" * 80)
    print("Overall RNG Consistency:")
    print("-" * 80)
    overall_rate = 100.0 * total_matches / total_checks if total_checks > 0 else 0
    print(f"Total matches: {total_matches}/{total_checks} ({overall_rate:.2f}%)")
    print(f"Total differences: {total_checks - total_matches}/{total_checks} ({100-overall_rate:.2f}%)")
    
    print("\n" + "=" * 80)
    print("Interpretation:")
    print("-" * 80)
    if overall_rate >= 99.0:
        print("✓ Excellent RNG consistency (≥99% match)")
        print("  Remaining differences likely due to:")
        print("    - Numerical precision in spline interpolation")
        print("    - Different k-vector sampling (P(k) evaluated at different k)")
        print("    - Minor rounding differences in power spectrum evaluation")
    elif overall_rate >= 95.0:
        print("⚠ Good RNG consistency (≥95% match)")
        print("  Some differences may indicate:")
        print("    - Spline interpolation accuracy")
        print("    - k-space sampling differences")
    else:
        print("✗ Poor RNG consistency (<95% match)")
        print("  Significant differences detected - may indicate RNG issues")
    
    return {
        'coordinates': coordinates,
        'overall_rate': overall_rate,
        'total_matches': total_matches,
        'total_checks': total_checks,
        'stats': stats
    }

if __name__ == '__main__':
    filename = 'compare_256_vs_512_res512.txt'
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    
    analyze_rng_consistency(filename)

