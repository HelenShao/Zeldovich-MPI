#!/usr/bin/env python3
"""
Compare D, F, G, H values between hermitian and zeldovich RNG debug logs.

Extracts D, F, G, H values (random Gaussian draws and PLT corrections) for each coordinate
and shows numerical differences when they don't match.

Tolerance notes:
- For single precision (float): default 1e-6 (machine epsilon ~1.19e-7)
- For double precision (double): use 1e-10 or tighter
- The default 1e-6 is appropriate for single precision comparisons
"""

import re
import sys
import argparse
from collections import defaultdict


def parse_log_line(line):
    """Parse a log line and extract coordinates, k values, and D, F, G, H values.
    
    Returns tuple: (Y, x, z, kx, ky, kz, D_real, D_imag, F_real, F_imag, G_real, G_imag, H_real, H_imag) or None if parsing fails.
    """
    # Match: Y=(\d+) (x,z)=\((\d+),(\d+)\): k=\(([^,]+),([^,]+),([^)]+)\) ... D=\(([^,]+),([^)]+)\) F=\(([^,]+),([^)]+)\) G=\(([^,]+),([^)]+)\) H=\(([^,]+),([^)]+)\)
    match = re.search(
        r'Y=(\d+).*\(x,z\)=\((\d+),(\d+)\): k=\(([^,]+),([^,]+),([^)]+)\).*D=\(([^,]+),([^)]+)\).*F=\(([^,]+),([^)]+)\).*G=\(([^,]+),([^)]+)\).*H=\(([^,]+),([^)]+)\)',
        line
    )
    if match:
        Y = int(match.group(1))
        x = int(match.group(2))
        z = int(match.group(3))
        kx = int(match.group(4))
        ky = int(match.group(5))
        kz = int(match.group(6))
        D_real = float(match.group(7))
        D_imag = float(match.group(8))
        F_real = float(match.group(9))
        F_imag = float(match.group(10))
        G_real = float(match.group(11))
        G_imag = float(match.group(12))
        H_real = float(match.group(13))
        H_imag = float(match.group(14))
        return (Y, x, z, kx, ky, kz, D_real, D_imag, F_real, F_imag, G_real, G_imag, H_real, H_imag)
    return None


def load_log_file(filename):
    """Load a log file and return a dictionary keyed by (Y, x, z) -> (kx, ky, kz, D, F, G, H values)."""
    data = {}
    try:
        with open(filename, 'r') as f:
            for line_num, line in enumerate(f, 1):
                parsed = parse_log_line(line)
                if parsed:
                    Y, x, z, kx, ky, kz, D_real, D_imag, F_real, F_imag, G_real, G_imag, H_real, H_imag = parsed
                    key = (Y, x, z)
                    if key in data:
                        print(f"Warning: Duplicate coordinate {key} at line {line_num} in {filename}", 
                              file=sys.stderr)
                    data[key] = (kx, ky, kz, 
                                 (D_real, D_imag), (F_real, F_imag), 
                                 (G_real, G_imag), (H_real, H_imag))
        return data
    except FileNotFoundError:
        print(f"Error: File not found: {filename}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"Error reading {filename}: {e}", file=sys.stderr)
        return None


def compare_values(hermitian_data, zeldovich_data, output_file=None, tolerance=1e-6, compare_which='all'):
    """Compare D, F, G, H values between the two datasets and report differences.
    
    compare_which: 'all', 'D', 'F', 'G', 'H', or comma-separated list like 'D,F,G'
    """
    all_coords = set(hermitian_data.keys()) | set(zeldovich_data.keys())
    
    # Parse compare_which
    if compare_which == 'all':
        compare_list = ['D', 'F', 'G', 'H']
    else:
        compare_list = [v.strip().upper() for v in compare_which.split(',')]
    
    differences = []
    missing_in_zeldovich = []
    missing_in_hermitian = []
    
    # Threshold for "effectively zero"
    zero_threshold = 1e-20
    
    for coord in sorted(all_coords):
        Y, x, z = coord
        
        if coord not in hermitian_data:
            missing_in_hermitian.append(coord)
            continue
        if coord not in zeldovich_data:
            missing_in_zeldovich.append(coord)
            continue
        
        kx_h, ky_h, kz_h, D_h, F_h, G_h, H_h = hermitian_data[coord]
        kx_z, ky_z, kz_z, D_z, F_z, G_z, H_z = zeldovich_data[coord]
        
        # Use k values from hermitian (should match, but use hermitian as reference)
        kx, ky, kz = kx_h, ky_h, kz_h
        
        # Compare each requested value type
        coord_diffs = {}
        for val_type, h_val, z_val in [('D', D_h, D_z), ('F', F_h, F_z), ('G', G_h, G_z), ('H', H_h, H_z)]:
            if val_type not in compare_list:
                continue
                
            h_real, h_imag = h_val
            z_real, z_imag = z_val
            
            diff_real = h_real - z_real
            diff_imag = h_imag - z_imag
            diff_magnitude = (diff_real**2 + diff_imag**2)**0.5
            
            # Check if one is zero and the other is not
            h_real_is_zero = abs(h_real) < zero_threshold
            h_imag_is_zero = abs(h_imag) < zero_threshold
            z_real_is_zero = abs(z_real) < zero_threshold
            z_imag_is_zero = abs(z_imag) < zero_threshold
            
            zero_mismatch_real = (h_real_is_zero and not z_real_is_zero) or (not h_real_is_zero and z_real_is_zero)
            zero_mismatch_imag = (h_imag_is_zero and not z_imag_is_zero) or (not h_imag_is_zero and z_imag_is_zero)
            
            # Check if different
            if (abs(diff_real) > tolerance or abs(diff_imag) > tolerance) or zero_mismatch_real or zero_mismatch_imag:
                coord_diffs[val_type] = {
                    'hermitian': (h_real, h_imag),
                    'zeldovich': (z_real, z_imag),
                    'diff': (diff_real, diff_imag),
                    'magnitude': diff_magnitude,
                    'zero_mismatch': zero_mismatch_real or zero_mismatch_imag
                }
        
        if coord_diffs:
            differences.append({
                'coord': coord,
                'k': (kx, ky, kz),
                'diffs': coord_diffs
            })
    
    # Prepare output
    output_lines = []
    output_lines.append("=" * 80)
    output_lines.append("RNG D, F, G, H Value Comparison Report")
    output_lines.append("=" * 80)
    output_lines.append(f"\nTotal coordinates in hermitian: {len(hermitian_data)}")
    output_lines.append(f"Total coordinates in zeldovich: {len(zeldovich_data)}")
    output_lines.append(f"Total unique coordinates: {len(all_coords)}")
    output_lines.append(f"Tolerance: {tolerance}")
    output_lines.append(f"Comparing: {', '.join(compare_list)}")
    output_lines.append("")
    
    if missing_in_hermitian:
        output_lines.append(f"WARNING: {len(missing_in_hermitian)} coordinates missing in hermitian:")
        for coord in sorted(missing_in_hermitian)[:10]:
            output_lines.append(f"  Y={coord[0]} (x,z)=({coord[1]},{coord[2]})")
        if len(missing_in_hermitian) > 10:
            output_lines.append(f"  ... and {len(missing_in_hermitian) - 10} more")
        output_lines.append("")
    
    if missing_in_zeldovich:
        output_lines.append(f"WARNING: {len(missing_in_zeldovich)} coordinates missing in zeldovich:")
        for coord in sorted(missing_in_zeldovich)[:10]:
            output_lines.append(f"  Y={coord[0]} (x,z)=({coord[1]},{coord[2]})")
        if len(missing_in_zeldovich) > 10:
            output_lines.append(f"  ... and {len(missing_in_zeldovich) - 10} more")
        output_lines.append("")
    
    if differences:
        # Count differences by type
        diff_counts = defaultdict(int)
        zero_mismatch_counts = defaultdict(int)
        max_magnitudes = defaultdict(float)
        
        for diff in differences:
            for val_type, val_diff in diff['diffs'].items():
                diff_counts[val_type] += 1
                if val_diff['zero_mismatch']:
                    zero_mismatch_counts[val_type] += 1
                max_magnitudes[val_type] = max(max_magnitudes[val_type], val_diff['magnitude'])
        
        output_lines.append(f"Found {len(differences)} coordinates with different values:")
        for val_type in compare_list:
            if val_type in diff_counts:
                output_lines.append(f"  {val_type}: {diff_counts[val_type]} differences (max magnitude: {max_magnitudes[val_type]:.10e})")
                if zero_mismatch_counts[val_type] > 0:
                    output_lines.append(f"    ({zero_mismatch_counts[val_type]} have zero/non-zero mismatches)")
        output_lines.append("")
        
        # Show detailed differences (first 50)
        output_lines.append(f"{'Y':<4} {'x':<4} {'z':<4} {'Type':<4} {'Hermitian (real, imag)':<40} {'Zeldovich (real, imag)':<40} {'Difference (real, imag)':<40} {'|Diff|':<15}")
        output_lines.append("-" * 240)
        
        # Sort by maximum difference magnitude
        differences.sort(key=lambda x: -max(v['magnitude'] for v in x['diffs'].values()))
        
        shown = 0
        for diff in differences:
            if shown >= 50:
                break
            Y, x, z = diff['coord']
            kx, ky, kz = diff['k']
            
            for val_type in sorted(diff['diffs'].keys()):
                if shown >= 50:
                    break
                val_diff = diff['diffs'][val_type]
                h_real, h_imag = val_diff['hermitian']
                z_real, z_imag = val_diff['zeldovich']
                diff_real, diff_imag = val_diff['diff']
                magnitude = val_diff['magnitude']
                is_zero_mismatch = val_diff['zero_mismatch']
                
                marker = " [ZERO-MISMATCH]" if is_zero_mismatch else ""
                
                output_lines.append(
                    f"{Y:<4} {x:<4} {z:<4} {val_type:<4} "
                    f"({h_real:15.10e}, {h_imag:15.10e})  "
                    f"({z_real:15.10e}, {z_imag:15.10e})  "
                    f"({diff_real:15.10e}, {diff_imag:15.10e})  "
                    f"{magnitude:15.10e}{marker}"
                )
                shown += 1
        
        if len(differences) > 50:
            output_lines.append(f"\n... and {len(differences) - 50} more coordinates with differences")
        
        output_lines.append("")
        output_lines.append(f"Summary: {len(differences)} coordinates have different values")
    else:
        output_lines.append("SUCCESS: All values match (within tolerance)!")
    
    output_lines.append("")
    output_lines.append("=" * 80)
    
    output_text = "\n".join(output_lines)
    
    # Write to file if specified, otherwise print to stdout
    if output_file:
        try:
            with open(output_file, 'w') as f:
                f.write(output_text)
            print(f"Comparison report saved to: {output_file}")
            print(f"Found {len(differences)} differences out of {len(all_coords)} coordinates")
        except Exception as e:
            print(f"Error writing to {output_file}: {e}", file=sys.stderr)
            print(output_text)
            return 1
    else:
        print(output_text)
    
    return 0 if not differences else 1


def main():
    parser = argparse.ArgumentParser(
        description='Compare D, F, G, H values between hermitian and zeldovich RNG debug logs'
    )
    parser.add_argument(
        '--hermitian', '-H',
        default='rng_logs/hermitian_rng_debug_sorted.txt',
        help='Hermitian sorted log file (default: rng_logs/hermitian_rng_debug_sorted.txt)'
    )
    parser.add_argument(
        '--zeldovich', '-z',
        default='rng_logs/zeldovich_rng_debug_sorted.txt',
        help='Zeldovich sorted log file (default: rng_logs/zeldovich_rng_debug_sorted.txt)'
    )
    parser.add_argument(
        '--output', '-o',
        default='rng_logs/all_values_comparison.txt',
        help='Output file for comparison report (default: rng_logs/all_values_comparison.txt)'
    )
    parser.add_argument(
        '--tolerance', '-t',
        type=float,
        default=1e-6,
        help='Tolerance for floating point comparison (default: 1e-6 for single precision)'
    )
    parser.add_argument(
        '--compare', '-c',
        default='all',
        help='Which values to compare: all, D, F, G, H, or comma-separated like D,F,G (default: all)'
    )
    
    args = parser.parse_args()
    
    # Load both files
    print(f"Loading hermitian log: {args.hermitian}")
    hermitian_data = load_log_file(args.hermitian)
    if hermitian_data is None:
        return 1
    
    print(f"Loading zeldovich log: {args.zeldovich}")
    zeldovich_data = load_log_file(args.zeldovich)
    if zeldovich_data is None:
        return 1
    
    print(f"Comparing values ({args.compare})...")
    return compare_values(hermitian_data, zeldovich_data, args.output, args.tolerance, args.compare)


if __name__ == '__main__':
    sys.exit(main())

