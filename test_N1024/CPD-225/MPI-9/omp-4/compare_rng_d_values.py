#!/usr/bin/env python3
"""
Compare D values between hermitian and zeldovich RNG debug logs.

Extracts D values (random Gaussian draws) for each coordinate and shows
numerical differences when they don't match.

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
    """Parse a log line and extract coordinates, k values, and D values.
    
    Returns tuple: (Y, x, z, kx, ky, kz, D_real, D_imag) or None if parsing fails.
    """
    # Match: Y=(\d+) (x,z)=\((\d+),(\d+)\): k=\(([^,]+),([^,]+),([^)]+)\) ... D=\(([^,]+),([^)]+)\)
    match = re.search(
        r'Y=(\d+).*\(x,z\)=\((\d+),(\d+)\): k=\(([^,]+),([^,]+),([^)]+)\).*D=\(([^,]+),([^)]+)\)',
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
        return (Y, x, z, kx, ky, kz, D_real, D_imag)
    return None


def load_log_file(filename):
    """Load a log file and return a dictionary keyed by (Y, x, z) -> (kx, ky, kz, D_real, D_imag)."""
    data = {}
    try:
        with open(filename, 'r') as f:
            for line_num, line in enumerate(f, 1):
                parsed = parse_log_line(line)
                if parsed:
                    Y, x, z, kx, ky, kz, D_real, D_imag = parsed
                    key = (Y, x, z)
                    if key in data:
                        print(f"Warning: Duplicate coordinate {key} at line {line_num} in {filename}", 
                              file=sys.stderr)
                    data[key] = (kx, ky, kz, D_real, D_imag)
        return data
    except FileNotFoundError:
        print(f"Error: File not found: {filename}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"Error reading {filename}: {e}", file=sys.stderr)
        return None


def compare_d_values(hermitian_data, zeldovich_data, output_file=None, tolerance=1e-6):
    """Compare D values between the two datasets and report differences."""
    all_coords = set(hermitian_data.keys()) | set(zeldovich_data.keys())
    
    differences = []
    missing_in_zeldovich = []
    missing_in_hermitian = []
    
    for coord in sorted(all_coords):
        Y, x, z = coord
        
        if coord not in hermitian_data:
            missing_in_hermitian.append(coord)
            continue
        if coord not in zeldovich_data:
            missing_in_zeldovich.append(coord)
            continue
        
        kx_h, ky_h, kz_h, D_h_real, D_h_imag = hermitian_data[coord]
        kx_z, ky_z, kz_z, D_z_real, D_z_imag = zeldovich_data[coord]
        
        # Use k values from hermitian (should match, but use hermitian as reference)
        kx, ky, kz = kx_h, ky_h, kz_h
        
        # Calculate differences
        diff_real = D_h_real - D_z_real
        diff_imag = D_h_imag - D_z_imag
        diff_magnitude = (diff_real**2 + diff_imag**2)**0.5
        
        # Threshold for "effectively zero" (very small to catch exact zeros)
        zero_threshold = 1e-20
        
        # Check if one is zero and the other is not (special case)
        h_real_is_zero = abs(D_h_real) < zero_threshold
        h_imag_is_zero = abs(D_h_imag) < zero_threshold
        z_real_is_zero = abs(D_z_real) < zero_threshold
        z_imag_is_zero = abs(D_z_imag) < zero_threshold
        
        # Flag if one is zero and the other is not
        zero_mismatch_real = (h_real_is_zero and not z_real_is_zero) or (not h_real_is_zero and z_real_is_zero)
        zero_mismatch_imag = (h_imag_is_zero and not z_imag_is_zero) or (not h_imag_is_zero and z_imag_is_zero)
        
        # Check if different (accounting for floating point tolerance OR zero mismatch)
        if (abs(diff_real) > tolerance or abs(diff_imag) > tolerance) or zero_mismatch_real or zero_mismatch_imag:
            differences.append({
                'coord': coord,
                'k': (kx, ky, kz),
                'hermitian': (D_h_real, D_h_imag),
                'zeldovich': (D_z_real, D_z_imag),
                'diff': (diff_real, diff_imag),
                'magnitude': diff_magnitude,
                'zero_mismatch': zero_mismatch_real or zero_mismatch_imag
            })
    
    # Prepare output
    output_lines = []
    output_lines.append("=" * 80)
    output_lines.append("RNG D Value Comparison Report")
    output_lines.append("=" * 80)
    output_lines.append(f"\nTotal coordinates in hermitian: {len(hermitian_data)}")
    output_lines.append(f"Total coordinates in zeldovich: {len(zeldovich_data)}")
    output_lines.append(f"Total unique coordinates: {len(all_coords)}")
    output_lines.append(f"Tolerance: {tolerance}")
    output_lines.append("")
    
    if missing_in_hermitian:
        output_lines.append(f"WARNING: {len(missing_in_hermitian)} coordinates missing in hermitian:")
        for coord in sorted(missing_in_hermitian)[:10]:  # Show first 10
            output_lines.append(f"  Y={coord[0]} (x,z)=({coord[1]},{coord[2]})")
        if len(missing_in_hermitian) > 10:
            output_lines.append(f"  ... and {len(missing_in_hermitian) - 10} more")
        output_lines.append("")
    
    if missing_in_zeldovich:
        output_lines.append(f"WARNING: {len(missing_in_zeldovich)} coordinates missing in zeldovich:")
        for coord in sorted(missing_in_zeldovich)[:10]:  # Show first 10
            output_lines.append(f"  Y={coord[0]} (x,z)=({coord[1]},{coord[2]})")
        if len(missing_in_zeldovich) > 10:
            output_lines.append(f"  ... and {len(missing_in_zeldovich) - 10} more")
        output_lines.append("")
    
    if differences:
        # Count zero mismatches
        zero_mismatches = sum(1 for d in differences if d.get('zero_mismatch', False))
        
        output_lines.append(f"Found {len(differences)} coordinates with different D values:")
        if zero_mismatches > 0:
            output_lines.append(f"  ({zero_mismatches} have zero/non-zero mismatches)")
        output_lines.append("")
        output_lines.append(f"{'Y':<4} {'x':<4} {'z':<4} {'kx':<6} {'ky':<6} {'kz':<6} {'Hermitian D (real, imag)':<40} {'Zeldovich D (real, imag)':<40} {'Difference (real, imag)':<40} {'|Diff|':<15}")
        output_lines.append("-" * 240)
        
        # Sort by difference magnitude (largest first), but put zero mismatches first
        differences.sort(key=lambda x: (not x.get('zero_mismatch', False), -x['magnitude']))
        
        for diff in differences:
            Y, x, z = diff['coord']
            kx, ky, kz = diff['k']
            D_h_real, D_h_imag = diff['hermitian']
            D_z_real, D_z_imag = diff['zeldovich']
            diff_real, diff_imag = diff['diff']
            magnitude = diff['magnitude']
            is_zero_mismatch = diff.get('zero_mismatch', False)
            
            # Add marker for zero mismatches
            marker = " [ZERO-MISMATCH]" if is_zero_mismatch else ""
            
            output_lines.append(
                f"{Y:<4} {x:<4} {z:<4} {kx:<6} {ky:<6} {kz:<6} "
                f"({D_h_real:15.10e}, {D_h_imag:15.10e})  "
                f"({D_z_real:15.10e}, {D_z_imag:15.10e})  "
                f"({diff_real:15.10e}, {diff_imag:15.10e})  "
                f"{magnitude:15.10e}{marker}"
            )
        
        output_lines.append("")
        output_lines.append(f"Summary: {len(differences)} coordinates have different D values")
        output_lines.append(f"Largest difference magnitude: {differences[0]['magnitude']:.10e}")
        output_lines.append(f"Smallest difference magnitude: {differences[-1]['magnitude']:.10e}")
    else:
        output_lines.append("SUCCESS: All D values match (within tolerance)!")
    
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
        description='Compare D values between hermitian and zeldovich RNG debug logs'
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
        default='rng_logs/d_value_comparison.txt',
        help='Output file for comparison report (default: rng_logs/d_value_comparison.txt)'
    )
    parser.add_argument(
        '--tolerance', '-t',
        type=float,
        default=1e-6,
        help='Tolerance for floating point comparison (default: 1e-6 for single precision)'
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
    
    print(f"Comparing D values...")
    return compare_d_values(hermitian_data, zeldovich_data, args.output, args.tolerance)


if __name__ == '__main__':
    sys.exit(main())

