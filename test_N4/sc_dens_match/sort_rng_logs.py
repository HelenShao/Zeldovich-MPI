#!/usr/bin/env python3
"""
Sort RNG debug log files by Y, x, z coordinates for comparison.

The hermitian and zeldovich codes process Y-slices in different orders,
so logs must be sorted before comparison.
"""

import re
import sys
import argparse
import subprocess


def extract_coords(line):
    """Extract Y, x, z coordinates from a log line.
    
    Returns tuple (Y, x, z) for sorting, or (999, 999, 999) if not found.
    """
    match = re.search(r'Y=(\d+).*\(x,z\)=\((\d+),(\d+)\)', line)
    if match:
        return (int(match.group(1)), int(match.group(2)), int(match.group(3)))
    return (999, 999, 999)  # Put unmatched lines at end


def sort_log_file(input_file, output_file):
    """Sort a log file by Y, x, z coordinates."""
    try:
        with open(input_file, 'r') as f:
            lines = f.readlines()
        
        lines.sort(key=extract_coords)
        
        with open(output_file, 'w') as f:
            f.writelines(lines)
        
        print(f"Sorted {len(lines)} lines: {input_file} -> {output_file}")
        return True
    except FileNotFoundError:
        print(f"Error: File not found: {input_file}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"Error processing {input_file}: {e}", file=sys.stderr)
        return False


def main():
    parser = argparse.ArgumentParser(
        description='Sort RNG debug log files by Y, x, z coordinates'
    )
    parser.add_argument(
        '--hermitian', '-H',
        default='rng_logs/hermitian_rng_debug_filtered.txt',
        help='Input hermitian log file (default: rng_logs/hermitian_rng_debug_filtered.txt)'
    )
    parser.add_argument(
        '--zeldovich', '-z',
        default='rng_logs/zeldovich_rng_debug_filtered.txt',
        help='Input zeldovich log file (default: rng_logs/zeldovich_rng_debug_filtered.txt)'
    )
    parser.add_argument(
        '--output-dir', '-o',
        default='rng_logs',
        help='Output directory for sorted files (default: rng_logs)'
    )
    parser.add_argument(
        '--diff-output', '-d',
        default=None,
        help='Save diff output to this file (default: no diff)'
    )
    
    args = parser.parse_args()
    
    # Determine output filenames
    import os
    hermitian_base = os.path.basename(args.hermitian)
    zeldovich_base = os.path.basename(args.zeldovich)
    
    hermitian_out = os.path.join(args.output_dir, hermitian_base.replace('_filtered.txt', '_sorted.txt'))
    zeldovich_out = os.path.join(args.output_dir, zeldovich_base.replace('_filtered.txt', '_sorted.txt'))
    
    # Sort both files
    success1 = sort_log_file(args.hermitian, hermitian_out)
    success2 = sort_log_file(args.zeldovich, zeldovich_out)
    
    if success1 and success2:
        print(f"\nSorted files ready for comparison:")
        print(f"  {hermitian_out}")
        print(f"  {zeldovich_out}")
        
        # Run diff if requested
        if args.diff_output:
            try:
                result = subprocess.run(
                    ['diff', hermitian_out, zeldovich_out],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True
                )
                with open(args.diff_output, 'w') as f:
                    f.write(result.stdout)
                    f.write(result.stderr)
                
                if result.returncode == 0:
                    print(f"\nNo differences found (diff saved to {args.diff_output})")
                else:
                    print(f"\nDifferences found (saved to {args.diff_output})")
                    print(f"  {len(result.stdout.splitlines())} differing lines")
            except Exception as e:
                print(f"Error running diff: {e}", file=sys.stderr)
                return 1
        else:
            print(f"\nCompare with: diff {hermitian_out} {zeldovich_out}")
        
        return 0
    else:
        return 1


if __name__ == '__main__':
    sys.exit(main())

