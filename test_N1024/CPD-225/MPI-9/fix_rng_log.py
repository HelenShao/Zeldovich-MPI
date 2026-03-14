#!/usr/bin/env python3
"""
Fix corrupted RNG debug log files from MPI runs.

When multiple MPI ranks write to stderr simultaneously, lines can get
concatenated together. This script:
1. Splits concatenated lines (e.g., [FACTOR-DEBUG]...[RNG-DEBUG])
2. Filters out [FACTOR-DEBUG] lines (not needed for RNG comparison)
3. Extracts only [RNG-DEBUG] lines
"""

import re
import sys

def is_valid_rng_debug_line(line):
    """Check if a line is a valid, complete RNG-DEBUG line.
    
    Returns (True, cleaned_line) if valid, (False, None) otherwise.
    Also cleans corrupted lines by extracting only the first valid RNG-DEBUG portion.
    """
    if not line.startswith('[RNG-DEBUG]'):
        return False, None
    
    # Check for corruption: line contains another [RNG-DEBUG] or [FACTOR-DEBUG]
    first_rng_idx = line.find('[RNG-DEBUG]')
    second_rng_idx = line.find('[RNG-DEBUG]', first_rng_idx + 1)
    factor_idx = line.find('[FACTOR-DEBUG]')
    
    # If line contains another tag, truncate at that point
    truncate_at = len(line)
    if second_rng_idx > 0:
        truncate_at = min(truncate_at, second_rng_idx)
    if factor_idx > 0:
        truncate_at = min(truncate_at, factor_idx)
    
    # Also check for pattern "k=(" appearing twice (sign of concatenation)
    first_k_idx = line.find('k=(')
    if first_k_idx > 0:
        second_k_idx = line.find('k=(', first_k_idx + 1)
        if second_k_idx > 0:
            truncate_at = min(truncate_at, second_k_idx)
    
    cleaned_line = line[:truncate_at].rstrip()
    
    # Check for all required components in the cleaned line
    required = ['Y=', '(x,z)=', 'k=(', 'D=(', 'F=(', 'G=(', 'H=(']
    if not all(comp in cleaned_line for comp in required):
        return False, None
    
    # Validate it can be parsed by the comparison script's regex
    import re
    pattern = r'Y=(\d+).*\(x,z\)=\((\d+),(\d+)\): k=\(([^,]+),([^,]+),([^)]+)\).*D=\(([^,]+),([^)]+)\)'
    match = re.search(pattern, cleaned_line)
    if not match:
        return False, None
    
    # Validate that all extracted values are valid
    try:
        Y = int(match.group(1))
        x = int(match.group(2))
        z = int(match.group(3))
        kx = int(match.group(4))
        ky = int(match.group(5))
        kz = int(match.group(6))
        D_real = float(match.group(7))
        D_imag = float(match.group(8))
        # Check they're reasonable (not NaN, not Inf, coordinates in valid range)
        if not (0 <= Y < 256 and 0 <= x < 256 and 0 <= z < 256):
            return False, None
        if not (-1e10 < D_real < 1e10 and -1e10 < D_imag < 1e10):
            return False, None
    except (ValueError, OverflowError):
        return False, None
    
    return True, cleaned_line


def fix_log_file(input_file, output_file):
    """Fix corrupted log file by splitting concatenated lines and filtering."""
    rng_debug_lines = []
    
    try:
        with open(input_file, 'r') as f:
            content = f.read()
        
        # Split by [RNG-DEBUG] pattern to extract all RNG-DEBUG lines
        # This handles cases where lines are concatenated
        parts = re.split(r'(\[RNG-DEBUG\])', content)
        
        # Reconstruct RNG-DEBUG lines
        i = 0
        while i < len(parts):
            if parts[i] == '[RNG-DEBUG]' and i + 1 < len(parts):
                # Found [RNG-DEBUG] tag, reconstruct the line
                line = '[RNG-DEBUG]' + parts[i + 1]
                # Extract up to the next [FACTOR-DEBUG] or [RNG-DEBUG] or end of line
                # Remove any trailing [FACTOR-DEBUG] content
                if '[FACTOR-DEBUG]' in line:
                    # Split at [FACTOR-DEBUG] and take only the RNG-DEBUG part
                    line = line.split('[FACTOR-DEBUG]')[0]
                # Remove newlines and clean up
                line = line.rstrip()
                # Only add if it's a valid, complete RNG-DEBUG line
                valid, cleaned = is_valid_rng_debug_line(line)
                if valid:
                    rng_debug_lines.append(cleaned)
                i += 2
            else:
                i += 1
        
        # Also check for any standalone [RNG-DEBUG] lines we might have missed
        # Process line by line for simpler cases
        with open(input_file, 'r') as f:
            for line in f:
                line = line.rstrip()
                # Skip [FACTOR-DEBUG] only lines
                if line.startswith('[FACTOR-DEBUG]') and '[RNG-DEBUG]' not in line:
                    continue
                # Handle corrupted lines with both tags
                if '[FACTOR-DEBUG]' in line and '[RNG-DEBUG]' in line:
                    # Extract RNG-DEBUG portion
                    rng_part = line.split('[RNG-DEBUG]')
                    if len(rng_part) > 1:
                        fixed = '[RNG-DEBUG]' + rng_part[1]
                        # Remove any trailing [FACTOR-DEBUG] content
                        if '[FACTOR-DEBUG]' in fixed:
                            fixed = fixed.split('[FACTOR-DEBUG]')[0]
                        fixed = fixed.rstrip()
                        valid, cleaned = is_valid_rng_debug_line(fixed)
                        if valid:
                            if cleaned not in rng_debug_lines:
                                rng_debug_lines.append(cleaned)
                elif line.startswith('[RNG-DEBUG]'):
                    # Normal RNG-DEBUG line - validate it's complete
                    valid, cleaned = is_valid_rng_debug_line(line)
                    if valid:
                        if cleaned not in rng_debug_lines:
                            rng_debug_lines.append(cleaned)
        
        # Remove duplicates while preserving order
        seen = set()
        unique_lines = []
        for line in rng_debug_lines:
            # Use first 100 chars as key to identify duplicates (should be unique per coordinate)
            key = line[:100] if len(line) > 100 else line
            if key not in seen:
                seen.add(key)
                unique_lines.append(line)
        
        # Write cleaned output
        with open(output_file, 'w') as f:
            for line in unique_lines:
                f.write(line + '\n')
        
        print(f"Fixed log file: {len(unique_lines)} unique RNG-DEBUG lines extracted")
        print(f"  (removed {len(rng_debug_lines) - len(unique_lines)} duplicates)")
        print(f"Output written to: {output_file}")
        return True
        
    except FileNotFoundError:
        print(f"Error: File not found: {input_file}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"Error processing {input_file}: {e}", file=sys.stderr)
        return False


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Fix corrupted RNG debug log files from MPI runs'
    )
    parser.add_argument(
        '--input', '-i',
        default='rng_logs/hermitian_rng_debug.log',
        help='Input log file (default: rng_logs/hermitian_rng_debug.log)'
    )
    parser.add_argument(
        '--output', '-o',
        default='rng_logs/hermitian_rng_debug_fixed.txt',
        help='Output fixed log file (default: rng_logs/hermitian_rng_debug_fixed.txt)'
    )
    
    args = parser.parse_args()
    
    if fix_log_file(args.input, args.output):
        sys.exit(0)
    else:
        sys.exit(1)
