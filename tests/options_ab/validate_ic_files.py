#!/usr/bin/env python3
"""
Validate all particle initial condition files.
Checks that i, j, k indices and displacements/velocities are reasonable.

Usage:
    ./validate_ic_files.py [output_directory] [N]

Examples:
    ./validate_ic_files.py output_option_a 256
    ./validate_ic_files.py output_option_b 256
    ./validate_ic_files.py output_option_a     # Will try to infer N from files
"""

import struct
import sys
import os
import glob
import math
from collections import defaultdict

def read_all_particles(filename):
    """Read all particles from a binary particle file."""
    try:
        with open(filename, 'rb') as f:
            data = f.read()
    except FileNotFoundError:
        return None, f"File not found: {filename}"
    except Exception as e:
        return None, f"Error reading file: {e}"
    
    # RVZelParticle format (single precision):
    #   unsigned short i, j, k (3 * 2 = 6 bytes) + 2 bytes padding = 8 bytes
    #   float displ[3] (3 * 4 = 12 bytes) 
    #   float vel[3] (3 * 4 = 12 bytes)
    #   Total: 32 bytes per particle
    particle_size = 32
    file_size = len(data)
    num_particles = file_size // particle_size
    
    if num_particles == 0:
        return None, "File too small or empty"
    
    if file_size % particle_size != 0:
        # Warn but continue - may have small padding at end
        pass
    
    particles = []
    for i in range(num_particles):
        offset = i * particle_size
        try:
            # Unpack: unsigned short i, j, k (6 bytes) + skip 2 padding bytes + 3 floats + 3 floats
            i_val, j_val, k_val = struct.unpack('HHH', data[offset:offset+6])
            # Skip 2 bytes padding (offset+6 to offset+8)
            displ0, displ1, displ2 = struct.unpack('fff', data[offset+8:offset+20])
            vel0, vel1, vel2 = struct.unpack('fff', data[offset+20:offset+32])
            particles.append({
                'i': i_val, 'j': j_val, 'k': k_val,
                'displ': [displ0, displ1, displ2],
                'vel': [vel0, vel1, vel2]
            })
        except struct.error as e:
            return None, f"Error unpacking particle {i}: {e}"
    
    return particles, None

def is_reasonable_float(value, max_magnitude=1e6):
    """Check if a float value is reasonable (not NaN, Inf, or too large)."""
    if math.isnan(value) or math.isinf(value):
        return False
    return abs(value) < max_magnitude

def is_reasonable_index(value, N):
    """Check if an index value is in reasonable range [0, N-1]."""
    return 0 <= value < N

def validate_particles(particles, N, filename):
    """Validate all particles in a file."""
    errors = []
    warnings = []
    stats = {
        'total': len(particles),
        'i_range': [float('inf'), float('-inf')],
        'j_range': [float('inf'), float('-inf')],
        'k_range': [float('inf'), float('-inf')],
        'displ_mag_max': 0.0,
        'displ_mag_mean': 0.0,
        'vel_mag_max': 0.0,
        'vel_mag_mean': 0.0,
        'invalid_displ': 0,
        'invalid_vel': 0,
        'out_of_range_i': 0,
        'out_of_range_j': 0,
        'out_of_range_k': 0,
    }
    
    if len(particles) == 0:
        errors.append("File contains no particles")
        return errors, warnings, stats
    
    displ_mags = []
    vel_mags = []
    
    for idx, p in enumerate(particles):
        # Check indices are in valid range [0, N-1]
        # Count all errors but only report first few per file to avoid spam
        i_ok = is_reasonable_index(p['i'], N)
        j_ok = is_reasonable_index(p['j'], N)
        k_ok = is_reasonable_index(p['k'], N)
        
        if not i_ok:
            if stats['out_of_range_i'] < 3:  # Only report first 3 per file
                errors.append(f"Particle {idx}: i={p['i']} out of range [0, {N-1}]")
            stats['out_of_range_i'] += 1
        
        if not j_ok:
            if stats['out_of_range_j'] < 3:
                errors.append(f"Particle {idx}: j={p['j']} out of range [0, {N-1}]")
            stats['out_of_range_j'] += 1
        
        if not k_ok:
            if stats['out_of_range_k'] < 3:
                errors.append(f"Particle {idx}: k={p['k']} out of range [0, {N-1}]")
            stats['out_of_range_k'] += 1
        
        # Update ranges
        stats['i_range'][0] = min(stats['i_range'][0], p['i'])
        stats['i_range'][1] = max(stats['i_range'][1], p['i'])
        stats['j_range'][0] = min(stats['j_range'][0], p['j'])
        stats['j_range'][1] = max(stats['j_range'][1], p['j'])
        stats['k_range'][0] = min(stats['k_range'][0], p['k'])
        stats['k_range'][1] = max(stats['k_range'][1], p['k'])
        
        # Check displacements
        for d in p['displ']:
            if not is_reasonable_float(d):
                errors.append(f"Particle {idx}: Invalid displacement component {d}")
                stats['invalid_displ'] += 1
                break
        
        # Check velocities
        for v in p['vel']:
            if not is_reasonable_float(v):
                errors.append(f"Particle {idx}: Invalid velocity component {v}")
                stats['invalid_vel'] += 1
                break
        
        # Calculate magnitudes (skip if any component is invalid)
        try:
            if all(is_reasonable_float(d) for d in p['displ']):
                displ_mag = math.sqrt(sum(d*d for d in p['displ']))
                if not (math.isnan(displ_mag) or math.isinf(displ_mag)):
                    displ_mags.append(displ_mag)
        except:
            pass
        
        try:
            if all(is_reasonable_float(v) for v in p['vel']):
                vel_mag = math.sqrt(sum(v*v for v in p['vel']))
                if not (math.isnan(vel_mag) or math.isinf(vel_mag)):
                    vel_mags.append(vel_mag)
        except:
            pass
    
    # Compute statistics
    if displ_mags:
        stats['displ_mag_max'] = max(displ_mags)
        stats['displ_mag_mean'] = sum(displ_mags) / len(displ_mags)
    
    if vel_mags:
        stats['vel_mag_max'] = max(vel_mags)
        stats['vel_mag_mean'] = sum(vel_mags) / len(vel_mags)
    
    # Check for suspicious values
    if stats['displ_mag_max'] > 100.0:
        warnings.append(f"Large displacement magnitude: {stats['displ_mag_max']:.6f}")
    
    if stats['vel_mag_max'] > 100.0:
        warnings.append(f"Large velocity magnitude: {stats['vel_mag_max']:.6f}")
    
    # Calculate percentage of valid indices
    total_particles = stats['total']
    invalid_particles = max(stats['out_of_range_i'], stats['out_of_range_j'], stats['out_of_range_k'])
    valid_pct = 100.0 * (total_particles - invalid_particles) / total_particles if total_particles > 0 else 0.0
    
    if valid_pct < 50.0 and total_particles > 10:
        errors.append(f"Only {valid_pct:.1f}% of particles have valid indices (i,j,k all in [0,{N-1}])")
    
    return errors, warnings, stats

def infer_N_from_files(file_list):
    """Try to infer N from file names or particle data."""
    # Try to extract from filename pattern ic_rank*_i*_k*_*
    # Or try to find maximum i, j, k values from reading files
    max_i = -1
    max_j = -1
    max_k = -1
    
    # Sample a few files to find max indices
    sample_files = file_list[:min(10, len(file_list))]
    for filename in sample_files:
        particles, err = read_all_particles(filename)
        if particles is None:
            continue
        for p in particles[:100]:  # Sample first 100 particles
            max_i = max(max_i, p['i'])
            max_j = max(max_j, p['j'])
            max_k = max(max_k, p['k'])
    
    if max_i >= 0:
        # N should be max_index + 1
        inferred_N = max(max_i, max_j, max_k) + 1
        return inferred_N
    
    return None

def main():
    if len(sys.argv) < 2:
        print("Usage: validate_ic_files.py <output_directory> [N]")
        print("Example: validate_ic_files.py output_option_a 256")
        sys.exit(1)
    
    output_dir = sys.argv[1]
    N = int(sys.argv[2]) if len(sys.argv) > 2 else None
    
    if not os.path.isdir(output_dir):
        print(f"Error: Directory not found: {output_dir}")
        sys.exit(1)
    
    # Find all IC files
    pattern = os.path.join(output_dir, "ic_*")
    files = sorted(glob.glob(pattern))
    
    if len(files) == 0:
        print(f"Error: No IC files found in {output_dir}")
        sys.exit(1)
    
    print("=" * 80)
    print(f"VALIDATING IC FILES: {output_dir}")
    print("=" * 80)
    print(f"Found {len(files)} files")
    
    # Infer N if not provided
    if N is None:
        print("N not provided, attempting to infer from files...")
        N = infer_N_from_files(files)
        if N is None:
            print("Error: Could not infer N from files. Please provide N explicitly.")
            sys.exit(1)
        print(f"Inferred N = {N}")
    
    print(f"Using N = {N}")
    print("=" * 80)
    print()
    
    # Statistics across all files
    total_files = 0
    total_particles = 0
    total_errors = 0
    total_warnings = 0
    file_errors = defaultdict(list)
    file_warnings = defaultdict(list)
    
    global_stats = {
        'i_range': [float('inf'), float('-inf')],
        'j_range': [float('inf'), float('-inf')],
        'k_range': [float('inf'), float('-inf')],
        'displ_mag_max': 0.0,
        'displ_mag_mean': 0.0,
        'vel_mag_max': 0.0,
        'vel_mag_mean': 0.0,
    }
    
    # Validate each file
    for filename in files:
        basename = os.path.basename(filename)
        particles, err = read_all_particles(filename)
        
        if particles is None:
            print(f"ERROR: {basename}: {err}")
            file_errors[basename].append(err)
            total_errors += 1
            continue
        
        errors, warnings, stats = validate_particles(particles, N, filename)
        
        total_files += 1
        total_particles += stats['total']
        
        if errors:
            file_errors[basename].extend(errors)
            total_errors += len(errors)
        
        if warnings:
            file_warnings[basename].extend(warnings)
            total_warnings += len(warnings)
        
        # Update global stats
        global_stats['i_range'][0] = min(global_stats['i_range'][0], stats['i_range'][0])
        global_stats['i_range'][1] = max(global_stats['i_range'][1], stats['i_range'][1])
        global_stats['j_range'][0] = min(global_stats['j_range'][0], stats['j_range'][0])
        global_stats['j_range'][1] = max(global_stats['j_range'][1], stats['j_range'][1])
        global_stats['k_range'][0] = min(global_stats['k_range'][0], stats['k_range'][0])
        global_stats['k_range'][1] = max(global_stats['k_range'][1], stats['k_range'][1])
        global_stats['displ_mag_max'] = max(global_stats['displ_mag_max'], stats['displ_mag_max'])
        global_stats['vel_mag_max'] = max(global_stats['vel_mag_max'], stats['vel_mag_max'])
    
    # Print summary
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total files processed: {total_files}")
    print(f"Total particles: {total_particles}")
    print(f"Total errors: {total_errors}")
    print(f"Total warnings: {total_warnings}")
    print()
    
    print("Global Statistics:")
    print(f"  i index range: [{int(global_stats['i_range'][0])}, {int(global_stats['i_range'][1])}] (expected [0, {N-1}])")
    print(f"  j index range: [{int(global_stats['j_range'][0])}, {int(global_stats['j_range'][1])}] (expected [0, {N-1}])")
    print(f"  k index range: [{int(global_stats['k_range'][0])}, {int(global_stats['k_range'][1])}] (expected [0, {N-1}])")
    print(f"  Max displacement magnitude: {global_stats['displ_mag_max']:.6f}")
    print(f"  Max velocity magnitude: {global_stats['vel_mag_max']:.6f}")
    print()
    
    # Print errors if any (limit output)
    if file_errors:
        print("=" * 80)
        print("ERRORS FOUND:")
        print("=" * 80)
        error_count = 0
        max_files_to_show = 10
        for filename, errors in sorted(file_errors.items()):
            if error_count >= max_files_to_show:
                remaining = len(file_errors) - max_files_to_show
                print(f"\n... and {remaining} more files with errors (total {total_errors} errors)")
                break
            print(f"\n{filename}:")
            for err in errors[:5]:  # Limit to first 5 errors per file
                print(f"  - {err}")
            if len(errors) > 5:
                print(f"  ... and {len(errors) - 5} more errors in this file")
            error_count += 1
        print()
    
    # Print warnings if any
    if file_warnings:
        print("=" * 80)
        print("WARNINGS:")
        print("=" * 80)
        for filename, warnings in sorted(file_warnings.items()):
            print(f"\n{filename}:")
            for warn in warnings:
                print(f"  - {warn}")
        print()
    
    # Final verdict
    print("=" * 80)
    if total_errors == 0:
        print("✓ VALIDATION PASSED: All files are valid")
        sys.exit(0)
    else:
        print(f"✗ VALIDATION FAILED: {total_errors} errors found")
        sys.exit(1)

if __name__ == '__main__':
    main()

