#!/usr/bin/env python3
"""
Compare two particle IC text files and compute differences in displacement and velocity fields.
Shows that differences are approximately zero for matching particles.
"""

import sys
import numpy as np

def parse_ic_file(filename):
    """
    Parse particle IC file and return dictionary keyed by (i, j, k).
    
    Returns:
        dict: {(i, j, k): {'displ': [Z, Y, X], 'vel': [Z, Y, X], 'idx': particle_idx}}
    """
    particles = {}
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip header lines
            if not line or line.startswith('#'):
                continue
            
            # Parse data line
            parts = line.split()
            if len(parts) < 10:
                continue
            
            try:
                idx = int(parts[0])
                i = int(parts[1])
                j = int(parts[2])
                k = int(parts[3])
                displ_z = float(parts[4])
                displ_y = float(parts[5])
                displ_x = float(parts[6])
                vel_z = float(parts[7])
                vel_y = float(parts[8])
                vel_x = float(parts[9])
                
                key = (i, j, k)
                particles[key] = {
                    'idx': idx,
                    'displ': [displ_z, displ_y, displ_x],  # Z, Y, X
                    'vel': [vel_z, vel_y, vel_x]  # Z, Y, X
                }
            except (ValueError, IndexError) as e:
                print(f"Warning: Could not parse line: {line}", file=sys.stderr)
                continue
    
    return particles

def compare_ic_files(file1, file2, output_file=None):
    """
    Compare two IC files and compute differences.
    """
    if output_file:
        fout = open(output_file, 'w')
        original_stdout = sys.stdout
        sys.stdout = fout
    
    try:
        print(f"Reading {file1}...")
        particles1 = parse_ic_file(file1)
        print(f"  Found {len(particles1)} particles")
        
        print(f"Reading {file2}...")
        particles2 = parse_ic_file(file2)
        print(f"  Found {len(particles2)} particles")
        
        # Find common particles
        common_keys = set(particles1.keys()) & set(particles2.keys())
        print(f"\nCommon particles: {len(common_keys)}")
        
        if len(common_keys) == 0:
            print("ERROR: No common particles found!", file=sys.stderr)
            return
        
        # Compute differences
        displ_diffs = {'Z': [], 'Y': [], 'X': []}
        vel_diffs = {'Z': [], 'Y': [], 'X': []}
        
        for key in sorted(common_keys):
            p1 = particles1[key]
            p2 = particles2[key]
            
            # Displacement differences (Z, Y, X)
            for comp_idx, comp_name in enumerate(['Z', 'Y', 'X']):
                diff = p1['displ'][comp_idx] - p2['displ'][comp_idx]
                displ_diffs[comp_name].append(diff)
                
                diff = p1['vel'][comp_idx] - p2['vel'][comp_idx]
                vel_diffs[comp_name].append(diff)
        
        # Convert to numpy arrays for statistics
        displ_diffs = {k: np.array(v) for k, v in displ_diffs.items()}
        vel_diffs = {k: np.array(v) for k, v in vel_diffs.items()}
        
        # Print statistics
        print("\n" + "="*80)
        print("DISPLACEMENT DIFFERENCES (Hermitian - Zeldovich)")
        print("="*80)
        print(f"{'Component':<12} {'Mean':<20} {'Std':<20} {'Max Abs':<20} {'Min':<20} {'Max':<20}")
        print("-"*80)
        
        for comp in ['Z', 'Y', 'X']:
            diff = displ_diffs[comp]
            mean_diff = np.mean(diff)
            std_diff = np.std(diff)
            max_abs_diff = np.max(np.abs(diff))
            min_diff = np.min(diff)
            max_diff = np.max(diff)
            
            print(f"{comp:<12} {mean_diff:20.12e} {std_diff:20.12e} {max_abs_diff:20.12e} {min_diff:20.12e} {max_diff:20.12e}")
        
        print("\n" + "="*80)
        print("VELOCITY DIFFERENCES (Hermitian - Zeldovich)")
        print("="*80)
        print(f"{'Component':<12} {'Mean':<20} {'Std':<20} {'Max Abs':<20} {'Min':<20} {'Max':<20}")
        print("-"*80)
        
        for comp in ['Z', 'Y', 'X']:
            diff = vel_diffs[comp]
            mean_diff = np.mean(diff)
            std_diff = np.std(diff)
            max_abs_diff = np.max(np.abs(diff))
            min_diff = np.min(diff)
            max_diff = np.max(diff)
            
            print(f"{comp:<12} {mean_diff:20.12e} {std_diff:20.12e} {max_abs_diff:20.12e} {min_diff:20.12e} {max_diff:20.12e}")
        
        # Print summary
        print("\n" + "="*80)
        print("SUMMARY")
        print("="*80)
        
        all_displ_diffs = np.concatenate([displ_diffs['Z'], displ_diffs['Y'], displ_diffs['X']])
        all_vel_diffs = np.concatenate([vel_diffs['Z'], vel_diffs['Y'], vel_diffs['X']])
        
        print(f"Total displacement differences: {len(all_displ_diffs)}")
        print(f"  Mean absolute difference: {np.mean(np.abs(all_displ_diffs)):.12e}")
        print(f"  Max absolute difference: {np.max(np.abs(all_displ_diffs)):.12e}")
        print(f"  Standard deviation: {np.std(all_displ_diffs):.12e}")
        
        print(f"\nTotal velocity differences: {len(all_vel_diffs)}")
        print(f"  Mean absolute difference: {np.mean(np.abs(all_vel_diffs)):.12e}")
        print(f"  Max absolute difference: {np.max(np.abs(all_vel_diffs)):.12e}")
        print(f"  Standard deviation: {np.std(all_vel_diffs):.12e}")
        
        # Check if differences are approximately zero
        # Use tolerance appropriate for single-precision float differences
        tolerance = 1e-5
        displ_max_abs = np.max(np.abs(all_displ_diffs))
        vel_max_abs = np.max(np.abs(all_vel_diffs))
        
        print(f"\n" + "="*80)
        print("VERIFICATION")
        print("="*80)
        print(f"Tolerance: {tolerance:.1e} (appropriate for single-precision float comparisons)")
        print(f"Max displacement difference: {displ_max_abs:.12e}")
        print(f"Max velocity difference: {vel_max_abs:.12e}")
        
        if displ_max_abs < tolerance and vel_max_abs < tolerance:
            print("✓ All differences are within tolerance (approximately zero)")
            print("  Differences are consistent with floating-point precision.")
        else:
            print("⚠ Some differences exceed tolerance")
            if displ_max_abs >= tolerance:
                print(f"  Displacement differences exceed tolerance by {displ_max_abs/tolerance:.2f}x")
            if vel_max_abs >= tolerance:
                print(f"  Velocity differences exceed tolerance by {vel_max_abs/tolerance:.2f}x")
            print("  Note: Differences on order of 1e-6 are typical for single-precision float comparisons.")
        
        # Show first few particle differences as examples
        print("\n" + "="*80)
        print("EXAMPLE PARTICLE DIFFERENCES (first 10)")
        print("="*80)
        print(f"{'i':<6} {'j':<6} {'k':<6} {'Δdispl_Z':<20} {'Δdispl_Y':<20} {'Δdispl_X':<20} {'Δvel_Z':<20} {'Δvel_Y':<20} {'Δvel_X':<20}")
        print("-"*80)
        
        count = 0
        for key in sorted(common_keys):
            if count >= 10:
                break
            i, j, k = key
            p1 = particles1[key]
            p2 = particles2[key]
            
            displ_diff = [p1['displ'][idx] - p2['displ'][idx] for idx in range(3)]
            vel_diff = [p1['vel'][idx] - p2['vel'][idx] for idx in range(3)]
            
            print(f"{i:<6} {j:<6} {k:<6} "
                  f"{displ_diff[0]:<20.12e} {displ_diff[1]:<20.12e} {displ_diff[2]:<20.12e} "
                  f"{vel_diff[0]:<20.12e} {vel_diff[1]:<20.12e} {vel_diff[2]:<20.12e}")
            count += 1
    
    finally:
        if output_file:
            sys.stdout = original_stdout
            fout.close()
            print(f"\nResults saved to {output_file}", file=sys.stderr)

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 compare_ic_differences.py <hermitian_file> <zeldovich_file> [output_file]")
        print("Example: python3 compare_ic_differences.py hermitian_ic_0_all_slabs.txt zeldovich_ic_0_all_slabs.txt")
        print("Example: python3 compare_ic_differences.py hermitian_ic_0_all_slabs.txt zeldovich_ic_0_all_slabs.txt ic_differences.txt")
        sys.exit(1)
    
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) > 3 else None
    
    compare_ic_files(file1, file2, output_file)

if __name__ == '__main__':
    main()

