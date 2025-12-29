#!/usr/bin/env python3
"""
Compare particle IC files from two different locations.
Reads all particles and compares by (i, j, k) indices.
"""

import struct
import sys
import os
from collections import defaultdict

def read_all_particles(filename):
    """Read all particles from a binary IC file."""
    try:
        with open(filename, 'rb') as f:
            data = f.read()
    except FileNotFoundError:
        return None, f"File not found: {filename}"
    except Exception as e:
        return None, f"Error reading file: {e}"
    
    particle_size = 32
    file_size = len(data)
    max_particles = file_size // particle_size
    
    if max_particles == 0:
        return None, "File too small or empty"
    
    particles = {}
    
    for i in range(max_particles):
        offset = i * particle_size
        if offset + particle_size > file_size:
            break
        
        try:
            i_val, j_val, k_val, padding, displ0, displ1, displ2, vel0, vel1, vel2 = struct.unpack(
                'HHH H fff fff', data[offset:offset+particle_size]
            )
            key = (i_val, j_val, k_val)
            particles[key] = {
                'i': i_val, 'j': j_val, 'k': k_val,
                'displ': [displ0, displ1, displ2],
                'vel': [vel0, vel1, vel2]
            }
        except struct.error as e:
            return None, f"Error unpacking particle {i}: {e}"
    
    return particles, max_particles

def compare_particles(particles1, particles2, label1="Set 1", label2="Set 2"):
    """Compare two sets of particles by (i, j, k) indices."""
    keys1 = set(particles1.keys())
    keys2 = set(particles2.keys())
    
    common_keys = keys1 & keys2
    only_in_1 = keys1 - keys2
    only_in_2 = keys2 - keys1
    
    print(f"\n{'='*80}")
    print(f"COMPARISON: {label1} vs {label2}")
    print(f"{'='*80}")
    print(f"Total particles in {label1}: {len(keys1)}")
    print(f"Total particles in {label2}: {len(keys2)}")
    print(f"Common particles (same i,j,k): {len(common_keys)}")
    print(f"Only in {label1}: {len(only_in_1)}")
    print(f"Only in {label2}: {len(only_in_2)}")
    
    if len(common_keys) == 0:
        print("\nWARNING: No particles with matching (i, j, k) indices found!")
        return
    
    # Compare values for common particles
    exact_matches = 0
    close_matches = 0  # Within 1e-5
    differences = []
    max_diff_displ = [0.0, 0.0, 0.0]
    max_diff_vel = [0.0, 0.0, 0.0]
    max_diff_key = None
    
    for key in sorted(common_keys)[:100]:  # Check first 100 common particles
        p1 = particles1[key]
        p2 = particles2[key]
        
        # Check if values match exactly
        if (p1['displ'] == p2['displ'] and p1['vel'] == p2['vel']):
            exact_matches += 1
        else:
            # Calculate differences
            diff_displ = [abs(p1['displ'][i] - p2['displ'][i]) for i in range(3)]
            diff_vel = [abs(p1['vel'][i] - p2['vel'][i]) for i in range(3)]
            
            max_diff_displ = [max(max_diff_displ[i], diff_displ[i]) for i in range(3)]
            max_diff_vel = [max(max_diff_vel[i], diff_vel[i]) for i in range(3)]
            
            if sum(diff_displ) + sum(diff_vel) > sum(max_diff_displ) + sum(max_diff_vel) - 1e-5:
                max_diff_key = key
            
            if all(d < 1e-5 for d in diff_displ + diff_vel):
                close_matches += 1
            else:
                differences.append({
                    'key': key,
                    'displ_diff': diff_displ,
                    'vel_diff': diff_vel,
                    'p1': p1,
                    'p2': p2
                })
    
    print(f"\nValue comparison (first 100 common particles):")
    print(f"  Exact matches: {exact_matches}")
    print(f"  Close matches (diff < 1e-5): {close_matches}")
    print(f"  Different values: {len(differences)}")
    
    if max_diff_displ[0] > 0 or max_diff_vel[0] > 0:
        print(f"\nMaximum differences found:")
        print(f"  Displacement: [{max_diff_displ[0]:.6f}, {max_diff_displ[1]:.6f}, {max_diff_displ[2]:.6f}]")
        print(f"  Velocity: [{max_diff_vel[0]:.6f}, {max_diff_vel[1]:.6f}, {max_diff_vel[2]:.6f}]")
        
        if max_diff_key:
            p1 = particles1[max_diff_key]
            p2 = particles2[max_diff_key]
            print(f"\nExample particle with large difference (i={max_diff_key[0]}, j={max_diff_key[1]}, k={max_diff_key[2]}):")
            print(f"  {label1}: displ={p1['displ']}, vel={p1['vel']}")
            print(f"  {label2}: displ={p2['displ']}, vel={p2['vel']}")
    
    # Check for normalization/scaling factors
    if len(differences) > 0:
        print(f"\nChecking for normalization/scaling factors...")
        sample_keys = list(common_keys)[:min(50, len(common_keys))]
        ratios_displ = []
        ratios_vel = []
        
        for key in sample_keys:
            p1 = particles1[key]
            p2 = particles2[key]
            
            # Calculate ratios (avoid division by zero)
            for i in range(3):
                if abs(p2['displ'][i]) > 1e-10:
                    ratios_displ.append(p1['displ'][i] / p2['displ'][i])
                if abs(p2['vel'][i]) > 1e-10:
                    ratios_vel.append(p1['vel'][i] / p2['vel'][i])
        
        if ratios_displ:
            avg_ratio_displ = sum(ratios_displ) / len(ratios_displ)
            print(f"  Average displacement ratio ({label1}/{label2}): {avg_ratio_displ:.6f}")
        
        if ratios_vel:
            avg_ratio_vel = sum(ratios_vel) / len(ratios_vel)
            print(f"  Average velocity ratio ({label1}/{label2}): {avg_ratio_vel:.6f}")
    
    # Show some example differences
    if len(differences) > 0:
        print(f"\nFirst 5 particles with different values:")
        for i, diff in enumerate(differences[:5]):
            key = diff['key']
            print(f"  Particle (i={key[0]}, j={key[1]}, k={key[2]}):")
            print(f"    {label1}: displ={diff['p1']['displ']}, vel={diff['p1']['vel']}")
            print(f"    {label2}: displ={diff['p2']['displ']}, vel={diff['p2']['vel']}")
            print(f"    Diff: displ={diff['displ_diff']}, vel={diff['vel_diff']}")

def main():
    # Accept command-line arguments if provided, otherwise use defaults
    if len(sys.argv) >= 3:
        dir1 = sys.argv[1]
        dir2 = sys.argv[2]
    else:
        # Default paths (PLT enabled)
        dir1 = "/home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_assembly/particle_ics"
        dir2 = "/home/helenshao/InitialConditions/zeldovich-PLT/output"
    
    print("Reading particles from both locations...")
    
    # Read all files from directory 1
    particles1_all = {}
    files1 = sorted([f for f in os.listdir(dir1) if f.startswith('ic_')])
    print(f"\nDirectory 1 ({dir1}):")
    for f in files1:
        filepath = os.path.join(dir1, f)
        particles, count = read_all_particles(filepath)
        if particles is None:
            print(f"  {f}: ERROR - {count}")
        else:
            print(f"  {f}: {count} particles")
            particles1_all.update(particles)
    
    # Read all files from directory 2
    particles2_all = {}
    files2 = sorted([f for f in os.listdir(dir2) if f.startswith('ic_')])
    print(f"\nDirectory 2 ({dir2}):")
    for f in files2:
        filepath = os.path.join(dir2, f)
        particles, count = read_all_particles(filepath)
        if particles is None:
            print(f"  {f}: ERROR - {count}")
        else:
            print(f"  {f}: {count} particles")
            particles2_all.update(particles)
    
    # Compare
    compare_particles(particles1_all, particles2_all, 
                     label1="hermitian_3d_matrix_production", 
                     label2="zeldovich-PLT")

if __name__ == '__main__':
    main()

