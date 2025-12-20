#!/usr/bin/env python3
"""
Compare final output matrices between N=256 and N=512 runs
Tests RNG consistency by comparing real-space output in overlapping region

This tests the complete pipeline: RNG → Fourier → IFFT → output
For overlapping region (x,y,z < 256), output should match if RNG is consistent

Usage:
    python compare_final_outputs.py --N256-dir /path/to/N256/output --N512-dir /path/to/N512/output
"""

import numpy as np
import glob
import os
import sys
import argparse

def read_zslab_binary(filename):
    """
    Read Z-slab binary file from hermitian_3d_matrix output
    
    Format: [Array][Y][X] as complex128
    Returns: array of shape (narray, N, x_count) or None if file doesn't exist
    """
    if not os.path.exists(filename):
        return None
    
    data = np.fromfile(filename, dtype=np.complex128)
    return data

def extract_full_z_slab(z, output_dir, N, narray=4):
    """
    Extract full Z-slab from all rank files by reassembling
    
    Strategy:
    1. Find all files for this Z-slab across all ranks
    2. Read each file and determine its x_count (from file size)
    3. Reassemble assuming contiguous X-distribution sorted by rank number
       (This works for small N where X is typically distributed contiguously)
    
    Returns: array of shape (narray, N, N) or None if files not found
    """
    pattern = os.path.join(output_dir, f"rank_*/z{z}_slab_N{N}.bin")
    files = sorted(glob.glob(pattern))
    
    if len(files) == 0:
        return None
    
    # Collect data from all ranks
    rank_data = []
    total_x_count = 0
    
    for filename in files:
        # Parse rank from path: rank_*/z*_slab_N*.bin
        rank_dir = os.path.dirname(filename)
        try:
            rank = int(os.path.basename(rank_dir).split('_')[1])
        except (ValueError, IndexError):
            continue
        
        # Read file
        data_raw = read_zslab_binary(filename)
        if data_raw is None:
            continue
        
        # Infer x_count from data size
        x_count = len(data_raw) // (narray * N)
        if len(data_raw) != narray * N * x_count:
            continue
        
        data = data_raw.reshape(narray, N, x_count)
        rank_data.append((rank, data, x_count))
        total_x_count += x_count
    
    if len(rank_data) == 0:
        return None
    
    # Sort by rank number (assuming contiguous X-distribution)
    rank_data.sort(key=lambda x: x[0])
    
    # Verify total X count matches N (or is close for verification)
    if total_x_count != N:
        # This might happen due to rounding, but warn if too far off
        if abs(total_x_count - N) > 2:
            print(f"Warning: Z={z} total x_count={total_x_count} != N={N}")
    
    # Reassemble full slab
    full_slab = np.zeros((narray, N, N), dtype=np.complex128)
    x_start = 0
    
    for rank, data, x_count in rank_data:
        x_end = min(x_start + x_count, N)
        actual_x_count = x_end - x_start
        
        if actual_x_count > 0 and x_start < N:
            full_slab[:, :, x_start:x_end] = data[:, :, :actual_x_count]
        x_start = x_end
    
    return full_slab

def compare_overlapping_region(data256, data512, N_overlap=256):
    """
    Compare overlapping region (x,y < N_overlap) between N=256 and N=512 Z-slabs
    
    Note: Each Z-slab is 2D (Y, X), not 3D
    We compare multiple Z-slabs (z = 0 to N_overlap-1) to get 3D comparison
    
    Returns: statistics about differences
    """
    if data256 is None or data512 is None:
        return None
    
    narray = data256.shape[0]
    N256 = data256.shape[1]  # Y dimension, should be 256
    N512 = data512.shape[1]   # Y dimension, should be 512
    
    # Extract overlapping region (Y, X < N_overlap)
    overlap256 = data256[:, :N_overlap, :N_overlap]
    overlap512 = data512[:, :N_overlap, :N_overlap]
    
    results = {}
    
    for array_idx in range(narray):
        arr256 = overlap256[array_idx]
        arr512 = overlap512[array_idx]
        
        # Compare real parts
        diff_real = np.abs(arr256.real - arr512.real)
        diff_imag = np.abs(arr256.imag - arr512.imag)
        
        results[array_idx] = {
            'max_diff_real': np.max(diff_real),
            'mean_diff_real': np.mean(diff_real),
            'max_diff_imag': np.max(diff_imag),
            'mean_diff_imag': np.mean(diff_imag),
            'rms_256_real': np.sqrt(np.mean(arr256.real**2)),
            'rms_512_real': np.sqrt(np.mean(arr512.real**2)),
            'rms_256_imag': np.sqrt(np.mean(arr256.imag**2)),
            'rms_512_imag': np.sqrt(np.mean(arr512.imag**2)),
            'match_count_real': np.sum(diff_real < 1e-10),
            'match_count_imag': np.sum(diff_imag < 1e-10),
            'total_points': arr256.size,
        }
    
    return results

def main():
    parser = argparse.ArgumentParser(description='Compare final output matrices between N=256 and N=512')
    parser.add_argument('--N256-dir', type=str, required=True,
                       help='Directory containing N=256 output files (rank_*/z*_slab_N256.bin)')
    parser.add_argument('--N512-dir', type=str, required=True,
                       help='Directory containing N=512 output files (rank_*/z*_slab_N512.bin)')
    parser.add_argument('--narray', type=int, default=4,
                       help='Number of arrays (default: 4)')
    parser.add_argument('--z-start', type=int, default=0,
                       help='First Z-slab to compare (default: 0)')
    parser.add_argument('--z-end', type=int, default=256,
                       help='Last Z-slab to compare (default: 256, compares z=0 to 255)')
    parser.add_argument('--N-overlap', type=int, default=256,
                       help='Overlapping region size in Y and X (default: 256)')
    
    args = parser.parse_args()
    
    print("=" * 80)
    print("Final Output Matrix Comparison (RNG Consistency Test)")
    print("=" * 80)
    print()
    print(f"Comparing Z-slabs {args.z_start} to {args.z_end-1} in overlapping region (y,x < {args.N_overlap})")
    print(f"N=256 directory: {args.N256_dir}")
    print(f"N=512 directory: {args.N512_dir}")
    print()
    
    # Read and compare multiple Z-slabs
    print("Reading and comparing Z-slabs...")
    all_results = []
    
    for z in range(args.z_start, args.z_end):
        data256 = extract_full_z_slab(z, args.N256_dir, 256, args.narray)
        data512 = extract_full_z_slab(z, args.N512_dir, 512, args.narray)
        
        if data256 is None or data512 is None:
            if z % 10 == 0:  # Print every 10th missing slab to avoid spam
                print(f"  Z={z}: Skipping (files not found)")
            continue
        
        # Compare overlapping region for this Z-slab
        results = compare_overlapping_region(data256, data512, args.N_overlap)
        if results:
            all_results.append((z, results))
            if z % 50 == 0:  # Print progress every 50 slabs
                print(f"  Z={z}: Compared successfully")
    
    if len(all_results) == 0:
        print("ERROR: No Z-slabs could be compared")
        sys.exit(1)
    
    print(f"\nAggregating results across {len(all_results)} Z-slabs...")
    
    # Aggregate results across all Z-slabs
    aggregated = {}
    for array_idx in range(args.narray):
        aggregated[array_idx] = {
            'max_diff_real': 0.0,
            'mean_diff_real': 0.0,
            'max_diff_imag': 0.0,
            'mean_diff_imag': 0.0,
            'match_count_real': 0,
            'match_count_imag': 0,
            'total_points': 0,
        }
    
    for z, results in all_results:
        for array_idx in range(args.narray):
            r = results[array_idx]
            agg = aggregated[array_idx]
            agg['max_diff_real'] = max(agg['max_diff_real'], r['max_diff_real'])
            agg['max_diff_imag'] = max(agg['max_diff_imag'], r['max_diff_imag'])
            agg['mean_diff_real'] += r['mean_diff_real'] * r['total_points']
            agg['mean_diff_imag'] += r['mean_diff_imag'] * r['total_points']
            agg['match_count_real'] += r['match_count_real']
            agg['match_count_imag'] += r['match_count_imag']
            agg['total_points'] += r['total_points']
    
    # Normalize means
    for array_idx in range(args.narray):
        agg = aggregated[array_idx]
        if agg['total_points'] > 0:
            agg['mean_diff_real'] /= agg['total_points']
            agg['mean_diff_imag'] /= agg['total_points']
    
    results = aggregated
    
    print("=" * 80)
    print("Results:")
    print("=" * 80)
    print()
    
    array_names = ['Density+X-disp', 'Y-disp+Z-disp', 'X-velocity', 'Y-velocity+Z-velocity']
    
    for array_idx in range(args.narray):
        r = results[array_idx]
        name = array_names[array_idx] if array_idx < len(array_names) else f"Array {array_idx}"
        
        print(f"Array {array_idx} ({name}):")
        print(f"  Real part:")
        print(f"    Max difference: {r['max_diff_real']:.6e}")
        print(f"    Mean difference: {r['mean_diff_real']:.6e}")
        print(f"    Matches (< 1e-10): {r['match_count_real']}/{r['total_points']} ({100*r['match_count_real']/r['total_points']:.2f}%)")
        print(f"  Imaginary part:")
        print(f"    Max difference: {r['max_diff_imag']:.6e}")
        print(f"    Mean difference: {r['mean_diff_imag']:.6e}")
        print(f"    Matches (< 1e-10): {r['match_count_imag']}/{r['total_points']} ({100*r['match_count_imag']/r['total_points']:.2f}%)")
        print()
    
    # Summary
    total_matches_real = sum(r['match_count_real'] for r in results.values())
    total_matches_imag = sum(r['match_count_imag'] for r in results.values())
    total_points = sum(r['total_points'] for r in results.values())
    
    print("=" * 80)
    print("Summary:")
    print(f"  Real part matches: {total_matches_real}/{total_points} ({100*total_matches_real/total_points:.2f}%)")
    print(f"  Imaginary part matches: {total_matches_imag}/{total_points} ({100*total_matches_imag/total_points:.2f}%)")
    print("=" * 80)
    
    if total_matches_real == total_points and total_matches_imag == total_points:
        print("SUCCESS: Perfect match! RNG consistency verified.")
    elif total_matches_real > 0.99 * total_points and total_matches_imag > 0.99 * total_points:
        print("SUCCESS: Near-perfect match (>99%). RNG consistency verified.")
    else:
        print("WARNING: Some differences found. Check max differences above.")

if __name__ == '__main__':
    main()

