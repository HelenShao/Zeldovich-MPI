#!/usr/bin/env python3
"""
Reassemble Z-slabs for N=16 with 4 MPI ranks

Domain decomposition:
- grid_x=2, grid_z=2 (2x2 grid)
- Rank 0: X=[0,8), Z=[0,8)   -> z=0-7
- Rank 1: X=[0,8), Z=[8,16)  -> z=8-15
- Rank 2: X=[8,16), Z=[0,8)  -> z=0-7
- Rank 3: X=[8,16), Z=[8,16) -> z=8-15

For each z:
- z=0-7:  Need rank_0 (X=0-7) + rank_2 (X=8-15)
- z=8-15: Need rank_1 (X=0-7) + rank_3 (X=8-15)
"""

import numpy as np
import os
import sys

def reassemble_z_slab(z, output_dir, N=16, narray=4):
    """
    Reassemble a single Z-slab from rank files
    
    Args:
        z: Z value to reassemble
        output_dir: Directory containing rank_*/ directories
        N: Grid size (default: 16)
        narray: Number of arrays (default: 4)
    
    Returns:
        Full Z-slab as numpy array: [narray][N][N] (complex128)
    """
    # Determine which ranks own this Z value
    grid_z = 2  # For 4 ranks
    z_block = z // (N // grid_z)  # Which Z block: 0 or 1
    
    if z_block == 0:
        # Z=0-7: Need ranks 0 and 2
        ranks = [0, 2]
        x_ranges = [(0, 8), (8, 16)]  # Rank 0: X=0-7, Rank 2: X=8-15
    else:
        # Z=8-15: Need ranks 1 and 3
        ranks = [1, 3]
        x_ranges = [(0, 8), (8, 16)]  # Rank 1: X=0-7, Rank 3: X=8-15
    
    # Allocate full slab
    full_slab = np.zeros((narray, N, N), dtype=np.complex128)
    
    # Read and assemble from each rank
    for rank, (x_start, x_end) in zip(ranks, x_ranges):
        filename = os.path.join(output_dir, f"rank_{rank}/z{z}_slab_N{N}.bin")
        
        if not os.path.exists(filename):
            print(f"ERROR: File not found: {filename}")
            return None
        
        # Read file
        data = np.fromfile(filename, dtype=np.complex128)
        file_size = len(data) * 16  # 16 bytes per complex128
        
        # Verify file size
        x_count = x_end - x_start
        expected_size = narray * N * x_count * 16
        if file_size != expected_size:
            print(f"WARNING: {filename} size mismatch: expected {expected_size}, got {file_size}")
            # Try to infer x_count from file size
            x_count = file_size // (narray * N * 16)
            x_end = x_start + x_count
        
        # Reshape: [narray][N][x_count]
        data = data.reshape(narray, N, x_count)
        
        # Copy into full slab
        full_slab[:, :, x_start:x_end] = data
        
        print(f"  Rank {rank}: X=[{x_start},{x_end}), read {x_count} X values")
    
    return full_slab


def reassemble_all_z_slabs(output_dir, N=16, narray=4, z_start=0, z_end=None):
    """
    Reassemble all Z-slabs and save them
    
    Args:
        output_dir: Directory containing rank_*/ directories
        N: Grid size
        narray: Number of arrays
        z_start: First Z to reassemble (default: 0)
        z_end: Last Z to reassemble (default: N)
    """
    if z_end is None:
        z_end = N
    
    print("=" * 80)
    print(f"REASSEMBLING Z-SLABS FOR N={N}")
    print("=" * 80)
    print(f"Output directory: {output_dir}")
    print(f"Z range: [{z_start}, {z_end})")
    print(f"Domain decomposition: 2x2 grid (4 ranks)")
    print("=" * 80)
    print()
    
    # Create output directory for reassembled slabs
    reassembled_dir = os.path.join(output_dir, "reassembled_z_slabs")
    os.makedirs(reassembled_dir, exist_ok=True)
    
    # Process each Z-slab
    for z in range(z_start, z_end):
        print(f"Reassembling Z={z}...")
        
        full_slab = reassemble_z_slab(z, output_dir, N, narray)
        
        if full_slab is None:
            print(f"  ERROR: Failed to reassemble Z={z}")
            continue
        
        # Verify shape
        if full_slab.shape != (narray, N, N):
            print(f"  ERROR: Wrong shape: {full_slab.shape}, expected ({narray}, {N}, {N})")
            continue
        
        # Save reassembled slab
        output_file = os.path.join(reassembled_dir, f"z{z}_slab_N{N}_full.bin")
        full_slab.tofile(output_file)
        
        # Print statistics
        print(f"   Reassembled: shape {full_slab.shape}")
        print(f"   Saved to: {output_file}")
        
        # Print some statistics for first slab
        if z == z_start:
            print(f"\n  Statistics for Z={z}:")
            for array_idx in range(narray):
                arr = full_slab[array_idx]
                print(f"    Array {array_idx}:")
                print(f"      Real: min={np.min(arr.real):.6e}, max={np.max(arr.real):.6e}, mean={np.mean(arr.real):.6e}")
                print(f"      Imag: min={np.min(arr.imag):.6e}, max={np.max(arr.imag):.6e}, mean={np.mean(arr.imag):.6e}")
        
        print()
    
    print("=" * 80)
    print(f"COMPLETE: Reassembled {z_end - z_start} Z-slabs")
    print(f"Output directory: {reassembled_dir}")
    print("=" * 80)


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Reassemble Z-slabs for N=16 with 4 MPI ranks')
    parser.add_argument('output_dir', type=str, default='.', nargs='?',
                       help='Directory containing rank_*/ directories (default: current directory)')
    parser.add_argument('--N', type=int, default=16,
                       help='Grid size (default: 16)')
    parser.add_argument('--narray', type=int, default=4,
                       help='Number of arrays (default: 4)')
    parser.add_argument('--z-start', type=int, default=0,
                       help='First Z to reassemble (default: 0)')
    parser.add_argument('--z-end', type=int, default=None,
                       help='Last Z to reassemble (default: N)')
    parser.add_argument('--z-only', type=int, default=None,
                       help='Reassemble only this Z value (for testing)')
    
    args = parser.parse_args()
    
    if args.z_only is not None:
        # Reassemble only one Z value (for manual testing)
        print(f"Reassembling Z={args.z_only} only...")
        full_slab = reassemble_z_slab(args.z_only, args.output_dir, args.N, args.narray)
        
        if full_slab is not None:
            print(f"\n Successfully reassembled Z={args.z_only}")
            print(f"  Shape: {full_slab.shape}")
            print(f"  Dtype: {full_slab.dtype}")
            
            # Save it
            output_file = f"z{args.z_only}_slab_N{args.N}_reassembled.bin"
            full_slab.tofile(output_file)
            print(f"  Saved to: {output_file}")
            
            # Print statistics
            print(f"\n  Statistics:")
            for array_idx in range(args.narray):
                arr = full_slab[array_idx]
                print(f"    Array {array_idx}:")
                print(f"      Real: min={np.min(arr.real):.6e}, max={np.max(arr.real):.6e}")
                print(f"      Imag: min={np.min(arr.imag):.6e}, max={np.max(arr.imag):.6e}")
        else:
            print(f"✗ Failed to reassemble Z={args.z_only}")
            sys.exit(1)
    else:
        # Reassemble all Z values
        reassemble_all_z_slabs(args.output_dir, args.N, args.narray, 
                              args.z_start, args.z_end)


if __name__ == '__main__':
    main()
