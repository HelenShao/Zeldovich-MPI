#!/usr/bin/env python3
"""
Reassemble Z-slabs from MPI-parallelized hermitian_3d_matrix_production output

This script reads partial Z-slabs written by each MPI rank and reassembles them
into complete Z-slabs containing all X and Y values for each Z coordinate.

Domain Decomposition:
- MPI ranks are organized in a grid_x × grid_z decomposition
- Each rank owns a rectangular region: X=[x_start, x_end), Z=[z_start, z_end)
- For a given Z value, multiple ranks contribute their X ranges

Input Format:
- Files: rank_{rank}/z{z}_slab_N{N}.bin
- Format: [Array][Y][X] (complex128)
- Each file contains: narray arrays × N Y-values × x_count X-values

Output Format:
- Files: reassembled_z_slabs/z{z}_slab_N{N}_full.bin
- Format: [Array][Y][X] (complex128)
- Each file contains: narray arrays × N Y-values × N X-values (complete)
"""

import numpy as np
import os
import sys
import math


def calculate_grid_decomposition(num_ranks):
    """
    Calculate grid_x and grid_z from total number of ranks.
    Tries to make grid_x and grid_z as close as possible (square-like).
    
    Args:
        num_ranks: Total number of MPI ranks
    
    Returns:
        (grid_x, grid_z): Grid dimensions
    """
    # Find factors closest to sqrt(num_ranks)
    grid_x = int(math.sqrt(num_ranks))
    while num_ranks % grid_x != 0:
        grid_x -= 1
    grid_z = num_ranks // grid_x
    
    # Prefer wider grid (grid_x >= grid_z) for better cache locality
    if grid_x < grid_z:
        grid_x, grid_z = grid_z, grid_x
    
    return grid_x, grid_z


def get_rank_bounds(rank, grid_x, grid_z, N):
    """
    Calculate the X and Z bounds for a given rank.
    
    Args:
        rank: MPI rank number
        grid_x: Number of ranks in X dimension
        grid_z: Number of ranks in Z dimension
        N: Grid size
    
    Returns:
        (x_start, x_end, z_start, z_end): Bounds for this rank
    """
    # Calculate rank's position in grid
    rank_x = rank % grid_x
    rank_z = rank // grid_x
    
    # Calculate bounds
    x_per_rank = N // grid_x
    z_per_rank = N // grid_z
    
    x_start = rank_x * x_per_rank
    x_end = (rank_x + 1) * x_per_rank if rank_x < grid_x - 1 else N
    
    z_start = rank_z * z_per_rank
    z_end = (rank_z + 1) * z_per_rank if rank_z < grid_z - 1 else N
    
    return x_start, x_end, z_start, z_end


def find_ranks_for_z(z, grid_x, grid_z, N):
    """
    Find which ranks contribute to a given Z value and their X ranges.
    
    Args:
        z: Z coordinate
        grid_x: Number of ranks in X dimension
        grid_z: Number of ranks in Z dimension
        N: Grid size
    
    Returns:
        List of (rank, x_start, x_end) tuples
    """
    z_per_rank = N // grid_z
    z_block = z // z_per_rank
    
    ranks_and_ranges = []
    for rank_x in range(grid_x):
        rank = z_block * grid_x + rank_x
        
        x_per_rank = N // grid_x
        x_start = rank_x * x_per_rank
        x_end = (rank_x + 1) * x_per_rank if rank_x < grid_x - 1 else N
        
        ranks_and_ranges.append((rank, x_start, x_end))
    
    return ranks_and_ranges


def reassemble_z_slab(z, output_dir, N, narray, grid_x, grid_z):
    """
    Reassemble a single Z-slab from rank files.
    
    Args:
        z: Z value to reassemble
        output_dir: Directory containing rank_*/ directories
        N: Grid size
        narray: Number of arrays
        grid_x: Number of ranks in X dimension
        grid_z: Number of ranks in Z dimension
    
    Returns:
        Full Z-slab as numpy array: [narray][N][N] (complex128)
    """
    # Find which ranks contribute to this Z value
    ranks_and_ranges = find_ranks_for_z(z, grid_x, grid_z, N)
    
    # Allocate full slab: [Array][Y][X]
    full_slab = np.zeros((narray, N, N), dtype=np.complex128)
    
    # Read and assemble from each rank
    for rank, x_start, x_end in ranks_and_ranges:
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
        
        # Reshape: [narray][N][x_count] (matches input format [Array][Y][X])
        data = data.reshape(narray, N, x_count)
        
        # Copy into full slab: [narray][N][x_start:x_end]
        full_slab[:, :, x_start:x_end] = data
        
        print(f"  Rank {rank}: X=[{x_start},{x_end}), read {x_count} X values")
    
    return full_slab


def reassemble_all_z_slabs(output_dir, N, narray, grid_x, grid_z, z_start=0, z_end=None):
    """
    Reassemble all Z-slabs and save them.
    
    Args:
        output_dir: Directory containing rank_*/ directories
        N: Grid size
        narray: Number of arrays
        grid_x: Number of ranks in X dimension
        grid_z: Number of ranks in Z dimension
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
    print(f"Domain decomposition: {grid_x}×{grid_z} grid ({grid_x * grid_z} ranks)")
    print(f"Data format: [Array][Y][X] (complex128)")
    print("=" * 80)
    print()
    
    # Create output directory for reassembled slabs
    reassembled_dir = os.path.join(output_dir, "reassembled_z_slabs")
    os.makedirs(reassembled_dir, exist_ok=True)
    
    # Process each Z-slab
    for z in range(z_start, z_end):
        print(f"Reassembling Z={z}...")
        
        full_slab = reassemble_z_slab(z, output_dir, N, narray, grid_x, grid_z)
        
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
    
    parser = argparse.ArgumentParser(
        description='Reassemble Z-slabs from MPI-parallelized hermitian_3d_matrix_production output',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-detect grid decomposition from number of rank directories
  python reassemble_z_slabs.py output_dir --N 1024 --narray 4
  
  # Specify grid decomposition explicitly
  python reassemble_z_slabs.py output_dir --N 1024 --narray 4 --grid-x 4 --grid-z 4
  
  # Reassemble only one Z value for testing
  python reassemble_z_slabs.py output_dir --N 1024 --narray 4 --z-only 0
        """
    )
    parser.add_argument('output_dir', type=str,
                       help='Directory containing rank_*/ directories')
    parser.add_argument('--N', type=int, required=True,
                       help='Grid size (required)')
    parser.add_argument('--narray', type=int, default=4,
                       help='Number of arrays (default: 4)')
    parser.add_argument('--grid-x', type=int, default=None,
                       help='Number of ranks in X dimension (auto-detect if not specified)')
    parser.add_argument('--grid-z', type=int, default=None,
                       help='Number of ranks in Z dimension (auto-detect if not specified)')
    parser.add_argument('--z-start', type=int, default=0,
                       help='First Z to reassemble (default: 0)')
    parser.add_argument('--z-end', type=int, default=None,
                       help='Last Z to reassemble (default: N)')
    parser.add_argument('--z-only', type=int, default=None,
                       help='Reassemble only this Z value (for testing)')
    
    args = parser.parse_args()
    
    # Auto-detect grid decomposition if not specified
    if args.grid_x is None or args.grid_z is None:
        # Count number of rank directories
        rank_dirs = [d for d in os.listdir(args.output_dir) 
                     if os.path.isdir(os.path.join(args.output_dir, d)) and d.startswith('rank_')]
        num_ranks = len(rank_dirs)
        
        if num_ranks == 0:
            print(f"ERROR: No rank_* directories found in {args.output_dir}")
            sys.exit(1)
        
        grid_x, grid_z = calculate_grid_decomposition(num_ranks)
        print(f"Auto-detected grid decomposition: {grid_x}×{grid_z} ({num_ranks} ranks)")
    else:
        grid_x = args.grid_x
        grid_z = args.grid_z
        num_ranks = grid_x * grid_z
        print(f"Using specified grid decomposition: {grid_x}×{grid_z} ({num_ranks} ranks)")
    
    if args.z_only is not None:
        # Reassemble only one Z value (for manual testing)
        print(f"Reassembling Z={args.z_only} only...")
        full_slab = reassemble_z_slab(args.z_only, args.output_dir, args.N, args.narray, grid_x, grid_z)
        
        if full_slab is not None:
            print(f"\n✓ Successfully reassembled Z={args.z_only}")
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
        reassemble_all_z_slabs(args.output_dir, args.N, args.narray, grid_x, grid_z,
                              args.z_start, args.z_end)


if __name__ == '__main__':
    main()
