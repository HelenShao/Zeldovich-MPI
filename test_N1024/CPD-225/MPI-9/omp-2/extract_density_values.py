#!/usr/bin/env python3
"""
Extract density field values from binary files and save to text files.
Each density file contains N^3 float32 values (for N=4: 64 values).
"""

import struct
import numpy as np
import sys
import os

def read_density_file(filename):
    """Read density file as float32 array."""
    try:
        with open(filename, 'rb') as f:
            data = f.read()
        density = np.frombuffer(data, dtype=np.float32)
        return density, None
    except FileNotFoundError:
        return None, f"File not found: {filename}"
    except Exception as e:
        return None, f"Error reading file: {e}"

def write_density_text(density, output_file, N):
    """Write density values to text file with index information."""
    with open(output_file, 'w') as f:
        f.write(f"# Density field values (N={N}, N³={N**3} values)\n")
        f.write("# Format: index i j k density_value\n")
        f.write("# Coordinates (i,j,k) correspond to grid cell (x,y,z)\n")
        f.write("#\n")
        
        idx = 0
        # C++ write order: WriteParticlesSlab_unified loops: for (j=Y) { for (k_local=X) }
        # For each Z-slab (i), writes: (Y=0, X=0), (Y=0, X=1), ..., (Y=0, X=N-1),
        #                              (Y=1, X=0), ..., (Y=N-1, X=N-1)
        # Z-slabs written sequentially, so overall file order is: (Y, X, Z)
        # This means: Z varies slowest (between Z-slabs), Y varies middle (within each Z-slab), X varies fastest (within each Y)
        # Python must match this order: Z (outer), Y (middle), X (inner)
        # Note: Loop order matches C++ write order: Z outer, Y middle, X inner
        for k in range(N):      # Z coordinate (outer, slowest)
            for j in range(N):  # Y coordinate (middle)
                for i in range(N):  # X coordinate (inner, fastest)
                    if idx < len(density):
                        f.write(f"{idx:4d}  {i:2d}  {j:2d}  {k:2d}  {density[idx]:15.10e}\n")
                        idx += 1

def main():
    if len(sys.argv) < 4:
        print("Usage: python3 extract_density_values.py <input_binary_file> <output_text_file> <N>")
        print("Example: python3 extract_density_values.py density4 density_values.txt 4")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    N = int(sys.argv[3])
    
    # Read density file
    density, error = read_density_file(input_file)
    if density is None:
        print(f"ERROR: {error}")
        sys.exit(1)
    
    expected_size = N**3
    if len(density) != expected_size:
        print(f"WARNING: Expected {expected_size} values, got {len(density)}")
    
    # Write to text file
    write_density_text(density, output_file, N)
    
    # Print statistics
    print(f"Extracted {len(density)} density values from {input_file}")
    print(f"Saved to {output_file}")
    print()
    print("Statistics:")
    print(f"  Min: {density.min():.10e}")
    print(f"  Max: {density.max():.10e}")
    print(f"  Mean: {density.mean():.10e}")
    print(f"  Std: {density.std():.10e}")
    print(f"  RMS: {np.sqrt(np.mean(density**2)):.10e}")

if __name__ == '__main__':
    main()

