#!/usr/bin/env python3
"""
Read zeldovich-PLT output files (ic_* and density files)

Usage:
    python3 read_particle_files.py output/ic_0
    python3 read_particle_files.py output/density16 --density
"""

import numpy as np
import struct
import sys
import os

def read_ic_file(filename, format_type="RVZel"):
    """
    Read particle file (ic_*)
    
    Args:
        filename: Path to ic_* file
        format_type: "RVZel" (float, 32 bytes) or "RVdoubleZel" (double, 54 bytes)
    
    Returns:
        Dictionary with particle data
    """
    if format_type == "RVZel":
        # RVZel format: uint16 i,j,k (6 bytes) + 2 padding + float displ[3] (12 bytes) + float vel[3] (12 bytes) = 32 bytes
        dtype = np.dtype([
            ('i', np.uint16),
            ('j', np.uint16),
            ('k', np.uint16),
            ('padding', np.uint16),  # 2 bytes padding for alignment
            ('displ', np.float32, 3),  # Single precision (4 bytes each)
            ('vel', np.float32, 3)     # Single precision (4 bytes each)
        ])
        particle_size = 32
    elif format_type == "RVdoubleZel":
        # RVdoubleZel format: uint16 i,j,k (6 bytes) + 2 padding + double displ[3] (24 bytes) + double vel[3] (24 bytes) = 54 bytes
        dtype = np.dtype([
            ('i', np.uint16),
            ('j', np.uint16),
            ('k', np.uint16),
            ('padding', np.uint16),
            ('displ', np.float64, 3),  # Double precision (8 bytes each)
            ('vel', np.float64, 3)      # Double precision (8 bytes each)
        ])
        particle_size = 54
    else:
        raise ValueError(f"Unknown format: {format_type}")
    
    file_size = os.path.getsize(filename)
    num_particles = file_size // particle_size
    
    print(f"Reading {filename}")
    print(f"  File size: {file_size} bytes")
    print(f"  Particle size: {particle_size} bytes ({format_type})")
    print(f"  Number of particles: {num_particles}")
    
    # Read file
    data = np.fromfile(filename, dtype=dtype)
    
    if len(data) != num_particles:
        print(f"  WARNING: Expected {num_particles} particles, got {len(data)}")
    
    return {
        'i': data['i'],
        'j': data['j'],
        'k': data['k'],
        'displ': data['displ'],  # Shape: (N, 3)
        'vel': data['vel'],       # Shape: (N, 3)
        'num_particles': num_particles
    }


def read_density_file(filename, ppd=16):
    """
    Read density file
    
    Args:
        filename: Path to density file
        ppd: Particles per dimension (to infer structure)
    
    Returns:
        Dictionary with density data
    """
    file_size = os.path.getsize(filename)
    
    # Density file contains float values (4 bytes each)
    # Total should be ppd^3 floats
    expected_size = ppd * ppd * ppd * 4
    num_floats = file_size // 4
    
    print(f"Reading {filename}")
    print(f"  File size: {file_size} bytes")
    print(f"  Float size: 4 bytes (single precision)")
    print(f"  Number of floats: {num_floats}")
    print(f"  Expected for N={ppd}: {expected_size} bytes ({ppd}^3 = {ppd*ppd*ppd} floats)")
    
    # Read as float32 (single precision)
    density = np.fromfile(filename, dtype=np.float32)
    
    # Reshape to 3D grid if possible
    if len(density) == ppd * ppd * ppd:
        density_3d = density.reshape(ppd, ppd, ppd)
        print(f"  Reshaped to 3D: ({ppd}, {ppd}, {ppd})")
    else:
        density_3d = None
        print(f"  WARNING: Cannot reshape to ({ppd}, {ppd}, {ppd})")
    
    return {
        'density': density,
        'density_3d': density_3d,
        'num_values': len(density)
    }


def print_particle_stats(particles):
    """Print statistics about particles"""
    print("\n" + "="*60)
    print("PARTICLE STATISTICS")
    print("="*60)
    
    print(f"\nNumber of particles: {particles['num_particles']}")
    
    print("\nGrid indices:")
    print(f"  i: min={np.min(particles['i'])}, max={np.max(particles['i'])}")
    print(f"  j: min={np.min(particles['j'])}, max={np.max(particles['j'])}")
    print(f"  k: min={np.min(particles['k'])}, max={np.max(particles['k'])}")
    
    print("\nDisplacements (single precision float = 4 bytes each):")
    for i, axis in enumerate(['x', 'y', 'z']):
        displ = particles['displ'][:, i]
        print(f"  {axis}: min={np.min(displ):.6e}, max={np.max(displ):.6e}, "
              f"mean={np.mean(displ):.6e}, std={np.std(displ):.6e}")
    
    print("\nVelocities (single precision float = 4 bytes each):")
    for i, axis in enumerate(['x', 'y', 'z']):
        vel = particles['vel'][:, i]
        print(f"  {axis}: min={np.min(vel):.6e}, max={np.max(vel):.6e}, "
              f"mean={np.mean(vel):.6e}, std={np.std(vel):.6e}")
    
    # Show first few particles
    print("\nFirst 5 particles:")
    print("  i   j   k   displ_x      displ_y      displ_z      vel_x        vel_y        vel_z")
    print("  " + "-"*80)
    for i in range(min(5, particles['num_particles'])):
        print(f"  {particles['i'][i]:3d} {particles['j'][i]:3d} {particles['k'][i]:3d} "
              f"{particles['displ'][i,0]:11.6e} {particles['displ'][i,1]:11.6e} {particles['displ'][i,2]:11.6e} "
              f"{particles['vel'][i,0]:11.6e} {particles['vel'][i,1]:11.6e} {particles['vel'][i,2]:11.6e}")


def print_density_stats(density_data):
    """Print statistics about density"""
    print("\n" + "="*60)
    print("DENSITY STATISTICS")
    print("="*60)
    
    print(f"\nNumber of density values: {density_data['num_values']}")
    print(f"  min={np.min(density_data['density']):.6e}")
    print(f"  max={np.max(density_data['density']):.6e}")
    print(f"  mean={np.mean(density_data['density']):.6e}")
    print(f"  std={np.std(density_data['density']):.6e}")
    
    if density_data['density_3d'] is not None:
        print(f"\n3D grid shape: {density_data['density_3d'].shape}")
        print(f"  Z-slab 0 (first 16x16): min={np.min(density_data['density_3d'][0]):.6e}, "
              f"max={np.max(density_data['density_3d'][0]):.6e}")


def main():
    if len(sys.argv) < 2:
        print("Usage:")
        print("  python3 read_particle_files.py <ic_file> [--format RVZel|RVdoubleZel]")
        print("  python3 read_particle_files.py <density_file> --density [--ppd N]")
        print("\nExamples:")
        print("  python3 read_particle_files.py output/ic_0")
        print("  python3 read_particle_files.py output/ic_0 --format RVdoubleZel")
        print("  python3 read_particle_files.py output/density16 --density --ppd 16")
        sys.exit(1)
    
    filename = sys.argv[1]
    
    if '--density' in sys.argv:
        # Read density file
        ppd = 16
        if '--ppd' in sys.argv:
            idx = sys.argv.index('--ppd')
            ppd = int(sys.argv[idx + 1])
        
        density_data = read_density_file(filename, ppd)
        print_density_stats(density_data)
    else:
        # Read particle file
        format_type = "RVZel"
        if '--format' in sys.argv:
            idx = sys.argv.index('--format')
            format_type = sys.argv[idx + 1]
        
        particles = read_ic_file(filename, format_type)
        print_particle_stats(particles)


if __name__ == '__main__':
    main()
