#!/usr/bin/env python3
"""
Read and display particle data from binary IC files.
Supports RVZelParticle format (36 bytes per particle).
"""

import struct
import sys

def read_particles(filename, num_particles=3):
    """Read first num_particles from a binary particle file."""
    try:
        with open(filename, 'rb') as f:
            data = f.read()
    except FileNotFoundError:
        return None, f"File not found: {filename}"
    except Exception as e:
        return None, f"Error reading file: {e}"
    
    # RVZelParticle format: int i,j,k (12 bytes) + float displ[3] (12 bytes) + float vel[3] (12 bytes) = 36 bytes
    particle_size = 36
    file_size = len(data)
    max_particles = file_size // particle_size
    
    if max_particles == 0:
        return None, "File too small or empty"
    
    show_count = min(num_particles, max_particles)
    particles = []
    
    for i in range(show_count):
        offset = i * particle_size
        if offset + particle_size > file_size:
            break
        
        # Unpack: 3 ints (i,j,k) + 6 floats (displ[3], vel[3])
        try:
            i_val, j_val, k_val, displ0, displ1, displ2, vel0, vel1, vel2 = struct.unpack(
                'iii fff fff', data[offset:offset+particle_size]
            )
            particles.append({
                'i': i_val, 'j': j_val, 'k': k_val,
                'displ': [displ0, displ1, displ2],
                'vel': [vel0, vel1, vel2]
            })
        except struct.error as e:
            return None, f"Error unpacking particle {i}: {e}"
    
    return particles, max_particles

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: read_particles.py <filename> [num_particles]")
        sys.exit(1)
    
    filename = sys.argv[1]
    num_particles = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    
    particles, max_particles = read_particles(filename, num_particles)
    
    if particles is None:
        print(f"Error: {max_particles}")
        sys.exit(1)
    
    print(f"File: {filename}")
    print(f"Total particles: {max_particles}")
    print(f"Showing first {len(particles)} particles:")
    print()
    print("  Particle |   i   j   k  |     displ[0]     displ[1]     displ[2]     |     vel[0]      vel[1]      vel[2]")
    print("  ---------|---------------|--------------------------------------------------|-----------------------------------")
    
    for idx, p in enumerate(particles):
        print(f"    {idx:8d} | {p['i']:4d} {p['j']:4d} {p['k']:4d} | "
              f"{p['displ'][0]:13.6f} {p['displ'][1]:13.6f} {p['displ'][2]:13.6f} | "
              f"{p['vel'][0]:13.6f} {p['vel'][1]:13.6f} {p['vel'][2]:13.6f}")

