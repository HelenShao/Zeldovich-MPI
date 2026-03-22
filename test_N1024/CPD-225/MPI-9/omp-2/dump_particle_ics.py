#!/usr/bin/env python3
"""
Brute force dump of all particle IC values for manual comparison.
Prints all particles in a readable format, one per line.
"""

import struct
import sys
import os

def dump_particle_ic_file(filename, N, slab_idx=None, output_file=None):
    """
    Dump all particle values from IC file.
    
    Args:
        filename: Path to particle IC file
        N: Grid size (each slab has N*N particles)
        slab_idx: Which slab to read (0-based). If None, read all slabs.
        output_file: Optional file to write output to. If None, prints to stdout.
    """
    try:
        filesize = os.path.getsize(filename)
        bytes_per_particle = 32
        particles_per_slab = N * N
        bytes_per_slab = particles_per_slab * bytes_per_particle
        
        if filesize % bytes_per_slab != 0:
            print(f"ERROR: File size {filesize} is not a multiple of {bytes_per_slab} bytes", file=sys.stderr)
            return False
        
        num_slabs = filesize // bytes_per_slab
        num_particles_in_file = filesize // bytes_per_particle
        
        if slab_idx is not None:
            if slab_idx >= num_slabs:
                print(f"ERROR: Slab index {slab_idx} >= number of slabs {num_slabs}", file=sys.stderr)
                return False
            num_particles = particles_per_slab
            start_offset = slab_idx * bytes_per_slab
            end_offset = start_offset + bytes_per_slab
            slab_info = f"slab_{slab_idx}"
        else:
            num_particles = num_particles_in_file
            start_offset = 0
            end_offset = filesize
            slab_info = f"all_{num_slabs}_slabs"
        
        with open(filename, 'rb') as f:
            f.seek(start_offset)
            data = f.read(end_offset - start_offset)
        
        pattern = '=HHH xx fff fff'
        
        if output_file:
            fout = open(output_file, 'w')
        else:
            fout = sys.stdout
        
        # Print header
        print(f"# Particle IC dump: {filename}", file=fout)
        print(f"# Slab info: {slab_info}", file=fout)
        print(f"# N={N}, particles_per_slab={particles_per_slab}, num_particles={num_particles}", file=fout)
        print(f"# Format: particle_idx i j k displ[0] displ[1] displ[2] vel[0] vel[1] vel[2]", file=fout)
        print(f"# displ[0]=Z, displ[1]=Y, displ[2]=X", file=fout)
        print(f"# vel[0]=Z, vel[1]=Y, vel[2]=X", file=fout)
        print("#", file=fout)
        
        # Print all particles
        for idx in range(num_particles):
            offset = idx * bytes_per_particle
            i, j, k, dx, dy, dz, vx, vy, vz = struct.unpack(pattern, data[offset:offset+bytes_per_particle])
            print(f"{idx:4d} {i:4d} {j:4d} {k:4d} "
                  f"{dx:15.10f} {dy:15.10f} {dz:15.10f} "
                  f"{vx:15.10f} {vy:15.10f} {vz:15.10f}", file=fout)
        
        if output_file:
            fout.close()
            print(f"Wrote {num_particles} particles to {output_file}", file=sys.stderr)
        else:
            fout.flush()
        
        return True
        
    except FileNotFoundError:
        print(f"ERROR: File not found: {filename}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return False

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 dump_particle_ics.py <ic_file> <N> [slab_idx] [output_file]")
        print("Example: python3 dump_particle_ics.py particle_ics/ic_0 4")
        print("Example: python3 dump_particle_ics.py particle_ics/ic_0 4 0 hermitian_ic_0_slab0.txt")
        print("  (slab_idx is optional: 0 for first slab, 1 for second slab, etc. If omitted, dumps all slabs)")
        print("  (output_file is optional: if omitted, prints to stdout)")
        sys.exit(1)
    
    filename = sys.argv[1]
    N = int(sys.argv[2])
    slab_idx = int(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3].isdigit() else None
    output_file = sys.argv[4] if len(sys.argv) > 4 else None
    
    # Handle case where output_file is provided but slab_idx is not
    if len(sys.argv) == 4 and not sys.argv[3].isdigit():
        output_file = sys.argv[3]
        slab_idx = None
    
    success = dump_particle_ic_file(filename, N, slab_idx, output_file)
    sys.exit(0 if success else 1)

if __name__ == '__main__':
    main()

