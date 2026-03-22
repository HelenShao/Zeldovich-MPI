#!/usr/bin/env python3
"""
Compare particle IC files between hermitian_3d_matrix_production and zeldovich-PLT.
Extracts and compares displacement and velocity fields.
"""

import struct
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MaxNLocator

# Load plotting parameters
sys.path.insert(0, os.path.expanduser('~/python_modules'))
try:
    from matplotlib_plotting_params import apply_plotting_params
    apply_plotting_params()
except ImportError:
    print("Warning: Could not import matplotlib_plotting_params, using default settings")

# RVZelParticle structure (actual format in zeldovich-PLT):
# Grid indices (8 bytes):
#   - i (uint16): Z-coordinate index (2 bytes)
#   - j (uint16): Y-coordinate index (2 bytes)
#   - k (uint16): X-coordinate index (2 bytes)
#   - Padding (2 bytes) for 4-byte alignment
# Displacement vector (12 bytes, 3 floats):
#   - displ[0]: Z-displacement (4 bytes)
#   - displ[1]: Y-displacement (4 bytes)
#   - displ[2]: X-displacement (4 bytes)
# Velocity vector (12 bytes, 3 floats):
#   - vel[0]: Z-velocity (4 bytes)
#   - vel[1]: Y-velocity (4 bytes)
#   - vel[2]: X-velocity (4 bytes)
# Total: 8 + 12 + 12 = 32 bytes per particle

def read_particle_ic_file(filename, N, slab_idx=None):
    """
    Read particle IC file in RVZel format.
    File may contain multiple Z-slabs (each with N*N particles).
    
    Args:
        filename: Path to particle IC file
        N: Grid size (each slab has N*N particles)
        slab_idx: Which slab to read (0-based). If None, read all slabs.
    
    Returns: dict with 'i', 'j', 'k', 'displ', 'vel' arrays
    """
    try:
        filesize = os.path.getsize(filename)
        bytes_per_particle = 32
        particles_per_slab = N * N
        bytes_per_slab = particles_per_slab * bytes_per_particle
        
        # File might contain multiple slabs, so check if size is a multiple of bytes_per_slab
        if filesize % bytes_per_slab != 0:
            return None, f"File size {filesize} is not a multiple of {bytes_per_slab} bytes (expected multiple of {particles_per_slab} particles per slab)"
        
        num_slabs = filesize // bytes_per_slab
        num_particles_in_file = filesize // bytes_per_particle
        
        if slab_idx is not None:
            if slab_idx >= num_slabs:
                return None, f"Slab index {slab_idx} >= number of slabs {num_slabs} in file"
            num_particles = particles_per_slab
            start_offset = slab_idx * bytes_per_slab
            end_offset = start_offset + bytes_per_slab
            print(f"Reading slab {slab_idx} from file (particles {slab_idx * particles_per_slab} to {(slab_idx + 1) * particles_per_slab - 1})")
        else:
            # Read all slabs
            num_particles = num_particles_in_file
            start_offset = 0
            end_offset = filesize
            if num_slabs > 1:
                print(f"File contains {num_slabs} Z-slabs, reading all {num_particles} particles")
        
        with open(filename, 'rb') as f:
            f.seek(start_offset)
            data = f.read(end_offset - start_offset)
        
        particles = {
            'i': np.zeros(num_particles, dtype=np.uint16),
            'j': np.zeros(num_particles, dtype=np.uint16),
            'k': np.zeros(num_particles, dtype=np.uint16),
            'displ': np.zeros((num_particles, 3), dtype=np.float32),
            'vel': np.zeros((num_particles, 3), dtype=np.float32)
        }
        
        # Structure layout (32 bytes per particle):
        #   HHH (6 bytes: i, j, k) + xx (2 bytes padding) = 8 bytes
        #   fff (12 bytes: displ[0], displ[1], displ[2]) = 12 bytes
        #   fff (12 bytes: vel[0], vel[1], vel[2]) = 12 bytes
        # Total: 8 + 12 + 12 = 32 bytes per particle
        pattern = '=HHH xx fff fff'
        
        for idx in range(num_particles):
            offset = idx * bytes_per_particle
            # Unpack: 3 uint16_t (i, j, k) + 2 bytes padding, 3 float (displ), 3 float (vel)
            i, j, k, dx, dy, dz, vx, vy, vz = struct.unpack(pattern, data[offset:offset+bytes_per_particle])
            particles['i'][idx] = i
            particles['j'][idx] = j
            particles['k'][idx] = k
            particles['displ'][idx, 0] = dx  # displ[0] = Z displacement
            particles['displ'][idx, 1] = dy  # displ[1] = Y displacement
            particles['displ'][idx, 2] = dz  # displ[2] = X displacement
            particles['vel'][idx, 0] = vx    # vel[0] = Z velocity
            particles['vel'][idx, 1] = vy    # vel[1] = Y velocity
            particles['vel'][idx, 2] = vz    # vel[2] = X velocity
        
        return particles, None
    except FileNotFoundError:
        return None, f"File not found: {filename}"
    except Exception as e:
        return None, f"Error reading file: {e}"

def compare_particle_ics(particles1, particles2, label1, label2, output_dir, N, output_suffix=None):
    """Compare two particle IC datasets and create visualizations.
    If output_suffix is set (e.g. 'ic_0'), save to displacement_comparison_ic_0.png etc.
    """
    os.makedirs(output_dir, exist_ok=True)
    disp_basename = 'displacement_comparison' + ('_' + output_suffix if output_suffix else '') + '.png'
    vel_basename = 'velocity_comparison' + ('_' + output_suffix if output_suffix else '') + '.png'
    
    # Sort by (i, j, k) to ensure consistent ordering
    def sort_key(idx):
        return (particles1['i'][idx], particles1['j'][idx], particles1['k'][idx])
    
    indices1 = np.arange(len(particles1['i']))
    indices1_sorted = sorted(indices1, key=sort_key)
    
    def sort_key2(idx):
        return (particles2['i'][idx], particles2['j'][idx], particles2['k'][idx])
    
    indices2 = np.arange(len(particles2['i']))
    indices2_sorted = sorted(indices2, key=sort_key2)
    
    # Extract sorted arrays
    i1 = particles1['i'][indices1_sorted]
    j1 = particles1['j'][indices1_sorted]
    k1 = particles1['k'][indices1_sorted]
    displ1 = particles1['displ'][indices1_sorted, :]
    vel1 = particles1['vel'][indices1_sorted, :]
    
    i2 = particles2['i'][indices2_sorted]
    j2 = particles2['j'][indices2_sorted]
    k2 = particles2['k'][indices2_sorted]
    displ2 = particles2['displ'][indices2_sorted, :]
    vel2 = particles2['vel'][indices2_sorted, :]
    
    # Verify coordinates match
    if not np.array_equal(i1, i2) or not np.array_equal(j1, j2) or not np.array_equal(k1, k2):
        print("WARNING: Coordinate arrays don't match exactly")
        print(f"  i match: {np.array_equal(i1, i2)}")
        print(f"  j match: {np.array_equal(j1, j2)}")
        print(f"  k match: {np.array_equal(k1, k2)}")
        # Find matching particles by coordinates
        coord_dict1 = {(i1[idx], j1[idx], k1[idx]): idx for idx in range(len(i1))}
        coord_dict2 = {(i2[idx], j2[idx], k2[idx]): idx for idx in range(len(i2))}
        common_coords = set(coord_dict1.keys()) & set(coord_dict2.keys())
        print(f"  Common coordinates: {len(common_coords)} / {len(coord_dict1)}")
        
        # Reorder to match
        common_coords_sorted = sorted(common_coords)
        matched_indices1 = [coord_dict1[coord] for coord in common_coords_sorted]
        matched_indices2 = [coord_dict2[coord] for coord in common_coords_sorted]
        
        displ1 = displ1[matched_indices1, :]
        displ2 = displ2[matched_indices2, :]
        vel1 = vel1[matched_indices1, :]
        vel2 = vel2[matched_indices2, :]
        num_compared = len(common_coords_sorted)
    else:
        num_compared = len(i1)
    
    print(f"\nComparing {num_compared} particles")
    
    # Collect status lines for writing to file when output_suffix is set
    status_lines = []
    
    # Compare displacements
    print("\nDisplacement Comparison:")
    print("-" * 80)
    for axis_idx, axis_name in enumerate(['Z', 'Y', 'X']):
        d1 = displ1[:, axis_idx]
        d2 = displ2[:, axis_idx]
        diff = d1 - d2
        max_diff = np.max(np.abs(diff))
        mean_diff = np.mean(diff)
        std_diff = np.std(diff)
        corr = np.corrcoef(d1, d2)[0, 1]
        
        # Calculate mean and std for each dataset
        mean_d1 = np.mean(d1)
        std_d1 = np.std(d1)
        mean_d2 = np.mean(d2)
        std_d2 = np.std(d2)
        
        print(f"  {axis_name}-displacement:")
        print(f"    {label1}: mean={mean_d1:.2e}, std={std_d1:.2e}")
        print(f"    {label2}: mean={mean_d2:.2e}, std={std_d2:.2e}")
        print(f"    Mean difference: {mean_diff:.2e}")
        print(f"    Std difference: {std_diff:.2e}")
        print(f"    Max |difference|: {max_diff:.2e}")
        print(f"    Correlation: {corr:.2e}")
        if max_diff < 1e-5:
            disp_status = "Status: MATCH (max diff < 1e-5)"
            print(f"    {disp_status}")
        else:
            disp_status = "Status: MISMATCH (max diff >= 1e-5)"
            print(f"    {disp_status}")
        status_lines.append(f"  {axis_name}-displacement: {disp_status}")
    
    # Compare velocities
    print("\nVelocity Comparison:")
    print("-" * 80)
    for axis_idx, axis_name in enumerate(['Z', 'Y', 'X']):
        v1 = vel1[:, axis_idx]
        v2 = vel2[:, axis_idx]
        diff = v1 - v2
        max_diff = np.max(np.abs(diff))
        mean_diff = np.mean(diff)
        std_diff = np.std(diff)
        corr = np.corrcoef(v1, v2)[0, 1]
        
        # Calculate mean and std for each dataset
        mean_v1 = np.mean(v1)
        std_v1 = np.std(v1)
        mean_v2 = np.mean(v2)
        std_v2 = np.std(v2)
        
        print(f"  {axis_name}-velocity:")
        print(f"    {label1}: mean={mean_v1:.2e}, std={std_v1:.2e}")
        print(f"    {label2}: mean={mean_v2:.2e}, std={std_v2:.2e}")
        print(f"    Mean difference: {mean_diff:.2e}")
        print(f"    Std difference: {std_diff:.2e}")
        print(f"    Max |difference|: {max_diff:.2e}")
        print(f"    Correlation: {corr:.2e}")
        if max_diff < 1e-5:
            vel_status = "Status: MATCH (max diff < 1e-5)"
            print(f"    {vel_status}")
        else:
            vel_status = "Status: MISMATCH (max diff >= 1e-5)"
            print(f"    {vel_status}")
        status_lines.append(f"  {axis_name}-velocity: {vel_status}")
    
    # Write status to text file when output_suffix is set (e.g. batch ic_0..ic_224)
    if output_suffix:
        status_file = os.path.join(output_dir, 'comparison_status.txt')
        with open(status_file, 'a') as f:
            f.write(f"=== {output_suffix} ===\n")
            f.write("\n".join(status_lines) + "\n\n")
        print(f"Status written to: {status_file}")
    
    # Create visualizations
    print("\nGenerating visualizations...")
    
    # Displacement comparison plots
    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.3)
    fig.suptitle('Displacement Field Comparison', fontsize=16, fontweight='bold')
    
    for axis_idx, axis_name in enumerate(['Z', 'Y', 'X']):
        d1 = displ1[:, axis_idx]
        d2 = displ2[:, axis_idx]
        
        # Histogram
        ax_hist = fig.add_subplot(gs[0, axis_idx])
        all_vals = np.concatenate([d1, d2])
        vmin, vmax = np.percentile(all_vals, [0.1, 99.9])
        bins = np.linspace(vmin, vmax, 50)
        ax_hist.hist(d1, bins=bins, alpha=0.6, label=label1, color='blue', density=True, edgecolor='black', linewidth=0.5)
        ax_hist.hist(d2, bins=bins, alpha=0.6, label=label2, color='red', density=True, edgecolor='black', linewidth=0.5)
        ax_hist.set_xlabel(f'{axis_name}-displacement', fontsize=10)
        ax_hist.set_ylabel('Density', fontsize=10)
        ax_hist.set_title(f'{axis_name}-displacement Distribution', fontsize=11)
        ax_hist.legend(fontsize=9)
        ax_hist.grid(True, alpha=0.3)
        # More x-axis ticks and labels
        ax_hist.xaxis.set_major_locator(MaxNLocator(nbins=10, prune=None))
        ax_hist.tick_params(axis='x', labelsize=8, rotation=45)
        
        # Scatter plot
        ax_scatter = fig.add_subplot(gs[1, axis_idx])
        if len(d1) > 50000:
            # Sample for large datasets
            rng = np.random.default_rng(0)
            sample_idx = rng.choice(len(d1), size=50000, replace=False)
            x = d1[sample_idx]
            y = d2[sample_idx]
            sampled = True
        else:
            x = d1
            y = d2
            sampled = False
        
        corr = np.corrcoef(x, y)[0, 1]
        ax_scatter.scatter(x, y, alpha=0.3, s=1)
        
        # Fit line
        if len(x) > 1:
            m, b = np.polyfit(x, y, 1)
            x_line = np.linspace(x.min(), x.max(), 100)
            ax_scatter.plot(x_line, m * x_line + b, 'r--', linewidth=2, label=f'fit: y={m:.6f}x+{b:.6f}')
        
        # Perfect agreement line
        lim_min, lim_max = min(x.min(), y.min()), max(x.max(), y.max())
        ax_scatter.plot([lim_min, lim_max], [lim_min, lim_max], 'k-', linewidth=1, alpha=0.6, label='y=x')
        
        ax_scatter.set_xlabel(f'{label1} {axis_name}-displacement', fontsize=10)
        ax_scatter.set_ylabel(f'{label2} {axis_name}-displacement', fontsize=10)
        title = f'{axis_name}-displacement Comparison (corr≈{corr:.10f}'
        if sampled:
            title += f', sampled'
        title += ')'
        ax_scatter.set_title(title, fontsize=11)
        ax_scatter.legend(fontsize=9)
        ax_scatter.grid(True, alpha=0.3)
        ax_scatter.set_aspect('equal', adjustable='box')
        # More x-axis ticks and labels
        ax_scatter.xaxis.set_major_locator(MaxNLocator(nbins=10, prune=None))
        ax_scatter.yaxis.set_major_locator(MaxNLocator(nbins=10, prune=None))
        ax_scatter.tick_params(axis='x', labelsize=8, rotation=45)
        ax_scatter.tick_params(axis='y', labelsize=8)
        
        # Difference histogram
        ax_diff = fig.add_subplot(gs[2, axis_idx])
        diff = d1 - d2
        ax_diff.hist(diff, bins=50, alpha=0.7, color='green', edgecolor='black', linewidth=0.5)
        mean_diff = np.mean(diff)
        std_diff = np.std(diff)
        ax_diff.axvline(mean_diff, color='red', linestyle='--', linewidth=2, label=f'mean: {mean_diff:.2e}')
        ax_diff.axvline(0, color='black', linestyle='-', linewidth=1, alpha=0.5)
        ax_diff.set_xlabel(f'{axis_name}-displacement difference', fontsize=10)
        ax_diff.set_ylabel('Count', fontsize=10)
        ax_diff.set_title(f'{axis_name}-displacement Difference (std={std_diff:.2e})', fontsize=11)
        ax_diff.legend(fontsize=9)
        ax_diff.grid(True, alpha=0.3)
        # More x-axis ticks and labels
        ax_diff.xaxis.set_major_locator(MaxNLocator(nbins=10, prune=None))
        ax_diff.tick_params(axis='x', labelsize=8, rotation=45)
    
    plt.savefig(os.path.join(output_dir, disp_basename), dpi=150, bbox_inches='tight')
    print(f"Saved: {os.path.join(output_dir, disp_basename)}")
    plt.close()
    
    # Velocity comparison plots
    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.3)
    fig.suptitle('Velocity Field Comparison', fontsize=16, fontweight='bold')
    
    for axis_idx, axis_name in enumerate(['Z', 'Y', 'X']):
        v1 = vel1[:, axis_idx]
        v2 = vel2[:, axis_idx]
        
        # Histogram
        ax_hist = fig.add_subplot(gs[0, axis_idx])
        all_vals = np.concatenate([v1, v2])
        vmin, vmax = np.percentile(all_vals, [0.1, 99.9])
        bins = np.linspace(vmin, vmax, 50)
        ax_hist.hist(v1, bins=bins, alpha=0.6, label=label1, color='blue', density=True, edgecolor='black', linewidth=0.5)
        ax_hist.hist(v2, bins=bins, alpha=0.6, label=label2, color='red', density=True, edgecolor='black', linewidth=0.5)
        ax_hist.set_xlabel(f'{axis_name}-velocity', fontsize=10)
        ax_hist.set_ylabel('Density', fontsize=10)
        ax_hist.set_title(f'{axis_name}-velocity Distribution', fontsize=11)
        ax_hist.legend(fontsize=9)
        ax_hist.grid(True, alpha=0.3)
        # More x-axis ticks and labels
        ax_hist.xaxis.set_major_locator(MaxNLocator(nbins=10, prune=None))
        ax_hist.tick_params(axis='x', labelsize=8, rotation=45)
        
        # Scatter plot
        ax_scatter = fig.add_subplot(gs[1, axis_idx])
        if len(v1) > 50000:
            # Sample for large datasets
            rng = np.random.default_rng(0)
            sample_idx = rng.choice(len(v1), size=50000, replace=False)
            x = v1[sample_idx]
            y = v2[sample_idx]
            sampled = True
        else:
            x = v1
            y = v2
            sampled = False
        
        corr = np.corrcoef(x, y)[0, 1]
        ax_scatter.scatter(x, y, alpha=0.3, s=1)
        
        # Fit line
        if len(x) > 1:
            m, b = np.polyfit(x, y, 1)
            x_line = np.linspace(x.min(), x.max(), 100)
            ax_scatter.plot(x_line, m * x_line + b, 'r--', linewidth=2, label=f'fit: y={m:.6f}x+{b:.6f}')
        
        # Perfect agreement line
        lim_min, lim_max = min(x.min(), y.min()), max(x.max(), y.max())
        ax_scatter.plot([lim_min, lim_max], [lim_min, lim_max], 'k-', linewidth=1, alpha=0.6, label='y=x')
        
        ax_scatter.set_xlabel(f'{label1} {axis_name}-velocity', fontsize=10)
        ax_scatter.set_ylabel(f'{label2} {axis_name}-velocity', fontsize=10)
        title = f'{axis_name}-velocity Comparison (corr≈{corr:.10f}'
        if sampled:
            title += f', sampled'
        title += ')'
        ax_scatter.set_title(title, fontsize=11)
        ax_scatter.legend(fontsize=9)
        ax_scatter.grid(True, alpha=0.3)
        ax_scatter.set_aspect('equal', adjustable='box')
        # More x-axis ticks and labels
        ax_scatter.xaxis.set_major_locator(MaxNLocator(nbins=10, prune=None))
        ax_scatter.yaxis.set_major_locator(MaxNLocator(nbins=10, prune=None))
        ax_scatter.tick_params(axis='x', labelsize=8, rotation=45)
        ax_scatter.tick_params(axis='y', labelsize=8)
        
        # Difference histogram
        ax_diff = fig.add_subplot(gs[2, axis_idx])
        diff = v1 - v2
        ax_diff.hist(diff, bins=50, alpha=0.7, color='green', edgecolor='black', linewidth=0.5)
        mean_diff = np.mean(diff)
        std_diff = np.std(diff)
        ax_diff.axvline(mean_diff, color='red', linestyle='--', linewidth=2, label=f'mean: {mean_diff:.2e}')
        ax_diff.axvline(0, color='black', linestyle='-', linewidth=1, alpha=0.5)
        ax_diff.set_xlabel(f'{axis_name}-velocity difference', fontsize=10)
        ax_diff.set_ylabel('Count', fontsize=10)
        ax_diff.set_title(f'{axis_name}-velocity Difference (std={std_diff:.2e})', fontsize=11)
        ax_diff.legend(fontsize=9)
        ax_diff.grid(True, alpha=0.3)
        # More x-axis ticks and labels
        ax_diff.xaxis.set_major_locator(MaxNLocator(nbins=10, prune=None))
        ax_diff.tick_params(axis='x', labelsize=8, rotation=45)
    
    plt.savefig(os.path.join(output_dir, vel_basename), dpi=150, bbox_inches='tight')
    print(f"Saved: {os.path.join(output_dir, vel_basename)}")
    plt.close()
    
    print("\nComparison complete!")

def main():
    if len(sys.argv) < 4:
        print("Usage: python3 compare_particle_ics.py <hermitian_ic_file> <zeldovich_ic_file> <N> [output_dir] [slab_idx] [output_suffix]")
        print("Example: python3 compare_particle_ics.py particle_ics/ic_0 /path/to/zeldovich/ic_0 4 visualizations/")
        print("Example: python3 compare_particle_ics.py particle_ics/ic_0 /path/to/zeldovich/ic_0 4 visualizations/ 0 ic_0")
        print("  (slab_idx optional: 0 for first slab, etc. If omitted, compares all particles)")
        print("  (output_suffix optional: e.g. ic_0 → displacement_comparison_ic_0.png, velocity_comparison_ic_0.png)")
        sys.exit(1)
    
    hermitian_file = sys.argv[1]
    zeldovich_file = sys.argv[2]
    N = int(sys.argv[3])
    output_dir = sys.argv[4] if len(sys.argv) > 4 else "visualizations"
    slab_idx = None
    if len(sys.argv) > 5 and sys.argv[5] not in ('', 'all'):
        try:
            slab_idx = int(sys.argv[5])
        except ValueError:
            pass
    output_suffix = sys.argv[6] if len(sys.argv) > 6 else None
    
    print("Reading particle IC files...")
    print(f"  Hermitian: {hermitian_file}")
    print(f"  Zeldovich: {zeldovich_file}")
    print(f"  N: {N}")
    if slab_idx is not None:
        print(f"  Slab index: {slab_idx}")
    else:
        print(f"  Slab index: All slabs")
    
    particles_hermitian, error1 = read_particle_ic_file(hermitian_file, N, slab_idx)
    if particles_hermitian is None:
        print(f"ERROR reading {hermitian_file}: {error1}")
        sys.exit(1)
    
    particles_zeldovich, error2 = read_particle_ic_file(zeldovich_file, N, slab_idx)
    if particles_zeldovich is None:
        print(f"ERROR reading {zeldovich_file}: {error2}")
        sys.exit(1)
    
    print(f"Read {len(particles_hermitian['i'])} particles from {hermitian_file}")
    print(f"Read {len(particles_zeldovich['i'])} particles from {zeldovich_file}")
    
    compare_particle_ics(
        particles_hermitian, particles_zeldovich,
        "my_mpi_code", "zeldovich-PLT",
        output_dir, N, output_suffix
    )

if __name__ == '__main__':
    main()

