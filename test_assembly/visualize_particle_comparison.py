import struct
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import pearsonr
from scipy.optimize import curve_fit

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

def load_all_particles_from_dir(directory):
    """Load all particles from all ic_* files in a directory."""
    particles_all = {}
    files = sorted([f for f in os.listdir(directory) if f.startswith('ic_')])
    
    for f in files:
        filepath = os.path.join(directory, f)
        particles, count = read_all_particles(filepath)
        if particles is not None:
            particles_all.update(particles)
    
    return particles_all

def extract_common_particles(particles1, particles2):
    """Extract arrays of values for common particles."""
    common_keys = set(particles1.keys()) & set(particles2.keys())
    
    if len(common_keys) == 0:
        return None
    
    # Extract displacement and velocity arrays
    displ1 = np.array([[particles1[k]['displ'][i] for k in common_keys] for i in range(3)])
    displ2 = np.array([[particles2[k]['displ'][i] for k in common_keys] for i in range(3)])
    vel1 = np.array([[particles1[k]['vel'][i] for k in common_keys] for i in range(3)])
    vel2 = np.array([[particles2[k]['vel'][i] for k in common_keys] for i in range(3)])
    
    return {
        'displ1': displ1,
        'displ2': displ2,
        'vel1': vel1,
        'vel2': vel2,
        'n_particles': len(common_keys)
    }

def plot_histograms_comparison(data, field_name, axes_labels, output_dir):
    """Create histogram comparison plots for each axis."""
    if field_name == 'displacement':
        field1 = data['displ1']
        field2 = data['displ2']
    else:
        field1 = data['vel1']
        field2 = data['vel2']
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f'{field_name.capitalize()} Distribution Comparison', fontsize=16, fontweight='bold')
    
    for axis_idx, axis_label in enumerate(axes_labels):
        ax = axes[axis_idx]
        values1 = field1[axis_idx]
        values2 = field2[axis_idx]
        
        # Determine bin range
        all_values = np.concatenate([values1, values2])
        vmin, vmax = np.percentile(all_values, [0.1, 99.9])
        bins = np.linspace(vmin, vmax, 50)
        
        # Plot histograms
        ax.hist(values1, bins=bins, alpha=0.6, label='hermitian_3d_matrix', 
                color='blue', density=True, edgecolor='black', linewidth=0.5)
        ax.hist(values2, bins=bins, alpha=0.6, label='zeldovich-PLT', 
                color='red', density=True, edgecolor='black', linewidth=0.5)
        
        # Calculate statistics
        mean1, std1 = np.mean(values1), np.std(values1)
        mean2, std2 = np.mean(values2), np.std(values2)
        
        ax.axvline(mean1, color='blue', linestyle='--', linewidth=2, alpha=0.7)
        ax.axvline(mean2, color='red', linestyle='--', linewidth=2, alpha=0.7)
        
        ax.set_xlabel(f'{axis_label} {field_name.capitalize()}', fontsize=12)
        ax.set_ylabel('Density', fontsize=12)
        ax.set_title(f'{axis_label}-axis', fontsize=12)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, f'{field_name}_histograms.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

def plot_scatter_comparison(data, field_name, axes_labels, output_dir):
    """Create scatter plots comparing values for each axis."""
    if field_name == 'displacement':
        field1 = data['displ1']
        field2 = data['displ2']
    else:
        field1 = data['vel1']
        field2 = data['vel2']
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f'{field_name.capitalize()} Scatter Comparison', fontsize=16, fontweight='bold')
    
    for axis_idx, axis_label in enumerate(axes_labels):
        ax = axes[axis_idx]
        values1 = field1[axis_idx]
        values2 = field2[axis_idx]
        
        # Subsample if too many points for visualization
        n_points = len(values1)
        if n_points > 10000:
            indices = np.random.choice(n_points, 10000, replace=False)
            values1_plot = values1[indices]
            values2_plot = values2[indices]
        else:
            values1_plot = values1
            values2_plot = values2
        
        # Scatter plot
        ax.scatter(values1_plot, values2_plot, alpha=0.3, s=1, color='blue')
        
        # Calculate correlation
        corr, p_value = pearsonr(values1, values2)
        
        # Fit a line
        if len(values1) > 1:
            # Fit y = a*x + b
            coeffs = np.polyfit(values1, values2, 1)
            x_line = np.linspace(values1.min(), values1.max(), 100)
            y_line = np.polyval(coeffs, x_line)
            ax.plot(x_line, y_line, 'r--', linewidth=2, 
                   label=f'y={coeffs[0]:.4f}x+{coeffs[1]:.4f}')
        
        # Perfect agreement line
        all_vals = np.concatenate([values1, values2])
        lim_min, lim_max = np.percentile(all_vals, [0.1, 99.9])
        ax.plot([lim_min, lim_max], [lim_min, lim_max], 'k-', linewidth=1, 
               alpha=0.5, label='y=x (perfect agreement)')
        
        ax.set_xlabel(f'hermitian_3d_matrix {axis_label}', fontsize=12)
        ax.set_ylabel(f'zeldovich-PLT {axis_label}', fontsize=12)
        ax.set_title(f'{axis_label}-axis\nCorrelation: {corr:.6f}\np-value: {p_value:.2e}', 
                    fontsize=10)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal', adjustable='box')
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, f'{field_name}_scatter.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

def plot_difference_histograms(data, field_name, axes_labels, output_dir):
    """Create histograms of differences between the two codes."""
    if field_name == 'displacement':
        field1 = data['displ1']
        field2 = data['displ2']
    else:
        field1 = data['vel1']
        field2 = data['vel2']
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f'{field_name.capitalize()} Difference Distribution', fontsize=16, fontweight='bold')
    
    for axis_idx, axis_label in enumerate(axes_labels):
        ax = axes[axis_idx]
        values1 = field1[axis_idx]
        values2 = field2[axis_idx]
        differences = values1 - values2
        
        # Plot histogram of differences
        ax.hist(differences, bins=50, alpha=0.7, color='purple', 
                edgecolor='black', linewidth=0.5)
        
        # Statistics
        mean_diff = np.mean(differences)
        std_diff = np.std(differences)
        rms_diff = np.sqrt(np.mean(differences**2))
        
        ax.axvline(mean_diff, color='red', linestyle='--', linewidth=2, 
                  label=f'Mean: {mean_diff:.6f}')
        ax.axvline(0, color='black', linestyle='-', linewidth=1, alpha=0.5)
        
        ax.set_xlabel(f'Difference ({axis_label})', fontsize=12)
        ax.set_ylabel('Count', fontsize=12)
        ax.set_title(f'{axis_label}-axis\nMean={mean_diff:.6f}, Std={std_diff:.6f}, RMS={rms_diff:.6f}', 
                    fontsize=10)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, f'{field_name}_differences.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

def plot_correlation_matrix(data, output_dir):
    """Create correlation matrices for both fields."""
    fig, axes = plt.subplots(2, 1, figsize=(10, 12))
    
    # Displacement correlations
    displ1 = data['displ1']
    displ2 = data['displ2']
    corr_displ = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            corr, _ = pearsonr(displ1[i], displ2[j])
            corr_displ[i, j] = corr
    
    im1 = axes[0].imshow(corr_displ, cmap='coolwarm', vmin=-1, vmax=1, aspect='auto')
    axes[0].set_xticks([0, 1, 2])
    axes[0].set_yticks([0, 1, 2])
    axes[0].set_xticklabels(['X', 'Y', 'Z'])
    axes[0].set_yticklabels(['X', 'Y', 'Z'])
    axes[0].set_title('Displacement Cross-Correlation Matrix\n(hermitian_3d_matrix vs zeldovich-PLT)', 
                     fontsize=12, fontweight='bold')
    axes[0].set_xlabel('zeldovich-PLT axis', fontsize=11)
    axes[0].set_ylabel('hermitian_3d_matrix axis', fontsize=11)
    
    # Add text annotations
    for i in range(3):
        for j in range(3):
            text = axes[0].text(j, i, f'{corr_displ[i, j]:.4f}',
                              ha="center", va="center", color="black", fontweight='bold')
    plt.colorbar(im1, ax=axes[0])
    
    # Velocity correlations
    vel1 = data['vel1']
    vel2 = data['vel2']
    corr_vel = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            corr, _ = pearsonr(vel1[i], vel2[j])
            corr_vel[i, j] = corr
    
    im2 = axes[1].imshow(corr_vel, cmap='coolwarm', vmin=-1, vmax=1, aspect='auto')
    axes[1].set_xticks([0, 1, 2])
    axes[1].set_yticks([0, 1, 2])
    axes[1].set_xticklabels(['X', 'Y', 'Z'])
    axes[1].set_yticklabels(['X', 'Y', 'Z'])
    axes[1].set_title('Velocity Cross-Correlation Matrix\n(hermitian_3d_matrix vs zeldovich-PLT)', 
                     fontsize=12, fontweight='bold')
    axes[1].set_xlabel('zeldovich-PLT axis', fontsize=11)
    axes[1].set_ylabel('hermitian_3d_matrix axis', fontsize=11)
    
    # Add text annotations
    for i in range(3):
        for j in range(3):
            text = axes[1].text(j, i, f'{corr_vel[i, j]:.4f}',
                              ha="center", va="center", color="black", fontweight='bold')
    plt.colorbar(im2, ax=axes[1])
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'correlation_matrices.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

def print_statistics(data):
    """Print summary statistics."""
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    print(f"Number of common particles: {data['n_particles']}")
    
    print("\nDisplacement Statistics:")
    print("-" * 80)
    axes_labels = ['X', 'Y', 'Z']
    for axis_idx, axis_label in enumerate(axes_labels):
        displ1 = data['displ1'][axis_idx]
        displ2 = data['displ2'][axis_idx]
        corr, p_val = pearsonr(displ1, displ2)
        mean_diff = np.mean(displ1 - displ2)
        std_diff = np.std(displ1 - displ2)
        rms_diff = np.sqrt(np.mean((displ1 - displ2)**2))
        print(f"{axis_label}-axis: Corr={corr:.6f}, Mean_diff={mean_diff:.6e}, "
              f"Std_diff={std_diff:.6e}, RMS_diff={rms_diff:.6e}")
    
    print("\nVelocity Statistics:")
    print("-" * 80)
    for axis_idx, axis_label in enumerate(axes_labels):
        vel1 = data['vel1'][axis_idx]
        vel2 = data['vel2'][axis_idx]
        corr, p_val = pearsonr(vel1, vel2)
        mean_diff = np.mean(vel1 - vel2)
        std_diff = np.std(vel1 - vel2)
        rms_diff = np.sqrt(np.mean((vel1 - vel2)**2))
        print(f"{axis_label}-axis: Corr={corr:.6f}, Mean_diff={mean_diff:.6e}, "
              f"Std_diff={std_diff:.6e}, RMS_diff={rms_diff:.6e}")

def main():
    # Accept command-line arguments if provided, otherwise use defaults
    if len(sys.argv) >= 3:
        dir1 = sys.argv[1]
        dir2 = sys.argv[2]
    else:
        # Default paths
        dir1 = "/home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_assembly/particle_ics"
        dir2 = "/home/helenshao/InitialConditions/zeldovich-PLT/output"
    
    # Output directory for plots (accept optional 3rd argument)
    if len(sys.argv) >= 4:
        output_dir = sys.argv[3]
    else:
        # Default output directory
        output_dir = "/home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_assembly/visualizations"
    os.makedirs(output_dir, exist_ok=True)
    
    print("Loading particles from both locations...")
    particles1 = load_all_particles_from_dir(dir1)
    particles2 = load_all_particles_from_dir(dir2)
    
    print(f"Loaded {len(particles1)} particles from hermitian_3d_matrix_production")
    print(f"Loaded {len(particles2)} particles from zeldovich-PLT")
    
    # Extract common particles
    data = extract_common_particles(particles1, particles2)
    if data is None:
        print("ERROR: No common particles found!")
        return
    
    print(f"Found {data['n_particles']} common particles")
    
    # Print statistics
    print_statistics(data)
    
    # Create visualizations
    axes_labels = ['X', 'Y', 'Z']
    
    print("\nGenerating visualizations...")
    
    # Histograms
    plot_histograms_comparison(data, 'displacement', axes_labels, output_dir)
    plot_histograms_comparison(data, 'velocity', axes_labels, output_dir)
    
    # Scatter plots
    plot_scatter_comparison(data, 'displacement', axes_labels, output_dir)
    plot_scatter_comparison(data, 'velocity', axes_labels, output_dir)
    
    # Difference histograms
    plot_difference_histograms(data, 'displacement', axes_labels, output_dir)
    plot_difference_histograms(data, 'velocity', axes_labels, output_dir)
    
    # Correlation matrices
    plot_correlation_matrix(data, output_dir)
    
    print(f"\nAll visualizations saved to: {output_dir}")

if __name__ == '__main__':
    main()

