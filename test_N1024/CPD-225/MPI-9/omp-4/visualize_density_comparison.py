#!/usr/bin/env python3
"""
Visualize density field comparison between hermitian_3d_matrix_production and zeldovich-PLT.
Similar to visualize_particle_comparison.py but for density fields.
"""

import struct
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import pearsonr

# Load plotting parameters
sys.path.insert(0, os.path.expanduser('~/python_modules'))
try:
    from matplotlib_plotting_params import apply_plotting_params
    apply_plotting_params()
except ImportError:
    print("Warning: Could not import matplotlib_plotting_params, using default settings")

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

def plot_histograms_comparison(density1, density2, label1, label2, output_dir):
    """Create histogram comparison plot."""
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    fig.suptitle('Density Field Distribution Comparison', fontweight='bold')
    
    # Determine bin range
    all_values = np.concatenate([density1, density2])
    vmin, vmax = np.percentile(all_values, [0.1, 99.9])
    bins = np.linspace(vmin, vmax, 50)
    
    # Plot histograms
    ax.hist(density1, bins=bins, alpha=0.6, label=label1, 
            color='blue', density=True, edgecolor='black', linewidth=0.5)
    ax.hist(density2, bins=bins, alpha=0.6, label=label2, 
            color='red', density=True, edgecolor='black', linewidth=0.5)
    
    # Calculate statistics
    mean1, std1 = np.mean(density1), np.std(density1)
    mean2, std2 = np.mean(density2), np.std(density2)
    
    ax.axvline(mean1, color='blue', linestyle='--', linewidth=2, alpha=0.7, 
              label=f'{label1} mean: {mean1:.2e}')
    ax.axvline(mean2, color='red', linestyle='--', linewidth=2, alpha=0.7,
              label=f'{label2} mean: {mean2:.2e}')
    
    ax.set_xlabel('Density Value')
    ax.set_ylabel('Histogram')
    ax.set_title(f'Density Distribution\n{label1}: mean={mean1:.2e}, std={std1:.2e}\n{label2}: mean={mean2:.2e}, std={std2:.2e}')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout(rect=[0, 0, 0.88, 1])
    output_file = os.path.join(output_dir, 'density_histograms.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

def plot_scatter_comparison(density1, density2, label1, label2, output_dir):
    """Create scatter plot comparing values."""
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    fig.suptitle('Density Field Scatter Comparison', fontweight='bold')
    
    # Scatter plot
    ax.scatter(density1, density2, alpha=0.5, s=20, color='blue')
    
    # Calculate correlation
    corr, p_value = pearsonr(density1, density2)
    
    # Fit a line
    if len(density1) > 1:
        coeffs = np.polyfit(density1, density2, 1)
        x_line = np.linspace(density1.min(), density1.max(), 100)
        y_line = np.polyval(coeffs, x_line)
        ax.plot(x_line, y_line, 'r--', linewidth=2, 
               label=f'y={coeffs[0]:.2e}x+{coeffs[1]:.2e}')
    
    # Perfect agreement line
    all_vals = np.concatenate([density1, density2])
    lim_min, lim_max = np.percentile(all_vals, [0.1, 99.9])
    ax.plot([lim_min, lim_max], [lim_min, lim_max], 'k-', linewidth=1, 
           alpha=0.5, label='y=x (perfect agreement)')
    
    ax.set_xlabel(f'{label1} Density')
    ax.set_ylabel(f'{label2} Density')
    ax.set_title(f'Density Comparison\nCorrelation: {corr:.2e}\np-value: {p_value:.2e}')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal', adjustable='box')
    
    plt.tight_layout(rect=[0, 0, 0.88, 1])
    output_file = os.path.join(output_dir, 'density_scatter.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

def plot_difference_histogram(density1, density2, label1, label2, output_dir):
    """Create histogram of differences."""
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    fig.suptitle('Density Field Difference Distribution', fontweight='bold')
    
    differences = density1 - density2
    
    # Plot histogram of differences
    ax.hist(differences, bins=50, alpha=0.7, color='purple', 
            edgecolor='black', linewidth=0.5)
    
    # Statistics
    mean_diff = np.mean(differences)
    std_diff = np.std(differences)
    rms_diff = np.sqrt(np.mean(differences**2))
    max_abs_diff = np.max(np.abs(differences))
    
    ax.axvline(mean_diff, color='red', linestyle='--', linewidth=2, 
              label=f'Mean: {mean_diff:.2e}')
    ax.axvline(0, color='black', linestyle='-', linewidth=1, alpha=0.5)
    
    ax.set_xlabel(f'Difference ({label1} - {label2})')
    ax.set_ylabel('Count')
    ax.set_title(f'Density Difference Distribution\nMean={mean_diff:.2e}, Std={std_diff:.2e}, RMS={rms_diff:.2e}, Max|diff|={max_abs_diff:.2e}')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout(rect=[0, 0, 0.88, 1])
    output_file = os.path.join(output_dir, 'density_differences.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

def print_statistics(density1, density2, label1, label2):
    """Print summary statistics."""
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    print(f"Number of density values: {len(density1)}")
    
    # Individual statistics
    print(f"\n{label1} Statistics:")
    print("-" * 80)
    print(f"  Min: {density1.min():.2e}")
    print(f"  Max: {density1.max():.2e}")
    print(f"  Mean: {density1.mean():.2e}")
    print(f"  Std: {density1.std():.2e}")
    print(f"  RMS: {np.sqrt(np.mean(density1**2)):.2e}")
    
    print(f"\n{label2} Statistics:")
    print("-" * 80)
    print(f"  Min: {density2.min():.2e}")
    print(f"  Max: {density2.max():.2e}")
    print(f"  Mean: {density2.mean():.2e}")
    print(f"  Std: {density2.std():.2e}")
    print(f"  RMS: {np.sqrt(np.mean(density2**2)):.2e}")
    
    # Comparison statistics
    differences = density1 - density2
    corr, p_value = pearsonr(density1, density2)
    
    print(f"\nComparison Statistics:")
    print("-" * 80)
    print(f"  Correlation: {corr:.2e}")
    print(f"  p-value: {p_value:.2e}")
    print(f"  Mean difference: {np.mean(differences):.2e}")
    print(f"  Std difference: {np.std(differences):.2e}")
    print(f"  RMS difference: {np.sqrt(np.mean(differences**2)):.2e}")
    print(f"  Max absolute difference: {np.max(np.abs(differences)):.2e}")
    
    # Count matches
    tol = 1e-6
    exact_matches = np.sum(differences == 0)
    near_matches = np.sum(np.abs(differences) < tol) - exact_matches
    differences_count = np.sum(np.abs(differences) >= tol)
    
    print(f"\nMatch Statistics (tolerance = {tol}):")
    print("-" * 80)
    print(f"  Exact matches (diff == 0): {exact_matches}")
    print(f"  Near matches (diff < {tol}): {near_matches}")
    print(f"  Differences (diff >= {tol}): {differences_count}")

def get_N_from_param_file(param_file):
    """Read NP from parameter file and return N = cube_root(NP)."""
    try:
        with open(param_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('NP') and '=' in line:
                    # Parse NP = value (handle comments)
                    parts = line.split('#')[0].split('=')
                    if len(parts) >= 2:
                        np_value = int(parts[1].strip())
                        N = round(np_value ** (1/3))
                        return N
    except Exception as e:
        print(f"Warning: Could not read N from {param_file}: {e}")
    return None

def main():
    # Use current working directory for relative paths (not script directory)
    # This allows scripts to work from any directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    cwd = os.getcwd()
    
    # Default paths (relative to current working directory or script directory)
    if len(sys.argv) >= 3:
        file1 = sys.argv[1]
        file2 = sys.argv[2]
    else:
        # Auto-detect N from parameter file or directory structure
        # Try to detect from current directory or script directory
        if os.path.basename(cwd) == os.path.basename(script_dir):
            test_dir = cwd
        else:
            test_dir = script_dir
        
        # Try to find parameter file and read N from it
        N = None
        import glob
        param_files = glob.glob(os.path.join(test_dir, "*.par"))
        if param_files:
            N = get_N_from_param_file(param_files[0])
            if N:
                print(f"Detected N={N} from parameter file: {param_files[0]}")
        
        # Fallback: try to detect from directory name
        if N is None:
            parent_dir = os.path.basename(os.path.dirname(test_dir))
            import re
            match = re.search(r'N(\d+)', parent_dir)
            if match:
                N = int(match.group(1))
                print(f"Detected N={N} from directory name: {parent_dir}")
        
        # Default fallback
        if N is None:
            N = 4
            print(f"Warning: Could not detect N, defaulting to N={N}")
        
        density_filename = f'density{N}'
        
        # Determine zeldovich output directory based on test type
        test_subdir = os.path.basename(test_dir)  # e.g., "MPI", "PLT", "sc_ic", "sc_dens", "MPI-16-nodes", "MPI-32-nodes"
        if test_subdir == 'MPI':
            zeldovich_output_dir = f'output_N{N}_MPI'
        elif test_subdir == 'MPI-16-nodes':
            zeldovich_output_dir = f'output_N{N}_MPI_16_nodes'
        elif test_subdir == 'MPI-32-nodes':
            zeldovich_output_dir = f'output_N{N}_MPI_32_nodes'
        elif test_subdir == 'PLT':
            zeldovich_output_dir = f'output_N{N}_PLT'
        elif 'sc_ic' in test_subdir:
            zeldovich_output_dir = f'output_N{N}_sc_ic'
        elif 'sc_dens' in test_subdir:
            zeldovich_output_dir = f'output_N{N}_sc_dens'
        else:
            # Default: try to match test subdirectory name
            zeldovich_output_dir = f'output_N{N}_{test_subdir}'
        
        # Use relative paths from current working directory
        file1 = os.path.join(cwd, "particle_ics", density_filename)
        if not os.path.exists(file1):
            # Fallback: try relative to script directory
            file1 = os.path.join(script_dir, "particle_ics", density_filename)
        
        # Zeldovich output: try common location (relative to user home)
        # User can override by providing as command-line argument
        zeldovich_base = os.path.expanduser("~/InitialConditions/zeldovich-PLT")
        file2 = os.path.join(zeldovich_base, zeldovich_output_dir, density_filename)
    
    # Output directory for plots (accept optional 3rd argument)
    if len(sys.argv) >= 4:
        output_dir = sys.argv[3]
    else:
        # Default to current directory's visualizations subdirectory
        output_dir = os.path.join(cwd, "visualizations")
        if not os.path.exists(output_dir):
            # Fallback: use script directory
            output_dir = os.path.join(script_dir, "visualizations")
    os.makedirs(output_dir, exist_ok=True)
    
    label1 = "my_mpi_code"
    label2 = "zeldovich-PLT"
    
    print("Reading density files...")
    density1, error1 = read_density_file(file1)
    density2, error2 = read_density_file(file2)
    
    if density1 is None:
        print(f"ERROR reading {file1}: {error1}")
        return
    if density2 is None:
        print(f"ERROR reading {file2}: {error2}")
        return
    
    print(f"Read {len(density1)} values from {file1}")
    print(f"Read {len(density2)} values from {file2}")
    
    if len(density1) != len(density2):
        print(f"WARNING: Different array sizes ({len(density1)} vs {len(density2)})")
        min_len = min(len(density1), len(density2))
        density1 = density1[:min_len]
        density2 = density2[:min_len]
        print(f"Using first {min_len} values for comparison")
    
    # Print statistics
    print_statistics(density1, density2, label1, label2)
    
    # Create visualizations
    print("\nGenerating visualizations...")
    
    plot_histograms_comparison(density1, density2, label1, label2, output_dir)
    plot_scatter_comparison(density1, density2, label1, label2, output_dir)
    plot_difference_histogram(density1, density2, label1, label2, output_dir)
    
    print(f"\nAll visualizations saved to: {output_dir}")

if __name__ == '__main__':
    main()

