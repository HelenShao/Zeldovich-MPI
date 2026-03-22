#!/usr/bin/env python3
"""
Compare density file histograms regardless of value order.
This helps determine if values are correct but just in wrong order,
or if the values themselves are wrong.
"""

import numpy as np
import sys
import matplotlib.pyplot as plt

def read_density_file(filename):
    """Read density file as float32 array (single precision), ignoring order."""
    try:
        with open(filename, 'rb') as f:
            data = f.read()
        density = np.frombuffer(data, dtype=np.float32)  # Single precision float32
        return density, None
    except FileNotFoundError:
        return None, f"File not found: {filename}"
    except Exception as e:
        return None, f"Error reading file: {e}"

def write_all_values(density, output_file):
    """Write all density values to text file (no coordinate information)."""
    with open(output_file, 'w') as f:
        f.write(f"# All density values (N³={len(density)} values, order ignored)\n")
        f.write("# Format: index value\n")
        f.write("# Data type: float32 (single precision)\n")
        f.write("#\n")
        for idx, val in enumerate(density):
            f.write(f"{idx:4d}  {val:15.10e}\n")

def compare_histograms(hermitian_file, zeldovich_file, output_dir="."):
    """Compare histograms of hermitian and zeldovich density values."""
    
    # Read both files
    hermitian, error = read_density_file(hermitian_file)
    if hermitian is None:
        print(f"ERROR reading hermitian file: {error}")
        return
    
    zeldovich, error = read_density_file(zeldovich_file)
    if zeldovich is None:
        print(f"ERROR reading zeldovich file: {error}")
        return
    
    # Check sizes match
    if len(hermitian) != len(zeldovich):
        print(f"WARNING: Size mismatch - hermitian: {len(hermitian)}, zeldovich: {len(zeldovich)}")
    
    # Write all values to text files
    write_all_values(hermitian, f"{output_dir}/hermitian_all_values.txt")
    write_all_values(zeldovich, f"{output_dir}/zeldovich_all_values.txt")
    print(f"Extracted all values:")
    print(f"  Hermitian: {len(hermitian)} values -> {output_dir}/hermitian_all_values.txt")
    print(f"  Zeldovich: {len(zeldovich)} values -> {output_dir}/zeldovich_all_values.txt")
    
    # Statistics
    print("\nStatistics (order ignored):")
    print(f"Hermitian: min={hermitian.min():.10e}, max={hermitian.max():.10e}, mean={hermitian.mean():.10e}, std={hermitian.std():.10e}")
    print(f"Zeldovich: min={zeldovich.min():.10e}, max={zeldovich.max():.10e}, mean={zeldovich.mean():.10e}, std={zeldovich.std():.10e}")
    
    # Sort values for comparison
    hermitian_sorted = np.sort(hermitian)
    zeldovich_sorted = np.sort(zeldovich)
    
    # Compare sorted values
    # Using single precision (float32) tolerances: ~7 decimal digits of precision
    max_diff = np.max(np.abs(hermitian_sorted - zeldovich_sorted))
    mean_diff = np.mean(np.abs(hermitian_sorted - zeldovich_sorted))
    # Relative tolerance for single precision: ~1e-6 to 1e-7 for values near 1.0
    # For density values, use absolute tolerance based on typical value range
    abs_tol = 1e-6  # Appropriate for single precision float32
    rel_tol = 1e-5  # Relative tolerance
    
    print(f"\nSorted value comparison (single precision float32):")
    print(f"  Max absolute difference: {max_diff:.10e}")
    print(f"  Mean absolute difference: {mean_diff:.10e}")
    print(f"  Tolerance: {abs_tol:.1e} (absolute) or {rel_tol:.1e} (relative)")
    
    # Check both absolute and relative differences
    rel_diff = max_diff / (np.max(np.abs(zeldovich_sorted)) + 1e-10)
    
    if max_diff < abs_tol or rel_diff < rel_tol:
        print("  ✓ Values match (within single precision tolerance) - likely just an ordering issue!")
    elif max_diff < 1e-4:
        print("  ~ Values close (within 1e-4) - may be ordering or small numerical differences")
    else:
        print("  ✗ Values differ significantly - not just an ordering issue")
    
    # Create histograms
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Histogram 1: Overlaid histograms
    ax = axes[0, 0]
    bins = 50
    ax.hist(hermitian, bins=bins, alpha=0.5, label='Hermitian', color='blue', density=True)
    ax.hist(zeldovich, bins=bins, alpha=0.5, label='Zeldovich', color='red', density=True)
    ax.set_xlabel('Density Value')
    ax.set_ylabel('Probability Density')
    ax.set_title('Density Value Histograms (Overlaid)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Histogram 2: Side-by-side
    ax = axes[0, 1]
    ax.hist(hermitian, bins=bins, alpha=0.7, label='Hermitian', color='blue', density=True)
    ax.set_xlabel('Density Value')
    ax.set_ylabel('Probability Density')
    ax.set_title('Hermitian Density Histogram')
    ax.grid(True, alpha=0.3)
    
    ax = axes[1, 0]
    ax.hist(zeldovich, bins=bins, alpha=0.7, label='Zeldovich', color='red', density=True)
    ax.set_xlabel('Density Value')
    ax.set_ylabel('Probability Density')
    ax.set_title('Zeldovich Density Histogram')
    ax.grid(True, alpha=0.3)
    
    # Scatter plot: sorted values
    ax = axes[1, 1]
    min_len = min(len(hermitian_sorted), len(zeldovich_sorted))
    ax.scatter(hermitian_sorted[:min_len], zeldovich_sorted[:min_len], alpha=0.5, s=1)
    ax.plot([hermitian_sorted.min(), hermitian_sorted.max()], 
            [hermitian_sorted.min(), hermitian_sorted.max()], 
            'r--', label='y=x')
    ax.set_xlabel('Hermitian (Sorted)')
    ax.set_ylabel('Zeldovich (Sorted)')
    ax.set_title('Sorted Value Comparison')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = f"{output_dir}/density_histogram_comparison.png"
    plt.savefig(output_file, dpi=150)
    print(f"\nHistogram comparison saved to: {output_file}")
    
    # Value-by-value comparison (if sizes match)
    if len(hermitian) == len(zeldovich):
        # Check if values match when sorted (using single precision tolerance)
        # float32 has ~7 decimal digits, so atol=1e-6 is appropriate
        matches_when_sorted = np.allclose(hermitian_sorted, zeldovich_sorted, atol=1e-6, rtol=1e-5)
        print(f"\nValues match when sorted (atol=1e-6, rtol=1e-5): {matches_when_sorted}")
        
        # Check if values match in original order
        matches_original = np.allclose(hermitian, zeldovich, atol=1e-6, rtol=1e-5)
        print(f"Values match in original order (atol=1e-6, rtol=1e-5): {matches_original}")
        
        if matches_when_sorted and not matches_original:
            print("  → This confirms: values are correct but in wrong order!")
        elif matches_when_sorted and matches_original:
            print("  → Values match in both sorted and original order!")
        else:
            print("  → Values differ even when sorted - not just an ordering issue")

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 compare_density_histograms.py <hermitian_density_file> <zeldovich_density_file> [output_dir]")
        print("Example: python3 compare_density_histograms.py particle_ics/density4 /path/to/zeldovich/output/density4 .")
        sys.exit(1)
    
    hermitian_file = sys.argv[1]
    zeldovich_file = sys.argv[2]
    output_dir = sys.argv[3] if len(sys.argv) >= 4 else "."
    
    compare_histograms(hermitian_file, zeldovich_file, output_dir)

if __name__ == '__main__':
    main()

