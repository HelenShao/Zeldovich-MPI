#!/usr/bin/env python3
"""
Compare output fields between Hermitian 3D and Zeldovich-PLT codes
Tests normalization by computing RMS ratio
"""

import numpy as np
import glob
import os
import sys

def read_hermitian_zslab(z, rank, N, narray, x_count):
    """
    Read one Z-slab from hermitian 3D output
    
    Returns: complex array of shape (narray, N, x_count)
    """
    filename = f"hermitian_3d_z{z}_rank{rank}.bin"
    
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File not found: {filename}")
    
    # Read binary complex data (each element: 2 doubles = 16 bytes)
    data = np.fromfile(filename, dtype=np.complex128)
    
    # Reshape to [Array][Y][X]
    expected_size = narray * N * x_count
    if len(data) != expected_size:
        raise ValueError(f"File size mismatch: expected {expected_size}, got {len(data)}")
    
    data = data.reshape(narray, N, x_count)
    
    return data


def extract_density_full_slab(z, N, narray=4, output_dir="./"):
    """
    Extract full density field for one Z-slab from all ranks
    
    Returns: density array of shape (N, N)
    """
    # Find all files for this Z
    pattern = f"{output_dir}/hermitian_3d_z{z}_rank*.bin"
    files = sorted(glob.glob(pattern))
    
    if len(files) == 0:
        raise FileNotFoundError(f"No files found matching: {pattern}")
    
    print(f"Found {len(files)} rank files for Z={z}")
    
    density_slices = []
    
    for filename in files:
        # Parse rank number from filename
        basename = os.path.basename(filename)
        rank = int(basename.split('rank')[1].split('.')[0])
        
        # Read entire file to infer x_count
        data_raw = np.fromfile(filename, dtype=np.complex128)
        x_count = len(data_raw) // (narray * N)
        
        print(f"  Rank {rank}: x_count={x_count}")
        
        # Reshape
        data = data_raw.reshape(narray, N, x_count)
        
        # Extract density (array 0, real part)
        density_slice = data[0, :, :].real
        
        # Check imaginary part (should be ~0 after FFT)
        imag_max = np.max(np.abs(data[0, :, :].imag))
        print(f"    Max imaginary part: {imag_max:.3e} (should be ~0)")
        
        density_slices.append(density_slice)
    
    # Concatenate along X dimension
    density_full = np.concatenate(density_slices, axis=1)
    
    if density_full.shape != (N, N):
        print(f"Warning: Expected shape ({N}, {N}), got {density_full.shape}")
    
    return density_full


def read_zeldovich_density(z, N, output_dir="./"):
    """
    Read density field from Zeldovich-PLT output
    
    Zeldovich writes density to: zeldovich.dens_{N}
    Format: binary float32 array [Z][Y][X]
    """
    dens_filename = f"{output_dir}/zeldovich.dens_{N}"
    
    if not os.path.exists(dens_filename):
        raise FileNotFoundError(
            f"Density file not found: {dens_filename}\n"
            f"Make sure Zeldovich-PLT was run with qdensity=1"
        )
    
    # Read full density cube
    dens_data = np.fromfile(dens_filename, dtype=np.float32)
    expected_size = N * N * N
    
    if len(dens_data) != expected_size:
        raise ValueError(f"File size mismatch: expected {expected_size}, got {len(dens_data)}")
    
    dens_data = dens_data.reshape(N, N, N)
    
    # Extract Z-slab
    zeldovich_dens = dens_data[z, :, :]
    
    return zeldovich_dens


def compare_density_fields(z, N, hermitian_dir="./", zeldovich_dir="./", narray=4):
    """
    Compare density fields from both codes for one Z-slab
    """
    print(f"\n{'='*60}")
    print(f"Comparing density fields for Z={z}, N={N}")
    print(f"{'='*60}")
    
    # Read hermitian 3D density
    print(f"\nReading Hermitian 3D output from: {hermitian_dir}")
    hermitian_dens = extract_density_full_slab(z, N, narray, hermitian_dir)
    
    # Read zeldovich density
    print(f"\nReading Zeldovich-PLT output from: {zeldovich_dir}")
    zeldovich_dens = read_zeldovich_density(z, N, zeldovich_dir)
    
    # Statistics
    print(f"\n{'='*60}")
    print(f"Field Statistics")
    print(f"{'='*60}")
    
    print(f"\nHermitian 3D:")
    print(f"  Shape:   {hermitian_dens.shape}")
    print(f"  Mean:    {np.mean(hermitian_dens):15.6e}")
    print(f"  Std:     {np.std(hermitian_dens):15.6e}")
    print(f"  Min:     {np.min(hermitian_dens):15.6e}")
    print(f"  Max:     {np.max(hermitian_dens):15.6e}")
    
    print(f"\nZeldovich-PLT:")
    print(f"  Shape:   {zeldovich_dens.shape}")
    print(f"  Mean:    {np.mean(zeldovich_dens):15.6e}")
    print(f"  Std:     {np.std(zeldovich_dens):15.6e}")
    print(f"  Min:     {np.min(zeldovich_dens):15.6e}")
    print(f"  Max:     {np.max(zeldovich_dens):15.6e}")
    
    # Normalization test
    print(f"\n{'='*60}")
    print(f"Normalization Test")
    print(f"{'='*60}")
    
    rms_hermitian = np.std(hermitian_dens)
    rms_zeldovich = np.std(zeldovich_dens)
    rms_ratio = rms_hermitian / rms_zeldovich
    
    print(f"\nRMS (Hermitian):      {rms_hermitian:.6e}")
    print(f"RMS (Zeldovich):      {rms_zeldovich:.6e}")
    print(f"RMS ratio (H/Z):      {rms_ratio:.6e}")
    
    # Check different normalization scenarios
    N3 = N**3
    print(f"\nExpected ratios:")
    print(f"  If normalized:      {1.0:.6e}")
    print(f"  If missing 1/N³:    {N3:.6e}  (N³ = {N3})")
    print(f"  If extra 1/N³:      {1.0/N3:.6e}")
    
    print(f"\nDiagnosis:")
    if abs(rms_ratio - 1.0) < 0.01:
        print("  ✓ Normalization MATCHES! No correction needed.")
        status = "MATCH"
    elif abs(rms_ratio - N3) < N3 * 0.01:
        print(f"  ✗ Hermitian is {N3:.0f}x too large - Missing 1/N³ normalization")
        print(f"    → Add: multiply output by 1/N³ = {1.0/N3:.6e}")
        status = "MISSING_N3"
    elif abs(rms_ratio - 1.0/N3) < 1.0/N3 * 0.01:
        print(f"  ✗ Hermitian is {N3:.0f}x too small - Extra 1/N³ normalization")
        print(f"    → Remove: divide output by 1/N³")
        status = "EXTRA_N3"
    else:
        print(f"  ? Unexpected ratio: {rms_ratio:.6e}")
        print(f"    Not matching standard normalization patterns")
        status = "UNKNOWN"
    
    # Mean difference (should be ~0 for density perturbations)
    mean_hermitian = np.mean(hermitian_dens)
    mean_zeldovich = np.mean(zeldovich_dens)
    print(f"\nMean values (should be ~0 for density perturbations):")
    print(f"  Hermitian: {mean_hermitian:.6e}")
    print(f"  Zeldovich: {mean_zeldovich:.6e}")
    
    return {
        'hermitian_dens': hermitian_dens,
        'zeldovich_dens': zeldovich_dens,
        'rms_ratio': rms_ratio,
        'status': status
    }


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Compare Hermitian 3D and Zeldovich-PLT outputs'
    )
    parser.add_argument('--z', type=int, default=0, 
                       help='Z-slab to compare (default: 0)')
    parser.add_argument('--N', type=int, required=True,
                       help='Grid size (ppd)')
    parser.add_argument('--hermitian-dir', type=str, default='./',
                       help='Directory with Hermitian 3D output')
    parser.add_argument('--zeldovich-dir', type=str, default='./',
                       help='Directory with Zeldovich-PLT output')
    parser.add_argument('--narray', type=int, default=4,
                       help='Number of arrays in Hermitian output (default: 4)')
    
    args = parser.parse_args()
    
    try:
        result = compare_density_fields(
            z=args.z,
            N=args.N,
            hermitian_dir=args.hermitian_dir,
            zeldovich_dir=args.zeldovich_dir,
            narray=args.narray
        )
        
        print(f"\n{'='*60}")
        print(f"Comparison complete!")
        print(f"{'='*60}\n")
        
        # Exit with status code based on result
        if result['status'] == 'MATCH':
            sys.exit(0)
        else:
            sys.exit(1)
            
    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(2)


if __name__ == '__main__':
    main()

