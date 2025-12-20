# Field Extraction Guide: Hermitian 3D vs Zeldovich-PLT

## How Fields are Packed in Arrays

### Zeldovich-PLT Array Layout (from zeldovich.cpp lines 453-458)

```cpp
AYZX(slab, 0, yres, z, x) = D + I * F;      // Array 0
AYZX(slab, 1, yres, z, x) = G + I * H;      // Array 1
AYZX(slab, 2, yres, z, x) = 0. + I * F * f; // Array 2 (PLT only)
AYZX(slab, 3, yres, z, x) = G * f + I * H * f; // Array 3 (PLT only)
```

Where:
- `D` = density field (real-valued in real space)
- `F` = x-displacement
- `G` = y-displacement
- `H` = z-displacement
- `f` = growth rate factor (PLT velocity correction)

### Field Mapping

| Field | Array | Component | Extraction |
|-------|-------|-----------|------------|
| **Density** | 0 | Real part | `real(array[0])` |
| **X-displacement** | 0 | Imag part | `imag(array[0])` |
| **Y-displacement** | 1 | Real part | `real(array[1])` |
| **Z-displacement** | 1 | Imag part | `imag(array[1])` |
| **X-velocity** | 2 | Imag part | `imag(array[2])` (PLT only) |
| **Y-velocity** | 3 | Real part | `real(array[3])` (PLT only) |
| **Z-velocity** | 3 | Imag part | `imag(array[3])` (PLT only) |

### WriteParticlesSlab Extraction (output.cpp lines 93-97)

```cpp
dens   = real(YX(slab1, y, x))  // = real(array[0]) = D
pos[0] = imag(YX(slab1, y, x))  // = imag(array[0]) = F  (x-displ)
pos[1] = real(YX(slab2, y, x))  // = real(array[1]) = G  (y-displ)
pos[2] = imag(YX(slab2, y, x))  // = imag(array[1]) = H  (z-displ)
```

## Your Hermitian 3D Output Format

### File Structure
- **Files per Z-slab**: `hermitian_3d_z{z}_rank{rank}.bin`
- **Layout**: `[Array][Y][X]` for each Z-slab
- **Data type**: `fftw_complex_t` (2 doubles per element)
- **Arrays**: `narray = 4` (or 2 if not doing velocities)

### Reading Your Binary Files

```python
import numpy as np

def read_hermitian_zslab(z, rank, N, narray, x_count):
    """
    Read one Z-slab from hermitian 3D output
    
    Parameters:
    -----------
    z : int
        Z-index of slab
    rank : int
        MPI rank that wrote this file
    N : int
        Grid size (ppd)
    narray : int
        Number of arrays (typically 4)
    x_count : int
        Number of X-values owned by this rank
    
    Returns:
    --------
    data : complex array of shape (narray, N, x_count)
    """
    filename = f"hermitian_3d_z{z}_rank{rank}.bin"
    
    # Read binary complex data (each element: 2 doubles)
    data = np.fromfile(filename, dtype=np.complex128)
    
    # Reshape to [Array][Y][X]
    data = data.reshape(narray, N, x_count)
    
    return data

# Usage:
z = 0
rank = 0
N = 256
narray = 4
x_count = 64  # Depends on decomposition

data = read_hermitian_zslab(z, rank, N, narray, x_count)

# Extract fields:
density       = data[0, :, :].real  # Array 0, real part
x_displacement = data[0, :, :].imag  # Array 0, imag part
y_displacement = data[1, :, :].real  # Array 1, real part
z_displacement = data[1, :, :].imag  # Array 1, imag part
```

## Important: After 3D FFT - Real Fields in Complex Storage

**After the complete 3D FFT (Fourier → Real space):**
- All fields are **mathematically real-valued** (not complex)
- But they are **stored** in complex arrays for memory efficiency

**Key distinction:**
- **Storage format:** `complex128` (2 doubles per element)
- **Mathematical content:** Two real fields packed together

**Example:**
```python
# Array[0] stores TWO real fields:
density       = data[0, :, :].real  # Real field #1 (float64)
x_displacement = data[0, :, :].imag  # Real field #2 (float64)

# Both are real-valued! The .imag is just the storage location.
assert density.dtype == np.float64
assert x_displacement.dtype == np.float64
```

**Why this works:**
- Hermitian symmetry in Fourier space guarantees real output
- Complex array holds two real fields: one in `.real`, one in `.imag`
- Both `.real` and `.imag` contain real numbers (not complex!)

**See `REAL_VS_COMPLEX_PACKING.md` for detailed explanation.**

### Example: Extract density field for z=0

```python
import numpy as np
import glob

def extract_density_full_slab(z, N, output_dir="./"):
    """
    Extract full density field for one Z-slab from all ranks
    
    Returns:
    --------
    density : array of shape (N, N)
        Full density field for this Z-slab
    """
    # Find all files for this Z
    files = sorted(glob.glob(f"{output_dir}/hermitian_3d_z{z}_rank*.bin"))
    
    density_slices = []
    
    for filename in files:
        # Parse rank number from filename
        rank = int(filename.split('rank')[1].split('.')[0])
        
        # Determine x_count for this rank (may vary)
        # For now, read entire file and infer shape
        data = np.fromfile(filename, dtype=np.complex128)
        
        # Assuming narray=4
        narray = 4
        x_count = len(data) // (narray * N)
        
        # Reshape
        data = data.reshape(narray, N, x_count)
        
        # Extract density (array 0, real part)
        density_slice = data[0, :, :].real
        
        density_slices.append(density_slice)
    
    # Concatenate along X dimension
    density_full = np.concatenate(density_slices, axis=1)
    
    return density_full

# Usage:
density_z0 = extract_density_full_slab(z=0, N=256)
print(f"Density field shape: {density_z0.shape}")  # Should be (256, 256)
print(f"Density RMS: {np.std(density_z0)}")
```

## Reading Zeldovich-PLT Particle Files

```python
import numpy as np
import struct

def read_zeldovich_slab(z, ppd, output_dir="./", format="RVdoubleZel"):
    """
    Read Zeldovich-PLT particle file
    
    Parameters:
    -----------
    z : int
        Z-index of slab
    ppd : int
        Particles per dimension
    output_dir : str
        Output directory
    format : str
        Particle format (RVdoubleZel, RVZel, Zeldovich, ZelSimple)
    
    Returns:
    --------
    density : array of shape (ppd, ppd)
    displacement : array of shape (ppd, ppd, 3)
    velocity : array of shape (ppd, ppd, 3) or None
    """
    filename = f"{output_dir}/ic_{z}"
    
    if format == "RVdoubleZel":
        # Struct: unsigned short i,j,k (6 bytes) + double displ[3] (24 bytes) + double vel[3] (24 bytes)
        dtype = np.dtype([
            ('i', np.uint16),
            ('j', np.uint16), 
            ('k', np.uint16),
            ('displ', np.float64, 3),
            ('vel', np.float64, 3)
        ])
    elif format == "RVZel":
        # Struct: unsigned short i,j,k (6 bytes) + float displ[3] (12 bytes) + float vel[3] (12 bytes)
        dtype = np.dtype([
            ('i', np.uint16),
            ('j', np.uint16),
            ('k', np.uint16),
            ('displ', np.float32, 3),
            ('vel', np.float32, 3)
        ])
    else:
        raise NotImplementedError(f"Format {format} not implemented yet")
    
    # Read all particles
    particles = np.fromfile(filename, dtype=dtype)
    
    # Reshape to (ppd, ppd)
    particles = particles.reshape(ppd, ppd)
    
    # Extract fields
    # Note: In Zeldovich output, the density field is NOT written directly
    # You need to compute it from the displacement field divergence
    # OR read from separate density file if qdensity=1
    
    displacement = particles['displ']  # Shape: (ppd, ppd, 3)
    velocity = particles['vel']        # Shape: (ppd, ppd, 3)
    
    return displacement, velocity

# Usage:
displ_z0, vel_z0 = read_zeldovich_slab(z=0, ppd=256, format="RVdoubleZel")
print(f"Displacement field shape: {displ_z0.shape}")  # (256, 256, 3)
```

## Complete Comparison Script

```python
import numpy as np
import glob

def compare_density_fields(z, N, hermitian_dir="./", zeldovich_dir="./"):
    """
    Compare density fields from both codes for one Z-slab
    """
    # Read hermitian 3D density
    hermitian_dens = extract_density_full_slab(z, N, hermitian_dir)
    
    # Read zeldovich density (from separate density file)
    # Note: Zeldovich writes density to separate file if qdensity=1
    # Format: binary float array [Z][Y][X]
    dens_filename = f"{zeldovich_dir}/zeldovich.dens_{N}"
    if os.path.exists(dens_filename):
        # Read full density cube
        dens_data = np.fromfile(dens_filename, dtype=np.float32)
        dens_data = dens_data.reshape(N, N, N)
        zeldovich_dens = dens_data[z, :, :]
    else:
        print(f"Warning: {dens_filename} not found")
        print("Zeldovich must be run with qdensity=1 to output density field")
        return None
    
    # Compare
    print(f"\n=== Comparison for Z={z} ===")
    print(f"Hermitian 3D:")
    print(f"  Shape: {hermitian_dens.shape}")
    print(f"  Mean: {np.mean(hermitian_dens):.6e}")
    print(f"  RMS: {np.std(hermitian_dens):.6e}")
    print(f"  Min/Max: {np.min(hermitian_dens):.6e} / {np.max(hermitian_dens):.6e}")
    
    print(f"\nZeldovich-PLT:")
    print(f"  Shape: {zeldovich_dens.shape}")
    print(f"  Mean: {np.mean(zeldovich_dens):.6e}")
    print(f"  RMS: {np.std(zeldovich_dens):.6e}")
    print(f"  Min/Max: {np.min(zeldovich_dens):.6e} / {np.max(zeldovich_dens):.6e}")
    
    # Compute ratio
    rms_ratio = np.std(hermitian_dens) / np.std(zeldovich_dens)
    print(f"\n=== Normalization Test ===")
    print(f"RMS ratio (Hermitian/Zeldovich): {rms_ratio:.6e}")
    
    if abs(rms_ratio - 1.0) < 0.01:
        print("✓ Normalization matches!")
    elif abs(rms_ratio - N**3) < 0.01:
        print(f"✗ Missing N³ normalization (N³ = {N**3})")
    elif abs(rms_ratio - 1.0/N**3) < 0.01:
        print(f"✗ Extra N³ normalization (1/N³ = {1.0/N**3})")
    else:
        print(f"? Unexpected ratio: {rms_ratio}")
    
    return hermitian_dens, zeldovich_dens, rms_ratio

# Run comparison
hermitian_dens, zeldovich_dens, ratio = compare_density_fields(
    z=0, 
    N=256,
    hermitian_dir="/path/to/hermitian/output",
    zeldovich_dir="/path/to/zeldovich/output"
)
```

## Summary

**To extract density for z=0 from your output:**
```python
data = read_hermitian_zslab(z=0, rank=0, N=256, narray=4, x_count=64)
density = data[0, :, :].real  # Array 0, real part
```

**Key points:**
1. Density is in **array[0], real part**
2. After 3D FFT, use `.real` even though fields are real-valued
3. May need to combine multiple ranks to get full slab
4. Zeldovich density comes from separate file (qdensity=1) or computed from displacement divergence

