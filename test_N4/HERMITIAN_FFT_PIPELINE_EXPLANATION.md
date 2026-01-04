# Hermitian Code FFT Pipeline: 2D FFT → 1D FFT

## Overview

The Hermitian code decomposes the 3D FFT as: **2D FFT first, then 1D FFT**. This document explains how this works in terms of the XY slices.

## Data Structure: Y-Slices

### What is a Y-Slice?

A **Y-slice** is a 2D plane at a fixed Y-coordinate:
- **Format**: `[Array][Z][X]` - narray arrays, N Z-values, N X-values
- **Memory order**: X is fastest varying (stride-1), then Z (stride N), then array_idx (stride N²)
- **Each Y-slice** represents all (X, Z) values at a specific Y
- **Total Y-slices**: N/2 + 1 pairs (primary + conjugate for Hermitian symmetry)

### Example for N=4:
- Y-slices: Y=0, Y=1, Y=2 (Y=3 is conjugate of Y=1, Y=4 is conjugate of Y=0)
- Each slice: 4 arrays × 4 X-values × 4 Z-values = 64 complex values per slice

## Stage 1: Generate Y-Slices + 2D FFT

### Step 1.1: Generate Y-Slices in Fourier Space

**Location**: `generate_hermitian_slice_pair_local()` in `hermitian_generation.c`

**Process**:
1. For each Y-slice (e.g., Y=0, Y=1, Y=2):
   - Generate D, F, G, H values in **Fourier space** (k-space) with power spectrum weighting
   - Store in `local_y_slices` buffer: `[Array][Z][X]` format (X is fastest varying)
   - **Initial state**: Data is in **Fourier space** in (kx, kz) dimensions
   - **Y coordinate maps to ky**: Each Y-slice at coordinate Y corresponds to Fourier mode ky
     - ky = (Y > N/2) ? Y - N : Y
     - Example: Y=0 → ky=0, Y=1 → ky=1, Y=N/2 → ky=N/2, Y=N/2+1 → ky=-N/2+1
   - So the data is: `f(kx, ky, kz)` where ky is determined by the Y coordinate

**Memory Layout**:
```
local_y_slices[slice_idx][array_idx][z][x]
  Memory order: [Slice][Array][Z][X] where X is stride-1 (fastest varying)
  - slice_idx: 0 = primary, 1 = conjugate
  - array_idx: 0-3 (D, F, G, H)
  - z: 0 to N-1 (Z coordinate, stride N)
  - x: 0 to N-1 (X coordinate, stride 1 - fastest varying)
  
  Index formula: x + N * (z + N * (array_idx + narray * slice_idx))
```

### Step 1.2: Apply 2D FFT (X, Z axes)

**Location**: `hermitian_generation.c` lines 1866, 1896

**Process**:
```c
FFTW_EXECUTE_DFT(plan_2d, prim_array_start, prim_array_start);
```

**What happens**:
- **Input**: Y-slice in Fourier space `f1(kx, ky, kz)` 
  - Generated in Fourier space (k-space) for kx, ky, and kz
  - ky is determined by the Y coordinate: ky = (Y > N/2) ? Y - N : Y
  - Data is organized by Y coordinate (real-space index), but represents Fourier mode ky
- **2D FFT (inverse FFT) along X and Z axes**: Transforms (kx, kz) from Fourier space to real space
  - `f1(kx, ky, kz)` → `f2(x, ky, z)`
  - ky remains in Fourier space (it's determined by which Y-slice we're on)
- **Output**: Y-slice in **real space** for X and Z dimensions
  - **ky is still in Fourier space** (determined by Y coordinate)
  - Y is a real-space coordinate index, but it represents Fourier mode ky

**After 2D FFT**:
- Data is now: `f2(x, ky, z)` - **real space in (x, z), ky is in Fourier space**
- This is stored in `local_y_slices` buffer
- Each Y-slice is now a 2D real-space plane at fixed Y coordinate (which corresponds to Fourier mode ky)

**Key Point**: After the 2D FFT, we have:
- **Real space**: X and Z dimensions (transformed from kx, kz)
- **Fourier space**: ky dimension (determined by Y coordinate: ky = (Y > N/2) ? Y - N : Y)
- Y is a real-space coordinate index, but it already represents Fourier mode ky
- This is the **intermediate state** before the 1D FFT along Y

## Stage 2: MPI Communication (All-to-All)

**Location**: `main.cpp` around line 1292

**Process**:
- Y-slices are packed and sent via `MPI_Ialltoallv`
- Data is reorganized from Y-distributed to X-distributed
- Each rank receives **pencils** (Y-lines at fixed X, Z) for its X-range

**After Communication**:
- Data arrives in `recv_buffer` organized by source rank
- Format: `[src][array][batch_slice][pencil]`
- Each pencil is a line along Y at fixed (X, Z)
  - `pencil_idx = x_idx * z_count + z_idx`
  - Contains all Y values (0 to N-1) for that (X, Z) point

## Stage 3: Unpack + 1D FFT

### Step 3.1: Unpack into Z-Slabs

**Location**: `z_streaming_unpack()` in `z_streaming.c` lines 60-135

**Process**:
- Extract one Z-slab at a time from `recv_buffer`
- Z-slab format: `[Array][X][Y]` - narray arrays, x_count X-values, N Y-values
- Store in `local_z_slab` buffer

**Memory Layout**:
```
local_z_slab[array_idx][x_idx][y]
  - array_idx: 0-3 (D, F, G, H)
  - x_idx: local X index (0 to x_count-1)
  - y: 0 to N-1 (Y coordinate)
```

**After Unpacking**:
- Data is in `local_z_slab` in `[Array][X][Y]` format
- **Real space**: X and Z dimensions (from 2D FFT)
- **Real-space coordinate**: Y dimension (will be transformed to ky in next step)
- This is the same intermediate state, just reorganized

### Step 3.2: Apply 1D FFT (Y-axis)

**Location**: `z_streaming.c` lines 215-223

**Process**:
```c
#pragma omp parallel for collapse(2)
for (int array_idx = 0; array_idx < narray; array_idx++) {
    for (int x_idx = 0; x_idx < x_count; x_idx++) {
        fftw_complex_t *y_data = &ZSLAB(array_idx, x_idx, 0, N, narray, x_count);
        FFTW_EXECUTE_DFT(plan_1d_y, y_data, y_data);
    }
}
```

**What happens**:
- For each (Array, X) pair, apply 1D FFT along Y-direction
- **Input**: `f2(x, ky, z)` - real space in (x, z), ky is in Fourier space (from Y coordinate)
- **1D FFT (inverse FFT) along Y-axis**: Transforms ky from Fourier space to real space
  - The data is organized by Y coordinate (which maps to ky)
  - The 1D FFT transforms ky to y in real space
  - `f2(x, ky, z)` → `F(x, y, z)`
- **Output**: Full 3D real space `F(x, y, z)`

**After 1D FFT**:
- Data is now in **full 3D real space**
- Format: `[Array][X][Y]` where X and Y are now in real space
- This is the **final result** after complete 3D FFT: `F(x, y, z)`

## Complete Pipeline Summary

```
1. Generate Y-slices in Fourier space
   Input:  f(kx, ky, kz)  [Generated in Fourier space, ky determined by Y coordinate]
   
2. Apply 2D FFT (inverse FFT along X, Z axes)
   Transform: Fourier space → Real space for (kx,kz) → (x,z)
   Output: f2(x, ky, z)   [Real space in x,z, ky is in Fourier space]
   ↓
   [MPI All-to-All: Reorganize from Y-distributed to X-distributed]
   ↓
3. Unpack into Z-slabs
   Format: [Array][X][Y] - same data, different organization
   
4. Apply 1D FFT (Y-axis)
   Output: F(x, y, z) [Full 3D real space]
```

## Key Insights

1. **Y-Slices are 2D planes**: Each Y-slice is a 2D array in (X, Z) at fixed Y
2. **2D FFT transforms (X, Z)**: Converts from Fourier space to real space in X and Z
3. **Data reorganization**: MPI all-to-all redistributes from Y-slices to Z-slabs
4. **1D FFT transforms Y**: Converts ky from Fourier space to real space (y)
5. **Final result**: Full 3D real space `F(x, y, z)`

## Comparison with Zeldovich-PLT

### Zeldovich: 1D FFT → 2D FFT
```
1. Generate data in Fourier space: f(kx, ky, kz)
2. 1D FFT (inverse) along Z: f(kx, ky, kz) → f1(kx, ky, z) [Transform kz to z in real space]
3. 2D FFT (inverse) along Y, X: f1(kx, ky, z) → F(x, y, z) [Transform (kx, ky) to (x, y) in real space]
```

### Hermitian: 2D FFT → 1D FFT
```
1. Generate Y-slices in Fourier space: f1(kx, ky, kz) [Fourier in kx,ky,kz, ky determined by Y coordinate]
2. 2D FFT (inverse) along X, Z: f1(kx, ky, kz) → f2(x, ky, z) [Transform (kx,kz) to (x,z) real space, ky remains in Fourier space]
3. [MPI reorganization]
4. 1D FFT (inverse) along Y: f2(x, ky, z) → F(x, y, z) [Transform ky to y in real space, complete 3D FFT]
```

**Critical Comparison Point**:
- After Zeldovich's 2D FFT = After Hermitian's 1D FFT
- Both should produce `F(x, y, z)` - the complete 3D FFT result in real space

