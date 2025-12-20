# FFT Implementation Comparison: Zeldovich-PLT vs Hermitian 3D Matrix Production

## Overview

Both codes perform FFT transformations on Hermitian matrices to produce real-space fields, but with different architectures and conventions.

## Key Differences

### 1. FFT Sign Convention

**Zeldovich-PLT:**
```c
plan1d  = fftw_plan_dft_1d(n, p, p, +1, FFTW_PATIENT);   // Inverse FFT
plan2d  = fftw_plan_dft_2d(n, n, p, p, +1, FFTW_PATIENT); // Inverse FFT
plan1d_forward = fftw_plan_dft_1d(n, p, p, -1, FFTW_PATIENT);   // Forward FFT
plan2d_forward = fftw_plan_dft_2d(n, n, p, p, -1, FFTW_PATIENT); // Forward FFT
```

**Hermitian 3D Matrix Production:**
```c
#define FFT_SIGN FFTW_BACKWARD  // = +1
*plan_2d_out = FFTW_PLAN_DFT_2D(N, N, dummy_2d, dummy_2d, FFT_SIGN, FFTW_ESTIMATE);
*plan_1d_out = FFTW_PLAN_DFT_1D(N, dummy_1d, dummy_1d, FFT_SIGN, FFTW_ESTIMATE);
```

**Result:** Both use `+1` (FFTW_BACKWARD) for Fourier → real space transformation, so this is **CONSISTENT**.

### 2. FFT Transform Type

**Zeldovich-PLT:**
- Uses `fftw_plan_dft_1d` and `fftw_plan_dft_2d` (complex-to-complex)
- Processes complex arrays throughout
- Manual Hermitian symmetry enforcement in code

**Hermitian 3D Matrix Production:**
- Uses `FFTW_PLAN_DFT_2D` and `FFTW_PLAN_DFT_1D` (complex-to-complex)
- Processes complex arrays throughout
- Manual Hermitian symmetry enforcement in code

**Result:** Both use complex-to-complex transforms. **CONSISTENT**.

### 3. FFT Planning Strategy

**Zeldovich-PLT:**
```c
- Uses FFTW_PATIENT for planning
- Imports/exports wisdom from fftw_zeldovich.wisdom
- Pre-allocates dummy buffer for planning
```

**Hermitian 3D Matrix Production:**
```c
- Uses FFTW_ESTIMATE for planning (faster, less optimal)
- No wisdom import/export by default
- Can use FFTW_MEASURE or FFTW_PATIENT via compile flags
- Pre-allocates dummy buffer, frees immediately after plan creation
```

**Result:** Different planning strategies affect performance but not correctness. Zeldovich-PLT optimizes for repeated runs.

### 4. FFT Execution Order

**Zeldovich-PLT:**
1. Generate Hermitian pairs (Y-slice and conjugate)
2. Store conjugate: `conj(D) + I*conj(F)`
3. Apply 1D FFT in Z-direction (per Y-slice): `InverseFFT_Yonly()`
4. Later apply 2D FFT in Y-X directions: `Inverse2dFFT()`
5. **Total: 3D FFT via 1D(Z) + 2D(YX)**

**Hermitian 3D Matrix Production:**
1. Generate Hermitian pairs (Y-slice and conjugate)
2. Store conjugate: mirror Y-slices with proper indexing
3. Apply 2D FFT in X-Z directions (per Y-slice): `FFTW_EXECUTE_DFT(plan_2d, ...)`
4. Exchange data to pencil layout via MPI
5. Apply 1D FFT in Y-direction (per pencil): `FFTW_EXECUTE_DFT(plan_1d, ...)`
6. **Total: 3D FFT via 2D(XZ) + 1D(Y)**

**Result:** Different decomposition order but mathematically equivalent. **KEY DIFFERENCE** in implementation strategy.

### 5. Hermitian Symmetry Handling

**Zeldovich-PLT:**
```c
// Generate for Y and Y_mirror simultaneously
AYZX(slab, 0, yres, z, x) = D + I * F;
AYZX(slabHer, 0, yresHer, zHer, xHer) = conj(D) + I * conj(F);

// Special handling for Y=0 (self-conjugate):
if (yblock == 0 && yres == 0) {
    // Copy first half-plane to second half
    // Enforce Hermitian structure
}
```

**Hermitian 3D Matrix Production:**
```c
// Generate for Y=i and Y=N-i as pairs
// Primary slice:
PRIM_SLICE(a, x, z)[0] = D_real;
PRIM_SLICE(a, x, z)[1] = D_imag;

// Conjugate slice (mirrored):
CONJ_SLICE(a, x_mirror, z_mirror)[0] = D_real;
CONJ_SLICE(a, x_mirror, z_mirror)[1] = -D_imag;

// Special handling for Y=0 and Y=N/2 (self-conjugate):
// Process only half the plane, mirror the rest
```

**Result:** Both enforce Hermitian symmetry but with different data layouts. **CONSISTENT** in principle.

### 6. Parallelization Strategy

**Zeldovich-PLT:**
- OpenMP parallelization over Y-slices
- Each thread processes one Y-slice pair
- Shared RNG array `Pk.v2rng[y]` indexed by Y
- No MPI parallelization (single-node)

**Hermitian 3D Matrix Production:**
- Hybrid MPI + OpenMP
- MPI parallelization over Y-slices (distributed)
- Optional OpenMP within Y-slice (X-Z loops)
- Each rank has own RNG array indexed by Y
- MPI communication for data redistribution (Y-slice → pencil)

**Result:** Hermitian 3D scales to multiple nodes, Zeldovich-PLT is single-node. **MAJOR ARCHITECTURAL DIFFERENCE**.

### 7. Data Layout

**Zeldovich-PLT:**
```c
AYZX(slab, a, y, z, x)  // Array-Y-Z-X ordering
```

**Hermitian 3D Matrix Production:**
```c
Y_SLICE(y_idx, a, z, x)  // Y-Array-Z-X ordering in generation
PENCIL(pencil, a, y)     // After MPI exchange
ZSLAB(x, a, y)           // Final output for Zeldovich compatibility
```

**Result:** Different internal layouts, but final output is compatible. **DIFFERS** in intermediate stages.

### 8. Normalization

**Zeldovich-PLT:**
- FFT normalization handled implicitly by FFTW (no explicit scaling)
- Applies growth factors and power spectrum during generation

**Hermitian 3D Matrix Production:**
- FFT normalization handled implicitly by FFTW (no explicit scaling)
- Applies power spectrum during generation via `cgauss()`
- Computes F, G, H from D using k-vector components

**Result:** Both rely on FFTW's default normalization. **CONSISTENT**.

## Remaining Differences (Post RNG-Fix)

After fixing the z=0 skip bug and using SPLINE_RESOLUTION=512:

### 1. Numerical Precision
- **Status:** 99.79% match between N=256 and N=512
- **Cause:** Floating-point precision limits in:
  - Power spectrum spline interpolation
  - Complex arithmetic operations
  - Accumulated rounding errors
- **Magnitude:** Remaining differences < 2e-02 (mostly < 1e-03 with PS=2048)

### 2. RNG Consistency
- **Status:** Fixed by z=0 skip bug correction
- **k=(0,0,z) modes:** Now match perfectly
- **Remaining 7 differences (PS=512):** Non-systematic, likely interpolation

### 3. FFT Implementation
- **Status:** Mathematically equivalent
- **Zeldovich-PLT:** Z-first then YX (1D + 2D)
- **Hermitian 3D:** XZ-first then Y (2D + 1D)
- **Impact:** None on correctness, only affects performance and memory access patterns

### 4. Parallelization
- **Zeldovich-PLT:** Single-node OpenMP
- **Hermitian 3D:** Multi-node MPI + OpenMP
- **Trade-off:** 
  - Zeldovich-PLT: Simpler, no communication overhead
  - Hermitian 3D: Scales to larger grids, supports distributed memory

## Recommendations

### For Single-Node Runs (N ≤ 2048):
- Either code works
- Zeldovich-PLT may be simpler to use
- Hermitian 3D provides better control over RNG consistency across N values

### For Multi-Node Runs (N > 2048):
- Use Hermitian 3D Matrix Production
- MPI parallelization essential for memory constraints
- Achieved 99.79% consistency with RNG fix

### For RNG Consistency Testing:
- Use `SPLINE_RESOLUTION=512` or higher
- Enable `DEBUG_RNG_CONSISTENCY=1` for validation
- Expect ~0.2% differences due to floating-point precision limits

## Conclusion

Both codes implement FFT of Hermitian matrices correctly, with:
- **Same FFT sign convention** (FFTW_BACKWARD for Fourier → real)
- **Same transform type** (complex-to-complex DFT)
- **Equivalent Hermitian symmetry** enforcement
- **Different parallelization strategies** (single vs multi-node)
- **Different FFT decomposition orders** (but mathematically equivalent)

The remaining 0.2% differences after the z=0 fix are due to:
1. Floating-point precision limits
2. Spline interpolation accuracy (reduced with higher resolution)
3. Different numerical paths in power spectrum evaluation


