# Comparison: `just_density` Mode in Zeldovich-PLT vs Hermitian Code

## Zeldovich-PLT Behavior (`qdensity == 2`)

### 1. **Initialization** (line 306):
```cpp
int just_density = param.qdensity == 2;  // Don't generate displacements
```

### 2. **Array Count** (lines 878-882):
```cpp
int narray;
if (param.qdensity == 2) {
    narray = 1;  // Only one array (density)
} else {
    narray = param.qPLT ? 4 : 2;  // 4 arrays if PLT, 2 otherwise
}
```
- **Only Array 0 is allocated** (narray = 1)
- Arrays 1, 2, 3 are not allocated

### 3. **F, G, H Computation** (lines 409-444):
**IMPORTANT**: F, G, H are **still computed** even in density-only mode!
```cpp
if (D != 0.) {
    eigenmode e = get_eigenmode(kx, ky, kz, ppd, param.qPLT);
    // ... compute rescale, f ...
    F = rescale * I * e.vec[0] * param.fundamental * ik2 * D;
    G = rescale * I * e.vec[1] * param.fundamental * ik2 * D;
    H = rescale * I * e.vec[2] * param.fundamental * ik2 * D;
} else {
    F = G = H = 0.0;
    f = 0.;
}
```
- F, G, H are computed but **not stored**
- This is different from the Hermitian code!

### 4. **Storage** (lines 473-476):
```cpp
if (!just_density) {
    // Store D+iF, G+iH, velocities (Arrays 0, 1, 2, 3)
    AYZX(slab, 0, yres, z, x) = D + I * F;
    AYZX(slab, 1, yres, z, x) = G + I * H;
    // ... velocities ...
} else {
    // Density-only: Only store D (Array 0)
    AYZX(slab, 0, yres, z, x) = D;
    AYZX(slabHer, 0, yresHer, zHer, xHer) = conj(D);
}
```
- **Only Array 0 stores D** (density field)
- Conjugate slice stores `conj(D)`
- Arrays 1, 2, 3 are not used (not even allocated)

### 5. **FFT Operations**:
- 1D FFT along Z: Applied to Array 0 only (line 515)
- 2D FFT along Y, X: Applied to Array 0 only (line 662)
- Only the density field is transformed

### 6. **Output** (line 1004):
```cpp
if (param.qdensity != 2) {
    // Print displacement statistics
    fmt::print(stderr, "The maximum component-wise displacements are ...");
}
```
- Displacement statistics are not printed in density-only mode

## Hermitian Code Behavior (`just_density = 1`)

### Key Differences:

1. **F, G, H Computation**:
   - **Hermitian**: F, G, H are **set to zero** and **not computed** (lines 399-401)
   - **Zeldovich**: F, G, H are **computed** but **not stored** (lines 438-440)

2. **Memory Allocation**:
   - **Hermitian**: All 4 arrays are allocated, but only Array 0 is used
   - **Zeldovich**: Only Array 0 is allocated (narray = 1)

3. **Storage**:
   - **Hermitian**: Only Array 0 stores D; Arrays 1-3 remain zero/unused
   - **Zeldovich**: Only Array 0 stores D; Arrays 1-3 don't exist

## Summary

**Zeldovich-PLT (`just_density = 1`):**
- Computes F, G, H from D (but doesn't store them)
- Allocates only Array 0 (narray = 1)
- Stores only D in Array 0
- Applies FFTs only to Array 0
- More memory-efficient (only 1 array instead of 4)

**Hermitian Code (`just_density = 1`):**
- Sets F, G, H to zero (doesn't compute them)
- Allocates all 4 arrays (but only uses Array 0)
- Stores only D in Array 0
- Applies FFTs only to Array 0
- Less memory-efficient (allocates unused arrays)

## Why the Difference?

The Zeldovich code computes F, G, H even in density-only mode because:
- The computation is relatively cheap (just a few multiplications)
- It maintains code consistency (same computation path)
- The results are simply discarded, not stored

The Hermitian code skips the computation entirely for efficiency, but still allocates all arrays (which could be optimized).

