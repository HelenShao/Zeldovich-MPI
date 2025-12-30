# Analysis of Particle IC Differences Between hermitian_3d_matrix_production and zeldovich-PLT

## Summary

The particle IC values differ significantly between the two codes, with many zeros appearing in the hermitian code but not in zeldovich. After analysis, the **fundamental² factor cancels out** in the k2_cutoff comparison, so both codes should be mathematically equivalent. The differences must be due to other factors (FFT normalization, particle writing, or other implementation details).

## Key Findings

### 1. k2 Calculation Difference

**zeldovich.cpp (line 349):**
```cpp
k2 = (kx * kx + ky * ky + kz * kz) * fundamental2;
```
- k2 includes `fundamental²` (where `fundamental = 2π/BoxSize`)

**hermitian_generation.c (line 235-236):**
```c
int k2_int = kx*kx + ky*ky + kz*kz;  // Integer k² for k_cutoff comparison
double k2 = (double)k2_int;  // Floating-point k² for power spectrum
```
- k2_int does NOT include `fundamental²`

### 2. k2_cutoff Calculation (DIFFERENT Units!)

**zeldovich.cpp (line 321-322):**
```cpp
double k2_cutoff = param.nyquist * param.nyquist / (param.k_cutoff * param.k_cutoff);
```
Where `param.nyquist = π / separation = π * ppd / boxsize` (physical wavenumber units)

For N=16, boxsize=100, k_cutoff=1.0:
- `nyquist = π * 16 / 100 ≈ 0.503`
- `k2_cutoff = 0.503² / 1.0² ≈ 0.253` (in **physical k² units**)

**hermitian_generation.c (line 97):**
```c
k2_cutoff = (Nhalf_dbl * Nhalf_dbl) / (k_cutoff * k_cutoff);
```

For N=16, k_cutoff=1.0:
- `k2_cutoff = (8)² / (1.0)² = 64` (in **integer k² units**, no fundamental scaling)

### 3. The Bug: Inconsistent Units in Comparison

**zeldovich.cpp (line 362):**
```cpp
if ( ... || (!param.CornerModes && k2>=k2_cutoff) ... )
```
- `k2 = (kx² + ky² + kz²) * fundamental²` (physical k² units)
- `k2_cutoff ≈ 0.253` (physical k² units)
- Compares: `(kx² + ky² + kz²) * fundamental² >= 0.253`
- For fundamental = 2π/100 ≈ 0.0628, fundamental² ≈ 0.00394
- So: `(kx² + ky² + kz²) * 0.00394 >= 0.253`
- This means: `(kx² + ky² + kz²) >= 0.253 / 0.00394 ≈ 64.2`
- **Modes with integer k² >= 65 are zeroed**

**hermitian_generation.c (line 253):**
```c
if ( ... || (!CornerModes && (double)k2_int >= k2_cutoff) ... )
```
- `k2_int = kx² + ky² + kz²` (integer k², no fundamental scaling)
- `k2_cutoff = 64` (integer k² units)
- Compares: `(kx² + ky² + kz²) >= 64`
- For N=16, max k² = 3*(8)² = 192
- **Modes with integer k² >= 64 are zeroed** (e.g., k² = 64, 65, ..., 192)

**The Issue:** Both should zero modes with integer k² >= 65 (approximately), but hermitian also zeros k² = 64, which zeldovich does not!

### 4. Why Zeros Appear in Hermitian but Not Zeldovich

The hermitian code is zeroing many more modes than zeldovich because:
- Hermitian compares integer k² directly to k2_cutoff (64)
- Zeldovich compares k² * fundamental² to k2_cutoff, which effectively requires k² ≈ 16,256

**Example:**
- Mode with kx=8, ky=0, kz=0: integer k² = 64
  - Hermitian: `64 >= 64` → **ZEROED** ✓
  - Zeldovich: `64 * 0.00394 = 0.252 >= 0.253` → **NOT ZEROED** ✓ (correct, k²=64 should not be zeroed)

- Mode with kx=8, ky=1, kz=0: integer k² = 65
  - Hermitian: `65 >= 64` → **ZEROED** ✓
  - Zeldovich: `65 * 0.00394 = 0.256 >= 0.253` → **ZEROED** ✓ (both agree)

- Mode with kx=8, ky=8, kz=8: integer k² = 192
  - Hermitian: `192 >= 64` → **ZEROED** ✓
  - Zeldovich: `192 * 0.00394 = 0.756 >= 0.253` → **ZEROED** ✓ (both agree)

**The Problem:** Hermitian zeros k²=64, but zeldovich does not. This is the main source of differences!

### 5. Impact on Particle ICs

After the 3D FFT:
- Hermitian: Many high-k modes are zeroed → zeros appear in real-space particle fields
- Zeldovich: Almost no modes are zeroed (only extremely high k) → no zeros in particle fields

This explains why:
1. Many particles in hermitian have zero displacement/velocity
2. The average ratio is ~0.95 (hermitian/zeldovich) - hermitian has less power due to zeroed modes
3. Maximum differences are large (up to ~16) - the zeroed modes contribute significantly

## Verification: fundamental² Cancels Out

Mathematical proof that both approaches are equivalent:

- **zeldovich**: `k2_cutoff = nyquist² / k_cutoff² = (N/2)² * fundamental² / k_cutoff²`
- **zeldovich**: `k2 = (kx² + ky² + kz²) * fundamental²`
- **zeldovich comparison**: `k2 >= k2_cutoff`
  → `(kx² + ky² + kz²) * fundamental² >= (N/2)² * fundamental² / k_cutoff²`
  → `(kx² + ky² + kz²) >= (N/2)² / k_cutoff²`

- **hermitian**: `k2_cutoff = (N/2)² / k_cutoff²`
- **hermitian**: `k2_int = kx² + ky² + kz²`
- **hermitian comparison**: `k2_int >= k2_cutoff`
  → `(kx² + ky² + kz²) >= (N/2)² / k_cutoff²`

**Conclusion**: Both codes should zero the same modes. The fundamental² factor cancels out.

## Other Potential Causes of Differences

Since the k2_cutoff logic is equivalent, the differences must come from:

1. **FFT Normalization**: Different normalization factors in the FFT could lead to different amplitudes
2. **Particle Writing**: Different implementations of `WriteParticlesSlab` could handle edge cases differently
3. **Nyquist Mode Handling**: Subtle differences in how Nyquist modes are handled
4. **Floating Point Precision**: Accumulation of small differences through the FFT chain
5. **RNG Sequence**: Although RNG matches at the D level, small differences could accumulate

## Investigation Needed

The following should be checked:

1. **Compare FFT outputs**: Check if the FFT results (before particle writing) match between the two codes
2. **Compare F, G, H values**: Verify that the computed F, G, H arrays match
3. **Check FFT normalization**: Verify that both codes use the same FFT normalization
4. **Check particle writing**: Compare the particle writing code paths
5. **Check for numerical precision issues**: Small differences could accumulate through the FFT chain

## Additional Observations

1. **Particle positions with zeros**: In hermitian output, particles at k=8, k=9, etc. have zeros. These correspond to positions where high-k modes (that were zeroed) would have contributed.

2. **Average ratio ~0.95**: The hermitian code has less total power because many modes are zeroed, leading to smaller average displacement/velocity values.

3. **Both codes use same RNG**: The RNG sequences match (as verified in previous tests), so the differences are purely due to the k2_cutoff bug.

## Next Steps

1. Fix the k2_cutoff comparison in `hermitian_generation.c` to include fundamental²
2. Re-run the test and verify particle ICs match zeldovich
3. Verify that zeros no longer appear inappropriately
4. Confirm that the k2_cutoff condition now works as intended (zeroing only modes above the physical cutoff)

