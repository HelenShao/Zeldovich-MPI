# F, G, H Value Discrepancy Diagnosis

## Problem (RESOLVED ✓)
D values matched perfectly between hermitian and zeldovich codes, but F, G, H values differed by a factor that varied with kz:
- For small kz (≤8): ratio ≈ 8.0 (constant)
- For large kz (>8): ratio increased with kz (9.0, 10.0, 20.0, etc.)

## Example Coordinate (Before Fix)
- Y=0, (x,z)=(0,1), k=(0,0,1), k2=1.0
- D (both codes): (-6.8841943804e-02, 2.5901576221e-01) ✓ MATCH
- H (hermitian): (-5.1527865809e-01, -1.3695222299e-01)
- H (zeldovich): (-4.1222292647e+00, -1.0956177839e+00)
- Ratio: 8.0 (exactly)

## Formulas

### Zeldovich (zeldovich.cpp lines 451-453):
```cpp
F = rescale * I * e.vec[0] * param.fundamental * ik2 * D;
G = rescale * I * e.vec[1] * param.fundamental * ik2 * D;
H = rescale * I * e.vec[2] * param.fundamental * ik2 * D;
```
Where:
- `k2 = (kx² + ky² + kz²) * fundamental²` (line 350)
- `ik2 = 1.0 / k2` (line 388)
- `fundamental = 2π / BoxSize = 2π / 100.0 ≈ 0.0628`

So:
- `F = rescale * I * e.vec[0] * fundamental / (k2_index * fundamental²) * D`
- `F = rescale * I * e.vec[0] / (k2_index * fundamental) * D`

### Hermitian (hermitian_generation.c lines 502, 518-525):
```c
double factor = rescale / (k2 * fundamental);
F[0] = -e.vec[0] * factor * D[1];  // Real part
F[1] =  e.vec[0] * factor * D[0];  // Imaginary part
```
Where:
- `k2 = kx² + ky² + kz²` (no fundamental²)
- `fundamental` retrieved from `zeldovich_params_get_fundamental(params_handle)` (line 444)

So:
- `F = rescale * I * e.vec[0] / (k2_index * fundamental) * D`

## Expected Behavior
The formulas should match! Both should give:
- `F = rescale * I * e.vec[0] / (k2_index * fundamental) * D`

## Initial Investigation

Initial hypotheses considered:
1. **`fundamental` not retrieved correctly**: If `params_handle` is NULL or `fundamental` defaults to 1.0, ratio would be ~15.9 (not 8.0)
2. **Missing factor in formula**: Could explain 8x if combined with other issues
3. **`rescale` computed differently**: Checked - values matched exactly
4. **Eigenvector normalization**: This was the actual root cause

## ROOT CAUSE IDENTIFIED ✓

The F, G, H differences were caused by **missing eigenvector normalization** in the hermitian code.

### The Problem

**Zeldovich** applies normalization after getting the eigenvector:
```cpp
double norm = k2 / (kx * ehat.vec[0] + ky * ehat.vec[1] + kz * ehat.vec[2]);
e.vec[i] = norm * ehat.vec[i];
```

**Hermitian** was using the eigenvector directly from `plt_get_eigenmode()` without applying this normalization.

### Evidence

For pure kz modes (kx=0, ky=0):
- **kz ≤ 8**: 
  - Hermitian (before fix): `e.vec[2] = kz/8` (scaled)
  - Zeldovich: `e.vec[2] = kz` (from norm formula)
  - Ratio = 8.0 (constant)
  
- **kz > 8**:
  - Hermitian (before fix): `e.vec[2] = 1.0` (capped)
  - Zeldovich: `e.vec[2] = kz` (continues increasing)
  - Ratio = kz (increases with kz)

This explained why:
- Small kz: ratio was constant 8.0
- Large kz: ratio increased (9.0, 10.0, 20.0, etc.)

## Solution Applied ✓

The fix was implemented in `hermitian_generation.c` by adding the same normalization step that zeldovich uses. After getting the eigenvector from `plt_get_eigenmode()`, the code now:

1. **Sets the sign of z component** (for negative kz, because real FFT only stores +kz half-space):
   ```c
   if (kz < 0) {
       e.vec[2] = -e.vec[2];
   }
   ```

2. **Normalizes eigenvector to unit length** (interpolation might not preserve |e| = 1):
   ```c
   double e_mag = sqrt(e.vec[0] * e.vec[0] + e.vec[1] * e.vec[1] + e.vec[2] * e.vec[2]);
   if (e_mag > 0.0) {
       e.vec[0] /= e_mag;
       e.vec[1] /= e_mag;
       e.vec[2] /= e_mag;
   }
   ```

3. **Applies normalization factor** (the critical missing step):
   ```c
   double k_dot_e = kx * e.vec[0] + ky * e.vec[1] + kz * e.vec[2];
   double norm = (k2 > 0.0 && k_dot_e != 0.0) ? k2 / k_dot_e : 0.0;
   if (!isfinite(norm)) norm = 0.0;
   
   e.vec[0] *= norm;
   e.vec[1] *= norm;
   e.vec[2] *= norm;
   ```

This normalization was applied to all three eigenmode lookup paths:
- Primary path (`e`)
- Self-conjugate path (`e_sc`)
- Another self-conjugate path (`e_sc2`)

## Verification ✓

After applying the fix and re-running the tests:

**Test Results:**
- **Total coordinates compared**: 3,520 in hermitian, 3,521 in zeldovich
- **SUCCESS: All D, F, G, H values match (within tolerance 1e-6)!**
- One coordinate missing in hermitian: Y=511 (x,z)=(511,511) - this is above k²_cutoff and all values are zero, so it's a minor boundary case

**Conclusion:** The eigenvector normalization fix successfully resolved the F, G, H discrepancies. Both codes now produce identical results for all PLT correction values.

