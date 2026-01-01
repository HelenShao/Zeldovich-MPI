# F, G, H Value Discrepancy Diagnosis

## Problem
D values match perfectly between hermitian and zeldovich codes, but F, G, H values differ by a factor of ~8.0.

## Example Coordinate
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

## Possible Causes

1. **`fundamental` not retrieved correctly**: If `params_handle` is NULL or `fundamental` defaults to 1.0, then:
   - Hermitian: `factor = rescale / (k2 * 1.0) = rescale / k2`
   - Zeldovich: `factor = rescale / (k2 * 0.0628)`
   - Ratio would be: `1.0 / 0.0628 ≈ 15.9` (NOT 8.0)

2. **Missing factor in formula**: If there's a missing `fundamental` in the denominator AND something else, could explain 8x.

3. **`rescale` computed differently**: Check if `a_NL`, `a0`, `target_f`, or `plt_f` differ between codes.

4. **Eigenvector normalization**: Check if `e.vec[2]` values differ between codes.

## ROOT CAUSE IDENTIFIED

The F, G, H differences are caused by **missing eigenvector normalization** in the hermitian code.

### The Problem

**Zeldovich** applies normalization after getting the eigenvector:
```cpp
double norm = k2 / (kx * ehat.vec[0] + ky * ehat.vec[1] + kz * ehat.vec[2]);
e.vec[i] = norm * ehat.vec[i];
```

**Hermitian** uses the eigenvector directly from `plt_get_eigenmode()` without applying this normalization.

### Evidence

For pure kz modes (kx=0, ky=0):
- **kz ≤ 8**: 
  - Hermitian: `e.vec[2] = kz/8` (scaled)
  - Zeldovich: `e.vec[2] = kz` (from norm formula)
  - Ratio = 8.0 (constant)
  
- **kz > 8**:
  - Hermitian: `e.vec[2] = 1.0` (capped)
  - Zeldovich: `e.vec[2] = kz` (continues increasing)
  - Ratio = kz (increases with kz)

This explains why:
- Small kz: ratio is constant 8.0
- Large kz: ratio increases (9.0, 10.0, 20.0, etc.)

## Possible Causes

1. **`fundamental` not retrieved**: If `params_handle` is NULL or `fundamental` defaults to 1.0, the ratio would be ~15.9 (not 8.0)
2. **Missing factor of 8**: Could be related to BoxSize, fundamental, or a normalization factor
3. **Different `fundamental` computation**: If hermitian uses `fundamental = 1.0/BoxSize` instead of `2π/BoxSize`, ratio would be ~6.3 (not 8.0)

## Critical Check Needed

**Verify that `fundamental` is being retrieved correctly in hermitian code when computing F, G, H.**

The code at line 444 should retrieve `fundamental` from `params_handle`, but if `params_handle` is NULL or the retrieval fails, it defaults to 1.0 (line 440).

## Next Steps

1. **Add debug output** to print `fundamental`, `rescale`, `e.vec[2]`, and `factor` for test coordinate Y=0, (x,z)=(0,1)
2. **Verify `params_handle` is not NULL** when computing F, G, H (check the call site)
3. **Compare `fundamental` values** - should be `2π/100 ≈ 0.0628` in both codes
4. **Check if there's a missing factor** - the 8x suggests a specific scaling issue

