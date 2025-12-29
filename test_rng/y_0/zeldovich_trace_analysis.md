# Brute Force Trace Analysis: zeldovich.cpp Zeroing Logic

## Summary

After tracing through the zeldovich.cpp code for coordinates where differences occur, the **root cause** has been identified: **The hermitian code is missing the `k_cutoff` zeroing condition**.

## Trace: (y=0, z=6, x=6)

### Step-by-Step Execution in zeldovich.cpp

**Step 1: k-vector calculation (lines 334-347)**
- `y = 0`, `z = 6`, `x = 6`
- `ky = 0` (since `0 <= ppdhalf=8`)
- `kz = 6` (since `6 <= ppdhalf=8`)
- `kx = 6` (since `6 <= ppdhalf=8`)
- Result: `k = (6, 0, 6)`

**Step 2: k² calculation (line 349)**
- `k2_int = 6² + 0² + 6² = 72`
- `k2_cutoff = (ppdhalf²) / (k_cutoff²) = 8² / 1.0² = 64.0`

**Step 3: Zeroing condition check (lines 358-367)**

The zeldovich code checks three conditions (line 360-364):
```cpp
if ( (abs(kx)==kmax || abs(kz)==kmax || abs(ky)==kmax)
     || (!param.CornerModes && k2>=k2_cutoff)
     || (param.qonemode && !(...)) ) {
    D = 0.0;
    if (ver == 2) nskip++;
}
```

For (6,6) at y=0:
- Condition 1 (Nyquist): `abs(kx)==8?` → **False** (6 ≠ 8)
- Condition 2 (Nyquist): `abs(ky)==8?` → **False** (0 ≠ 8)
- Condition 3 (Nyquist): `abs(kz)==8?` → **False** (6 ≠ 8)
- Condition 4 (k_cutoff): `!CornerModes && k2_int >= k2_cutoff?`
  - `!False && 72 >= 64.0?`
  - `True && True?`
  - → **True**

**Result**: `D = 0.0` (zeroed by **k_cutoff condition**)

**Step 4: Y=0 mirroring (lines 555-573)**
- Condition `yblock == 0 && yres == 0` is **True**
- Mirroring loop processes `z = 0` to `ppdhalf-1` (0 to 7)
- For `z > 0`: `xmax = ppd` (16)
- (6,6) is **in the mirroring loop**
- Copies from `slabHer[yresHer=7, zHer=10, xHer=10]` to `slab[yres=0, zHer=10, xHer=10]`
- **Important**: The mirroring loop does NOT modify the original (6,6) location
- Since D was zeroed during main processing, it remains zero

### Comparison with Hermitian Code

**Hermitian code zeroing conditions (lines 228-262)**:
1. DC mode: `if (k2 == 0.0)` → **False** (k2 = 72)
2. Nyquist mode: `else if (is_nyquist)` where `is_nyquist = (abs_kx == Nhalf || abs_ky == Nhalf || abs_kz == Nhalf)`
   - `abs_kx = 6`, `abs_ky = 0`, `abs_kz = 6`, `Nhalf = 8`
   - `is_nyquist = (6==8 || 0==8 || 6==8) = False`
3. **k_cutoff condition: NOT IMPLEMENTED**

**Result**: `D = cgauss()` (non-zero)

## Root Cause

**The hermitian code does NOT implement the `k_cutoff` zeroing condition!**

### Zeldovich Code (line 362):
```cpp
|| (!param.CornerModes && k2>=k2_cutoff)
```

This condition zeros all modes where:
- `CornerModes = false` (default)
- `k2_int >= k2_cutoff = (N/2)² / k_cutoff²`

For N=16, k_cutoff=1.0:
- `k2_cutoff = 8² / 1.0² = 64`
- All modes with `k2_int >= 64` are zeroed

### Hermitian Code:
- Only zeros DC modes (`k2 == 0`)
- Only zeros Nyquist modes (`abs_kx == Nhalf || abs_ky == Nhalf || abs_kz == Nhalf`)
- **Does NOT check k_cutoff condition**

## Impact on Differences

### Coordinates Affected

For N=16 with k_cutoff=1.0, modes with `k2 >= 64` are zeroed in zeldovich but not in hermitian.

Examples:
- (6,6) at y=0: `k2 = 6² + 0² + 6² = 72 >= 64` → zeroed in zeldovich
- (7,4) at y=0: `k2 = 7² + 0² + 4² = 65 >= 64` → zeroed in zeldovich
- (4,7) at y=0: `k2 = 4² + 0² + 7² = 65 >= 64` → zeroed in zeldovich
- (7,7) at y=0: `k2 = 7² + 0² + 7² = 98 >= 64` → zeroed in zeldovich

### Pattern Analysis

The 8 non-mirrored Y=0 differences:
- (4,7): `k2 = 4² + 0² + 7² = 65 >= 64` ✓
- (5,7): `k2 = 5² + 0² + 7² = 74 >= 64` ✓
- (6,6): `k2 = 6² + 0² + 6² = 72 >= 64` ✓
- (6,7): `k2 = 6² + 0² + 7² = 85 >= 64` ✓
- (7,4): `k2 = 7² + 0² + 4² = 65 >= 64` ✓
- (7,5): `k2 = 7² + 0² + 5² = 74 >= 64` ✓
- (7,6): `k2 = 7² + 0² + 6² = 85 >= 64` ✓
- (7,7): `k2 = 7² + 0² + 7² = 98 >= 64` ✓

**All 8 coordinates have `k2 >= 64`!**

## Solution

The hermitian code needs to implement the `k_cutoff` zeroing condition to match zeldovich.cpp.

### Required Changes

1. **Read k_cutoff parameter** from parameter file (default: 1.0)
2. **Calculate k2_cutoff**: `k2_cutoff = (Nhalf²) / (k_cutoff²)`
3. **Add zeroing condition**: After Nyquist check, before cgauss() call:
   ```c
   else if (!params->CornerModes && k2_int >= k2_cutoff) {
       D[0] = D[1] = 0.0;
       nskip++;
   }
   ```

### Implementation Notes

- The condition should be checked **after** DC and Nyquist checks
- When D=0, increment `nskip` for RNG consistency
- This applies to **all Y slices**, not just Y=0
- The `CornerModes` parameter should be read from the parameter file (default: false)

## Verification

After implementing the fix:
1. All coordinates with `k2 >= k2_cutoff` should have D=0 in both codes
2. The 8 non-mirrored Y=0 differences should disappear
3. Similar differences at other Y slices should also be resolved
4. RNG alignment should remain correct (nskip accounts for zeroed modes)

## Additional Notes

- The `k_cutoff` parameter controls the maximum wavenumber: `k_max = k_Nyquist / k_cutoff`
- For `k_cutoff = 1.0`, `k_max = 8` (Nyquist frequency)
- For `k_cutoff = 2.0`, `k_max = 4` (half Nyquist)
- The condition `k2 >= k2_cutoff` zeros all modes with `|k| >= k_max`

