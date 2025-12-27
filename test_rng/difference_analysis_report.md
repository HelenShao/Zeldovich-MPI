# In-Depth Analysis of D Value Differences Between Hermitian and Zeldovich Codes

## Executive Summary

Analysis of 419 coordinate differences reveals several critical patterns:

1. **363 out of 419 differences (87%) are zero/non-zero mismatches** - one code generates D=0 while the other generates a non-zero value
2. **Y=0 slice has 70 differences, with 56 being actual numerical differences** (not just zero mismatches)
3. **Most differences occur at high k values** (|k| >= 6), particularly at kx=±7, ±6
4. **Missing coordinates**: 24 missing in hermitian (all at Y=0), 97 missing in zeldovich (all at Y=8)

## Key Findings

### 1. Y=0 Slice Special Handling

**Critical Difference**: The two codes handle Y=0 (self-conjugate slice) differently:

#### Zeldovich Code (`zeldovich.cpp`):
- Line 325: `y = yres + yblock * array.block;` - computes actual y coordinate
- Line 334: `ky = y > ppdhalf ? y - ppd : y;` - computes ky from y
- Line 555-573: **Special post-processing for Y=0**:
  - After generating D values, it **copies the first half plane onto the second half**
  - Sets origin (0,0,0) to zero
  - This means Y=0 is processed as a full plane, then symmetrized

#### Hermitian Code (`hermitian_generation.c`):
- Line 141: `if (y_mirror != global_y)` - processes as conjugate pair
- Line 1177: `#else` branch - processes as self-conjugate case (when `y_mirror == global_y`)
- Line 1164: `if (global_y == 0)` - sets origin to zero
- **No copying/symmetrization step** - processes all (x,z) coordinates independently

**Impact**: This fundamental difference in Y=0 handling explains why Y=0 has 56 actual numerical differences (not just zero mismatches).

### 2. Zero/Non-Zero Mismatches (363 cases)

**Pattern**: One code generates D=0 while the other generates a non-zero value.

**Likely Causes**:
1. **Different Nyquist frequency detection logic**:
   - Hermitian: Line 225 - `is_nyquist = (abs_kx == Nhalf || abs_ky == Nhalf || abs_kz == Nhalf)`
   - Zeldovich: Line 360 - `abs(kx)==kmax || abs(kz)==kmax || abs(ky)==kmax` where `kmax = (double) ppdhalf * ik_cutoff + .5`
   - The `ik_cutoff` factor in zeldovich may cause different modes to be zeroed

2. **Different k_cutoff handling**:
   - Zeldovich: Line 362 - `!param.CornerModes && k2>=k2_cutoff` - zeros modes above cutoff
   - Hermitian: No equivalent k_cutoff check in the D generation logic

3. **RNG skip differences**: When D=0, RNG calls are skipped, but the skip accumulation/application may differ between codes

### 3. High k Value Concentration

**Statistics**:
- 385 out of 419 differences (92%) occur at high k (|k| >= 6)
- Top kx values: ±7 (66 each), ±6 (55 each)
- Top ky values: 7 (85), 6 (72), 0 (70)
- Top kz values: 7 (69), -7 (59), 6 (58)

**Analysis**: High k values are more sensitive to:
- Nyquist frequency detection differences
- k_cutoff filtering differences
- RNG skip accumulation differences (missing grid points at high k)

### 4. Missing Coordinates

**Hermitian missing (24 coordinates, all Y=0)**:
- Pattern: Missing at (x,z) where x <= 4 and z >= 9
- Example: Y=0 (x,z)=(0,9), (0,10), (1,9), (1,10), etc.

**Zeldovich missing (97 coordinates, all Y=8)**:
- Pattern: Missing at Y=8 for all (x,z) combinations
- Y=8 corresponds to ky = 8 - 16 = -8 (negative k region)

**Analysis**: 
- Hermitian may be skipping certain coordinates in Y=0 processing
- Zeldovich may not be processing Y=8 at all (possibly treated as conjugate of Y=0?)

### 5. RNG Index Differences

**Potential Issue**: The RNG indexing may differ:
- Hermitian: Uses `global_y` directly as RNG index (line 275, 1105, 1254)
- Zeldovich: Uses `y` (computed from `yres + yblock * array.block`) as RNG index (line 332, 373)

**For Y=0**:
- Hermitian: `rng_index = global_y = 0`
- Zeldovich: `y = yres + yblock * array.block` - depends on block structure

This could cause RNG state misalignment if the block structure differs.

## Root Cause Analysis

### Primary Issue: Y=0 Processing Difference

The most significant difference is in how Y=0 (self-conjugate slice) is handled:

1. **Zeldovich**: Processes full plane, then symmetrizes by copying first half to second half
2. **Hermitian**: Processes all coordinates independently, no symmetrization step

This explains:
- Why Y=0 has 56 actual numerical differences (not just zero mismatches)
- Why coordinates are missing in different patterns
- Why RNG calls may be out of sync

### Secondary Issues:

1. **Nyquist frequency detection**: Different logic for determining which modes to zero
2. **k_cutoff filtering**: Zeldovich has additional filtering that Hermitian lacks
3. **RNG skip accumulation**: Different timing of when skips are accumulated and applied

## Recommendations

1. **Align Y=0 processing**: Make Hermitian code match Zeldovich's approach:
   - Process full plane for Y=0
   - Add symmetrization step (copy first half to second half)
   - Ensure origin is set to zero

2. **Unify Nyquist detection**: Ensure both codes use identical logic for detecting Nyquist frequencies

3. **Add k_cutoff check**: If Zeldovich's k_cutoff filtering is desired, add equivalent logic to Hermitian

4. **Verify RNG indexing**: Ensure RNG indices match between codes, especially for Y=0

5. **Debug RNG skip logic**: Add detailed logging to compare skip accumulation and application between codes

## Code Locations for Investigation

### Hermitian Code:
- `hermitian_generation.c` line 141: Conjugate pair vs self-conjugate branch
- `hermitian_generation.c` line 1177: Self-conjugate processing (Y=0 case)
- `hermitian_generation.c` line 225: Nyquist frequency detection
- `hermitian_generation.c` line 275: RNG index for Y slices

### Zeldovich Code:
- `zeldovich.cpp` line 325: y coordinate computation
- `zeldovich.cpp` line 334: ky computation
- `zeldovich.cpp` line 360: Nyquist frequency detection
- `zeldovich.cpp` line 555: Y=0 special handling
- `zeldovich.cpp` line 332: RNG checkpoint for y

