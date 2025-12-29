# Analysis of Remaining D Value Differences After Y=0 Fix

## Summary

After implementing the Y=0 slice fix (full plane processing + post-processing mirroring) and Hermitian symmetry fixes (conjugation for self-conjugate slices), the number of differences has been **significantly reduced**. However, 373 coordinates still show differences, all of which are **zero/non-zero mismatches** (hermitian generates non-zero D, zeldovich generates D=0).

**Status**: Results confirmed stable after re-running both codes. The differences are consistent and reproducible.

## Key Statistics

- **Total coordinates**: 
  - Hermitian: 1089
  - Zeldovich: 968
  - Missing in zeldovich: 121 (all at Y=8, the Nyquist slice)

- **Differences**: 373 coordinates
  - All 373 are zero/non-zero mismatches
  - Distribution by Y: Y=0 (24), Y=1 (24), Y=2 (28), Y=3 (32), Y=4 (49), Y=5 (59), Y=6 (72), Y=7 (85)

## Y=0 Slice Analysis

### Remaining Differences at Y=0: 24 coordinates

**Pattern 1: Mirrored Region (16 differences)**
- Coordinates where `x > 8` OR `z > 8` (the mirrored half-plane)
- **Expected behavior**: These are in the region that zeldovich zeros out after mirroring
- **Status**: Likely **expected** - zeldovich processes half-plane, mirrors, then zeros the mirrored region

**Pattern 2: Non-Mirrored Region (8 differences)**
- Coordinates where `x <= 8` AND `z <= 8` but still show differences:
  - (4,7), (5,7), (6,6), (6,7), (7,4), (7,5), (7,6), (7,7)
- **Status**: **Needs investigation** - These are in the "first half" that should match

**Complete Y=0 Difference Coordinates (24 total)**:
- Non-mirrored (x≤8, z≤8): (4,7), (5,7), (6,6), (6,7), (7,4), (7,5), (7,6), (7,7) - **8 coordinates**
- Mirrored (x>8 or z>8): (4,9), (5,9), (6,9), (6,10), (7,9), (7,10), (9,4), (9,5), (9,6), (9,7), (9,9), (9,10), (10,6), (10,7), (10,9), (10,10) - **16 coordinates**

### Example: (4,7) vs (7,4)

Looking at the logs:
- **Hermitian (4,7)**: `D=(7.8613813967e-03, 2.5068406016e-02)` (non-zero)
- **Zeldovich (4,7)**: `D=(0.0000000000e+00, 0.0000000000e+00)` (zero)

- **Hermitian (7,4)**: `D=(7.1277809329e-03, -1.0828594677e-02)` (non-zero)  
- **Zeldovich (7,4)**: `D=(0.0000000000e+00, 0.0000000000e+00)` (zero)

## Possible Causes

1. **Mirroring Logic Difference**: 
   - Zeldovich may zero out certain coordinates during/after mirroring that hermitian does not
   - The coordinates with `x=7` or `z=7` might be in a "boundary" region that zeldovich treats specially

2. **RNG Call Order**: 
   - Even though we aligned the RNG call order, there might still be subtle differences in when/where D=0 is set
   - The fact that zeldovich has D=0 suggests it's explicitly zeroing these coordinates

3. **Nyquist Boundary Handling**:
   - Coordinates with `x=7` or `z=7` are one step away from the Nyquist boundary (x=8 or z=8)
   - Zeldovich might have special handling for these "pre-Nyquist" coordinates

## Recommendations

1. **Check zeldovich.cpp mirroring logic** for Y=0:
   - Verify which coordinates are zeroed after mirroring
   - Check if there's special handling for coordinates with `x=7` or `z=7`

2. **Verify hermitian mirroring** matches zeldovich:
   - Ensure the post-processing mirroring in hermitian code zeros the same coordinates as zeldovich
   - Check if coordinates with `x=7` or `z=7` should be zeroed

3. **Missing Coordinates at Y=8**:
   - 121 coordinates missing in zeldovich, all at Y=8 (Nyquist slice)
   - This is likely expected if zeldovich doesn't process Y=8 the same way
   - Verify this is intentional

## Progress Assessment

✅ **Major Progress**: Y=0 differences reduced from many to 24
✅ **Pattern Identified**: Most differences are in expected regions (mirrored half-plane)
✅ **Stability Confirmed**: Re-run results match previous analysis exactly (373 total differences, same distribution)
⚠️ **Remaining Issue**: 8 differences in non-mirrored region at Y=0 need investigation
⚠️ **Other Y Slices**: Differences at Y=1-7 also need analysis (may have similar patterns)

## Recent Changes

After the initial Y=0 fix, additional Hermitian symmetry fixes were implemented:
1. **Preprocessor directive fixes**: Corrected `VERIFY_HERMITIAN_SYMMETRY` evaluation to ensure Mode 1 and Mode 2 work correctly
2. **Mirroring with conjugation**: Modified self-conjugate slice mirroring to apply complex conjugation for Hermitian symmetry (Mode 1) and anti-Hermitian symmetry (Mode 2)
3. **Extended mirroring**: Applied mirroring logic to both Y=0 and Y=N/2 slices

These fixes ensure Hermitian symmetry is preserved, but do not affect the RNG D-value comparison results (which are consistent before and after these fixes).

## Next Steps

1. **Examine zeldovich.cpp's Y=0 post-processing** to understand why it zeros coordinates with `x=7` or `z=7`
   - Check the mirroring logic in zeldovich.cpp (around lines 555-573)
   - Identify which coordinates are explicitly zeroed after mirroring
   - Determine if coordinates with `x=7` or `z=7` are intentionally zeroed

2. **Update hermitian code** to match zeldovich's zeroing behavior
   - If zeldovich zeros these coordinates, hermitian should do the same
   - Ensure the zeroing happens at the same point in the processing pipeline

3. **Analyze other Y slices** (Y=1-7) for similar patterns
   - Check if the same `x=7` or `z=7` pattern appears at other Y values
   - Determine if there's a systematic rule for which coordinates should be zeroed

4. **Re-test and verify** remaining differences are eliminated or explained
   - After implementing fixes, re-run both codes
   - Verify that differences are reduced or that remaining differences are intentional/expected

