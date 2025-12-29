# Detailed Implementation Plan: Fix Y=0 Slice Processing to Match Zeldovich Code

## Objective
Align Y=0 (self-conjugate slice) processing in hermitian code with zeldovich code behavior to eliminate the 56 numerical differences and 14 zero-mismatches observed at Y=0.

## Current State Analysis

### Zeldovich Code Behavior (`zeldovich.cpp` lines 555-573):
1. **Processes full Y=0 plane** (all x, z coordinates from 0 to ppd-1)
2. **After processing**, performs post-processing step:
   - If `yblock == 0 && yres == 0`:
     - Loops z = 0 to `ppdhalf`
     - For each z, loops x = 0 to `xmax` where:
       - `xmax = ppdhalf` if z == 0 (half line)
       - `xmax = ppd` if z > 0 (full line)
     - Copies from `slabHer` to `slab`: 
       - `AYZX(slab, a, yres, zHer, xHer) = AYZX(slabHer, a, yresHer, zHer, xHer)`
     - Sets origin (0,0,0) to zero

### Hermitian Code Current Behavior (`hermitian_generation.c`):
1. **Y=0 is correctly identified as self-conjugate** (line 238 in main.cpp: `y_mirror = y_primary` for Y=0)
2. **Enters self-conjugate branch** (line 639: `#else` when `y_mirror == global_y`)
3. **Uses zeldovich_method** (line 678: `#if USE_ZELDOVICH_METHOD`)
4. **Current zeldovich_method implementation**:
   - Processes **half plane only**: z = 0 to `Nhalf`, x varies by z
   - For z=0: x = 0 to `Nhalf+1`
   - For z>0: x = 0 to `N-1`
   - **Mirrors during processing** (lines 1052-1081): stores mirrored values immediately
   - Sets origin to zero (line 1164)

### Key Differences Identified:
1. **Processing scope**: Zeldovich processes full plane, Hermitian processes half plane
2. **Mirroring timing**: Zeldovich mirrors AFTER processing, Hermitian mirrors DURING processing
3. **RNG call order**: Different coordinate traversal order may cause RNG state misalignment

## Implementation Plan

### Phase 1: Verify Current Configuration
**Goal**: Ensure Y=0 is using zeldovich_method

**Steps**:
1. Verify `USE_ZELDOVICH_METHOD` is defined and set to 1 in `config.h` (line 32)
2. Add debug print to confirm Y=0 enters zeldovich_method branch
3. Document current behavior with test run

**Files to modify**:
- `src/config.h`: Verify `USE_ZELDOVICH_METHOD` definition
- `src/generation/hermitian_generation.c`: Add debug print at line 678

**Expected outcome**: Confirmation that Y=0 uses zeldovich_method branch

---

### Phase 2: Modify Zeldovich Method to Process Full Plane
**Goal**: Change zeldovich_method to process full Y=0 plane (like zeldovich.cpp), then mirror

**Current loop structure** (lines 685-1083):
```c
for (int z = 0; z <= Nhalf; z++) {
    int x_max = (z == 0 ? Nhalf + 1 : N);
    for (int x = 0; x < x_max; x++) {
        // Generate D, compute F,G,H
        // Store in primary slice
        // Mirror immediately if x != x_mirror || z != z_mirror
    }
}
```

**New loop structure** (match zeldovich.cpp):
```c
// Process FULL plane first (z = 0 to N-1, x = 0 to N-1)
for (int z = 0; z < N; z++) {
    for (int x = 0; x < N; x++) {
        // Generate D, compute F,G,H
        // Store ONLY in primary slice (no mirroring yet)
    }
}

// Post-processing: Mirror first half to second half (match zeldovich.cpp lines 555-573)
if (global_y == 0) {
    for (int z = 0; z <= Nhalf; z++) {
        int z_mirror = (z == 0) ? 0 : N - z;
        int x_max = (z == 0 ? Nhalf : N);  // Match zeldovich: ppdhalf for z=0, ppd for z>0
        for (int x = 0; x < x_max; x++) {
            int x_mirror = (x == 0) ? 0 : N - x;
            // Copy from primary to mirrored location
            for (int a = 0; a < narray; a++) {
                PRIM_SLICE(a, x_mirror, z_mirror) = PRIM_SLICE(a, x, z);
            }
        }
    }
    // Set origin to zero
    for (int a = 0; a < narray; a++) {
        PRIM_SLICE(a, 0, 0)[0] = 0.0;
        PRIM_SLICE(a, 0, 0)[1] = 0.0;
    }
}
```

**Files to modify**:
- `src/generation/hermitian_generation.c`:
  - Lines 685-1083: Replace half-plane loop with full-plane loop
  - Remove immediate mirroring logic (lines 1052-1081)
  - Add post-processing mirroring step after loop (new code block)

**Key changes**:
1. Change loop bounds: `z <= Nhalf` → `z < N`, `x_max` logic → `x < N`
2. Remove mirroring inside loop
3. Add post-processing mirroring block (only for Y=0)
4. Ensure origin is set to zero after mirroring

**RNG considerations**:
- Full plane processing will generate D values for all (x,z) coordinates
- RNG calls must match zeldovich.cpp order exactly
- Verify RNG skip logic still works correctly with full plane

---

### Phase 3: Align RNG Call Order
**Goal**: Ensure RNG calls match zeldovich.cpp exactly

**Analysis needed**:
1. Compare RNG call order between:
   - Zeldovich: Full plane (z=0 to ppd-1, x=0 to ppd-1)
   - Hermitian (new): Full plane (z=0 to N-1, x=0 to N-1)
2. Verify RNG skip accumulation matches:
   - Zeldovich: Lines 338, 344 (skip at Nyquist boundaries)
   - Hermitian: Lines 690-706 (current skip logic)

**Files to modify**:
- `src/generation/hermitian_generation.c`:
  - Update RNG skip logic to match full-plane traversal
  - Ensure skip accumulation happens at same points as zeldovich

**Verification**:
- Run RNG debug comparison after changes
- Verify RNG call order matches zeldovich exactly

---

### Phase 4: Handle Missing Coordinates
**Goal**: Ensure all coordinates are processed (no missing entries)

**Current issue**:
- Hermitian missing: 24 coordinates at Y=0 (x <= 4, z >= 9)
- Zeldovich missing: 97 coordinates at Y=8

**Analysis**:
- Missing coordinates suggest different loop bounds or filtering
- Full plane processing should eliminate missing coordinates in Y=0

**Files to modify**:
- `src/generation/hermitian_generation.c`:
  - Ensure full plane loop processes all coordinates
  - Verify no early exits or filtering that skips coordinates

**Verification**:
- Run comparison after changes
- Verify no missing coordinates at Y=0

---

### Phase 5: Testing and Validation
**Goal**: Verify changes eliminate Y=0 differences

**Test procedure**:
1. Rebuild hermitian code with changes
2. Run N=16 test case
3. Generate RNG debug logs
4. Run comparison: `python3 compare_rng_d_values.py`
5. Verify:
   - Y=0 differences reduced to 0 (or minimal, within tolerance)
   - No missing coordinates at Y=0
   - RNG call order matches zeldovich

**Success criteria**:
- Y=0 differences: 0 (or all within floating-point tolerance)
- Missing coordinates at Y=0: 0
- RNG debug logs show matching call order

---

## Implementation Order

1. **Phase 1** (Verification): Quick check, no code changes
2. **Phase 2** (Full plane processing): Core change, highest impact
3. **Phase 3** (RNG alignment): Critical for correctness
4. **Phase 4** (Missing coordinates): Should be resolved by Phase 2
5. **Phase 5** (Testing): Validate all changes

## Risk Assessment

**Low risk**:
- Phase 1: Read-only verification
- Phase 5: Testing only

**Medium risk**:
- Phase 2: Structural change but well-defined
- Phase 4: Should be automatic with Phase 2

**High risk**:
- Phase 3: RNG misalignment could cause cascading differences

## Rollback Plan

If issues arise:
1. Git commit before each phase
2. Can revert to previous commit if phase fails
3. Test after each phase before proceeding

## Notes

- The zeldovich.cpp code processes full plane, then mirrors. This is the key difference.
- Current hermitian code mirrors during processing, which may cause RNG state differences.
- Full plane processing should naturally eliminate missing coordinates.
- Post-processing mirroring matches zeldovich.cpp exactly (lines 555-573).

## Questions for Review

1. Should we maintain the `USE_ZELDOVICH_METHOD` flag, or always use this method for Y=0?
- Answer: Always use this method for self-conjugate slices, including Y=0 and nyquist slice. for now, keep the flag.
2. Should Y=N/2 (also self-conjugate) use the same approach?
- Answer: Yes, use the same approach for Y=N/2.
3. Do we need to handle any edge cases in the mirroring logic?

@InitialConditions/hermitian_3d_matrix_production/test_rng/Y0_FIX_PLAN.md i am confused why processing half the plane vs processing the full plane would cause rng misalignment for the first half of the plane since both codes should call the same rng values for the first half of the plane. later, these values are mirrored to the second half of the plane so they should still be identical for both codes
 
