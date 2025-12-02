# RNG Skip Logic Comparison: hermitian_generation.c vs zeldovich.cpp

## Overview
Both codes skip RNG calls for grid points that don't exist when N < MAX_PPD, ensuring RNG consistency across different matrix sizes. However, they implement this differently.

## Key Differences

### 1. **When Skip is Detected**

#### zeldovich.cpp (LoadPlane function):
- Detects skip **at the boundary** when crossing Nyquist:
  - `z == ppdhalf + 1`: Skip `(MAX_PPD - ppd) * MAX_PPD` (missing z-rows)
  - `x == ppdhalf + 1`: Skip `MAX_PPD - ppd` (missing x-values in current z-row)
- Also increments `nskip++` when `D = 0.0` (forced to zero)
- **Advance happens inside the loop** when `nskip > 0` before calling `cgauss()`

#### hermitian_generation.c:
- **Two strategies** depending on `PARALLELIZE_XZ_WITHIN_SLICE`:
  
  **Parallel mode** (`PARALLELIZE_XZ_WITHIN_SLICE=1`):
  - Calculates total skip **upfront** before loops
  - Conjugate pairs: `nskip = (MAX_PPD - N) * MAX_PPD + N * (MAX_PPD - N)`
  - Self-conjugate: `nskip = (MAX_PPD - N) * (1 + N/2)`
  
  **Sequential mode** (`PARALLELIZE_XZ_WITHIN_SLICE=0`):
  - Tracks skip **incrementally** like zeldovich.cpp
  - Detects at boundaries: `x == N - 1` and `z == N - 1`
  - Advances **immediately after** finishing each z-row or at end of all z-rows

### 2. **Skip Calculation Details**

#### zeldovich.cpp:
```cpp
// Line 338: When crossing z-wrap (z == ppdhalf + 1)
if (z == ppdhalf + 1 && ver == 2) nskip += (MAX_PPD - ppd) * MAX_PPD;

// Line 344: When crossing x-wrap (x == ppdhalf + 1)  
if (x == ppdhalf + 1 && ver == 2) nskip += MAX_PPD - ppd;

// Line 361: When D is forced to zero
if (ver == 2) nskip++;

// Line 363-366: Advance before cgauss() if nskip > 0
if (nskip) {
    Pk.v2rng[y].advance(2 * nskip);
    nskip = 0;
}
```

#### hermitian_generation.c (Sequential mode):
```c
// Conjugate pairs: After finishing x-loop (x == N - 1)
if (x == N - 1 && N < MAX_PPD) {
    nskip += MAX_PPD - N;  // Skip missing x-values
    // Advance immediately
    advance_rng(...);
    nskip = 0;
}

// Conjugate pairs: After finishing z-loop (z == N - 1)
if (z == N - 1 && N < MAX_PPD) {
    nskip += (MAX_PPD - N) * MAX_PPD;  // Skip missing z-rows
    // Advance immediately
    advance_rng(...);
    nskip = 0;
}
```

### 3. **Boundary Detection**

#### zeldovich.cpp:
- Uses `ppdhalf + 1` (where `ppdhalf = ppd / 2`)
- Detects **before** processing the first element after crossing Nyquist
- This is the point where we transition from positive to negative k-space

#### hermitian_generation.c:
- Uses `N - 1` (last valid index)
- Detects **after** processing the last element in a row/plane
- Advances immediately after finishing the row/plane

### 4. **Advance Timing**

#### zeldovich.cpp:
- Advances **before** calling `cgauss()` when `nskip > 0`
- Multiple advances possible within the loop
- Final advance at end: `Pk.v2rng[y].advance(2 * nskip)` (line 483)
- **Assertion**: `assert(Pk.v2rng[y] - checkpoint == 2 * MAX_PPD * MAX_PPD)`

#### hermitian_generation.c:
- Sequential mode: Advances **immediately after** detecting skip boundaries
- Parallel mode: Advances **before** entering loops (upfront calculation)
- No final assertion (but could add one for verification)

### 5. **D=0 Handling**

#### zeldovich.cpp:
- When `D = 0.0` (forced to zero), increments `nskip++` (line 361)
- This skip is accumulated and advanced later

#### hermitian_generation.c:
- When `D = 0.0` (DC mode or Nyquist), advances **immediately**:
  - `zeldovich_ps_advance_rng(..., 1)` or `advance_pcg_global(..., 2)`
- Does NOT accumulate D=0 skips into `nskip` variable

### 6. **Total Skip Amount**

Both should skip the same total amount for consistency:

**For conjugate pairs (full N×N plane):**
- Missing z-rows: `(MAX_PPD - N) * MAX_PPD`
- Missing x-values: `N * (MAX_PPD - N)` (one per z-row)
- **Total**: `(MAX_PPD - N) * (MAX_PPD + N)`

**For self-conjugate (half plane, z = 0 to N/2):**
- Missing x-values: `(MAX_PPD - N) * (1 + N/2)`
  - One skip for z=0 row (x goes 0 to N/2+1)
  - N/2 skips for z > 0 rows (x goes 0 to N-1)

### 7. **Key Architectural Difference**

#### zeldovich.cpp:
- **Single-pass accumulation**: Accumulates all skips, advances when needed
- **Lazy advance**: Only advances when about to use RNG
- Simpler logic, but requires checking `nskip` before each RNG call

#### hermitian_generation.c:
- **Immediate advance**: Advances as soon as skip is detected
- **Eager advance**: Resets `nskip` to 0 after each advance
- More complex logic with multiple advance points, but clearer separation

## Potential Issues

### 1. **D=0 Skip Handling**
- **zeldovich.cpp**: Accumulates D=0 skips into `nskip`, advances later
- **hermitian_generation.c**: Advances immediately for D=0, doesn't accumulate
- **Impact**: If there are many D=0 cases, the timing of advances differs, but total should be same

### 2. **Boundary Condition**
- **zeldovich.cpp**: Detects at `ppdhalf + 1` (first element after Nyquist)
- **hermitian_generation.c**: Detects at `N - 1` (last element before missing region)
- **Impact**: Should be equivalent, but different detection points

### 3. **Final Verification**
- **zeldovich.cpp**: Has assertion `assert(Pk.v2rng[y] - checkpoint == 2 * MAX_PPD * MAX_PPD)`
- **hermitian_generation.c**: No such assertion
- **Recommendation**: Add similar assertion to verify total RNG calls match expected

## Implementation Status

### Completed
1. **Added RNG verification** in hermitian_generation.c (similar to zeldovich.cpp assertion)
   - Added `VERIFY_RNG_CALLS` config flag (default: 0, disabled)
   - Tracks `total_rng_calls`, `total_rng_skips`, and `total_d_zero_skips`
   - Verifies: `2 * total_rng_calls + 2 * total_rng_skips + 2 * total_d_zero_skips == 2 * MAX_PPD * MAX_PPD`
   - Prints verification status for first few slices when enabled

### Next steps

1. **Test verification** - Enable `VERIFY_RNG_CALLS=1` and run tests to verify RNG consistency
2. **Consider unifying D=0 handling** - either accumulate or advance immediately, but be consistent
3. **Document the equivalence** of the two approaches mathematically
4. **Add debug output** to compare total skip amounts between implementations

