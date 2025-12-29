# Zeldovich Y=0 Mirroring Logic Analysis

## Summary

After analyzing `zeldovich.cpp` lines 555-573, here's what we found about the Y=0 post-processing mirroring logic.

## Mirroring Loop Structure (Lines 555-573)

```cpp
if (yblock == 0 && yres == 0) {
    // Copy the first half plane onto the second
    for (z = 0; z < ppdhalf; z++) {  // z = 0 to 7 for N=16
        zHer = ppd - z;
        if (z == 0) zHer = 0;
        // Treat y=z=0 as a half line
        int xmax = (z == 0 ? ppdhalf : ppd);  // xmax = 8 for z=0, 16 for z>0
        for (x = 0; x < xmax; x++) {
            xHer = ppd - x;
            if (x == 0) xHer = 0;
            for (a = 0; a < array.narray; a++) {
                AYZX(slab, a, yres, zHer, xHer) =
                   (AYZX(slabHer, a, yresHer, zHer, xHer));
            }
        }
    }
    // And the origin must be zero
    for (a = 0; a < array.narray; a++) AYZX(slab, a, 0, 0, 0) = 0.0;
}
```

## Key Observations

### 1. Mirroring Loop Coverage

For N=16 (ppd=16, ppdhalf=8):
- **z = 0**: Processes x = 0 to 7 (xmax = 8)
- **z = 1 to 7**: Processes x = 0 to 15 (xmax = 16)
- **Total coordinates in mirroring loop**: 120 coordinates

### 2. What the Mirroring Does

- **Copies from `slabHer` to `slab`** at mirrored locations
- **Does NOT modify** the original coordinates (x,z) in the first half
- **Only explicitly zeros** the origin (0,0,0)

### 3. Coordinates NOT in Mirroring Loop

Coordinates NOT processed by the mirroring loop:
- **z = 0**: x = 8 to 15 (8 coordinates)
- **z = 8 to 15**: x = 0 to 15 (128 coordinates)
- **Total**: 136 coordinates

These coordinates are in the "second half" and should be zeroed (they're in the mirrored region).

### 4. Coordinates with x=7 or z=7

**All coordinates with x=7 or z=7 that are in the first half ARE in the mirroring loop:**
- (7,0), (7,1), ..., (7,7) - all IN loop
- (0,7), (1,7), ..., (15,7) - all IN loop (for z=7, xmax=16)

**These coordinates should NOT be zeroed** because:
1. They're in the first half (z < 8)
2. They're processed in the mirroring loop
3. They're not explicitly zeroed in the code
4. They don't meet Nyquist zeroing conditions (|kx|≠8, |kz|≠8)

### 5. Nyquist Zeroing Conditions

Coordinates are zeroed during main processing if:
- `|kx| == kmax` OR `|kz| == kmax` OR `|ky| == kmax` (where kmax = ppdhalf * ik_cutoff)
- For N=16 with ik_cutoff=1.0: kmax = 8
- Coordinates with x=7 or z=7 have |kx|=7 or |kz|=7, so they are **NOT zeroed** by Nyquist conditions

### 6. Understanding k_cutoff Parameter

**k_cutoff is NOT the Nyquist wavenumber itself.** Instead:

- **Definition**: `k_max = k_Nyquist / k_cutoff`
  - If `k_cutoff = 1.0`: `k_max = k_Nyquist` (all modes up to Nyquist)
  - If `k_cutoff = 2.0`: `k_max = k_Nyquist/2` (only modes up to half-Nyquist)
  - If `k_cutoff = 1.5`: `k_max = k_Nyquist/1.5` (modes up to 2/3 of Nyquist)

- **In code** (line 305): `ik_cutoff = 1.0 / param.k_cutoff`
- **kmax calculation** (line 359): `kmax = (double) ppdhalf * ik_cutoff + .5`
  - This gives: `kmax = ppdhalf / k_cutoff`
  - For N=16 (ppdhalf=8) with k_cutoff=1.0: `kmax = 8` (the Nyquist)
  - For N=16 with k_cutoff=2.0: `kmax = 4` (half-Nyquist)

- **Default value**: `k_cutoff = 1.0` (from `parameters.cpp` line 37)

- **Effect on coordinates with x=7 or z=7**:
  - For k_cutoff=1.0: kmax=8, so |kx|=7 and |kz|=7 are **NOT zeroed**
  - For k_cutoff=1.5: kmax=5, so |kx|=7 and |kz|=7 **ARE zeroed**
  - For k_cutoff=2.0: kmax=4, so |kx|=7 and |kz|=7 **ARE zeroed**

**Important**: If the test uses k_cutoff > 1.0, then coordinates with |kx|=7 or |kz|=7 would be zeroed by the Nyquist condition check, which could explain the differences!

**Note**: Confirmed by user: the tests use the default `k_cutoff = 1.0` (the parameter file `param_N16_rng.par` does not specify `ZD_k_cutoff`). This means `kmax = 8` for N=16, so coordinates with |kx|=7 or |kz|=7 are **NOT** zeroed by the Nyquist condition. Therefore, **k_cutoff is NOT the explanation** for why coordinates with x=7 or z=7 are being zeroed.

### 7. The Shift and Y=ppd/2 Zeroing (Lines 710-711)

In the `ZeldovichXY` function (which does Y & X inverse FFTs):

```cpp
// The Nyquist frequency y=array.ppd/2 must now be set to 0
// because we shifted the data by one location.
// FLAW: this assumes PPD is even.
int y = array.ppd / 2;
for (int zres = 0; zres < array.block; zres++) {
    for (int a = 0; a < array.narray; a++) {
        for (int x = 0; x < array.ppd; x++) AZYX(slab, a, zres, y, x) = 0.0;
    }
}
```

**What this does:**
- **Zeros the entire Y=ppd/2 slice** (the Nyquist frequency slice in Y direction)
- This happens **AFTER** the Y=0 mirroring processing
- This is a **separate operation** from the Y=0 mirroring we're investigating
- The comment mentions "we shifted the data by one location" - this is related to the FFT implementation details

**Impact on our investigation:**
- This zeroing affects Y=8 (for N=16), not Y=0
- This is **NOT** the cause of the Y=0 differences we're seeing
- **However, this does NOT explain the missing Y=8 coordinates in RNG debug output**

## Why Only `yblock < numblock/2` is Processed (Hermitian Symmetry)

The main loop at line 628 only processes `yblock = 0` to `numblock/2 - 1` because of **Hermitian symmetry**. Here's how it works:

**For N=16 with numblock=2:**
- `block = ppd / numblock = 8`
- `yblock = 0`: Y = 0 to 7
- `yblock = 1`: Y = 8 to 15

**Processing yblock=0:**
1. Generates primary values in `slab` for Y=0 to 7
2. Simultaneously generates conjugated/reflected values in `slabHer` (indexed by `yresHer = block - 1 - yres = 7 - yres`)
3. At line 656:
   ```cpp
   array.StoreBlock(yblock, zblock, slab);                    // Y=0 to 7
   array.StoreBlock(numblock - 1 - yblock, zblock, slabHer); // Y=8 to 15
   ```
4. `slabHer` (with `yresHer = 0` to `7`) is stored into `yblock = 1`, which maps to global Y = `yresHer + 1 * block = 8` to `15`

**The key insight**: Due to Hermitian symmetry `f(-k) = conj(f(k))`, you only need to process the first half of Y values. The second half (Y=8 to 15) is automatically filled by storing the conjugated values from `slabHer`.

**Y=8 (Nyquist slice)**: Y=8 is `ppd/2`, the Nyquist frequency. It's filled from `slabHer[yresHer=0]`, which contains the conjugate of values at Y=0. However:
- Y=8 is **never directly processed** in the main loop (no RNG calls, no D value generation)
- It's filled by **conjugation from Y=0** when storing `slabHer` into `yblock=1`
- The RNG debug output (line 503) is printed **during** the main processing loop, so Y=8 never appears in it
- Y=8 is later zeroed in `ZeldovichXY` (line 718) due to the FFT shift, but this happens after the RNG debug output

## Why Y=8 coordinates are missing in RNG debug output:

**Confirmed by data**: The zeldovich RNG debug output contains exactly 121 coordinates for each Y value from Y=0 to Y=7. Y=8 is completely absent.

**Root cause**: The RNG debug output is printed during the main `ZeldovichZ` function (line 503), which processes Y values in blocks. Looking at line 628:

```cpp
for (int yblock = 0; yblock < array.numblock / 2; yblock++)
```

For N=16 with `ZD_NumBlock = 2`:
- `array.numblock = 2`, `array.block = 8`
- Loop processes: `yblock = 0` only (since `array.numblock / 2 = 1`)
- For `yblock = 0`: processes `yres = 0` to `7`
- Global Y = `yres + yblock * block = 0` to `7`

**Y=8 is never processed in the main loop**, so it never appears in the RNG debug output.

**Note on Y=8 zeroing**: The Y=8 zeroing in `ZeldovichXY` (line 718) is a **separate, later operation** that happens during the Y & X FFTs, **after** the RNG debug output has already been printed. It does not affect the RNG debug output. The missing Y=8 coordinates in the RNG debug output are **solely** because that Y slice is not processed in the main generation loop (`ZeldovichZ`), not because of the later zeroing operation in `ZeldovichXY`.

## Does Y=8 Zeroing Affect RNG Consistency or Alignment?

**Short answer: NO, we do not need to worry about it.**

**Detailed analysis:**

### RNG Calls for Y=8

**Hermitian code:**
- Processes Y=8: **YES** (121 coordinates logged)
- All Y=8 coordinates have `D=0` (Nyquist modes with `|ky|=8`)
- When `D=0`: `nskip++` (RNG call is skipped)
- **Total RNG calls for Y=8: 0** (all skipped)

**Zeldovich code:**
- Processes Y=8 in main loop: **NO**
- Y=8 is filled by conjugation from Y=0 (from `slabHer`)
- **Total RNG calls for Y=8: 0** (never processed)

### RNG Alignment Impact

Both codes end up with:
- **0 RNG calls for Y=8**
- **Same RNG state after Y=7**
- **Same RNG state before next Y value** (if there were more Y values to process)

**Conclusion**: Y=8 zeroing/processing does **NOT affect RNG alignment** because both codes effectively skip RNG calls for Y=8, just in different ways:
- **Hermitian**: Processes Y=8 but all D=0, so explicitly skips RNG calls
- **Zeldovich**: Doesn't process Y=8 at all, so no RNG calls are made

### Why This Doesn't Matter

1. **Y=8 is the Nyquist slice**: All modes have `|ky|=8`, so `D=0` by definition (Nyquist modes are zeroed)
2. **No RNG calls in either code**: Both codes end up with 0 RNG calls for Y=8
3. **RNG state remains aligned**: The RNG state after Y=7 is the same in both codes, and both skip Y=8, so the state remains aligned
4. **Missing Y=8 in comparison is expected**: Since zeldovich doesn't process Y=8 in the main loop, it won't appear in RNG debug output. This is expected and not a problem for RNG consistency.

**Bottom line**: The Y=8 zeroing in zeldovich is a post-processing step for FFT correctness and does not affect RNG consistency or alignment. The missing Y=8 coordinates in the comparison are expected and do not indicate an RNG alignment problem.

## The Mystery: Why Are (4,7), (5,7), (6,6), (6,7), (7,4), (7,5), (7,6), (7,7) Zeroed?

Based on the code analysis:
- These coordinates are **IN the mirroring loop** (first half)
- They are **NOT zeroed by Nyquist conditions** (assuming k_cutoff=1.0)
- They are **NOT explicitly zeroed** in the mirroring code
- They should have **non-zero D values** from main processing

**Possible explanations:**

1. ~~**k_cutoff > 1.0**~~: **RULED OUT** - Tests use k_cutoff=1.0, so this is not the cause.

2. **Implicit zeroing after mirroring**: There might be code after the mirroring that zeros coordinates in the first half, but it's not visible in lines 555-573. Need to check what happens after the mirroring loop.

3. **slabHer initialization**: If `slabHer` is not properly initialized for these coordinates, copying from `slabHer` might result in zero values.

4. **Memory initialization**: If `slab` is zero-initialized and these coordinates are never written to during main processing, they would remain zero.

5. **Different processing path**: These coordinates might follow a different code path that sets D=0.

## Next Steps

1. **Check the k_cutoff value** in the parameter file used for the test - if k_cutoff > 1.0, this would explain why coordinates with |kx|=7 or |kz|=7 are zeroed
2. **Check if slab/slabHer are zero-initialized** before main processing
3. **Check what happens AFTER the mirroring loop** - is there any code that zeros coordinates?
4. **Check the main processing loop** to see if coordinates with x=7 or z=7 are handled differently
5. **Compare with hermitian code** to see how it handles these coordinates and what k_cutoff value it uses

## Coordinates of Interest

The 8 differences in non-mirrored region at Y=0:
- (4,7), (5,7), (6,6), (6,7), (7,4), (7,5), (7,6), (7,7)

All of these:
- Are in the mirroring loop (first half)
- Should have non-zero D values
- Are being zeroed by zeldovich (but we don't see why in the code)

