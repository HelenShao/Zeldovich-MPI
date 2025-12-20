# Normalization Analysis: Zeldovich-PLT vs Hermitian 3D

## Why is `inv_ppd3` only in f_NL case?

## FFTW Behavior

**FFTW does NOT normalize:**
- Forward DFT: `F(k) = Σ f(x) exp(-2πikx/N)`  (no 1/N factor)
- Inverse DFT: `f(x) = Σ F(k) exp(+2πikx/N)`  (no 1/N factor)

**Round-trip:** `IFFT(FFT(f)) = N * f` (or `N³` for 3D)

## Standard Case (f_NL = 0) - ZeldovichXY()

### FFT Sequence:
1. **ZeldovichZ()**: `InverseFFT_Yonly()` - 1D iFFT along Y → multiplies by N
2. **ZeldovichXY()**: `Inverse2dFFT()` - 2D iFFT along Y,X → multiplies by N²
3. **Total scaling**: N × N² = **N³**

### Power Spectrum Normalization (line 226):
```cpp
normalization /= param.boxsize * param.boxsize * param.boxsize;
```

**This accounts for:**
1. Phase space volume: `(2π/L)³` per k-cell
2. **Implicitly accounts for FFTW N³ scaling** (so output has correct units)

### Output Function (output.cpp lines 60-65):
```cpp
// 2) The inverse FFT carries a pre-factor of 1/N
// 3) Our power spectrum convention requires a prefactor of sqrt(Volume)
//    We fixed these in the power spectrum.
```

**Note:** Comment says "1/N" but FFTW actually multiplies by N. The normalization in power spectrum accounts for this.

### Conclusion for Standard Case:
- Power spectrum normalization (BoxSize³) is for **phase space volume** `(2π/L)³`
- **Does NOT explicitly account for FFTW N³ scaling**
- **Question:** Is there a missing N³ normalization, or is it handled implicitly?

## Critical Observation:

The comment in power_spectrum.cpp (line 222-225) says:
```cpp
// Might need to normalize to the box volume.  This is appropriate
// if the iFFT is like FFTW, i.e., not dividing by N.
// This is because integrals over k bring in a phase space volume
// per cell of (2*pi/L)^3
normalization /= param.boxsize * param.boxsize * param.boxsize;
```

**This only mentions phase space volume, NOT FFTW N³ scaling!**

**Possible explanations:**
1. The BoxSize³ normalization is empirically chosen to give correct results (accounts for both implicitly)
2. There IS a missing N³ normalization, but outputs are still "correct" due to convention
3. The normalization is correct because BoxSize and N are related (BoxSize = N × grid_spacing)

**We need to test this empirically!**

## f_NL Case (f_NL ≠ 0) - ZeldovichXY_Phi()

### FFT Sequence:
1. **ZeldovichZ()**: Generate phi field → 1D iFFT along Y → N scaling
2. **ZeldovichXY_Phi()**: 
   - `Inverse2dFFT()` → N² scaling (total N³ so far)
   - Apply f_NL: `phi + f_NL * phi²`
   - **Apply `inv_ppd3`** → divide by N³
   - `Forward2dFFT()` → multiply by N²
   - **Net after forward FFT**: (N³ / N³) × N² = **N²**

3. **Why the extra normalization?**
   - Phi field goes through **extra round-trip** (iFFT → f_NL → FFT)
   - The forward FFT multiplies by N² again
   - Need to cancel the N³ from initial iFFTs before forward FFT
   - `inv_ppd3` corrects for the N³ scaling from the initial inverse FFTs

### Conclusion for f_NL Case:
- **Explicit `inv_ppd3` required** because of extra FFT round-trip
- Without it, phi field would be N³ times too large before forward FFT
- After forward FFT, would be N⁵ times too large

## Why This Matters

### Standard Case (f_NL = 0):
- Power spectrum normalization handles everything
- **No explicit N³ normalization needed in output**
- Both codes should match if using same PowerSpectrum class

### f_NL Case (f_NL ≠ 0):
- Requires explicit `inv_ppd3` for phi field
- **my code doesn't support f_NL yet**, so this doesn't apply

## Verification

**For my param_template.par (f_NL = 0):**
- v Both use same PowerSpectrum class
- v Both have BoxSize³ normalization in P(k)
- v Both use FFTW (no normalization)
- v **No explicit inv_ppd3 needed**

**Expected:** Outputs should match (within RNG consistency ~99.79%)

## Remaining Question

**Why does the comment in output.cpp say "The inverse FFT carries a pre-factor of 1/N"?**

This comment is **incorrect** - FFTW multiplies by N, not divides. But the key point is:
- The normalization is "fixed in the power spectrum"
- The BoxSize³ division accounts for both phase space volume AND FFTW scaling
- So the output has correct physical units without explicit N³ normalization

## To-do

1. **[?] Uncertain**: Whether standard case needs explicit N³ normalization
2. **[v] Confirmed**: f_NL case needs inv_ppd3 (but my code doesn't support f_NL)
3. **[ ] Run comparison test** to verify outputs match
4. **[ ] If outputs differ by N³ factor**: Add explicit `inv_ppd3` normalization
5. **[ ] Check field RMS ratio**: If ratio ≈ N³, normalization is missing

## Recommended Test

Run both codes with same parameters and compare:
```python
ratio = rms_hermitian / rms_zeldovich
if abs(ratio - 1.0) < 0.01:
    print("Normalization matches!")
elif abs(ratio - N**3) < 0.01:
    print("Missing N³ normalization - add inv_ppd3")
else:
    print(f"Unexpected ratio: {ratio}")
```

