# FFT Difference Analysis Summary

## Key Findings

### 1. **Differences Amplify After FFT**

The most critical finding is that **small differences before FFT become large differences after FFT**. 

**Example (Y=1, X=5, Z=2, Array=0):**
- **Before FFT**: Hermitian has Re=0.107, Im=-0.205; Zeldovich has Re=0.0, Im=0.0
- **Difference before**: Re=0.107, Im=0.205
- **After FFT**: Hermitian has Re=22.96, Im=0.373; Zeldovich has Re=-0.143, Im=-2.524
- **Difference after**: Re=23.11, Im=2.90
- **Amplification factor**: ~216x for real part, ~14x for imaginary part

This pattern is consistent across many locations - small pre-FFT differences (often where Hermitian has small non-zero values but Zeldovich has zero) explode after the FFT.

### 2. **No Systematic Scaling Factor**

- Mean ratio (hermitian/zeldovich): ~0.53 (real), ~0.78 (imag)
- Median ratio: ~0.17 (real), ~0.06 (imag)
- **High standard deviation** (80.9 for real, 17.5 for imag) indicates **no consistent scaling relationship**

This suggests the differences are not due to a simple normalization factor that could be corrected.

### 3. **Hotspot Location: Y=1, X=5, Z=2**

This coordinate appears in the top differences for **both arrays**:
- Array 0: Max real difference = 23.11
- Array 1: Max real difference = 22.14

This suggests a systematic issue at this specific location, possibly related to:
- FFT indexing/ordering differences
- Boundary condition handling
- Mode ordering differences

### 4. **Array-Specific Patterns**

- **Array 0**: Mean diff Re=5.36, Im=4.46
- **Array 1**: Mean diff Re=6.30, Im=6.19

Array 1 shows consistently larger differences, suggesting the issue may be more pronounced for one of the arrays (possibly G or H vs F).

### 5. **Statistical Summary**

- **Total elements**: 1024 (all match in coordinate space)
- **Mean difference**: Re=5.83, Im=5.33
- **Max difference**: Re=23.11, Im=21.56
- **Standard deviation**: Re=4.35, Im=3.87

The large standard deviation relative to the mean indicates **highly variable differences** across the grid, not a uniform offset.

## Root Cause Hypotheses

### Hypothesis 1: FFT Normalization Difference
**Most Likely**: The two codes may use different FFT normalization factors:
- FFTW typically uses unnormalized forward FFTs (multiply by N³)
- Some implementations normalize by 1/N³
- This would cause systematic scaling differences

**Test**: Check FFT normalization in both codes.

### Hypothesis 2: Pre-FFT Value Differences Propagate
**Observed**: Small differences before FFT (where Hermitian has small values but Zeldovich has zero) amplify dramatically after FFT.

**Possible causes**:
- Different `k2_cutoff` implementation (though both should zero modes with k² >= k2_cutoff)
- Different handling of edge cases
- Different precision in zero comparisons

**Test**: Verify that both codes are zeroing the same modes before FFT.

### Hypothesis 3: FFT Library/Implementation Differences
**Possible**: Different FFT libraries or FFT execution order could cause differences:
- FFTW vs other libraries
- Different FFT plans
- Different memory layouts affecting FFT

**Test**: Compare FFT library versions and FFT plan creation.

### Hypothesis 4: Indexing/Ordering Differences
**Observed**: Hotspot at Y=1, X=5, Z=2 suggests possible indexing mismatch.

**Possible causes**:
- Different array ordering (row-major vs column-major)
- Different FFT axis ordering
- Different coordinate system conventions

**Test**: Verify coordinate mapping between codes.

## Recommendations

1. **Check FFT Normalization**: Compare how both codes normalize their FFTs (before/after scaling factors)

2. **Verify k2_cutoff Implementation**: Ensure both codes are zeroing exactly the same modes before FFT

3. **Compare FFT Libraries**: Check FFTW versions and FFT plan parameters

4. **Investigate Hotspot**: Deep dive into location Y=1, X=5, Z=2 to understand why it's particularly problematic

5. **Check Array Ordering**: Verify that both codes use the same array ordering convention (F, G, H mapping)

6. **Compare Before FFT Values**: Since differences start before FFT, focus on understanding why Hermitian has small non-zero values where Zeldovich has zero

## Next Steps

1. Compare FFT normalization factors in both codes
2. Check if there's a systematic relationship between before-FFT differences and after-FFT differences
3. Investigate the specific hotspot location (Y=1, X=5, Z=2)
4. Verify that both codes are using identical FFT parameters and libraries

