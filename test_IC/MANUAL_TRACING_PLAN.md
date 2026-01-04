# Manual FFT Tracing Plan for N=4

## Goal
Manually trace where differences first appear in the FFT pipeline by comparing intermediate results at each step between Hermitian and Zeldovich codes.

## FFT Pipeline Comparison

### Zeldovich-PLT Pipeline:
1. Generate D, F, G, H (verified: same as Hermitian)
2. **1D FFT** (along Z-axis) → Intermediate state A
3. **2D FFT** (along Y, X axes) → Final result

### Hermitian Pipeline:
1. Generate D, F, G, H (verified: same as Zeldovich)
2. **2D FFT** (along Y, X axes) → Intermediate state B
3. **1D FFT** (along Z-axis) → Final result

## Key Comparison Points

### Critical Comparison: After Zeldovich's 2D FFT vs After Hermitian's 1D FFT

**These should be EQUIVALENT** if the FFTs are mathematically equivalent (just different order).

If these match:
- The FFTs are equivalent
- Differences must be in data reorganization/coordinate mapping between steps

If these differ:
- This is where the problem originates
- Need to check FFT implementation, normalization, or indexing

### Secondary Comparison: Intermediate States

After Zeldovich's 1D FFT vs After Hermitian's 2D FFT:
- These are DIFFERENT intermediate states (different FFT orders)
- They may not match directly, but should be related by coordinate transformation
- Useful for understanding data organization differences

## Recommended Approach

### Step 1: Set up N=4 Test Case

1. Create parameter file `param_N4_trace.par`:
   - CPD = 4 (N=4)
   - NP = 64 (4³)
   - Same RNG seed as current test
   - Same other parameters

2. Run both codes with N=4:
   ```bash
   # Hermitian
   ./hermitian_3d_matrix 4 test_IC/param_N4_trace.par
   
   # Zeldovich
   cd /home/helenshao/InitialConditions/zeldovich-PLT
   ./zeldovich 4 ../hermitian_3d_matrix_production/test_IC/param_N4_trace.par
   ```

### Step 2: Add Intermediate Matrix Dumps

**For Hermitian code (`src/main.cpp`):**
- Add dump after 2D FFT (before 1D FFT): `matrix_after_2d_fft.txt`
- Keep existing dump after 1D FFT (final): `matrix_after_fft.txt`

**For Zeldovich code:**
- Add dump after 1D FFT: `matrix_after_1d_fft.txt`
- Add dump after 2D FFT (before final): `matrix_after_2d_fft.txt`
- Keep existing dump after final FFT: `matrix_after_fft.txt`

### Step 3: Systematic Comparison

1. **Compare final results** (already done for N=8, verify for N=4)
2. **Compare intermediate states**:
   - `zeldovich_matrix_after_2d_fft.txt` vs `hermitian_matrix_after_1d_fft.txt`
   - These should match if FFTs are equivalent
3. **Check coordinate mapping**:
   - Verify same (Z, Y, X) coordinates have same values
   - If not, check for coordinate transformation
4. **Manual inspection**:
   - Pick 5-10 key coordinates: (0,0,0), (1,1,1), (2,2,2), (0,1,2), (3,3,3)
   - Trace each coordinate through both pipelines
   - Document where values first diverge

### Step 4: What to Look For

1. **FFT Normalization**: Check if values differ by a constant factor (N³, 1/N³, etc.)
2. **Coordinate Mapping**: Check if values match at shifted coordinates
3. **Array Indexing**: Verify [array][Y][X] vs [array][X][Y] ordering
4. **Hermitian Symmetry**: Check if symmetry constraints are applied differently
5. **Data Reorganization**: Check how data is reordered between FFT steps

## Expected Outcomes

### Scenario A: Intermediate states match
- FFTs are equivalent
- Problem is in coordinate mapping or data reorganization
- Focus on indexing/reordering code

### Scenario B: Intermediate states differ
- FFT implementations differ
- Check FFT normalization, plan creation, or library differences
- May need to align FFT implementations

### Scenario C: Partial matches
- Some coordinates match, others don't
- Suggests coordinate mapping issue
- Check how coordinates are transformed between steps

## Implementation Notes

### Where to Add Dumps in Hermitian Code

After 2D FFT (around line 1400-1500 in `main.cpp`):
- After `fftw_execute_dft_r2c` for 2D FFT
- Before all-to-all communication
- Dump in same format as `matrix_after_fft.txt`

### Where to Add Dumps in Zeldovich Code

Need to check zeldovich-PLT code structure:
- After 1D FFT completion
- After 2D FFT completion (before final processing)

## Next Steps

1. Create N=4 parameter file
2. Add intermediate dump code to both codes
3. Run both codes with N=4
4. Compare intermediate dumps systematically
5. Document findings

