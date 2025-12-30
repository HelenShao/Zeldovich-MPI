# Hermitian Symmetry Tests - Post k2_cutoff Implementation

This directory contains test scripts to verify that the k2_cutoff condition fix preserves Hermitian symmetry.

## Background

After implementing the k2_cutoff condition fix (zeroing modes with k² >= k2_cutoff when CornerModes=0), we need to verify that:
1. Hermitian symmetry is still preserved
2. The k2_cutoff zeroing doesn't introduce any symmetry violations
3. Both test modes (Mode 1 and Mode 2) pass
4. The zeroing condition is correctly applied while maintaining symmetry

## Test Scripts

### `test_hermitian_mode1_plt.pbs`
- **Mode**: VERIFY_HERMITIAN_SYMMETRY=1
- **Test**: Sets F=0, H=0
- **Expected Result**: After 3D FFT, final matrix should be purely real (imaginary parts ≈ 0)
- **Purpose**: Tests that with F=0 and H=0, only D and G remain, and the result is purely real

### `test_hermitian_mode2_plt.pbs`
- **Mode**: VERIFY_HERMITIAN_SYMMETRY=2
- **Test**: Sets D=0, G=0
- **Expected Result**: After 3D FFT, final matrix should be purely imaginary (real parts ≈ 0)
- **Purpose**: Tests that with D=0 and G=0, F and H are computed correctly, and the result is purely imaginary

## Compilation Flags

Both scripts compile with:
- `-DUSE_DOUBLE_PRECISION`: Double precision floating point
- `-DVERIFY_HERMITIAN_SYMMETRY=1` or `=2`: Enables verification mode
- `-DPRINT_DETAILED_SLICES=1`: Prints detailed slice information
- `-DPRINT_Z_SLABS=1`: Prints selected Z-slabs for verification
- `-DSKIP_VERIFICATION=0`: Enables all verification checks
- `-DSPLINE_RESOLUTION=512`: High accuracy spline interpolation

## Test Parameters

- **Grid Size**: N=1024
- **Parameter File**: `examples/param_N1024_plt.par`
- **MPI Configuration**: 2 nodes, 16 ranks per node, 2 threads per rank (32 total ranks)
- **PLT**: ENABLED (using eigenmodes from zeldovich-PLT)
- **PLT Rescaling**: ENABLED (target_z = 0.0)

## Running the Tests

### Submit Mode 1 Test:
```bash
cd /home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_hermitian/k2_cutoff
qsub test_hermitian_mode1_plt.pbs
```

### Submit Mode 2 Test:
```bash
cd /home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_hermitian/k2_cutoff
qsub test_hermitian_mode2_plt.pbs
```

## Checking Results

After the jobs complete, check the output files:

### Log Files:
- `test_hermitian_mode1_plt.log`: Full output from Mode 1 test
- `test_hermitian_mode2_plt.log`: Full output from Mode 2 test

### Output Files:
- `test_hermitian_mode1_plt.out`: PBS job output (stdout/stderr)
- `test_hermitian_mode2_plt.out`: PBS job output (stdout/stderr)

### Key Things to Check:

1. **No Hermitian Symmetry Errors**:
   ```bash
   grep -i "ERROR.*hermitian\|ERROR.*symmetry" test_hermitian_mode1_plt.log
   ```
   Should return no results (or count should be 0)

2. **Verification Messages**:
   ```bash
   grep -i "hermitian\|symmetry\|imaginary\|real" test_hermitian_mode1_plt.log | tail -20
   ```

3. **Final Real-Space Verification**:
   ```bash
   grep -i "max_real\|max_imag\|rms_imag\|purely imaginary\|purely real\|NOT OK" test_hermitian_mode1_plt.log | tail -10
   ```

## Expected Results

### Mode 1 (F=0, H=0):
- ✅ No Hermitian symmetry errors
- ✅ Final matrix should be purely real (imaginary parts ≈ 0)
- ✅ Verification checks should all pass

### Mode 2 (D=0, G=0):
- ✅ No Hermitian symmetry errors
- ✅ Final matrix should be purely imaginary (real parts ≈ 0)
- ✅ Verification checks should all pass

## What Changed in k2_cutoff Implementation

The k2_cutoff condition fix added:
- **k2_cutoff Zeroing**: Modes with k² >= k2_cutoff are zeroed when CornerModes=0 (matching zeldovich-PLT behavior)
- **Combined Zeroing Conditions**: DC modes (k²=0), Nyquist frequencies, and k2_cutoff modes are all handled in a unified condition
- **RNG Skip Logic**: RNG state is correctly advanced when modes are zeroed, ensuring alignment with zeldovich-PLT

These changes should **preserve** Hermitian symmetry while matching zeldovich-PLT's k_cutoff behavior. The zeroing of high-k modes should not affect the Hermitian symmetry of the remaining modes.

## k2_cutoff Implementation Details

The k2_cutoff condition is implemented as:
```c
if ((k2 == 0.0) || (is_nyquist) || (!CornerModes && (double)k2_int >= k2_cutoff)) {
    D[0] = D[1] = 0.0;
    // RNG skip logic...
}
```

Where:
- `k2_cutoff = (N/2)² / (k_cutoff)²`
- For N=1024, k_cutoff=1.0: `k2_cutoff = 512² = 262144`
- Modes with k² >= 262144 are zeroed (matching zeldovich-PLT behavior)

## Comparison with Previous Tests

Compare results with:
- `/home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_hermitian/post_y0_fix/test_hermitian_mode1_plt.out`
- `/home/helenshao/InitialConditions/hermitian_3d_matrix_production/test_hermitian/post_y0_fix/test_hermitian_mode2_plt.out`

The results should be **identical** or **very similar** (within numerical precision), confirming that the k2_cutoff fix doesn't break Hermitian symmetry. The only difference should be that high-k modes (k² >= k2_cutoff) are now zeroed, which is expected and correct behavior.

