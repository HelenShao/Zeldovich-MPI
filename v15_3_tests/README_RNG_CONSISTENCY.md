# RNG Consistency Verification

This document explains how to verify that overlapping grid points use the same random numbers across different matrix sizes (N values).

## Overview

When you run the code with different `N` values (e.g., N=4 and N=6), the RNG consistency mechanism ensures that:
- The same grid point (x,y,z) uses the same random number sequence
- However, the final output values will differ because:
  - Different k-vectors (wavenumbers) for the same coordinate
  - Different power spectrum evaluations P(k)
  - Different FFT operations

## How to Enable Debug Output

1. **Compile with debug flag:**
   ```bash
   make clean
   make CFLAGS="-DUSE_DOUBLE_PRECISION -DDEBUG_RNG_CONSISTENCY=1"
   ```

2. **Run your tests:**
   ```bash
   # Run with N=4
   mpirun -np 2 ./hermitian_3d_matrix -N 4 examples/param_template.par > output_N4.out 2>&1
   
   # Run with N=6
   mpirun -np 2 ./hermitian_3d_matrix -N 6 examples/param_template.par > output_N6.out 2>&1
   ```

3. **Compare the debug output:**
   ```bash
   python compare_rng_consistency.py output_N4.out output_N6.out
   ```

## What Gets Printed

The debug output prints `[RNG-DEBUG]` lines for test coordinates:
- (x,y,z) = (0,0,0), (1,0,0), (2,0,0), (3,0,0)
- For each coordinate, it prints:
  - N, Y, (x,z) coordinates
  - k-vector (kx, ky, kz) and k²
  - Raw values: D, F, G, H (complex numbers)

## Understanding the Results

### Expected Behavior

For overlapping coordinates (e.g., (0,0,0), (1,0,0), (2,0,0), (3,0,0)):

1. **k-vectors will differ**: Different N means different k-space sampling
   - N=4: kx ∈ {-2, -1, 0, 1, 2}
   - N=6: kx ∈ {-3, -2, -1, 0, 1, 2, 3}

2. **D, F, G, H may differ**: Even though the same random numbers are used, the power spectrum P(k) is evaluated at different k values, so the final values differ.

3. **RNG sequence is consistent**: The random number generator advances in the same way for the same (x,y,z) coordinate, regardless of N.

### What to Look For

The comparison script will show:
- Whether D, F, G, H match (within tolerance)
- If they differ, the difference values
- The k-vectors and k² values (which should differ)

**Note**: If D, F, G, H differ, this is expected because:
- The power spectrum P(k) depends on k, and k differs for different N
- The final output is P(k) × random_number, so different k → different output

## Example Output

```
[RNG-DEBUG] N=4 Y=0 (x,z)=(0,0): k=(0,0,0) k2=0.000000 | D=(1.285e-01,7.973e-02) F=(0.000e+00,0.000e+00) G=(0.000e+00,0.000e+00) H=(0.000e+00,0.000e+00)
[RNG-DEBUG] N=6 Y=0 (x,z)=(0,0): k=(0,0,0) k2=0.000000 | D=(1.285e-01,7.973e-02) F=(0.000e+00,0.000e+00) G=(0.000e+00,0.000e+00) H=(0.000e+00,0.000e+00)
```

For k²=0 (DC mode), D should match exactly because P(k=0) is the same.

For k²≠0, D may differ because P(k) is evaluated at different k values.

## Troubleshooting

If you don't see `[RNG-DEBUG]` lines:
1. Make sure you compiled with `-DDEBUG_RNG_CONSISTENCY=1`
2. Check that the test coordinates (x=0,1,2,3, y=0, z=0) are being processed
3. Check stderr output (debug lines go to stderr)

