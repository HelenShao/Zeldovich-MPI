# Step-by-Step Guide: Testing RNG Consistency Across Different N Values

This guide explains how to verify that the code generates the same random values for shared modes when run with different N values.

## Prerequisites

1. **Debug mode is already enabled**: `DEBUG_RNG_CONSISTENCY=1` is set in `config.h` (line 101)
2. **Test coordinates**: The code prints debug info for coordinates where `x, y, z <= MAX_DEBUG_COORD` (default: 10)
3. **Comparison script**: `compare_rng_consistency.py` is available in this directory

## Step 1: Compile the Code

Make sure the code is compiled (debug mode is already enabled by default):

```bash
cd /home/helenshao/InitialConditions/hermitian_3d_matrix_production
make clean
make
```

Or if you want to explicitly enable it:

```bash
make clean
make CFLAGS="-DUSE_DOUBLE_PRECISION -DDEBUG_RNG_CONSISTENCY=1"
```

## Step 2: Run Tests with Different N Values

Run the code with two different N values. The debug output goes to **stderr**, so redirect both stdout and stderr:

### Example: Test N=256 vs N=512

```bash
# Run with N=256
mpirun -np 4 ./hermitian_3d_matrix 256 examples/param_template.par > output_N256.out 2>&1

# Run with N=512  
mpirun -np 4 ./hermitian_3d_matrix 512 examples/param_template.par > output_N512.out 2>&1
```

### Example: Test N=4 vs N=6 (smaller, faster test)

```bash
# Run with N=4
mpirun -np 2 ./hermitian_3d_matrix 4 examples/param_template.par > output_N4.out 2>&1

# Run with N=6
mpirun -np 2 ./hermitian_3d_matrix 6 examples/param_template.par > output_N6.out 2>&1
```

**Note**: The command-line format is: `./hermitian_3d_matrix N [param_file.par]`
- `N`: Grid size (must be positive even integer)
- `param_file.par`: Optional parameter file for power spectrum mode

## Step 3: Compare the Results

Use the comparison script to analyze overlapping coordinates:

```bash
cd v15_3_tests
python compare_rng_consistency.py output_N256.out output_N512.out
```

Or for the smaller test:

```bash
python compare_rng_consistency.py output_N4.out output_N6.out
```

## Step 4: Interpret the Results

### What Gets Compared

The script compares **D, F, G, H** values (complex numbers) for overlapping coordinates:
- **D**: Density field
- **F**: X-displacement field  
- **G**: Y-displacement field
- **H**: Z-displacement field

### Expected Behavior

1. **For k² = 0 (DC mode)**: D, F, G, H should **match exactly** because:
   - Same k-vector: (0,0,0)
   - Same power spectrum: P(k=0)
   - Same random numbers used

2. **For k² ≠ 0**: D, F, G, H may **differ** because:
   - Different k-vectors: Different N means different k-space sampling
   - Different power spectrum: P(k) evaluated at different k values
   - **However**: The underlying random numbers are still consistent (same RNG sequence)

### Understanding the Output

The comparison script will show:
- Number of overlapping coordinates found
- For each coordinate:
  - k-vectors (which will differ for different N)
  - k² values (which will differ)
  - D, F, G, H values and whether they match (within tolerance 1e-10)

**Key Point**: Even if D, F, G, H differ for k² ≠ 0, this is **expected** because the power spectrum P(k) depends on k, and k differs for different N. The RNG consistency ensures the same random numbers are used, but the final values differ due to k-dependence.

### Example Output

Here's what you should see when running the comparison script:

```
================================================================================
RNG Consistency Comparison
================================================================================

File 1 (output_N512.out): 3459 debug entries
File 2 (output_N256.out): 3431 debug entries
Overlapping coordinates: 3377

Comparing overlapping coordinates:
--------------------------------------------------------------------------------

Coordinate (x,y,z) = (0, 0, 0):
  N1=512, N2=256
  k1=(0,0,0), k2=(0,0,0)
  k2_1=0.000000, k2_2=0.000000
  D: MATCH: re_diff=0.00e+00 im_diff=0.00e+00
    D1 (N=512) = (0, 0)
    D2 (N=256) = (0, 0)
  F: MATCH: re_diff=0.00e+00 im_diff=0.00e+00
    F1 (N=512) = (0, 0)
    F2 (N=256) = (0, 0)
  G: MATCH: re_diff=0.00e+00 im_diff=0.00e+00
    G1 (N=512) = (0, 0)
    G2 (N=256) = (0, 0)
  H: MATCH: re_diff=0.00e+00 im_diff=0.00e+00
    H1 (N=512) = (0, 0)
    H2 (N=256) = (0, 0)

Note: The DC mode (k=0) is explicitly set to zero in the code, so all fields (D, F, G, H) are zero for this coordinate.

Coordinate (x,y,z) = (0, 0, 1):
  N1=512, N2=256
  k1=(0,0,1), k2=(0,0,1)
  k2_1=1.000000, k2_2=1.000000
  D: MATCH: re_diff=0.00e+00 im_diff=0.00e+00
    D1 (N=512) = (-0.0215, -0.0891)
    D2 (N=256) = (-0.0215, -0.0891)
  F: MATCH: re_diff=0.00e+00 im_diff=0.00e+00
    F1 (N=512) = (0, 0)
    F2 (N=256) = (0, 0)
  G: MATCH: re_diff=0.00e+00 im_diff=0.00e+00
    G1 (N=512) = (0, 0)
    G2 (N=256) = (0, 0)
  H: MATCH: re_diff=0.00e+00 im_diff=0.00e+00
    H1 (N=512) = (0.0891, -0.0215)
    H2 (N=256) = (0.0891, -0.0215)

Coordinate (x,y,z) = (0, 0, 2):
  N1=512, N2=256
  k1=(0,0,2), k2=(0,0,2)
  k2_1=4.000000, k2_2=4.000000
  D: MATCH: re_diff=0.00e+00 im_diff=0.00e+00
    D1 (N=512) = (-0.0829, -0.0903)
    D2 (N=256) = (-0.0829, -0.0903)
  F: MATCH: re_diff=0.00e+00 im_diff=0.00e+00
    F1 (N=512) = (0, 0)
    F2 (N=256) = (0, 0)
  G: MATCH: re_diff=0.00e+00 im_diff=0.00e+00
    G1 (N=512) = (0, 0)
    G2 (N=256) = (0, 0)
  H: MATCH: re_diff=0.00e+00 im_diff=0.00e+00
    H1 (N=512) = (0.0452, -0.0414)
    H2 (N=256) = (0.0452, -0.0414)

Coordinate (x,y,z) = (1, 1, 0):
  N1=512, N2=256
  k1=(1,1,0), k2=(1,1,0)
  k2_1=2.000000, k2_2=2.000000
  D: MATCH: re_diff=0.00e+00 im_diff=0.00e+00
    D1 (N=512) = (0.0422, 0.0407)
    D2 (N=256) = (0.0422, 0.0407)
  F: MATCH: re_diff=0.00e+00 im_diff=0.00e+00
    F1 (N=512) = (0, 0)
    F2 (N=256) = (0, 0)
  G: MATCH: re_diff=0.00e+00 im_diff=0.00e+00
    G1 (N=512) = (0, 0)
    G2 (N=256) = (0, 0)
  H: MATCH: re_diff=0.00e+00 im_diff=0.00e+00
    H1 (N=512) = (-0.0136, 0.0141)
    H2 (N=256) = (-0.0136, 0.0141)

[... more coordinates ...]

================================================================================
RESULT: Some values differ. This may be expected if:
  - k-vectors differ (different N means different k-space sampling)
  - Power spectrum P(k) is evaluated at different k values
  - RNG consistency ensures same random numbers, but final values differ due to k-dependence
================================================================================
```

**Key observations from the example:**

1. **Header**: Shows number of debug entries found in each file and how many overlapping coordinates were compared.

2. **For each coordinate**:
   - Shows the coordinate `(x,y,z)` being compared
   - Shows `N1` and `N2` values
   - Shows k-vectors `k1` and `k2` (may be the same or different)
   - Shows `k²` values (may be the same or different)
   - For each field (D, F, G, H):
     - **MATCH**: Values match within tolerance (1e-10) - shows `re_diff=0.00e+00 im_diff=0.00e+00`
     - **DIFFER**: Values differ - shows the difference magnitude
     - Shows the actual values from both runs

3. **Expected results**:
   - **k² = 0 (DC mode)**: All fields should **MATCH** exactly
   - **Same k-vector, same k²**: Fields should **MATCH** (same mode, same power spectrum)
   - **Different k-vector or k²**: Fields may **DIFFER** (different mode, different power spectrum evaluation)

4. **Final result**: The script summarizes whether all values matched or some differed, with an explanation of why differences may be expected.

## Step 5: Verify Debug Output is Generated

If you don't see `[RNG-DEBUG]` lines in the output:

1. **Check compilation**: Make sure `DEBUG_RNG_CONSISTENCY=1` is set
   ```bash
   grep "DEBUG_RNG_CONSISTENCY" src/config.h
   ```

2. **Check stderr**: Debug lines go to stderr, make sure you're capturing it:
   ```bash
   # Correct: captures both stdout and stderr
   command > output.out 2>&1
   
   # Wrong: only captures stdout
   command > output.out
   ```

3. **Check coordinates**: The code only prints for coordinates where `x, y, z <= MAX_DEBUG_COORD` (default: 10). For very small N, make sure your test coordinates are within this range.

## Example: Quick Test Script

Here's a complete test script you can run:

```bash
#!/bin/bash
# test_rng_consistency.sh

cd /home/helenshao/InitialConditions/hermitian_3d_matrix_production

# Compile
make clean
make

# Run tests
echo "Running N=4 test..."
mpirun -np 2 ./hermitian_3d_matrix 4 examples/param_template.par > v15_3_tests/output_N4.out 2>&1

echo "Running N=6 test..."
mpirun -np 2 ./hermitian_3d_matrix 6 examples/param_template.par > v15_3_tests/output_N6.out 2>&1

# Compare
echo "Comparing results..."
cd v15_3_tests
python compare_rng_consistency.py output_N4.out output_N6.out
```

## Troubleshooting

### No overlapping coordinates found

- **Cause**: The two N values don't share any coordinates in the test range
- **Solution**: Use N values that are both >= 2*MAX_DEBUG_COORD (default: 20). For example, N=4 and N=6 both have coordinates (0,0,0), (1,0,0), etc.

### All values differ

- **Cause**: This is expected for k² ≠ 0 modes (see "Expected Behavior" above)
- **Solution**: Focus on k² = 0 (DC mode) - these should match exactly. Also check that the k-vectors differ as expected.

### Parse errors in comparison script

- **Cause**: Output format may have changed or lines got wrapped
- **Solution**: Check the raw output files for `[RNG-DEBUG]` lines and verify the format matches the expected pattern

## Advanced: Testing with Larger N

For larger N values (e.g., N=256, N=512), the code uses **sampling** to reduce output volume:
- Always prints coordinates where `x, y, z <= MAX_DEBUG_COORD` (small cube)
- Samples other coordinates with stride `DEBUG_SAMPLE_STRIDE`
- Always prints boundary coordinates at `(N/2)-1`

This ensures you still get overlapping coordinates to compare, even for large N values.

