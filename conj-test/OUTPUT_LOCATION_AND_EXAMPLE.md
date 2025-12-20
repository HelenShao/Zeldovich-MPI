# Output Location and Example Format

## File Writes: DISABLED ✓

**Confirmed:** File writes are **disabled** for testing.

- `SKIP_FILE_WRITE` is defined to `1` in `config.h` (line 303)
- Binary files (`rank_*/z*_slab_N*.bin`) will **NOT** be written
- The code section that would write files is wrapped in `#ifndef SKIP_FILE_WRITE` and will be skipped

## Where Outputs Are Saved

### 1. Print Statements (stdout)
All `printf()` statements (including `print_z_slab()`) go to **stdout**, which is captured by:

- **PBS scripts**: Saved to the `.out` file specified in `#PBS -o`
- **Direct run**: Printed to terminal

### 2. Example PBS Output File Location

If using a PBS script like `test_v15_3_powerspectrum.pbs`:
```bash
# Output will be in:
/home/helenshao/InitialConditions/hermitian_3d_matrix_production/v15_3_tests/test_v15_3_powerspectrum_N64.out
```

The exact filename depends on your PBS script's `#PBS -o` directive.

## Example Output Format

### For N=64, the output will look like:

```
[RANK 0] Z-slab Z=0 (X=[0,64), after 3D FFT):
  Array 0:
    Summary: re=[-2.345678e+01, 2.123456e+01], im=[-1.234567e+00, 1.876543e+00], non-zero imag: 4096/4096
    First 4 rows (Y=0-3):
      Y=0: X0= 12.34+ 0.56i X1= -5.67+ 1.23i X2=  8.90- 0.12i X3= -3.45+ 2.34i X4=  6.78+ 0.90i X5= -1.23- 0.45i X6=  9.87+ 1.11i X7= -4.56+ 0.78i ...
      Y=1: X0= -8.90+ 0.34i X1=  3.21- 1.56i X2= -7.65+ 2.10i X3=  1.23+ 0.45i X4= -9.87- 0.67i X5=  5.43+ 1.89i X6= -2.34+ 0.12i X7=  6.78- 1.23i ...
      Y=2: X0=  4.56- 0.78i X1= -1.23+ 2.34i X2=  7.89+ 0.56i X3= -3.45- 1.12i X4=  9.01+ 0.34i X5= -6.78+ 1.67i X6=  2.34- 0.45i X7= -8.90+ 0.89i ...
      Y=3: X0= -7.65+ 1.23i X1=  2.34- 0.56i X2= -9.87+ 0.78i X3=  4.56+ 1.45i X4= -1.23- 0.67i X5=  8.90+ 2.10i X6= -5.67+ 0.12i X7=  3.21- 1.34i ...
    ... (showing first 4 rows, first 8 X values)

  Array 1:
    Summary: re=[-1.234567e+00, 1.876543e+00], im=[-0.567890e+00, 0.987654e+00], non-zero imag: 4096/4096
    First 4 rows (Y=0-3):
      Y=0: X0=  0.56+ 0.12i X1= -0.34+ 0.45i X2=  0.78- 0.23i X3= -0.12+ 0.67i X4=  0.45+ 0.34i X5= -0.89- 0.12i X6=  0.23+ 0.56i X7= -0.67+ 0.78i ...
      Y=1: X0= -0.78+ 0.34i X1=  0.12- 0.56i X2= -0.45+ 0.89i X3=  0.67+ 0.23i X4= -0.34- 0.45i X5=  0.89+ 0.12i X6= -0.56+ 0.78i X7=  0.23- 0.34i ...
      Y=2: X0=  0.45- 0.78i X1= -0.23+ 0.67i X2=  0.89+ 0.34i X3= -0.12- 0.56i X4=  0.67+ 0.23i X5= -0.34+ 0.89i X6=  0.78- 0.12i X7= -0.45+ 0.56i ...
      Y=3: X0= -0.67+ 0.56i X1=  0.34- 0.23i X2= -0.89+ 0.78i X3=  0.12+ 0.45i X4= -0.56- 0.34i X5=  0.78+ 0.89i X6= -0.23+ 0.12i X7=  0.45- 0.67i ...
    ... (showing first 4 rows, first 8 X values)

  Array 2:
    Summary: re=[-0.123456e+00, 0.234567e+00], im=[-0.056789e+00, 0.098765e+00], non-zero imag: 4096/4096
    First 4 rows (Y=0-3):
      Y=0: X0=  0.12+ 0.06i X1= -0.03+ 0.04i X2=  0.08- 0.02i X3= -0.01+ 0.07i X4=  0.05+ 0.03i X5= -0.09- 0.01i X6=  0.02+ 0.06i X7= -0.07+ 0.08i ...
      ...

  Array 3:
    Summary: re=[-0.098765e+00, 0.187654e+00], im=[-0.045678e+00, 0.087654e+00], non-zero imag: 4096/4096
    First 4 rows (Y=0-3):
      Y=0: X0=  0.05+ 0.01i X1= -0.03+ 0.05i X2=  0.08- 0.02i X3= -0.01+ 0.07i X4=  0.04+ 0.03i X5= -0.09- 0.01i X6=  0.02+ 0.06i X7= -0.07+ 0.08i ...
      ...


[RANK 0] Z-slab Z=1 (X=[0,64), after 3D FFT):
  Array 0:
    Summary: re=[...], im=[...], non-zero imag: 4096/4096
    ...
```

## Key Things to Look For

### 1. **Non-Zero Imaginary Parts** (Critical!)
```
non-zero imag: 4096/4096  ← Should be close to 100% (all values have non-zero imag)
```

**Before (main branch):** `non-zero imag: 0/4096` (all imaginary parts ~0)  
**After (zeldovich-conj branch):** `non-zero imag: 4096/4096` (all imaginary parts have data)

### 2. **Array 0 (Density + X-displacement)**
- `re` = density field
- `im` = x-displacement field
- Both should have significant values

### 3. **Array 1 (Y-displacement + Z-displacement)**
- `re` = y-displacement field
- `im` = z-displacement field
- Both should have significant values

### 4. **Summary Statistics**
- `re=[min, max]`: Range of real parts
- `im=[min, max]`: Range of imaginary parts
- `non-zero imag: count/total`: How many values have non-zero imaginary parts

## How to View Output

### If using PBS:
```bash
# After job completes, view output file:
cat /home/helenshao/InitialConditions/hermitian_3d_matrix_production/v15_3_tests/test_v15_3_powerspectrum_N64.out

# Or search for Z-slab output:
grep "Z-slab" /path/to/output.out

# Or search for summary statistics:
grep "Summary:" /path/to/output.out
```

### If running directly:
Output will print to terminal in real-time.

## Verification Checklist

- [ ] No binary files created (check: `ls rank_*/z*_slab_N*.bin` should fail or show no files)
- [ ] Output file contains Z-slab print statements
- [ ] `non-zero imag` count is high (close to total, e.g., 4096/4096 for N=64)
- [ ] Array 0.imag has significant values (not ~0)
- [ ] Array 1.real and Array 1.imag have significant values

