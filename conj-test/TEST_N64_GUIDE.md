# Testing N=64 on zeldovich-conj Branch

## Current Status

- **Branch**: `zeldovich-conj`
- **Changes**: Zeldovich packing scheme implemented
- **Print function**: Updated to support N=64

## Quick Test Steps

### 1. Verify You're on the Branch
```bash
cd /home/helenshao/InitialConditions/hermitian_3d_matrix_production
git branch --show-current  # Should show: zeldovich-conj
```

### 2. Compile with N=64 Support
```bash
make clean
make CFLAGS="-DUSE_DOUBLE_PRECISION -DDEBUG_RNG_CONSISTENCY=1 -DSPLINE_RESOLUTION=512"
```

### 3. Run with N=64

**Option A: Use existing PBS script (modify N)**
```bash
# Edit your PBS script to set N=64
# Then submit:
qsub your_script.pbs
```

**Option B: Run directly (if on compute node)**
```bash
./your_executable --N 64 --param-file examples/param_template.par
```

### 4. What to Look For in Output

The print function will show:

1. **For each Z-slab:**
   - Summary statistics: min/max real and imaginary values
   - Count of non-zero imaginary parts
   - Sample rows (first 4 rows, first 8 X values)

2. **Key Verification:**
   - **Imaginary parts should be NON-ZERO** (unlike before)
   - Array 0: `re` = density, `im` = x-displacement
   - Array 1: `re` = y-displacement, `im` = z-displacement

3. **Expected Output Format:**
```
[RANK 0] Z-slab Z=0 (X=[0,64), after 3D FFT):
  Array 0:
    Summary: re=[min, max], im=[min, max], non-zero imag: count/total
    First 4 rows (Y=0-3):
      Y=0: X0=... X1=... ...
      ...
  Array 1:
    ...
```

## What Changed

### Before (main branch):
- Array 0: `D[0] + i*F[0]` (only real parts)
- Imaginary parts → ~0 after IFFT

### After (zeldovich-conj branch):
- Array 0: `(D[0] - F[1]) + i*(D[1] + F[0])` (full complex packing)
- Imaginary parts → **NON-ZERO** after IFFT (contains displacement data)

## Verification Checklist

- [ ] Compiles successfully
- [ ] Runs without errors
- [ ] Imaginary parts are non-zero (check summary statistics)
- [ ] Array 0.imag contains x-displacement data
- [ ] Array 1.real contains y-displacement data
- [ ] Array 1.imag contains z-displacement data

## Next Steps After Testing

If tests pass:
1. Compare with Zeldovich-PLT output (same parameters)
2. Verify normalization matches
3. Merge to main: `git checkout main && git merge zeldovich-conj`

If tests fail:
1. Check error messages
2. Verify compilation flags
3. Check parameter file settings

