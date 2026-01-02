# Implementation Notes for Matrix Dump Feature

## Overview

To enable brute force tracing, both `hermitian_3d_matrix_production` and `zeldovich-PLT` need to be modified to dump matrix values at two critical points:

1. **Before FFT**: After generation in Fourier space (D, F, G, H arrays)
2. **After FFT**: After 3D FFT in real space (final displacement/velocity fields)

## Required Code Modifications

### For hermitian_3d_matrix_production

#### 1. Add compile flags
In the workflow, compile with:
```bash
make CFLAGS="-DUSE_DOUBLE_PRECISION -DDUMP_MATRIX_BEFORE_FFT=1 -DDUMP_MATRIX_AFTER_FFT=1"
```

#### 2. Add matrix dump code

**Location: `src/generation/hermitian_generation.c`**

After generating D, F, G, H for each Y-slice (before 2D FFT), add:
```c
#ifdef DUMP_MATRIX_BEFORE_FFT
    if (N <= 16) {  // Only for small N
        FILE *dump_fp = fopen("matrix_before_fft.txt", "a");
        if (dump_fp) {
            for (int x = 0; x < N; x++) {
                for (int z = 0; z < N; z++) {
                    for (int array_idx = 0; array_idx < narray; array_idx++) {
                        fftw_complex_t *array_ptr = GET_ARRAY_PTR(primary_slices, array_idx, x, z, N);
                        fprintf(dump_fp, "Y=%d X=%d Z=%d Array=%d Re=%.15e Im=%.15e\n",
                                global_y, x, z, array_idx, array_ptr[0], array_ptr[1]);
                    }
                }
            }
            fclose(dump_fp);
        }
    }
#endif
```

**Location: `src/main.cpp`**

After 3D FFT (in the Z-slab processing loop), add:
```c
#ifdef DUMP_MATRIX_AFTER_FFT
    if (N <= 16) {  // Only for small N
        FILE *dump_fp = fopen("matrix_after_fft.txt", "a");
        if (dump_fp) {
            for (int y = 0; y < N; y++) {
                for (int x_idx = 0; x_idx < x_count; x_idx++) {
                    int x_global = x_start + x_idx;
                    for (int array_idx = 0; array_idx < narray; array_idx++) {
                        int64_t idx = (int64_t)y + (int64_t)N * ((int64_t)array_idx + (int64_t)narray * (int64_t)x_idx);
                        fprintf(dump_fp, "Z=%d Y=%d X=%d Array=%d Re=%.15e Im=%.15e\n",
                                z, y, x_global, array_idx, local_z_slab[idx][0], local_z_slab[idx][1]);
                    }
                }
            }
            fclose(dump_fp);
        }
    }
#endif
```

### For zeldovich-PLT

Similar modifications need to be made to `zeldovich-PLT/src/zeldovich.cpp`:

1. After generating D, F, G, H (before FFT), dump to `matrix_before_fft.txt`
2. After 3D FFT, dump to `matrix_after_fft.txt`

The format should match the hermitian code format for easy comparison.

## File Format

### Before FFT (Fourier space)
```
# Matrix values before FFT (Fourier space)
# Format: Y=<y> X=<x> Z=<z> Array=<array_idx> Re=<real> Im=<imag>
Y=0 X=0 Z=0 Array=0 Re=1.234567e+00 Im=0.000000e+00
Y=0 X=0 Z=0 Array=1 Re=2.345678e+00 Im=1.111111e-01
...
```

### After FFT (Real space)
```
# Matrix values after FFT (Real space)
# Format: Z=<z> Y=<y> X=<x> Array=<array_idx> Re=<real> Im=<imag>
Z=0 Y=0 X=0 Array=0 Re=1.234567e+00 Im=0.000000e+00
Z=0 Y=0 X=0 Array=1 Re=2.345678e+00 Im=1.111111e-01
...
```

## Notes

- Matrix dumps should only be enabled for small N (N <= 16) to avoid huge files
- Use append mode (`"a"`) when opening files to handle multiple ranks/slices
- Ensure files are opened/closed properly to avoid data loss
- Consider using rank-specific filenames if running with MPI (e.g., `matrix_before_fft_rank0.txt`)
- For serial runs (1 rank), a single file is sufficient

