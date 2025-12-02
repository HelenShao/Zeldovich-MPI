// ====================================================================================
// HERMITIAN 3D MATRIX MPI - DEBUG PRINTING UTILITIES
// ====================================================================================

#include "utils/printing.h"
#include <stdio.h>

// ====================================================================================
// FUNCTION IMPLEMENTATIONS
// ====================================================================================

void print_y_slice_fourier(int rank, int y_global, fftw_complex_t *slice, int N, int array_idx, const char *label) {
    #if PRINT_DETAILED_SLICES
    if (N > 16) return;  // Only print for small N
    
    printf("\n[RANK %d] %s Y=%d, Array=%d (Fourier space, after 2D FFT):\n", rank, label, y_global, array_idx);
    for (int x = 0; x < N; x++) {
        for (int z = 0; z < N; z++) {
            double re = slice[x * N + z][0];
            double im = slice[x * N + z][1];
            if (fabs_t(re) > 1e-10 || fabs_t(im) > 1e-10) {  // Only print non-zero
                printf("  F(X=%d,Y=%d,Z=%d,Array=%d): %+.6e %+.6ei\n", x, y_global, z, array_idx, re, im);
            }
        }
    }
    #else
    (void)rank; (void)y_global; (void)slice; (void)N; (void)array_idx; (void)label;
    #endif
}

void print_3d_matrix_visual(int N, fftw_complex_t *global_matrix, const char* title) {
    printf("\n%s\n", title);
    for (int y = 0; y < N; y++) {
        printf("Y=%d\n", y);
        for (int x = 0; x < N; x++) {
            for (int z = 0; z < N; z++) {
                size_t idx = (size_t)y * (size_t)N * (size_t)N + (size_t)x * (size_t)N + (size_t)z;
                double re = global_matrix[idx][0];
                double im = global_matrix[idx][1];
                if (fabs_t(im) < 1e-10) {
                    printf("%7.3f ", re);
                } else {
                    printf("%6.2f%+5.2fi ", re, im);
                }
            }
            printf("\n");
        }
        printf("\n");
    }
}

void print_z_slab(int rank, int z, fftw_complex_t *local_z_slab, int x_count, int N, int narray, int x_start) {
    #if PRINT_Z_SLABS
    if (N > 16) return;  // Only print for small N
    
    printf("\n[RANK %d] Z-slab Z=%d (X=[%d,%d), after 3D FFT):\n", rank, z, x_start, x_start + x_count);
    // Print only first array (array_idx=0)
    int array_idx = 0;
    printf("  Array %d:\n", array_idx);
    for (int y = 0; y < N; y++) {
        printf("    Y=%d: ", y);
        for (int x_idx = 0; x_idx < x_count; x_idx++) {
            int x_global = x_start + x_idx;
            // Access using ZSLAB macro formula: y + N * (array_idx + narray * x_idx)
            int64_t idx = (int64_t)y + (int64_t)N * ((int64_t)array_idx + (int64_t)narray * (int64_t)x_idx);
            double re = local_z_slab[idx][0];
            double im = local_z_slab[idx][1];
            if (fabs_t(im) < 1e-10) {
                printf("X%d=%7.3f ", x_global, re);
            } else {
                printf("X%d=%6.2f%+5.2fi ", x_global, re, im);
            }
        }
        printf("\n");
    }
    printf("\n");
    #else
    (void)rank; (void)z; (void)local_z_slab; (void)x_count; (void)N; (void)narray; (void)x_start;
    #endif
}

