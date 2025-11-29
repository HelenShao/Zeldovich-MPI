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

