#ifndef HERMITIAN_PRINTING_H
#define HERMITIAN_PRINTING_H

// ====================================================================================
// HERMITIAN 3D MATRIX MPI - DEBUG PRINTING UTILITIES
// ====================================================================================
// This file contains functions for printing debug information about slices and matrices.
// All functions respect PRINT_DETAILED_SLICES and PRINT_MATRICES configuration flags.
//
// Depends on: config.h, precision.h, types.h
// ====================================================================================

#include "config.h"
#include "precision.h"
#include "types.h"

// ====================================================================================
// FUNCTION DECLARATIONS
// ====================================================================================

// Print a Y-slice in Fourier space (after 2D FFT)
// Only prints for small N (N <= 16) and if PRINT_DETAILED_SLICES is enabled
// label: Descriptive label for the output (e.g., "Primary slice", "Conjugate slice")
void print_y_slice_fourier(int rank, int y_global, fftw_complex_t *slice, int N, int array_idx, const char *label);

// Print 3D matrix for visual inspection (small N only)
// Only prints if PRINT_MATRICES is enabled
// title: Title to print before the matrix
void print_3d_matrix_visual(int N, fftw_complex_t *global_matrix, const char* title);

#endif // HERMITIAN_PRINTING_H

