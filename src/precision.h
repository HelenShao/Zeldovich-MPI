#ifndef HERMITIAN_PRECISION_H
#define HERMITIAN_PRECISION_H

// ====================================================================================
// HERMITIAN 3D MATRIX MPI - PRECISION SELECTION
// ====================================================================================
// This file handles compile-time precision selection between single and double.
// It provides a unified interface regardless of precision choice.
//
// Usage:
//   Default: Single precision (float)
//   Double:  Compile with -DUSE_DOUBLE_PRECISION
//
// Memory impact:
//   Single: 8 bytes per complex number (2× float)
//   Double: 16 bytes per complex number (2× double)
//   For N=32K, narray=4: Single=~160 GB, Double=~320 GB
//
// Precision trade-off:
//   Float:  ~7 decimal digits
//   Double: ~15 decimal digits
// ====================================================================================

#include <fftw3.h>
#include <mpi.h>
#include <math.h>

#ifdef USE_DOUBLE_PRECISION
    // ====================================================================================
    // DOUBLE PRECISION MODE
    // ====================================================================================
    
    // Basic types
    typedef double real_t;
    typedef fftw_complex fftw_complex_t;
    typedef fftw_plan fftw_plan_t;
    
    // FFTW function mappings
    #define FFTW_PLAN_DFT_2D fftw_plan_dft_2d
    #define FFTW_PLAN_DFT_1D fftw_plan_dft_1d
    #define FFTW_EXECUTE_DFT fftw_execute_dft
    #define FFTW_DESTROY_PLAN fftw_destroy_plan
    
    // MPI datatypes
    #define MPI_COMPLEX_TYPE MPI_C_DOUBLE_COMPLEX
    #define MPI_REAL_TYPE MPI_DOUBLE
    
    // Math functions
    #define fabs_t fabs
    #define fmax_t fmax
    #define sqrt_t sqrt
    
    // Constants
    #define PRECISION_NAME "Double"
    #define BYTES_PER_COMPLEX 16
    
#else
    // ====================================================================================
    // SINGLE PRECISION MODE (DEFAULT)
    // ====================================================================================
    
    // Basic types
    typedef float real_t;
    typedef fftwf_complex fftw_complex_t;
    typedef fftwf_plan fftw_plan_t;
    
    // FFTW function mappings
    #define FFTW_PLAN_DFT_2D fftwf_plan_dft_2d
    #define FFTW_PLAN_DFT_1D fftwf_plan_dft_1d
    #define FFTW_EXECUTE_DFT fftwf_execute_dft
    #define FFTW_DESTROY_PLAN fftwf_destroy_plan
    
    // MPI datatypes
    #define MPI_COMPLEX_TYPE MPI_C_FLOAT_COMPLEX
    #define MPI_REAL_TYPE MPI_FLOAT
    
    // Math functions
    #define fabs_t fabsf
    #define fmax_t fmaxf
    #define sqrt_t sqrtf
    
    // Constants
    #define PRECISION_NAME "Single (float)"
    #define BYTES_PER_COMPLEX 8
    
#endif

// ====================================================================================
// PRECISION-INDEPENDENT UTILITIES
// ====================================================================================

// Get real and imaginary parts (works for both precisions)
#define CREAL(c) ((c)[0])
#define CIMAG(c) ((c)[1])

// Set complex number (works for both precisions)
#define CSET(c, re, im) do { (c)[0] = (re); (c)[1] = (im); } while(0)

// Complex conjugate (works for both precisions)
#define CCONJ(dest, src) do { (dest)[0] = (src)[0]; (dest)[1] = -(src)[1]; } while(0)

// Complex magnitude squared (works for both precisions)
#define CMAG2(c) ((c)[0] * (c)[0] + (c)[1] * (c)[1])

// Complex magnitude (works for both precisions)
#define CMAG(c) sqrt_t(CMAG2(c))

// ====================================================================================
// PRECISION INFORMATION (for debugging and reporting)
// ====================================================================================

// Print precision information
static inline void print_precision_info(int rank) {
    if (rank == 0) {
        printf("Precision: %s (%d bytes per complex)\n", 
               PRECISION_NAME, BYTES_PER_COMPLEX);
    }
}

// Get memory size for N³ matrix with narray arrays
static inline size_t get_matrix_memory_bytes(int N, int narray) {
    return (size_t)N * N * N * narray * BYTES_PER_COMPLEX;
}

// Get memory size in human-readable format
static inline void print_matrix_memory(int N, int narray, int rank) {
    if (rank == 0) {
        size_t bytes = get_matrix_memory_bytes(N, narray);
        double gb = bytes / (1024.0 * 1024.0 * 1024.0);
        printf("Full matrix memory: %.2f GB (%s precision)\n", gb, PRECISION_NAME);
    }
}

#endif // HERMITIAN_PRECISION_H

