#ifndef FFT_SETUP_H
#define FFT_SETUP_H

// ====================================================================================
// HERMITIAN 3D MATRIX MPI - FFT SETUP MODULE
// ====================================================================================
// This module handles FFTW plan creation for 2D and 1D FFT operations.
//
// Depends on: config.h, precision.h, types.h
// External: FFTW3
// ====================================================================================

#include "../config.h"
#include "../precision.h"
#include "../types.h"

#ifdef __cplusplus
extern "C" {
#endif

// ====================================================================================
// FUNCTION DECLARATIONS
// ====================================================================================

// Create FFTW plans for 2D and 1D transforms
// - plan_2d_out: For Y-slice processing (X-Z plane, size N×N)
// - plan_1d_out: For Y-direction FFT on pencils (size N)
// Both plans use FFTW_ESTIMATE for fast planning
void setup_fftw_plans_full(int N, fftw_plan_t *plan_2d_out, fftw_plan_t *plan_1d_out);

#ifdef __cplusplus
}
#endif

#endif // FFT_SETUP_H

