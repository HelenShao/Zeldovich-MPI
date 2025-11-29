#ifndef HERMITIAN_GENERATION_H
#define HERMITIAN_GENERATION_H

// ====================================================================================
// HERMITIAN 3D MATRIX MPI - HERMITIAN GENERATION MODULE
// ====================================================================================
// This module handles the generation of Hermitian Y-slice pairs with proper symmetry
// constraints for real-valued inverse FFT.
//
// Depends on: config.h, precision.h, types.h, fft/fft_setup.h, utils/verification.h
// External: PCG RNG, OpenMP, FFTW3
// ====================================================================================

#include "../config.h"
#include "../precision.h"
#include "../types.h"

// External RNG functions (from utils/rng.h)
#include "../utils/rng.h"
// Power spectrum functions (from utils/power_spectrum.h)
#include "../utils/power_spectrum.h"
// zeldovich-PLT wrapper (for v15.2)
#include "../utils/zeldovich_wrapper.h"

#ifdef __cplusplus
extern "C" {
#endif

// ====================================================================================
// FUNCTION DECLARATIONS
// ====================================================================================

// Generate one pair of Hermitian Y-slices (primary + conjugate) with 2D FFT
// - N: Grid size
// - global_y: Primary Y-index
// - y_mirror: Mirror Y-index (N-y for conjugate pairs, y for self-conjugate)
// - primary_slices: Output buffer for primary slice [narray][N][N]
// - conjugate_slices: Output buffer for conjugate slice [narray][N][N] (same as primary if self-conjugate)
// - narray: Number of arrays (1, 2, or 4)
// - plan_2d: Precomputed 2D FFT plan
// - rank: MPI rank (for debug output only)
//
// Uses OpenMP to parallelize X-Z loops within the rank
// Uses global_y for deterministic RNG seeding (thread-safe)
// Applies 2D FFT to transform from Fourier space to real space (X,Z)
void generate_hermitian_slice_pair_local(
    int N,
    int global_y,
    int y_mirror,
    fftw_complex_t *primary_slices,   // Flat array (all narray arrays for primary slice)
    fftw_complex_t *conjugate_slices, // Flat array (all narray arrays for conjugate slice)
    int narray,                       // Number of arrays per slice
    fftw_plan_t plan_2d,              // 2D FFT plan
    int rank,                         // MPI rank (for debug output)
    const power_spectrum_params_t *ps_params,  // Legacy power spectrum parameters (NULL = use uniform RNG)
    PowerSpectrumHandle ps_handle,   // v15.2: zeldovich-PLT PowerSpectrum handle (NULL = use legacy or uniform RNG)
    ParametersHandle params_handle);  // v15.2: zeldovich-PLT Parameters handle (needed for fundamental wavenumber)

#ifdef __cplusplus
}
#endif

#endif // HERMITIAN_GENERATION_H

