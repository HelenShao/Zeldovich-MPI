#ifndef ZELDOVICH_WRAPPER_H
#define ZELDOVICH_WRAPPER_H

// ====================================================================================
// ZELDOVICH-PLT C WRAPPER - HEADER
// ====================================================================================
// C interface for zeldovich-PLT PowerSpectrum class
// ====================================================================================

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// Opaque pointer type for PowerSpectrum object
typedef void* PowerSpectrumHandle;

// Opaque pointer type for Parameters object
typedef void* ParametersHandle;

// ====================================================================================
// PARAMETERS INTERFACE
// ====================================================================================

// Create Parameters object from parameter file
// Returns NULL on error
ParametersHandle zeldovich_params_create(const char* param_file);

// Destroy Parameters object
void zeldovich_params_destroy(ParametersHandle params);

// Get fundamental wavenumber (2pi/boxsize)
double zeldovich_params_get_fundamental(ParametersHandle params);

// Get boxsize
double zeldovich_params_get_boxsize(ParametersHandle params);

// Get Pk_scale
double zeldovich_params_get_Pk_scale(ParametersHandle params);

// Get ppd (grid size)
int64_t zeldovich_params_get_ppd(ParametersHandle params);

// Get seed
int zeldovich_params_get_seed(ParametersHandle params);

// Get Pk_powerlaw_index
double zeldovich_params_get_Pk_powerlaw_index(ParametersHandle params);

// Get fundamental wavenumber (2π/BoxSize)
double zeldovich_params_get_fundamental(ParametersHandle params);

// ====================================================================================
// POWER SPECTRUM INTERFACE
// ====================================================================================

// Create PowerSpectrum object
// n: Spline resolution (typically 128)
// params: Parameters object (must remain valid for lifetime of PowerSpectrum)
// Returns NULL on error
PowerSpectrumHandle zeldovich_ps_create(int n, ParametersHandle params);

// Destroy PowerSpectrum object
void zeldovich_ps_destroy(PowerSpectrumHandle ps);

// Initialize from power law
// Returns 0 on success, non-zero on error
int zeldovich_ps_init_powerlaw(PowerSpectrumHandle ps, double powerlaw_index, ParametersHandle params);

// Initialize from file
// Returns 0 on success, non-zero on error
int zeldovich_ps_init_file(PowerSpectrumHandle ps, const char* filename, ParametersHandle params);

// Evaluate power spectrum at given wavenumber
double zeldovich_ps_power(PowerSpectrumHandle ps, double wavenumber);

// Generate random number using zeldovich-PLT's v2rng
// rng_index: Y-slice index for RNG (0 to ppd/2)
// Returns: Random double in (0, 1] (matches one_rand<2>)
// Note: Uses zeldovich-PLT's PCG64 generator array (v2rng)
double zeldovich_ps_one_rand(PowerSpectrumHandle ps, int64_t rng_index);

// Generate complex Gaussian random number with power spectrum variance
// wavenumber: Physical wavenumber (k)
// rng_index: Y-slice index for RNG (0 to ppd/2)
// real: Output real part (double precision from zeldovich-PLT)
// imag: Output imaginary part (double precision from zeldovich-PLT)
// Note: zeldovich-PLT uses double precision. Convert to float if needed.
void zeldovich_ps_cgauss(PowerSpectrumHandle ps, double wavenumber, int64_t rng_index, double* real, double* imag);

// Advance RNG for given Y-slice index
// ps: PowerSpectrum handle
// params: Parameters handle (needed to get ppd for bounds checking)
// rng_index: Y-slice index for RNG (0 to ppd/2)
// nskip: Number of random numbers to skip (will be multiplied by 2 internally, matching zeldovich.cpp)
// Note: This is needed to maintain RNG consistency when N < MAX_PPD
void zeldovich_ps_advance_rng(PowerSpectrumHandle ps, ParametersHandle params, int64_t rng_index, int64_t nskip);

// Get normalization
double zeldovich_ps_get_normalization(PowerSpectrumHandle ps);

// Get powerlaw_index
double zeldovich_ps_get_powerlaw_index(PowerSpectrumHandle ps);

// Get Pk_smooth2
double zeldovich_ps_get_Pk_smooth2(PowerSpectrumHandle ps);

// Get fixed_power flag
int zeldovich_ps_get_fixed_power(PowerSpectrumHandle ps);

#ifdef __cplusplus
}
#endif

#endif

