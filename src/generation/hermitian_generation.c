// ====================================================================================
// HERMITIAN GENERATION MODULE
// ====================================================================================

#include "hermitian_generation.h"
#include "../utils/verification.h"
#include "../utils/zeldovich_wrapper.h" 
#include "../utils/plt_eigenmodes.h"  // For PLT eigenmode support
#include "../config.h"  // For DEBUG_PRINTS, SKIP_VERIFICATION, MAX_PPD
#include "../precision.h"  // For real_t, fabs_t, fmax_t
#include <stdio.h>
#include <stdlib.h>  // For malloc/free
#include <stdint.h> 
#include <math.h>    
#include <omp.h>

// ====================================================================================
// Generates one pair of Hermitian Y-slices (primary + conjugate) in Fourier
// Gaussian with power spectrum weighting (via zeldovich-PLT or legacy code)
// 2D FFT: Fourier --> real space (X,Z)
// Handles RNG nskip tracking
// ====================================================================================

void generate_hermitian_slice_pair_local(
    int N,
    int global_y,
    int y_mirror,
    fftw_complex_t *primary_slices,   
    fftw_complex_t *conjugate_slices, 
    int narray,                       // 4 arrays per slice, 7 C numbers
    fftw_plan_t plan_2d,
    int rank,
    const power_spectrum_params_t *ps_params,  // Legacy power spectrum parameters (NULL = use uniform RNG)
    PowerSpectrumHandle ps_handle,   // v15.2: zeldovich-PLT PowerSpectrum handle (NULL = use legacy or uniform RNG)
    ParametersHandle params_handle)  // v15.2: zeldovich-PLT Parameters handle (needed for fundamental wavenumber)
{
    // Debug: Log entry for seg fault error
    // #if DEBUG_PRINTS
    // fprintf(stderr,
    //         "[Rank %d] ENTER generate_hermitian_slice_pair_local: "
    //         "Y_primary=%d, Y_mirror=%d, ps_handle=%p, params_handle=%p\n",
    //         rank, global_y, y_mirror, (void*)ps_handle, (void*)params_handle);
    // fflush(stderr);
    // #endif

    // Precompute N/2 to avoid repeated division in loops
    int Nhalf = N / 2;
    
    // ========== Check for density-only mode (qdensity == 2) ==========
    // In density-only mode, we skip computing F, G, H (displacements) and only store D (density)
    int just_density = 0;
    if (params_handle != NULL) {
        int qdensity = zeldovich_params_get_qdensity(params_handle);
        just_density = (qdensity == 2);
    }
    
    // ========== PLT Rescaling factors (computed once per function call) ==========
    // These are used to rescale F, G, H when qPLTrescale is enabled
    // target_f: continuum linear theory growth rate
    // a_NL: scale factor at target redshift
    // a0: scale factor at initial redshift
    double target_f = 1.0;
    double a_NL = 1.0;
    double a0 = 1.0;
    int qPLTrescale = 0;
    
    if (params_handle != NULL && !just_density) {
        // Only compute rescaling factors if we're computing displacements
        double f_cluster = zeldovich_params_get_f_cluster(params_handle);
        // target_f is the continuum linear theory growth rate (fluid limit)
        target_f = (sqrt(1. + 24. * f_cluster) - 1.) * 0.25;
        
        qPLTrescale = zeldovich_params_get_qPLTrescale(params_handle);
        if (qPLTrescale) {
            double z_initial = zeldovich_params_get_z_initial(params_handle);
            double PLT_target_z = zeldovich_params_get_PLT_target_z(params_handle);
            a_NL = 1. / (1. + PLT_target_z);  // Scale factor at target redshift
            a0 = 1. / (1. + z_initial);       // Scale factor at initial redshift
        }
    }

    // ========== k_cutoff filtering parameters (computed once per function call) ==========
    // Calculate k2_cutoff for filtering high-wavenumber modes
    // This matches zeldovich.cpp line 321-322: k2_cutoff = nyquist^2 / (k_cutoff^2)
    double k_cutoff = 1.0;  // Default value
    int CornerModes = 0;    // Default value
    double k2_cutoff = 0.0; // Will be calculated
    
    if (params_handle != NULL) {
        k_cutoff = zeldovich_params_get_k_cutoff(params_handle);
        CornerModes = zeldovich_params_get_CornerModes(params_handle);
    }
    
    // Calculate k2_cutoff: k2_cutoff = (N/2)^2 / (k_cutoff^2)
    // For N=16, k_cutoff=1.0: k2_cutoff = 8^2 / 1.0^2 = 64
    double Nhalf_dbl = (double)Nhalf;
    k2_cutoff = (Nhalf_dbl * Nhalf_dbl) / (k_cutoff * k_cutoff);

    // Create local aliases for macro compatibility
    // primary_slices points to slice 0, conjugate_slices points to slice 1
    // Need to access them as if they're part of a larger buffer

    // Helper function to access primary slice (slice_idx=0)
    #define PRIM_SLICE(array_idx, x, z) \
        primary_slices[(int64_t)(x) + (N) * ((z) + (N) * (array_idx))]
    
    // Helper function to access conjugate slice (slice_idx=1)  
    #define CONJ_SLICE(array_idx, x, z) \
        conjugate_slices[(int64_t)(x) + (N) * ((z) + (N) * (array_idx))]
    
    // RNG verification: Track total RNG calls and skips
    // Similar to zeldovich.cpp assertion: assert(Pk.v2rng[y] - checkpoint == 2 * MAX_PPD * MAX_PPD)
    #if VERIFY_RNG_CALLS
    int64_t total_rng_calls = 0;  // Count of cgauss() calls (each uses 2 random numbers)
    int64_t total_rng_skips = 0;  // Count of skipped random numbers (from nskip advances, including D=0 skips)
    #endif
    
    // RNG consistency: Skip random numbers for missing grid points when N < MAX_PPD
    // so we can reproduce sim boxes w/ up to MAX_PPD x MAX_PPD grid generators
    int64_t nskip = 0;
    if (N < MAX_PPD) {
        // Calculate skip for missing grid points
        // Missing z-rows: z = N, N+1, ..., MAX_PPD-1 (skip when z == N-1, after processing all z-rows we iterate)
        // Missing x-values: x = N, N+1, ..., MAX_PPD-1 in each z-row (skip when x == N/2 + 1, crossing Nyquist boundary)
        // For parallel execution, calculate upfront; for sequential, track incrementally
        #if PARALLELIZE_XZ_WITHIN_SLICE
        // Parallel case: calculate total skip upfront (deterministic)
        // Missing z-rows: (MAX_PPD - N) * MAX_PPD (all missing z-rows, each with MAX_PPD x-values)
        // Missing x-values: N * (MAX_PPD - N) (for each z from 0 to N-1, skip missing x-values when crossing Nyquist)
        // Total: (MAX_PPD - N) * MAX_PPD + N * (MAX_PPD - N) = (MAX_PPD - N) * (MAX_PPD + N)
        // Use explicit int64_t casts to avoid integer overflow in intermediate calculations
        nskip = (int64_t)(MAX_PPD - N) * (int64_t)MAX_PPD + (int64_t)N * (int64_t)(MAX_PPD - N);
        #else
        // Sequential case: track incrementally during iteration (like zeldovich.cpp)
        nskip = 0;  // Will be accumulated during loops
        #endif
    }
    
    // In parallel mode, advance upfront skip before starting parallel loops
    #if PARALLELIZE_XZ_WITHIN_SLICE
    if (N < MAX_PPD && nskip > 0) {
        // Advance upfront skip before parallel loops to maintain RNG consistency
        if (ps_handle != NULL && params_handle != NULL) {
            int64_t rng_index = global_y;
            zeldovich_ps_advance_rng(ps_handle, params_handle, rng_index, nskip);
            #if VERIFY_RNG_CALLS
            total_rng_skips += nskip;
            #endif
        } else if (params_handle == NULL) {
            // Cast to uint64_t explicitly to avoid overflow in multiplication
            uint64_t advance_amount = (uint64_t)2 * (uint64_t)nskip;
            advance_pcg_global(global_y, advance_amount);
            #if VERIFY_RNG_CALLS
            total_rng_skips += nskip;
            #endif
        }
        nskip = 0;  // Reset after advancing
    }
    #endif
    
    if (y_mirror != global_y) {
        // ========== CONJUGATE PAIR: Y=i and Y=N-i ==========
        // PARALLELIZATION STRATEGY:
        // - PARALLELIZE_XZ_WITHIN_SLICE=1: Parallel (x,z) loops with locks (faster, needs locks)
        // - PARALLELIZE_XZ_WITHIN_SLICE=0: Sequential (x,z) loops (no locks, like zeldovich.cpp)
        // CRITICAL: Only generate D randomly, compute F, G, H deterministically from D
        #if PARALLELIZE_XZ_WITHIN_SLICE
        // Parallel (x,z) loops: requires locks for thread-safe RNG access
        #pragma omp parallel for collapse(2)
        #else
        // Sequential (x,z) loops: no locks needed, each thread processes different Y-slice
        #endif
        for (int z = 0; z < N; z++) {
            // RNG consistency: When crossing Nyquist boundary (z == Nhalf + 1),
            // skip ALL missing z-rows (z = N to MAX_PPD-1, each containing MAX_PPD x-values)
            // This matches zeldovich.cpp: skip at Nyquist boundary before processing negative kz region
            // The missing frequencies are in the MIDDLE of the MAX_PPD array (high positive and negative k),
            // due to FFT ordering: [0, 1, ..., N/2, -N/2+1, ..., -1]
            // High frequencies (both positive and negative) are physically located in the middle,
            // so we skip them all at once when we first enter the mirrored region
            #if !PARALLELIZE_XZ_WITHIN_SLICE
            if (z == Nhalf + 1 && N < MAX_PPD) {
                // Skip ALL missing z-rows: from z=N to z=MAX_PPD-1
                // This accounts for the "gap" in frequency space (high frequencies in the middle of array)
                int64_t skip_amount = (int64_t)(MAX_PPD - N) * (int64_t)MAX_PPD;
                nskip += skip_amount; // skip missing z-rows
                
                #if DEBUG_RNG_SKIP
                // Debug: Log skip accumulation for test coordinates
                int log_skip = (global_y <= MAX_DEBUG_COORD) || (global_y == Nhalf - 1);
                if (log_skip) {
                    fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d z=%d: ACCUMULATE skip at Nyquist boundary: +%lld (missing z-rows), total nskip=%lld\n",
                            N, global_y, z, (long long)skip_amount, (long long)nskip);
                    fflush(stderr);
                }
                #endif
            }
            #endif
            
            for (int x = 0; x < N; x++) {
                int x_mirror = (x == 0) ? 0 : N - x;
                int z_mirror = (z == 0) ? 0 : N - z;
                
                // ========== STEP 1: Calculate k-vector components ==========
                // RNG consistency: When crossing Nyquist boundary (x == Nhalf + 1),
                // skip ALL missing x-values (x = N to MAX_PPD-1) in current z-row
                // This matches zeldovich.cpp: skip at Nyquist boundary before processing negative kx region
                // The missing frequencies are in the MIDDLE of the MAX_PPD array (high positive and negative k),
                // due to FFT ordering: [0, 1, ..., N/2, -N/2+1, ..., -1]
                // High frequencies (both positive and negative) are physically located in the middle,
                // so we skip them all at once when we first enter the mirrored region
                #if !PARALLELIZE_XZ_WITHIN_SLICE
                if (x == Nhalf + 1 && N < MAX_PPD) {
                    // Skip ALL missing x-values: from x=N to x=MAX_PPD-1
                    // This accounts for the "gap" in frequency space (high frequencies in the middle of array)
                    int64_t skip_amount = MAX_PPD - N;
                    nskip += skip_amount; // skip missing x-values in this z-row
                    
                    #if DEBUG_RNG_SKIP
                    // Debug: Log skip accumulation for test coordinates
                    int log_skip = (z <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD) ||
                                   (z == Nhalf - 1 && global_y == Nhalf - 1);
                    if (log_skip) {
                        fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d z=%d x=%d: ACCUMULATE skip at Nyquist boundary: +%lld (missing x-values), total nskip=%lld\n",
                                N, global_y, z, x, (long long)skip_amount, (long long)nskip);
                        fflush(stderr);
                    }
                    #endif
                }
                #endif
                
                int kx = (x > Nhalf) ? x - N : x;
                int ky = (global_y > Nhalf) ? global_y - N : global_y;
                int kz = (z > Nhalf) ? z - N : z;
                int k2_int = kx*kx + ky*ky + kz*kz;  // Integer k² for k_cutoff comparison
                double k2 = (double)k2_int;  // Floating-point k² for power spectrum
                
                // ========== STEP 2: Generate D using RNG or cgauss() ==========
                // Nyquist frequency zeroing: Force Nyquist elements to zero for all three axes
                // This matches zeldovich.cpp line 354: abs(kx)==kmax || abs(kz)==kmax || abs(ky)==kmax
                // The Nyquist frequency (k = N/2) is self-conjugate and doesn't have a separate partner.
                // Due to the Y-shift in the reflected shell, we need to zero these to align mirroring expectations.
                int abs_kx = (kx < 0) ? -kx : kx;
                int abs_ky = (ky < 0) ? -ky : ky;
                int abs_kz = (kz < 0) ? -kz : kz;
                int is_nyquist = (abs_kx == Nhalf || abs_ky == Nhalf || abs_kz == Nhalf);
                
                fftw_complex D;
                // Zero D for: DC mode, Nyquist frequency, or k_cutoff filtering
                if ((k2 == 0.0)
                     || (is_nyquist)
                    // Force all elements with wavenumber above k_cutoff (nominally k_Nyquist) to zero
                     || (!CornerModes && (double)k2_int >= k2_cutoff)) {
                    D[0] = D[1] = 0.0;
                    // RNG consistency: When D=0, we skip the RNG call, so accumulate skip
                    // This matches zeldovich.cpp line 361: nskip++ (accumulate, don't advance immediately)
                    #if !PARALLELIZE_XZ_WITHIN_SLICE
                    nskip++;  // Accumulate skip, will be applied before next cgauss() call
                    #if DEBUG_RNG_SKIP
                    int log_skip = (x <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD && z <= MAX_DEBUG_COORD) ||
                                   (x == Nhalf - 1 && global_y == Nhalf - 1 && z == Nhalf - 1);
                    if (log_skip) {
                        if (k2 == 0.0) {
                            fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d (x,z)=(%d,%d): D=0 (DC mode), ACCUMULATE nskip++ (total=%lld)\n",
                                    N, global_y, x, z, (long long)nskip);
                        } else if (is_nyquist) {
                            fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d (x,z)=(%d,%d): D=0 (Nyquist: kx=%d ky=%d kz=%d), ACCUMULATE nskip++ (total=%lld)\n",
                                    N, global_y, x, z, kx, ky, kz, (long long)nskip);
                        } else {
                            fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d (x,z)=(%d,%d): D=0 (k_cutoff: k2=%d >= %.1f), ACCUMULATE nskip++ (total=%lld)\n",
                                    N, global_y, x, z, k2_int, k2_cutoff, (long long)nskip);
                        }
                        fflush(stderr);
                    }
                    #endif
                    #endif
                } else if (ps_handle != NULL && params_handle != NULL) {
                    // v15.2: Use zeldovich-PLT power spectrum-weighted Gaussian
                    // If ps_handle is available, we ALWAYS use power spectrum mode (cgauss)
                    
                    // Convert k indices to physical wavenumber: k_phys = k_index * fundamental
                    double fundamental = zeldovich_params_get_fundamental(params_handle);
                    double k2_phys = k2 * fundamental * fundamental;
                    double kmag = sqrt(k2_phys);
                    
                    // zeldovich_ps_cgauss returns double precision, convert to real_t
                    // zeldovich-PLT's v2rng array is sized to ppd/2, so valid indices are 0 to (N/2 - 1)
                    // In conjugate pair branch, global_y is in range [1, N/2-1] or [N/2+1, N-1]
                    int64_t rng_index = global_y;
                    
                    // RNG consistency: Advance zeldovich-PLT's RNG when crossing Nyquist boundaries
                    // This matches zeldovich.cpp behavior: advance before calling cgauss
                    #if !PARALLELIZE_XZ_WITHIN_SLICE
                    #if DEBUG_RNG_SKIP
                    // Debug: Log skip application for test coordinates
                    int log_skip = (x <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD && z <= MAX_DEBUG_COORD) ||
                                   (x == Nhalf - 1 && global_y == Nhalf - 1 && z == Nhalf - 1);
                    if (log_skip && nskip > 0) {
                        fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d (x,z)=(%d,%d): APPLYING skip=%lld BEFORE cgauss\n",
                                N, global_y, x, z, (long long)nskip);
                        fflush(stderr);
                    }
                    #endif
                    if (nskip > 0) {
                        zeldovich_ps_advance_rng(ps_handle, params_handle, rng_index, nskip);
                        #if VERIFY_RNG_CALLS
                        total_rng_skips += nskip;  // Track accumulated skip value
                        #endif
                        #if DEBUG_RNG_SKIP
                        if (log_skip) {
                            fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d (x,z)=(%d,%d): RNG advanced, nskip reset to 0\n",
                                    N, global_y, x, z);
                            fflush(stderr);
                        }
                        #endif
                        nskip = 0;  // Reset after advancing
                    }
                    #endif
                    
                    double D_real, D_imag;
                    #if DEBUG_RNG_SKIP
                    int log_cgauss = (x <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD && z <= MAX_DEBUG_COORD) ||
                                     (x == Nhalf - 1 && global_y == Nhalf - 1 && z == Nhalf - 1);
                    if (log_cgauss) {
                        fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d (x,z)=(%d,%d): CALLING cgauss (nskip=%lld, kmag=%.6e)\n",
                                N, global_y, x, z, (long long)nskip, kmag);
                        fflush(stderr);
                    }
                    #endif
                    zeldovich_ps_cgauss(ps_handle, kmag, rng_index, &D_real, &D_imag);
                    #if VERIFY_RNG_CALLS
                    total_rng_calls++;  // Each cgauss() call uses 2 random numbers
                    #endif
                    D[0] = (real_t)D_real;
                    D[1] = (real_t)D_imag;
                    // Note: For VERIFY_HERMITIAN_SYMMETRY==2, we set D=0 AFTER computing F and H
                } else if (ps_params != NULL) {
                    // Legacy: Use my own power spectrum-weighted Gaussian (cgauss)
                    // Convert k indices to physical wavenumber: k_phys = k_index * fundamental
                    // For now, fundamental = 1.0 (can be added to ps_params later)
                    double fundamental = 1.0;  // TODO: Add to power_spectrum_params_t
                    double k2_phys = k2 * fundamental * fundamental;
                    double kmag = sqrt(k2_phys);
                    
                    // RNG consistency: Advance our own RNG when crossing Nyquist boundaries
                    // This matches the zeldovich-PLT branch behavior
                    #if !PARALLELIZE_XZ_WITHIN_SLICE
                    if (nskip > 0 && N < MAX_PPD) {
                        advance_pcg_global(global_y, 2 * nskip);
                        #if VERIFY_RNG_CALLS
                        total_rng_skips += nskip;  // Track accumulated skip value
                        #endif
                        nskip = 0;  // Reset after advancing
                    }
                    #endif
                    
                    cgauss(ps_params, kmag, global_y, &D);
                    #if VERIFY_RNG_CALLS
                    total_rng_calls++;  // Each cgauss() call uses 2 random numbers
                    #endif
                    // Note: For VERIFY_HERMITIAN_SYMMETRY==2, we set D=0 AFTER computing F and H
                } else {
                    // Fallback: use uniform random numbers (white noise mode, no power spectrum)
                    // This should only happen when ps_handle is NULL (no parameter file provided)
                    // Note: If ps_handle is available, we should have used cgauss() above
                    // This branch is only for backward compatibility (no power spectrum mode)
                    
                    // RNG consistency: Advance our own RNG when crossing Nyquist boundaries
                    // This matches the other branches' behavior
                    #if !PARALLELIZE_XZ_WITHIN_SLICE
                    if (nskip > 0 && N < MAX_PPD) {
                        advance_pcg_global(global_y, 2 * nskip);
                        #if VERIFY_RNG_CALLS
                        total_rng_skips += nskip;  // Track accumulated skip value
                        #endif
                        nskip = 0;  // Reset after advancing
                    }
                    #endif
                    
                    double D_re = random_real_pcg_global(global_y);
                    double D_im = random_real_pcg_global(global_y);
                    #if VERIFY_RNG_CALLS
                    total_rng_calls++;  // Two random_real calls = 1 complex number = 2 random numbers
                    #endif
                    D[0] = D_re;
                    D[1] = D_im;
                    // Note: For VERIFY_HERMITIAN_SYMMETRY==2, we set D=0 AFTER computing F and H
                }
                
                // ========== STEP 3: Compute F, G, H from D ==========
                // Skip F, G, H computation in density-only mode (qdensity == 2)
                fftw_complex F, G, H;
                int use_plt = 0;  // Declare outside to use in f computation
                eigenmode e;      // Declare outside to use in f computation
                
                // Declare f before if/else blocks (needed for velocity arrays in all modes)
                double f = 1.0;
                
                // Match zeldovich.cpp: set f=0 when D==0 (line 443)
                if (D[0] == 0.0 && D[1] == 0.0) {
                    f = 0.0;
                }
                
                if (just_density) {
                    // Density-only mode: Set F, G, H to zero (not used)
                    F[0] = F[1] = G[0] = G[1] = H[0] = H[1] = 0.0;

                    // Density-only mode: Compute F from D (for D+i*F test)
                    // Match zeldovich behavior: only compute F when D != 0
                    // if (D[0] != 0.0 || D[1] != 0.0) {
                    //     // Avoid division by zero: if k2==0, set F=0 (D is already 0 for k2==0)
                    //     if (k2 == 0.0) {
                    //         F[0] = F[1] = 0.0;
                    //     } else {
                    //         double fundamental = 1.0;
                    //         if (params_handle != NULL) {
                    //             fundamental = zeldovich_params_get_fundamental(params_handle);
                    //         }
                    //         // For VERIFY_HERMITIAN_SYMMETRY=1, we don't use PLT rescaling, so rescale=1.0
                    //         double rescale = 1.0;
                    //         // Compute factor the same way as normal mode
                    //         double factor = rescale / (k2 * fundamental);

                    //         F[0] = -kx * factor * D[1];
                    //         F[1] =  kx * factor * D[0];
                    //     }
                    // } else {
                    //     // D == 0, so F = 0 (match zeldovich line 442)
                    //     F[0] = F[1] = 0.0;
                    // }

                    // // D[0] = D[0] - F[1]; // Test D+iF calculation
                    // // D[1] = D[1] + F[0]; // Test D+iF calculation
                    // // F[0] = F[1] // Test D+iF calculation
                    // G[0] = G[1] = H[0] = H[1] = 0.0;
                } else {
                double ik2 = 1.0 / k2; // later: define this at top of code and use it instead of dividing
                
                #if defined(VERIFY_HERMITIAN_SYMMETRY) && VERIFY_HERMITIAN_SYMMETRY == 1
                // printf("Hermitian mode 1: F=0 and H=0 to test Hermitian symmetry\n");
                // Verification mode 1: Set F=0 and H=0 to test Hermitian symmetry
                // With F=0 and H=0, conjugate slices should be true conjugates of primary slices
                // After 3D FFT, the result should be purely real
                // IMPORTANT: Must use the same factor computation as normal mode to ensure
                // G is computed identically. Otherwise the test is invalid.
                // Get fundamental wavenumber (needed for factor computation)
                double fundamental = 1.0;
                if (params_handle != NULL) {
                    fundamental = zeldovich_params_get_fundamental(params_handle);
                }
                // For VERIFY_HERMITIAN_SYMMETRY=1, we don't use PLT rescaling, so rescale=1.0
                double rescale = 1.0;
                // Compute factor the same way as normal mode
                double factor = rescale / (k2 * fundamental);
                
                F[0] = 0.0;
                F[1] = 0.0;
                
                // Compute G using the same formula as normal mode
                G[0] = -ky * factor * D[1];
                G[1] =  ky * factor * D[0];
                
                H[0] = 0.0;
                H[1] = 0.0;
                #elif defined(VERIFY_HERMITIAN_SYMMETRY) && VERIFY_HERMITIAN_SYMMETRY == 2
                // Verification mode 2: Compute F and H from D, then set D=0 and G=0
                // This makes Array 0 and Array 1 purely imaginary (real parts = 0)
                // After 3D FFT, the result should be purely imaginary (real parts = 0)
                // IMPORTANT: Must use the same factor computation as normal mode to ensure
                // F is computed identically. Otherwise the test is invalid.
                // Get fundamental wavenumber (needed for factor computation)
                double fundamental = 1.0;
                if (params_handle != NULL) {
                    fundamental = zeldovich_params_get_fundamental(params_handle);
                }
                // For VERIFY_HERMITIAN_SYMMETRY=2, we don't use PLT rescaling, so rescale=1.0
                double rescale = 1.0;
                // Compute factor the same way as normal mode
                double factor = rescale / (k2 * fundamental);
                
                // Now compute F and H using the same formula as normal mode
                F[0] = -kx * factor * D[1];
                F[1] =  kx * factor * D[0];
                
                G[0] = 0.0;  // Set G=0
                G[1] = 0.0;
                
                H[0] = -kz * factor * D[1];
                H[1] =  kz * factor * D[0];
                
                // Now set D=0 (after F and H are computed)
                D[0] = 0.0;
                D[1] = 0.0;
                #else
                // Normal operation: Compute F, G, H from D
                // Check if PLT is enabled
                int qPLT = 0;
                // Default fundamental = 1.0 
                double fundamental = 1.0;
                
                // Get fundamental wavenumber (needed for both PLT and non-PLT cases)
                if (params_handle != NULL) {
                    fundamental = zeldovich_params_get_fundamental(params_handle);
                    qPLT = zeldovich_params_get_qPLT(params_handle);
                    if (qPLT) {
                        // Get PLT eigenmode for this k-vector
                        // Convert kx, ky, kz to array indices for plt_get_eigenmode
                        int ikx = (kx < 0) ? N + kx : kx;
                        int iky = (ky < 0) ? N + ky : ky;
                        int ikz = (kz < 0) ? N + kz : kz;
                        // Handle z index (only positive half-space stored)
                        if (ikz > N / 2) ikz = N - ikz;
                        
                        if (plt_get_eigenmode(ikx, iky, ikz, (int64_t)N, &e) == 0) {
                            use_plt = 1;
                            
                            // Apply the same normalization as zeldovich.cpp
                            // 1. Set sign of z component (real FFT only gives +kz half-space)
                            // 2. Ensure |e| = 1 (normalize after interpolation, like zeldovich)
                            // 3. Apply norm = k2 / (k * e) to get the final eigenvector
                            
                            // Set the sign of the z component (because the real FFT only gives the +kz half-space)
                            // see zeldovich.cpp line 255: ehat.vec[2] *= copysign(1, kz);
                            if (kz < 0) {
                                e.vec[2] = -e.vec[2];
                            }
                            
                            // Normalize eigenvector (interpolation might not preserve |e| = 1)
                            // see zeldovich.cpp lines 257-263
                            double e_mag = sqrt(e.vec[0] * e.vec[0] + e.vec[1] * e.vec[1] + e.vec[2] * e.vec[2]);
                            if (e_mag > 0.0) {
                                e.vec[0] /= e_mag;
                                e.vec[1] /= e_mag;
                                e.vec[2] /= e_mag;
                            }
                            
                            // Apply normalization: norm = k2 / (k * e)
                            // This upweights each mode by 1/(khat*ehat),see zeldovich.cpp line 266
                            double k_dot_e = kx * e.vec[0] + ky * e.vec[1] + kz * e.vec[2];
                            double norm = (k2 > 0.0 && k_dot_e != 0.0) ? k2 / k_dot_e : 0.0;
                            if (!isfinite(norm)) norm = 0.0;
                            
                            // Scale eigenvector by norm (see zeldovich.cpp lines 268-270)
                            e.vec[0] *= norm;
                            e.vec[1] *= norm;
                            e.vec[2] *= norm;
                        } else {
                            // If eigenmode lookup fails, fall back to normal computation
                            if (rank == 0 && x == 0 && z == 0) {
                                fprintf(stderr, "[WARNING] Failed to get PLT eigenmode for (kx=%d, ky=%d, kz=%d), using normal computation\n",
                                        kx, ky, kz);
                            }
                        }
                    }
                }
                
                // ========== STEP 3.5: Compute PLT growth rate f and rescale (before computing F, G, H) ==========
                // f is the logarithmic derivative of the growth factor that scales velocities
                // When PLT is enabled: f = (sqrt(1. + 24 * e.val * f_cluster) - 1) / 4.
                // When PLT is not enabled: f = 1.0 (default)
                // Skip in density-only mode (qdensity == 2)
                // NOTE: Use outer f variable (declared on 397), not a new inner f!
                // Otherwise the PLT f value goes out of scope before velocity storage.
                double rescale = 1.0;
                
                if (!just_density && use_plt && params_handle != NULL) {
                    double f_cluster = zeldovich_params_get_f_cluster(params_handle);
                    // PLT growth rate: f = (sqrt(1. + 24 * e.val * f_cluster) - 1) / 4.
                    // This is f_growth, the logarithmic derivative of the growth factor
                    // that scales the velocities. The corrections are sourced from:
                    // 1) PLT growth rate 2) Addition of a smooth, non-clustering component
                    // to the background (<= NOT A PLT EFFECT)
                    // If PLT is turned on, we have to combine the effects here.
                    // If not, we apply f_cluster during output.
                    f = (sqrt(1. + 24. * e.val * f_cluster) - 1.) * 0.25;
                    
                    // Compute rescaling if qPLTrescale is enabled
                    // rescale = pow(a_NL/a0, target_f - plt_f)
                    // where plt_f = f (PLT growth rate) and target_f is the continuum growth rate
                    if (qPLTrescale) {
                        double plt_f = f;  // PLT growth rate for this mode
                        rescale = pow(a_NL / a0, target_f - plt_f);
                    }
                }
                
                // Compute factor (used identically in both PLT and non-PLT cases)
                // In zeldovich.cpp: k2 includes fundamental^2, so ik2 = 1/(k2_index * fundamental^2)
                // F = rescale * I * vec * fundamental * ik2 * D
                //   = rescale * I * vec * fundamental / (k2_index * fundamental^2) * D
                //   = rescale * I * vec / (k2_index * fundamental) * D
                // where vec is either e.vec[i] (PLT) or k[i] (non-PLT)
                // Note: our ik2 = 1/k2_index (no fundamental), so we use 1/(k2*fundamental)
                double factor = rescale / (k2 * fundamental);
                
                if (use_plt) {
                    // PLT mode: Use eigenvector instead of k-vector
                    #if DEBUG_EIGENVECTOR
                    // Debug: Print eigenvector components for test coordinates
                    int debug_eigen = (rank == 0 && x <= 2 && global_y <= 2 && z <= 2);
                    if (debug_eigen) {
                        fprintf(stderr, "[EIGEN-DEBUG] N=%d Y=%d (x,z)=(%d,%d) k=(%d,%d,%d): e.vec=[%.6f, %.6f, %.6f] e.val=%.6f\n",
                                N, global_y, x, z, kx, ky, kz, e.vec[0], e.vec[1], e.vec[2], e.val);
                        fflush(stderr);
                    }
                    #endif
                    // F = rescale * i * e.vec[0] * fundamental * ik2 * D
                    //   = rescale * i * e.vec[0] * fundamental * ik2 * (D_re + i*D_im)
                    //   = rescale * (-e.vec[0] * fundamental * ik2 * D_im + i * e.vec[0] * fundamental * ik2 * D_re)
                    F[0] = -e.vec[0] * factor * D[1];  // Real part
                    F[1] =  e.vec[0] * factor * D[0];  // Imaginary part
                    
                    G[0] = -e.vec[1] * factor * D[1];
                    G[1] =  e.vec[1] * factor * D[0];
                    
                    H[0] = -e.vec[2] * factor * D[1];
                    H[1] =  e.vec[2] * factor * D[0];
                } else {
                    // Normal operation: Compute F, G, H from D using k-vector
                    // In zeldovich.cpp: k2 includes fundamental^2, so F = I * kx * fundamental * ik2 * D
                    // In our code: k2 doesn't include fundamental^2, so we need to multiply by fundamental
                    // F = rescale * i * kx * fundamental * ik2 * D
                    //   = rescale * i * kx * fundamental * ik2 * (D_re + i*D_im)
                    //   = rescale * (-kx * fundamental * ik2 * D_im + i * kx * fundamental * ik2 * D_re)
                    F[0] = -kx * factor * D[1];  // Real part
                    F[1] =  kx * factor * D[0];  // Imaginary part
                    
                    G[0] = -ky * factor * D[1];
                    G[1] =  ky * factor * D[0];
                    
                    H[0] = -kz * factor * D[1];
                    H[1] =  kz * factor * D[0];
                }
                
                // Print factor debug info while variables are in scope
                #if DEBUG_RNG_CONSISTENCY
                // Check if we should print debug info for this coordinate
                int boundary_coord = Nhalf - 1;
                int test_boundary = (boundary_coord >= 0 && 
                                     x == boundary_coord && global_y == boundary_coord && z == boundary_coord);
                int effective_boundary = (boundary_coord <= MAX_DEBUG_BOUNDARY_COORD) ? boundary_coord : MAX_DEBUG_BOUNDARY_COORD;
                int max_test_coord = (MAX_DEBUG_COORD > effective_boundary) ? MAX_DEBUG_COORD : effective_boundary;
                int in_test_range = (x <= max_test_coord && global_y <= max_test_coord && z <= max_test_coord);
                int should_print_factor = 0;
                if (test_boundary) {
                    should_print_factor = 1;
                } else if (N <= DEBUG_FULL_PRINT_MAX_N) {
                    should_print_factor = in_test_range;
                } else {
                    if (x <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD && z <= MAX_DEBUG_COORD) {
                        should_print_factor = 1;
                    } else if (in_test_range) {
                        should_print_factor = (x % DEBUG_SAMPLE_STRIDE == 0 &&
                                               global_y % DEBUG_SAMPLE_STRIDE == 0 &&
                                               z % DEBUG_SAMPLE_STRIDE == 0);
                    }
                }
                if (should_print_factor) {
                    if (use_plt && params_handle != NULL) {
                        fprintf(stderr, "[FACTOR-DEBUG] N=%d Y=%d (x,z)=(%d,%d): fundamental=%.10e rescale=%.10e "
                                "e.vec=[%.10e,%.10e,%.10e] e.val=%.10e factor=%.10e\n",
                                N, global_y, x, z, fundamental, rescale,
                                e.vec[0], e.vec[1], e.vec[2], e.val, factor);
                    } else {
                        fprintf(stderr, "[FACTOR-DEBUG] N=%d Y=%d (x,z)=(%d,%d): fundamental=%.10e rescale=%.10e "
                                "k=[%d,%d,%d] factor=%.10e\n",
                                N, global_y, x, z, fundamental, rescale, kx, ky, kz, factor);
                    }
                    fflush(stderr);
                }
                #endif
                
                // ========== DEBUG: Verify F Hermitian symmetry ==========
                #if DEBUG_RNG_CONSISTENCY && !SKIP_VERIFICATION
                // For small N, verify that F has proper Hermitian symmetry
                // Check: F(-k) should equal conj(F(k))
                if (N <= DEBUG_FULL_PRINT_MAX_N && y_mirror != global_y) {
                    // Only check for conjugate pairs (not self-conjugate)
                    int x_mirror_check = (x == 0) ? 0 : N - x;
                    int z_mirror_check = (z == 0) ? 0 : N - z;
                    int kx_mirror = (x_mirror_check > Nhalf) ? x_mirror_check - N : x_mirror_check;
                    int kz_mirror = (z_mirror_check > Nhalf) ? z_mirror_check - N : z_mirror_check;
                    
                    // Compute what F(-k) should be: F(-k) = i * (-kx) * factor * D(-k)
                    // But we're storing conj(F) at -k, so we need to verify that matches
                    // For now, just print F values for comparison
                    if (x <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD && z <= MAX_DEBUG_COORD) {
                        fprintf(stderr, "[F-HERMITIAN-DEBUG] N=%d Y=%d (x,z)=(%d,%d) k=(%d,%d,%d): "
                                "F=(%.10e,%.10e) factor=%.10e kx=%d\n",
                                N, global_y, x, z, kx, ky, kz, F[0], F[1], factor, kx);
                        if (x_mirror_check <= MAX_DEBUG_COORD && z_mirror_check <= MAX_DEBUG_COORD) {
                            fprintf(stderr, "[F-HERMITIAN-DEBUG] N=%d Y=%d (x,z)=(%d,%d) k_mirror=(%d,%d,%d): "
                                    "Will store conj(F)=(%.10e,%.10e) at mirror position\n",
                                    N, y_mirror, x_mirror_check, z_mirror_check, kx_mirror, 
                                    (global_y > Nhalf) ? global_y - N : N - global_y, kz_mirror,
                                    F[0], -F[1]);
                        }
                        fflush(stderr);
                    }
                }
                #endif
                #endif
                }  // End of else block for !just_density
                
                // ========== DEBUG: Print RNG values for consistency checking ==========
                #if DEBUG_RNG_CONSISTENCY
                // Print D, F, G, H for test coordinates to verify RNG consistency across N
                // Tests all coordinates (x,y,z) where x,y,z <= MAX_DEBUG_COORD
                // Also tests coordinates at (N/2)-1 to ensure boundary of overlapping region is tested
                // This ensures overlapping regions are tested for any pair of N values
                // where min(N1, N2) >= 2*MAX_DEBUG_COORD
                int boundary_coord = Nhalf - 1;  // (N/2)-1 for current N
                int test_boundary = (boundary_coord >= 0 && 
                                     x == boundary_coord && global_y == boundary_coord && z == boundary_coord);
                int effective_boundary = (boundary_coord <= MAX_DEBUG_BOUNDARY_COORD) ? boundary_coord : MAX_DEBUG_BOUNDARY_COORD;
                int max_test_coord = (MAX_DEBUG_COORD > effective_boundary) ? MAX_DEBUG_COORD : effective_boundary;
                // Check if coordinate is in test range
                int in_test_range = (x <= max_test_coord && global_y <= max_test_coord && z <= max_test_coord);
                // For large N, use sampling to reduce output volume
                // Always print boundary coordinates and coordinates in small N range
                int should_print = 0;
                if (test_boundary) {
                    // Always print boundary coordinate (N/2)-1, (N/2)-1, (N/2)-1
                    should_print = 1;
                } else if (N <= DEBUG_FULL_PRINT_MAX_N) {
                    // For small N, print all coordinates in test range
                    should_print = in_test_range;
                } else {
                    // For large N, print only sampled coordinates 
                    if (x <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD && z <= MAX_DEBUG_COORD) {
                        should_print = 1;  // Print small coordinate cube
                    } else if (in_test_range) {
                        // Sample: only print if x, y, z are multiples of stride
                        should_print = (x % DEBUG_SAMPLE_STRIDE == 0 &&
                                       global_y % DEBUG_SAMPLE_STRIDE == 0 &&
                                       z % DEBUG_SAMPLE_STRIDE == 0);
                    }
                }
                if (should_print) {
                    fprintf(stderr, "[RNG-DEBUG] N=%d Y=%d (x,z)=(%d,%d): k=(%d,%d,%d) k2=%.6f | "
                            "D=(%.10e,%.10e) F=(%.10e,%.10e) G=(%.10e,%.10e) H=(%.10e,%.10e)\n",
                            N, global_y, x, z, kx, ky, kz, k2,
                            D[0], D[1], F[0], F[1], G[0], G[1], H[0], H[1]);
                    fflush(stderr);
                }
                #endif
                
                // ========== STEP 4: Store in arrays ==========
                if (just_density) {
                    // Density-only mode: Only store D (density) in Array 0
                    // Array 0: D (density only, no displacement)
                    PRIM_SLICE(0, x, z)[0] = D[0];  // Real = D_re
                    PRIM_SLICE(0, x, z)[1] = D[1];  // Imag = D_im

                    // Test D+iF calculation
                    // PRIM_SLICE(0, x, z)[0] = D[0] - F[1];  // Real = D_re - F_im
                    // PRIM_SLICE(0, x, z)[1] = D[1] + F[0];  // Imag = D_im + F_re
                } else {
                    // Normal mode: Store D+iF, G+iH, and optionally velocities
                    // Array 0: D + i*F (density + X-displacement)
                    // D + i*F = (D[0] + i*D[1]) + i*(F[0] + i*F[1]) = (D[0] - F[1]) + i*(D[1] + F[0])
                    PRIM_SLICE(0, x, z)[0] = D[0] - F[1];  // Real = D_re - F_im
                    PRIM_SLICE(0, x, z)[1] = D[1] + F[0];  // Imag = D_im + F_re
                    
                    // ========== DEBUG: Verify stored D+i*F values ==========
                    #if DEBUG_RNG_CONSISTENCY && !SKIP_VERIFICATION
                    if (N <= DEBUG_FULL_PRINT_MAX_N && x <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD && z <= MAX_DEBUG_COORD) {
                        double stored_re = PRIM_SLICE(0, x, z)[0];
                        double stored_im = PRIM_SLICE(0, x, z)[1];
                        double expected_re = D[0] - F[1];
                        double expected_im = D[1] + F[0];
                        double diff_re = fabs_t(stored_re - expected_re);
                        double diff_im = fabs_t(stored_im - expected_im);
                        if (diff_re > 1e-10 || diff_im > 1e-10) {
                            fprintf(stderr, "[STORAGE-DEBUG] N=%d Y=%d (x,z)=(%d,%d): Storage mismatch! "
                                    "stored=(%.10e,%.10e) expected=(%.10e,%.10e) diff=(%.10e,%.10e)\n",
                                    N, global_y, x, z, stored_re, stored_im, expected_re, expected_im, diff_re, diff_im);
                            fflush(stderr);
                        }
                        fprintf(stderr, "[STORAGE-DEBUG] N=%d Y=%d (x,z)=(%d,%d): D=(%.10e,%.10e) F=(%.10e,%.10e) "
                                "stored D+i*F=(%.10e,%.10e)\n",
                                N, global_y, x, z, D[0], D[1], F[0], F[1], stored_re, stored_im);
                        fflush(stderr);
                    }
                    #endif
                    
                    // Array 1: G + i*H (Y-displacement + Z-displacement)
                    // G + i*H = (G[0] + i*G[1]) + i*(H[0] + i*H[1]) = (G[0] - H[1]) + i*(G[1] + H[0])
                    PRIM_SLICE(1, x, z)[0] = G[0] - H[1];  // Real = G_re - H_im
                    PRIM_SLICE(1, x, z)[1] = G[1] + H[0];  // Imag = G_im + H_re
                    
                    if (narray >= 4) { // qPLT == true
                        // Array 2: 0 + i*F*f (X-velocity)
                        // 0 + i*(F*f) = 0 + i*((F[0] + i*F[1])*f) = -F[1]*f + i*(F[0]*f)
                        // f is computed above (PLT growth rate if PLT enabled, else 1.0)
                        PRIM_SLICE(2, x, z)[0] = -F[1] * f;  // Real = -F_im * f
                        PRIM_SLICE(2, x, z)[1] = F[0] * f;   // Imag = F_re * f
                        
                        // Array 3: G*f + i*H*f (Y-velocity + Z-velocity)
                        // (G*f) + i*(H*f) = (G[0] - H[1])*f + i*((G[1] + H[0])*f)
                        PRIM_SLICE(3, x, z)[0] = (G[0] - H[1]) * f;  // Real = (G_re - H_im) * f
                        PRIM_SLICE(3, x, z)[1] = (G[1] + H[0]) * f;  // Imag = (G_im + H_re) * f
                    }
                }
                
                // ========== STEP 5: Store conjugates (Zeldovich scheme: conj(D) + i*conj(F)) ==========
                if (just_density) {
                    // Density-only mode: Only store D (density) in Array 0
                    // For conjugate, store conj(D) = (D[0], -D[1])
                    CONJ_SLICE(0, x_mirror, z_mirror)[0] = D[0];   // Real = D_re
                    CONJ_SLICE(0, x_mirror, z_mirror)[1] = -D[1];  // Imag = -D_im (conjugate)
                    
                    // Test D+iF calculation
                    // CONJ_SLICE(0, x_mirror, z_mirror)[0] = D[0] + F[1];  // Real = D_re + F_im
                    // CONJ_SLICE(0, x_mirror, z_mirror)[1] = F[0] - D[1];  // Imag = F_re - D_im
                } else {
                    // Normal mode: Store conj(D)+i*conj(F), conj(G)+i*conj(H), and optionally velocities
                    // For mode -k, store conj(D) + i*conj(F), NOT the conjugate of (D + i*F)!
                    // conj(D) + i*conj(F) = (D[0] - i*D[1]) + i*(F[0] - i*F[1]) = (D[0] + F[1]) + i*(F[0] - D[1])
                    CONJ_SLICE(0, x_mirror, z_mirror)[0] = D[0] + F[1];  // Real = D_re + F_im
                    CONJ_SLICE(0, x_mirror, z_mirror)[1] = F[0] - D[1];  // Imag = F_re - D_im
                    
                    // conj(G) + i*conj(H) = (G[0] - i*G[1]) + i*(H[0] - i*H[1]) = (G[0] + H[1]) + i*(H[0] - G[1])
                    CONJ_SLICE(1, x_mirror, z_mirror)[0] = G[0] + H[1];  // Real = G_re + H_im
                    CONJ_SLICE(1, x_mirror, z_mirror)[1] = H[0] - G[1];  // Imag = H_re - G_im
                    
                    if (narray >= 4) {
                        // Array 2: 0 + i*conj(F*f) = i*conj(F*f)
                        // conj(F*f) = F_re*f - i*F_im*f
                        // i*conj(F*f) = i*(F_re*f - i*F_im*f) = F_im*f + i*F_re*f
                        // f is computed above (PLT growth rate if PLT enabled, else 1.0)
                        CONJ_SLICE(2, x_mirror, z_mirror)[0] = F[1] * f;   // Real = F_im * f
                        CONJ_SLICE(2, x_mirror, z_mirror)[1] = F[0] * f;   // Imag = F_re * f (FIXED: was -F[0])
                        
                        // Array 3: conj(G*f) + i*conj(H*f) = (G[0] + H[1])*f + i*((H[0] - G[1])*f)
                        CONJ_SLICE(3, x_mirror, z_mirror)[0] = (G[0] + H[1]) * f;  // Real = (G_re + H_im) * f
                        CONJ_SLICE(3, x_mirror, z_mirror)[1] = (H[0] - G[1]) * f;  // Imag = (H_re - G_im) * f
                    }
                }
                
            }
            
        }
    } else {
        // ========== SELF-CONJUGATE: Y=0 or Y=N/2 ==========
        // Same approach as conjugate pair: generate D, compute F,G,H
        // Reset nskip for self-conjugate case
        // For self-conjugate slices, we process full N×N plane (like zeldovich.cpp processes ppd×ppd)
        // So we DO need to skip for missing MAX_PPD grid points to match zeldovich behavior
        if (N < MAX_PPD) {
            #if PARALLELIZE_XZ_WITHIN_SLICE
            // Parallel case: calculate total skip upfront
            // For self-conjugate, we iterate z = 0 to N-1, x = 0 to N-1 (full N×N plane)
            // We cross z boundary at z == Nhalf + 1: skip (MAX_PPD - N) * MAX_PPD
            // We cross x boundary at x == Nhalf + 1 for each z: skip (MAX_PPD - N) per z
            // Total: (MAX_PPD - N) * MAX_PPD + N * (MAX_PPD - N) = (MAX_PPD - N) * (MAX_PPD + N)
            nskip = (int64_t)(MAX_PPD - N) * (int64_t)MAX_PPD + (int64_t)N * (int64_t)(MAX_PPD - N);
            #else
            // Sequential case: track incrementally during iteration (like zeldovich.cpp)
            nskip = 0;  // Will be accumulated during loops
            #endif
        } else {
            nskip = 0;
        }
        
        // In parallel mode, advance upfront skip before starting parallel loops
        #if PARALLELIZE_XZ_WITHIN_SLICE
        if (N < MAX_PPD && nskip > 0) {
            // Advance upfront skip before parallel loops to maintain RNG consistency
            if (ps_handle != NULL && params_handle != NULL) {
                int64_t rng_index = global_y;
                zeldovich_ps_advance_rng(ps_handle, params_handle, rng_index, nskip);
                #if VERIFY_RNG_CALLS
                total_rng_skips += nskip;
                #endif
            } else if (params_handle == NULL) {
                advance_pcg_global(global_y, 2 * nskip);
                #if VERIFY_RNG_CALLS
                total_rng_skips += nskip;
                #endif
            }
            nskip = 0;  // Reset after advancing
        }
        #endif
        
        #if USE_ZELDOVICH_METHOD
        // Zeldovich method: Fill half the plane, mirror the rest
        #if DEBUG_PRINTS
        if (global_y == 0) {
            fprintf(stderr, "[Y0-DEBUG] Y=0 detected as self-conjugate, using zeldovich_method (USE_ZELDOVICH_METHOD=%d)\n", USE_ZELDOVICH_METHOD);
            fflush(stderr);
        }
        #endif
        #if PARALLELIZE_XZ_WITHIN_SLICE
        // Parallel z-loop: requires locks for thread-safe RNG access
        #pragma omp parallel for
        #else
        #endif
        // Process FULL plane first (like zeldovich.cpp): z = 0 to N-1, x = 0 to N-1
        for (int z = 0; z < N; z++) {
            // RNG skipping: match zeldovich.cpp - skip at Nyquist boundary (z == Nhalf + 1)
            // This applies to ALL Y slices including self-conjugate (Y=0 and Y=N/2)
            // This matches zeldovich.cpp line 338: skip at z == ppdhalf + 1
            #if !PARALLELIZE_XZ_WITHIN_SLICE
            if (z == Nhalf + 1 && N < MAX_PPD) {
                // Skip ALL missing z-rows: from z=N to z=MAX_PPD-1
                // Use explicit int64_t casts to avoid integer overflow
                int64_t skip_amount = (int64_t)(MAX_PPD - N) * (int64_t)MAX_PPD;
                nskip += skip_amount; // skip missing z-rows
                
                #if DEBUG_RNG_SKIP
                int log_skip = (global_y <= MAX_DEBUG_COORD) || (global_y == Nhalf - 1);
                if (log_skip) {
                    fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d z=%d: ACCUMULATE skip at Nyquist boundary (self-conj, full plane): +%lld (missing z-rows), total nskip=%lld\n",
                            N, global_y, z, (long long)skip_amount, (long long)nskip);
                    fflush(stderr);
                }
                #endif
            }
            #endif
            
            for (int x = 0; x < N; x++) {
                // RNG skipping: match zeldovich.cpp - skip at Nyquist boundary (x == Nhalf + 1)
                // This applies to ALL Y slices including self-conjugate (Y=0 and Y=N/2)
                // This matches zeldovich.cpp line 344: skip at x == ppdhalf + 1
                #if !PARALLELIZE_XZ_WITHIN_SLICE
                if (x == Nhalf + 1 && N < MAX_PPD) {
                    // Skip ALL missing x-values: from x=N to x=MAX_PPD-1
                    int64_t skip_amount = (int64_t)(MAX_PPD - N);
                    nskip += skip_amount; // skip missing x-values in this z-row
                    
                    #if DEBUG_RNG_SKIP
                    int log_skip = (z <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD) ||
                                   (z == Nhalf - 1 && global_y == Nhalf - 1);
                    if (log_skip) {
                        fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d z=%d x=%d: ACCUMULATE skip at Nyquist boundary (self-conj, full plane): +%lld (missing x-values), total nskip=%lld\n",
                                N, global_y, z, x, (long long)skip_amount, (long long)nskip);
                        fflush(stderr);
                    }
                    #endif
                }
                #endif
                
                // Calculate k-vector components first (needed for both cgauss and uniform RNG)
                int kx = (x > Nhalf) ? x - N : x;
                int ky = (global_y > Nhalf) ? global_y - N : global_y;
                int kz = (z > Nhalf) ? z - N : z;
                int k2_int = kx*kx + ky*ky + kz*kz;  // Integer k² for k_cutoff comparison
                double k2 = (double)k2_int;  // Floating-point k² for power spectrum
                
                // Nyquist frequency zeroing: Force Nyquist elements to zero for all three axes
                // This matches zeldovich.cpp line 354: abs(kx)==kmax || abs(kz)==kmax || abs(ky)==kmax
                // The Nyquist frequency (k = N/2) is self-conjugate and doesn't have a separate partner.
                // Due to the Y-shift in the reflected shell, we need to zero these to align mirroring expectations.
                int abs_kx = (kx < 0) ? -kx : kx;
                int abs_ky = (ky < 0) ? -ky : ky;
                int abs_kz = (kz < 0) ? -kz : kz;
                int is_nyquist = (abs_kx == Nhalf || abs_ky == Nhalf || abs_kz == Nhalf);
                
                // Generate D using RNG or cgauss()
                // For self-conjugate slices, we still need to check for power spectrum mode
                fftw_complex D;
                if ((k2 == 0.0) || (is_nyquist) || (!CornerModes && (double)k2_int >= k2_cutoff)) {
                    // Zero D for: DC mode, Nyquist frequency, or k_cutoff filtering
                    // This matches zeldovich.cpp line 360-364: zeroing conditions
                    D[0] = D[1] = 0.0;
                    // RNG consistency: When D=0, we skip the RNG call, so accumulate skip
                    // This matches zeldovich.cpp line 361: nskip++ (accumulate, don't advance immediately)
                    #if !PARALLELIZE_XZ_WITHIN_SLICE
                    nskip++;  // Accumulate skip, will be applied before next cgauss() call
                    #if DEBUG_RNG_SKIP
                    int log_skip = (x <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD && z <= MAX_DEBUG_COORD) ||
                                   (x == Nhalf - 1 && global_y == Nhalf - 1 && z == Nhalf - 1);
                    if (log_skip) {
                        if (k2 == 0.0) {
                            fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d (x,z)=(%d,%d): D=0 (DC mode, self-conj), ACCUMULATE nskip++ (total=%lld)\n",
                                    N, global_y, x, z, (long long)nskip);
                        } else if (is_nyquist) {
                            fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d (x,z)=(%d,%d): D=0 (Nyquist: kx=%d ky=%d kz=%d, self-conj), ACCUMULATE nskip++ (total=%lld)\n",
                                    N, global_y, x, z, kx, ky, kz, (long long)nskip);
                        } else {
                            fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d (x,z)=(%d,%d): D=0 (k_cutoff: k2=%d >= %.1f, self-conj), ACCUMULATE nskip++ (total=%lld)\n",
                                    N, global_y, x, z, k2_int, k2_cutoff, (long long)nskip);
                        }
                        fflush(stderr);
                    }
                    #endif
                    #endif
                } else if (ps_handle != NULL && params_handle != NULL) {
                    // v15.2: Use zeldovich-PLT power spectrum-weighted Gaussian
                    // If ps_handle is available, we ALWAYS use power spectrum mode (cgauss)
                    // Convert k indices to physical wavenumber: k_phys = k_index * fundamental
                    double fundamental = zeldovich_params_get_fundamental(params_handle);
                    double k2_phys = k2 * fundamental * fundamental;
                    double kmag = sqrt(k2_phys);
                    
                    // zeldovich_ps_cgauss returns double precision, convert to real_t
                    // zeldovich-PLT's v2rng array is sized to ppd/2, so valid indices are 0 to (N/2 - 1)
                    // Since we've already handled global_y == N/2 above, global_y is now < N/2
                    int64_t rng_index = global_y;
                    
                    // Advance zeldovich-PLT's RNG when crossing Nyquist boundaries
                    #if !PARALLELIZE_XZ_WITHIN_SLICE
                    if (nskip > 0) {
                        zeldovich_ps_advance_rng(ps_handle, params_handle, rng_index, nskip);
                        #if VERIFY_RNG_CALLS
                        total_rng_skips += nskip;
                        #endif
                        nskip = 0;  // Reset after advancing
                    }
                    #endif
                    
                    double D_real, D_imag;
                    zeldovich_ps_cgauss(ps_handle, kmag, rng_index, &D_real, &D_imag);
                    #if VERIFY_RNG_CALLS
                    total_rng_calls++;  // Each cgauss() call uses 2 random numbers
                    #endif
                    D[0] = (real_t)D_real;
                    D[1] = (real_t)D_imag;
                    // Note: For VERIFY_HERMITIAN_SYMMETRY==2, we set D=0 AFTER computing F and H
                } else if (ps_params != NULL) {
                    // Legacy: Use standalone power spectrum-weighted Gaussian (cgauss)
                    double fundamental = 1.0;  // TODO: Add to power_spectrum_params_t
                    double k2_phys = k2 * fundamental * fundamental;
                    double kmag = sqrt(k2_phys);
                    
                    // RNG consistency: Advance our own RNG when crossing Nyquist boundaries
                    // This matches the zeldovich-PLT branch behavior
                    #if !PARALLELIZE_XZ_WITHIN_SLICE
                    if (nskip > 0 && N < MAX_PPD) {
                        advance_pcg_global(global_y, 2 * nskip);
                        #if VERIFY_RNG_CALLS
                        total_rng_skips += nskip;
                        #endif  // Each complex number uses 2 random numbers
                        nskip = 0;  // Reset after advancing
                    }
                    #endif
                    
                    cgauss(ps_params, kmag, global_y, &D);
                    #if VERIFY_RNG_CALLS
                    total_rng_calls++;  // Each cgauss() call uses 2 random numbers
                    #endif
                    // Note: For VERIFY_HERMITIAN_SYMMETRY==2, we set D=0 AFTER computing F and H
                } else {
                    // Fallback: use uniform random numbers (white noise mode, no power spectrum)
                    // This should only happen when ps_handle is NULL (no parameter file provided)
                    // Note: If ps_handle is available, we should have used cgauss() above
                    
                    // RNG consistency: Advance our own RNG when crossing Nyquist boundaries
                    #if !PARALLELIZE_XZ_WITHIN_SLICE
                    if (nskip > 0 && N < MAX_PPD) {
                        advance_pcg_global(global_y, 2 * nskip);
                        #if VERIFY_RNG_CALLS
                        total_rng_skips += nskip;
                        #endif  // Each complex number uses 2 random numbers
                        nskip = 0;  // Reset after advancing
                    }
                    #endif
                    
                    double D_re = random_real_pcg_global(global_y);
                    double D_im = random_real_pcg_global(global_y);
                    #if VERIFY_RNG_CALLS
                    total_rng_calls++;  // Two random_real calls = 1 complex number = 2 random numbers
                    #endif
                    D[0] = D_re;
                    D[1] = D_im;
                    // Note: For VERIFY_HERMITIAN_SYMMETRY==2, we set D=0 AFTER computing F and H
                }
                
                // Compute F, G, H from D
                // Skip F, G, H computation in density-only mode (qdensity == 2)
                fftw_complex F, G, H;
                
                // Declare f_sc before if/else blocks (needed for velocity arrays in all modes)
                double f_sc = 1.0;
                
                // Match zeldovich.cpp: set f=0 when D==0 (line 443)
                if (D[0] == 0.0 && D[1] == 0.0) {
                    f_sc = 0.0;
                }
                
                if (just_density) {
                    // Density-only mode: Set F, G, H to zero (not used)
                    F[0] = F[1] = G[0] = G[1] = H[0] = H[1] = 0.0;

                    // Density-only mode: Compute F from D (for D+i*F test)
                    // Match zeldovich behavior: only compute F when D != 0
                    // if (D[0] != 0.0 || D[1] != 0.0) {
                    //     // Avoid division by zero: if k2==0, set F=0 (D is already 0 for k2==0)
                    //     if (k2 == 0.0) {
                    //         F[0] = F[1] = 0.0;
                    //     } else {
                    //         double fundamental_sc = 1.0;
                    //         if (params_handle != NULL) {
                    //             fundamental_sc = zeldovich_params_get_fundamental(params_handle);
                    //         }
                    //         // For VERIFY_HERMITIAN_SYMMETRY=2, we don't use PLT rescaling, so rescale=1.0
                    //         double rescale_sc = 1.0;
                    //         // Compute factor the same way as normal mode
                    //         double factor_sc = rescale_sc / (k2 * fundamental_sc);
                            
                    //         // Now compute F and H using the same formula as normal mode
                    //         F[0] = -kx * factor_sc * D[1];
                    //         F[1] =  kx * factor_sc * D[0];
                    //     }
                    // } else {
                    //     // D == 0, so F = 0 (match zeldovich line 442)
                    //     F[0] = F[1] = 0.0;
                    // }

                    // // D+i*F test: Store D + i*F (D should NOT be overwritten)
                    // // D[0] = F[0]; // REMOVED: This was overwriting D with F (wrong!)
                    // // D[1] = F[1]; // REMOVED: This was overwriting D with F (wrong!)
                    // // F[0] = F[1] = G[0] = G[1] = H[0] = H[1] = 0.0; // Test D+iF calculation
                    // G[0] = G[1] = H[0] = H[1] = 0.0;
                } else if (k2 == 0.0) {
                    F[0] = F[1] = G[0] = G[1] = H[0] = H[1] = 0.0;
                    f_sc = 0.0;  // Match zeldovich.cpp: set f=0 when k2==0 (D==0)
                } else {
                    double ik2 = 1.0 / k2;
                    
                    #if defined(VERIFY_HERMITIAN_SYMMETRY) && VERIFY_HERMITIAN_SYMMETRY == 1
                    // Verification mode 1: Set F=0 and H=0 to test Hermitian symmetry
                    // With F=0 and H=0, only D and G remain, and result should be purely real
                    F[0] = 0.0;
                    F[1] = 0.0;
                    G[0] = -ky * ik2 * D[1];
                    G[1] =  ky * ik2 * D[0];
                    H[0] = 0.0;
                    H[1] = 0.0;
                    #elif defined(VERIFY_HERMITIAN_SYMMETRY) && VERIFY_HERMITIAN_SYMMETRY == 2
                    // Verification mode 2: Compute F and H from D, then set D=0 and G=0
                    // This makes Array 0 and Array 1 purely imaginary (real parts = 0)
                    // After 3D FFT, the result should be purely imaginary (real parts = 0)
                    // IMPORTANT: Must use the same factor computation as normal mode to ensure
                    // F is computed identically. Otherwise the test is invalid.
                    // Get fundamental wavenumber (needed for factor computation)
                    double fundamental_sc = 1.0;
                    if (params_handle != NULL) {
                        fundamental_sc = zeldovich_params_get_fundamental(params_handle);
                    }
                    // For VERIFY_HERMITIAN_SYMMETRY=2, we don't use PLT rescaling, so rescale=1.0
                    double rescale_sc = 1.0;
                    // Compute factor the same way as normal mode
                    double factor_sc = rescale_sc / (k2 * fundamental_sc);
                    
                    // Now compute F and H using the same formula as normal mode
                    F[0] = -kx * factor_sc * D[1];
                    F[1] =  kx * factor_sc * D[0];
                    
                    G[0] = 0.0;  // Set G=0
                    G[1] = 0.0;
                    
                    H[0] = -kz * factor_sc * D[1];
                    H[1] =  kz * factor_sc * D[0];
                    
                    // Now set D=0 (after F and H are computed)
                    D[0] = 0.0;
                    D[1] = 0.0;
                    #else
                    // Normal operation: Compute F, G, H from D
                    // double ik2 = 1.0 / k2; --> UNCOMMENT LATER! AFTER D_0 TEST
                    // Check if PLT is enabled
                    int qPLT_sc = 0;
                    // Default fundamental = 1.0 
                    double fundamental_sc = 1.0;
                    eigenmode e_sc;
                    int use_plt_sc = 0;
                    
                    // ========== STEP 3.5: Compute PLT growth rate f and rescale (before computing F, G, H) ==========
                    // f is the logarithmic derivative of the growth factor that scales velocities
                    // When PLT is enabled: f = (sqrt(1. + 24 * e.val * f_cluster) - 1) / 4.
                    // When PLT is not enabled: f = 1.0 (default)
                    // Skip in density-only mode (qdensity == 2)
                    // f_sc already declared above, update it here if PLT is enabled
                    double rescale_sc = 1.0;
                    
                    // Get fundamental wavenumber (needed for both PLT and non-PLT cases)
                    if (!just_density && params_handle != NULL) {
                        fundamental_sc = zeldovich_params_get_fundamental(params_handle);
                        qPLT_sc = zeldovich_params_get_qPLT(params_handle);
                        if (qPLT_sc) {
                            // Get PLT eigenmode for this k-vector
                            // Convert kx, ky, kz to array indices for plt_get_eigenmode
                            int ikx = (kx < 0) ? N + kx : kx;
                            int iky = (ky < 0) ? N + ky : ky;
                            int ikz = (kz < 0) ? N + kz : kz;
                            // Handle z index (only positive half-space stored)
                            if (ikz > N / 2) ikz = N - ikz;
                            
                            if (plt_get_eigenmode(ikx, iky, ikz, (int64_t)N, &e_sc) == 0) {
                                use_plt_sc = 1;
        
                                // Eigenval normalization
                                if (kz < 0) {
                                    e_sc.vec[2] = -e_sc.vec[2];
                                }
                                
                                // Normalize eigenvector (interpolation might not preserve |e_sc| = 1)
                                double e_sc_mag = sqrt(e_sc.vec[0] * e_sc.vec[0] + e_sc.vec[1] * e_sc.vec[1] + e_sc.vec[2] * e_sc.vec[2]);
                                if (e_sc_mag > 0.0) {
                                    e_sc.vec[0] /= e_sc_mag;
                                    e_sc.vec[1] /= e_sc_mag;
                                    e_sc.vec[2] /= e_sc_mag;
                                }
                                
                                // Apply normalization: norm = k2 / (k * e_sc)
                                double k_dot_e_sc = kx * e_sc.vec[0] + ky * e_sc.vec[1] + kz * e_sc.vec[2];
                                double norm_sc = (k2 > 0.0 && k_dot_e_sc != 0.0) ? k2 / k_dot_e_sc : 0.0;
                                if (!isfinite(norm_sc)) norm_sc = 0.0;
                                
                                // Scale eigenvector by norm
                                e_sc.vec[0] *= norm_sc;
                                e_sc.vec[1] *= norm_sc;
                                e_sc.vec[2] *= norm_sc;
                                
                                // Compute f and rescale before computing F, G, H
                                double f_cluster = zeldovich_params_get_f_cluster(params_handle);
                                // PLT growth rate: f = (sqrt(1. + 24 * e.val * f_cluster) - 1) / 4.
                                f_sc = (sqrt(1. + 24. * e_sc.val * f_cluster) - 1.) * 0.25;
                                
                                // Compute rescaling if qPLTrescale is enabled
                                if (qPLTrescale) {
                                    double plt_f = f_sc;  // PLT growth rate for this mode
                                    rescale_sc = pow(a_NL / a0, target_f - plt_f);
                                }
                            } else {
                                // If eigenmode lookup fails, fall back to normal computation
                                if (rank == 0 && x == 0 && z == 0) {
                                    fprintf(stderr, "[WARNING] Failed to get PLT eigenmode for (kx=%d, ky=%d, kz=%d), using normal computation\n",
                                            kx, ky, kz);
                                }
                            }
                        }
                    }
                    
                    // Compute factor (used identically in both PLT and non-PLT cases)
                    // In zeldovich.cpp: k2 includes fundamental^2, so factor = rescale / (k2 * fundamental)
                    // where vec is either e.vec[i] (PLT) or k[i] (non-PLT)
                    double factor = rescale_sc / (k2 * fundamental_sc);
                    
                    if (use_plt_sc) {
                        // PLT mode: Use eigenvector instead of k-vector
                        // F = rescale * i * e.vec[0] * fundamental * ik2 * D
                        F[0] = -e_sc.vec[0] * factor * D[1];
                        F[1] =  e_sc.vec[0] * factor * D[0];
                        
                        G[0] = -e_sc.vec[1] * factor * D[1];
                        G[1] =  e_sc.vec[1] * factor * D[0];
                        
                        H[0] = -e_sc.vec[2] * factor * D[1];
                        H[1] =  e_sc.vec[2] * factor * D[0];
                    } else {
                        // Normal operation: Compute F, G, H from D using k-vector
                        // In zeldovich.cpp: k2 includes fundamental^2, so F = I * kx * fundamental * ik2 * D
                        // In our code: k2 doesn't include fundamental^2, so we need to multiply by fundamental
                        F[0] = -kx * factor * D[1];
                        F[1] =  kx * factor * D[0];
                        G[0] = -ky * factor * D[1];
                        G[1] =  ky * factor * D[0];
                        H[0] = -kz * factor * D[1];
                        H[1] =  kz * factor * D[0];
                    }
                    
                    // Print factor debug info while variables are in scope (self-conjugate path)
                    #if DEBUG_RNG_CONSISTENCY
                    int boundary_coord_sc = Nhalf - 1;
                    int test_boundary_sc = (boundary_coord_sc >= 0 && 
                                             x == boundary_coord_sc && global_y == boundary_coord_sc && z == boundary_coord_sc);
                    int effective_boundary_sc = (boundary_coord_sc <= MAX_DEBUG_BOUNDARY_COORD) ? boundary_coord_sc : MAX_DEBUG_BOUNDARY_COORD;
                    int max_test_coord_sc = (MAX_DEBUG_COORD > effective_boundary_sc) ? MAX_DEBUG_COORD : effective_boundary_sc;
                    int in_test_range_sc = (x <= max_test_coord_sc && global_y <= max_test_coord_sc && z <= max_test_coord_sc);
                    int should_print_factor_sc = 0;
                    if (test_boundary_sc) {
                        should_print_factor_sc = 1;
                    } else if (N <= DEBUG_FULL_PRINT_MAX_N) {
                        should_print_factor_sc = in_test_range_sc;
                    } else {
                        if (x <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD && z <= MAX_DEBUG_COORD) {
                            should_print_factor_sc = 1;
                        } else if (in_test_range_sc) {
                            should_print_factor_sc = (x % DEBUG_SAMPLE_STRIDE == 0 &&
                                                       global_y % DEBUG_SAMPLE_STRIDE == 0 &&
                                                       z % DEBUG_SAMPLE_STRIDE == 0);
                        }
                    }
                    if (should_print_factor_sc) {
                        if (use_plt_sc && params_handle != NULL) {
                            fprintf(stderr, "[FACTOR-DEBUG] N=%d Y=%d (x,z)=(%d,%d): fundamental=%.10e rescale=%.10e "
                                    "e.vec=[%.10e,%.10e,%.10e] e.val=%.10e factor=%.10e\n",
                                    N, global_y, x, z, fundamental_sc, rescale_sc,
                                    e_sc.vec[0], e_sc.vec[1], e_sc.vec[2], e_sc.val, factor);
                        } else {
                            fprintf(stderr, "[FACTOR-DEBUG] N=%d Y=%d (x,z)=(%d,%d): fundamental=%.10e rescale=%.10e "
                                    "k=[%d,%d,%d] factor=%.10e\n",
                                    N, global_y, x, z, fundamental_sc, rescale_sc, kx, ky, kz, factor);
                        }
                        fflush(stderr);
                    }
                    #endif
                    #endif
                }
                
                // ========== DEBUG: Print RNG values for consistency checking ==========
                #if DEBUG_RNG_CONSISTENCY
                // Print D, F, G, H for test coordinates to verify RNG consistency across N
                // Tests all coordinates (x,y,z) where x,y,z <= MAX_DEBUG_COORD
                // Also tests coordinates at (N/2)-1 to ensure boundary of overlapping region is tested
                // This ensures overlapping regions are tested for any pair of N values
                // where min(N1, N2) >= 2*MAX_DEBUG_COORD
                // For N=4 and N=6: overlapping region is x,y,z <= 2 (Nhalf for N=4)
                //   (N/2)-1 = 1 for N=4, so (1,1,1) should match between N=4 and N=6
                // For N=8 and N=10: overlapping region is x,y,z <= 4 (Nhalf for N=8)
                //   (N/2)-1 = 3 for N=8, so (3,3,3) should match between N=8 and N=10
                int boundary_coord = Nhalf - 1;  // (N/2)-1 for current N
                // Test coordinate (N/2)-1, (N/2)-1, (N/2)-1) to verify boundary of overlapping region
                // This ensures that when comparing N1 and N2, the boundary coordinate of the
                // smaller N (which is in the overlapping region) is tested for both N values
                int test_boundary = (boundary_coord >= 0 && 
                                     x == boundary_coord && global_y == boundary_coord && z == boundary_coord);
                // Test coordinates up to max(MAX_DEBUG_COORD, min((N/2)-1, MAX_DEBUG_BOUNDARY_COORD))
                int effective_boundary = (boundary_coord <= MAX_DEBUG_BOUNDARY_COORD) ? boundary_coord : MAX_DEBUG_BOUNDARY_COORD;
                int max_test_coord = (MAX_DEBUG_COORD > effective_boundary) ? MAX_DEBUG_COORD : effective_boundary;
                // Check if coordinate is in test range
                int in_test_range = (x <= max_test_coord && global_y <= max_test_coord && z <= max_test_coord);
                // For large N, use sampling to reduce output volume
                // Always print boundary coordinates and coordinates in small N range
                int should_print = 0;
                if (test_boundary) {
                    // Always print boundary coordinate (N/2)-1, (N/2)-1, (N/2)-1
                    should_print = 1;
                } else if (N <= DEBUG_FULL_PRINT_MAX_N) {
                    // For small N, print all coordinates in test range
                    should_print = in_test_range;
                } else {
                    // For large N, print only sampled coordinates (multiples of stride)
                    // Also always print coordinates <= MAX_DEBUG_COORD (small coordinate cube)
                    if (x <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD && z <= MAX_DEBUG_COORD) {
                        should_print = 1;  // print small coordinate cube
                    } else if (in_test_range) {
                        // Sample: only print if x, y, z are multiples of stride
                        should_print = (x % DEBUG_SAMPLE_STRIDE == 0 &&
                                       global_y % DEBUG_SAMPLE_STRIDE == 0 &&
                                       z % DEBUG_SAMPLE_STRIDE == 0);
                    }
                }
                if (should_print) {
                    fprintf(stderr, "[RNG-DEBUG] N=%d Y=%d (x,z)=(%d,%d): k=(%d,%d,%d) k2=%.6f | "
                            "D=(%.10e,%.10e) F=(%.10e,%.10e) G=(%.10e,%.10e) H=(%.10e,%.10e)\n",
                            N, global_y, x, z, kx, ky, kz, k2,
                            D[0], D[1], F[0], F[1], G[0], G[1], H[0], H[1]);
                    fflush(stderr);
                }
                #endif
                
                // Compute mirror indices for self-conjugate slices (like zeldovich's zHer, xHer)
                int x_mirror = (x == 0) ? 0 : N - x;
                int z_mirror = (z == 0) ? 0 : N - z;
                
                // Store in arrays (Zeldovich packing scheme)
                // For self-conjugate slices: store D+i*F in primary slice,
                // and conj(D)+i*conj(F) in conjugate slice at mirrored positions
                // (like zeldovich stores in slab and slabHer)
                if (just_density) {
                    // Density-only mode: Only store D (density) in Array 0
                    PRIM_SLICE(0, x, z)[0] = D[0];  // Real = D_re
                    PRIM_SLICE(0, x, z)[1] = D[1];  // Imag = D_im

                    CONJ_SLICE(0, x_mirror, z_mirror)[0] = D[0];  // Real = D_re
                    CONJ_SLICE(0, x_mirror, z_mirror)[1] = -D[1];  // Imag = -D_im (conjugate)

                    // Test D+iF calculation
                    // PRIM_SLICE(0, x, z)[0] = D[0] - F[1];  // Real = D_re - F_im
                    // PRIM_SLICE(0, x, z)[1] = D[1] + F[0];  // Imag = D_im + F_re
                    
                    // Conjugate slice (at mirrored position): conj(D) + i*conj(F) = (D[0] + F[1]) + i*(F[0] - D[1])
                    // CONJ_SLICE(0, x_mirror, z_mirror)[0] = D[0] + F[1];  // Real = D_re + F_im
                    // CONJ_SLICE(0, x_mirror, z_mirror)[1] = F[0] - D[1];  // Imag = F_re - D_im
                } else {
                    // Normal mode: Store D+iF, G+iH, and optionally velocities
                    // Primary slice: Array 0: D + i*F = (D[0] - F[1]) + i*(D[1] + F[0])
                    PRIM_SLICE(0, x, z)[0] = D[0] - F[1];
                    PRIM_SLICE(0, x, z)[1] = D[1] + F[0];
                    
                    // Array 1: G + i*H = (G[0] - H[1]) + i*(G[1] + H[0])
                    PRIM_SLICE(1, x, z)[0] = G[0] - H[1];
                    PRIM_SLICE(1, x, z)[1] = G[1] + H[0];
                    
                    // Conjugate slice (at mirrored position): conj(D)+i*conj(F), conj(G)+i*conj(H)
                    // conj(D) + i*conj(F) = (D[0] + F[1]) + i*(F[0] - D[1])
                    CONJ_SLICE(0, x_mirror, z_mirror)[0] = D[0] + F[1];  // Real = D_re + F_im
                    CONJ_SLICE(0, x_mirror, z_mirror)[1] = F[0] - D[1];  // Imag = F_re - D_im
                    
                    // conj(G) + i*conj(H) = (G[0] + H[1]) + i*(H[0] - G[1])
                    CONJ_SLICE(1, x_mirror, z_mirror)[0] = G[0] + H[1];  // Real = G_re + H_im
                    CONJ_SLICE(1, x_mirror, z_mirror)[1] = H[0] - G[1];  // Imag = H_re - G_im
                    
                    if (narray >= 4) {
                        // Primary slice: Array 2: 0 + i*F*f = -F[1]*f + i*(F[0]*f)
                        // f_sc is computed above (PLT growth rate if PLT enabled, else 1.0)
                        PRIM_SLICE(2, x, z)[0] = -F[1] * f_sc;
                        PRIM_SLICE(2, x, z)[1] = F[0] * f_sc;
                        // Array 3: G*f + i*H*f = (G[0] - H[1])*f + i*((G[1] + H[0])*f)
                        PRIM_SLICE(3, x, z)[0] = (G[0] - H[1]) * f_sc;
                        PRIM_SLICE(3, x, z)[1] = (G[1] + H[0]) * f_sc;
                        
                        // Conjugate slice (at mirrored position): i*conj(F*f), conj(G*f)+i*conj(H*f)
                        // i*conj(F*f) = i*(F_re*f - i*F_im*f) = F_im*f + i*F_re*f
                        CONJ_SLICE(2, x_mirror, z_mirror)[0] = F[1] * f_sc;   // Real = F_im * f
                        CONJ_SLICE(2, x_mirror, z_mirror)[1] = F[0] * f_sc;   // Imag = F_re * f
                        
                        // conj(G*f) + i*conj(H*f) = (G[0] + H[1])*f + i*((H[0] - G[1])*f)
                        CONJ_SLICE(3, x_mirror, z_mirror)[0] = (G[0] + H[1]) * f_sc;  // Real = (G_re + H_im) * f
                        CONJ_SLICE(3, x_mirror, z_mirror)[1] = (H[0] - G[1]) * f_sc;  // Imag = (H_re - G_im) * f
                    }
                }
                
            }
        }
        
        // Post-processing: Mirror first half to second half (match zeldovich.cpp lines 555-573)
        // This is done AFTER processing the full plane, matching zeldovich.cpp behavior
        // zeldovich.cpp: for (z = 0; z < ppdhalf; z++) where ppdhalf = N/2
        // Apply to ALL self-conjugate slices (Y=0 and Y=N/2) since they use the same processing method
        if (global_y == 0 || global_y == Nhalf) {
            for (int z = 0; z < Nhalf; z++) {
                int z_mirror = (z == 0) ? 0 : N - z;
                // Match zeldovich.cpp: xmax = ppdhalf for z=0, ppd for z>0
                // zeldovich.cpp: int xmax = (z == 0 ? ppdhalf : ppd);
                int x_max = (z == 0 ? Nhalf : N);
                for (int x = 0; x < x_max; x++) {
                    int x_mirror = (x == 0) ? 0 : N - x;
                    // Mirror by taking complex conjugate (match zeldovich.cpp line 566-567)
                    // zeldovich.cpp copies from slabHer (which contains conjugates) to slab
                    // For self-conjugate slices, we need to conjugate when mirroring to preserve Hermitian symmetry
                    // f(kx, 0, kz) = conj(f(-kx, 0, -kz)) for self-conjugate slice at Y=0 (normal operation and Mode 1)
                    // For Mode 2 (purely imaginary result), we need anti-Hermitian: f(kx, 0, kz) = -conj(f(-kx, 0, -kz))
                    #if defined(VERIFY_HERMITIAN_SYMMETRY) && VERIFY_HERMITIAN_SYMMETRY == 2
                    // Mode 2: Anti-Hermitian symmetry for purely imaginary result
                    // f(-k) = -conj(f(k)) means we negate the conjugate
                    if (x == 0 && z == 0) {
                        fprintf(stderr, "[MIRROR-DEBUG] Y=%d: Mode 2 branch executing (VERIFY_HERMITIAN_SYMMETRY=2)\n", global_y);
                        fflush(stderr);
                    }
                    for (int a = 0; a < narray; a++) {
                        // Negated complex conjugate: -(a + i*b)* = -(a - i*b) = -a + i*b
                        PRIM_SLICE(a, x_mirror, z_mirror)[0] = CONJ_SLICE(a, x_mirror, z_mirror)[0];  // -PRIM_SLICE(a, x, z)[0];  // Real part (negated)
                        PRIM_SLICE(a, x_mirror, z_mirror)[1] = CONJ_SLICE(a, x_mirror, z_mirror)[1]; // PRIM_SLICE(a, x, z)[1];   // Imaginary part (same)
                    }
                    #else
                    // Normal operation and Mode 1: Hermitian symmetry for purely real result
                    // Match zeldovich's scheme: copy conj(D)+i*conj(F) from conjugate_slices to primary_slices
                    // This copies the conjugate values that were stored in conjugate_slices during main loop
                    // (like zeldovich copies from slabHer to slab at lines 560-561)
                    if (x == 0 && z == 0) {
                        #if defined(VERIFY_HERMITIAN_SYMMETRY) && VERIFY_HERMITIAN_SYMMETRY == 1
                        fprintf(stderr, "[MIRROR-DEBUG] Y=%d: Mode 1 branch executing (VERIFY_HERMITIAN_SYMMETRY=1)\n", global_y);
                        #else
                        fprintf(stderr, "[MIRROR-DEBUG] Y=%d: Normal operation branch executing (VERIFY_HERMITIAN_SYMMETRY not defined or != 1,2)\n", global_y);
                        #endif
                        fflush(stderr);
                    }
                    for (int a = 0; a < narray; a++) {
                        // Copy from conjugate slice (which contains conj(D)+i*conj(F) at mirror positions)
                        // to primary slice at mirror positions
                        PRIM_SLICE(a, x_mirror, z_mirror)[0] = CONJ_SLICE(a, x_mirror, z_mirror)[0];  // Real part
                        PRIM_SLICE(a, x_mirror, z_mirror)[1] = CONJ_SLICE(a, x_mirror, z_mirror)[1];  // Imaginary part
                    }
                    #endif
                }
            }
            // Set origin to zero (match zeldovich.cpp line 572)
            // This applies to both Y=0 and Y=N/2 self-conjugate slices
            for (int a = 0; a < narray; a++) {
                PRIM_SLICE(a, 0, 0)[0] = 0.0;
                PRIM_SLICE(a, 0, 0)[1] = 0.0;
            }
        }
        // Set imaginary parts to 0 at special points
        for (int a = 0; a < narray; a++) {
            PRIM_SLICE(a, 0, 0)[1] = 0.0;
            PRIM_SLICE(a, Nhalf, 0)[1] = 0.0;
            PRIM_SLICE(a, 0, Nhalf)[1] = 0.0;
            PRIM_SLICE(a, Nhalf, Nhalf)[1] = 0.0;
        }
        // OBSOLETE CODE REMOVED: The non-zeldovich_method path has been removed
        // because zeldovich_method is the only supported path. The obsolete code
        // (previously at lines 1325-1734) was buggy and is no longer needed.
        #endif  // End of USE_ZELDOVICH_METHOD conditional
        // ========== DEBUG: Check nskip value immediately after #endif ==========
        #if DEBUG_RNG_SKIP
        if (global_y <= MAX_DEBUG_COORD || global_y == Nhalf - 1) {
            fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d: IMMEDIATELY AFTER #endif (line 1527) - nskip=%lld (0x%llx)\n",
                    N, global_y, (long long)nskip, (unsigned long long)nskip);
            fflush(stderr);
        }
        #endif
        // ========== DEBUG: Check nskip value right before loop ends ==========
        #if DEBUG_RNG_SKIP
        if (global_y <= MAX_DEBUG_COORD || global_y == Nhalf - 1) {
            fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d: RIGHT BEFORE LOOP ENDS (inside #endif) - nskip=%lld (0x%llx)\n",
                    N, global_y, (long long)nskip, (unsigned long long)nskip);
            fflush(stderr);
        }
        #endif
    }
    
    // ========== DEBUG: Check nskip value right after loop ends ==========
    #if DEBUG_RNG_SKIP
    if (global_y <= MAX_DEBUG_COORD || global_y == Nhalf - 1) {
        fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d: AFTER LOOP ENDS - nskip=%lld (0x%llx)\n",
                N, global_y, (long long)nskip, (unsigned long long)nskip);
        fflush(stderr);
    }
    #endif
    
    // ========== DEBUG: Check values BEFORE 2D FFT (Fourier space - expect non-zero imag) ==========
    #if DEBUG_PRINTS && !SKIP_VERIFICATION
    if (rank == 0 && global_y <= 2) {
        real_t debug_max_real_before[4] = {0, 0, 0, 0};
        real_t debug_max_imag_before[4] = {0, 0, 0, 0};
        for (int a = 0; a < narray; a++) {
            for (int x = 0; x < N; x++) {
                for (int z = 0; z < N; z++) {
                    double re = fabs_t(PRIM_SLICE(a, x, z)[0]);
                    double im = fabs_t(PRIM_SLICE(a, x, z)[1]);
                    debug_max_real_before[a] = fmax_t(debug_max_real_before[a], re);
                    debug_max_imag_before[a] = fmax_t(debug_max_imag_before[a], im);
                }
            }
        }
        printf("[DEBUG Y=%d] Before 2D FFT (Fourier space) - PRIMARY slice - Max real/imag per array: ", global_y);
        for (int a = 0; a < narray; a++) {
            printf("A%d: re=%.3e im=%.3e ", a, debug_max_real_before[a], debug_max_imag_before[a]);
        }
        printf("(NOTE: Non-zero imag expected in Fourier space)\n");
        
        // Also check conjugate slice for conjugate pairs
        if (y_mirror != global_y) {
            real_t debug_max_real_conj_before[4] = {0, 0, 0, 0};
            real_t debug_max_imag_conj_before[4] = {0, 0, 0, 0};
            for (int a = 0; a < narray; a++) {
                for (int x = 0; x < N; x++) {
                    for (int z = 0; z < N; z++) {
                        double re = fabs_t(CONJ_SLICE(a, x, z)[0]);
                        double im = fabs_t(CONJ_SLICE(a, x, z)[1]);
                        debug_max_real_conj_before[a] = fmax_t(debug_max_real_conj_before[a], re);
                        debug_max_imag_conj_before[a] = fmax_t(debug_max_imag_conj_before[a], im);
                    }
                }
            }
            printf("[DEBUG Y=%d] Before 2D FFT (Fourier space) - CONJUGATE slice - Max real/imag per array: ", global_y);
            for (int a = 0; a < narray; a++) {
                printf("A%d: re=%.3e im=%.3e ", a, debug_max_real_conj_before[a], debug_max_imag_conj_before[a]);
            }
            printf("(NOTE: Should be conjugate of primary!)\n");
            
            // Check a few sample points to verify conjugate relationship
            printf("[DEBUG Y=%d] Sample conjugate check (first 3 points): ", global_y);
            for (int a = 0; a < narray && a < 2; a++) {
                for (int x = 0; x < 2 && x < N; x++) {
                    for (int z = 0; z < 2 && z < N; z++) {
                        double prim_re = PRIM_SLICE(a, x, z)[0];
                        double prim_im = PRIM_SLICE(a, x, z)[1];
                        int x_mirror = (x == 0) ? 0 : N - x;
                        int z_mirror = (z == 0) ? 0 : N - z;
                        double conj_re = CONJ_SLICE(a, x_mirror, z_mirror)[0];
                        double conj_im = CONJ_SLICE(a, x_mirror, z_mirror)[1];
                        printf("A%d[%d,%d]: prim=(%.3e,%.3e) conj[%d,%d]=(%.3e,%.3e) ", 
                               a, x, z, prim_re, prim_im, x_mirror, z_mirror, conj_re, conj_im);
                    }
                }
            }
            printf("\n");
        }
    }
    #endif
    
    // Verify Hermitian symmetry BEFORE 2D FFT
    // STAGE 7: Updated to check all arrays independently
    #if VERIFY_HERMITIAN_SYMMETRY
    for (int a = 0; a < narray; a++) {
        fftw_complex_t *prim_array = &PRIM_SLICE(a, 0, 0);
        fftw_complex_t *conj_array = (y_mirror != global_y) ? &CONJ_SLICE(a, 0, 0) : NULL;
        verify_initial_fourier_hermitian_symmetry(N, prim_array, global_y, y_mirror, conj_array, a);
    }
    #endif
    
    // DEBUG: Check if conjugate pair is correct BEFORE 2D FFT
    // STAGE 7: Updated to check all arrays independently
    #if PRINT_DETAILED_SLICES
    if (y_mirror != global_y && N <= 16) {
        printf("\n[RANK %d - BEFORE FFT] Checking Y=%d and Y=%d conjugate relationship (all arrays):\n", rank, global_y, y_mirror);
        for (int a = 0; a < narray; a++) {
            fftw_complex_t *prim_array = &PRIM_SLICE(a, 0, 0);
            fftw_complex_t *conj_array = &CONJ_SLICE(a, 0, 0);
            verify_hermitian_pair(rank, global_y, y_mirror, prim_array, conj_array, N, a, false);
        }
    }
    #endif
    
    // ========== DUMP MATRIX BEFORE FFT (Fourier space, after generation) ==========
    #ifdef DUMP_MATRIX_BEFORE_FFT
    if (N <= 16) {  // Only for small N
        char dump_filename[256];
        snprintf(dump_filename, sizeof(dump_filename), "matrix_before_fft.txt");
        FILE *dump_fp = fopen(dump_filename, "a");
        if (dump_fp) {
            // Write header on first Y-slice
            if (global_y == 0 && rank == 0) {
                fprintf(dump_fp, "# Matrix values before FFT (Fourier space)\n");
                fprintf(dump_fp, "# Format: Y=<y> X=<x> Z=<z> Array=<array_idx> Re=<real> Im=<imag>\n");
                fflush(dump_fp);
            }
            // Dump primary slice
            int count = 0;
            for (int x = 0; x < N; x++) {
                for (int z = 0; z < N; z++) {
                    for (int array_idx = 0; array_idx < narray; array_idx++) {
                        int ret = fprintf(dump_fp, "Y=%d X=%d Z=%d Array=%d Re=%.15e Im=%.15e\n",
                                global_y, x, z, array_idx, PRIM_SLICE(array_idx, x, z)[0], PRIM_SLICE(array_idx, x, z)[1]);
                        if (ret < 0 && count == 0) {
                            fprintf(stderr, "[DUMP-DEBUG] fprintf failed for Y=%d X=%d Z=%d\n", global_y, x, z);
                        }
                        count++;
                    }
                }
            }
            // Dump conjugate slice if it exists
            if (y_mirror != global_y) {
                for (int x = 0; x < N; x++) {
                    for (int z = 0; z < N; z++) {
                        for (int array_idx = 0; array_idx < narray; array_idx++) {
                            fprintf(dump_fp, "Y=%d X=%d Z=%d Array=%d Re=%.15e Im=%.15e\n",
                                    y_mirror, x, z, array_idx, CONJ_SLICE(array_idx, x, z)[0], CONJ_SLICE(array_idx, x, z)[1]);
                            count++;
                        }
                    }
                }
            }
            fflush(dump_fp);  // Ensure data is written before closing
            fclose(dump_fp);
            if (rank == 0 && global_y == 0) {
                fprintf(stderr, "[DUMP-DEBUG] Wrote %d matrix elements for Y=%d to matrix_before_fft.txt\n", count, global_y);
            }
        } else {
            // Debug: Print error if file couldn't be opened
            if (rank == 0 && global_y == 0) {
                fprintf(stderr, "[DUMP-DEBUG] Failed to open matrix_before_fft.txt for writing\n");
            }
        }
    }
    #endif
    
    // Apply 2D FFT to all arrays independently
    // Loop over all arrays and apply FFT to each
    for (int a = 0; a < narray; a++) {
        // Primary slice: array 'a'
        fftw_complex_t *prim_array_start = &PRIM_SLICE(a, 0, 0);
        
        // DEBUG: Verify buffers are separate
        #if DEBUG_PRINTS && !SKIP_VERIFICATION
        if (rank == 0 && global_y <= 2 && y_mirror != global_y && a == 0) {
            fftw_complex_t *conj_array_start = &CONJ_SLICE(a, 0, 0);
            printf("[DEBUG Y=%d] Before 2D FFT - Buffer addresses: prim=%p conj=%p (diff=%ld)\n", 
                   global_y, (void*)prim_array_start, (void*)conj_array_start, 
                   (long)(conj_array_start - prim_array_start));
            printf("[DEBUG Y=%d] Before 2D FFT - Sample values: prim[0,0]=(%.3e,%.3e) conj[0,0]=(%.3e,%.3e)\n",
                   global_y, prim_array_start[0][0], prim_array_start[0][1],
                   conj_array_start[0][0], conj_array_start[0][1]);
        }
        #endif
        
        FFTW_EXECUTE_DFT(plan_2d, prim_array_start, prim_array_start);
        
        // DEBUG: Memory guard after primary FFT
        #if DEBUG_PRINTS
        if (rank < 4 && global_y <= 2) {
            // Verify buffer is still valid after FFT
            size_t prim_array_size = (size_t)N * N * sizeof(fftw_complex_t);
            volatile char checksum = 0;
            for (size_t i = 0; i < prim_array_size && i < 1024*1024; i += 4096) {
                checksum ^= ((char*)prim_array_start)[i];
            }
            fprintf(stderr, "[Rank %d] After primary FFT (array %d, Y=%d): buffer OK (checksum=%d)\n",
                   rank, a, global_y, (int)checksum);
            fflush(stderr);
        }
        #endif
        
        // Conjugate slice (if not self-conjugate)
        if (y_mirror != global_y) {
            fftw_complex_t *conj_array_start = &CONJ_SLICE(a, 0, 0);
            
            // DEBUG: Check values right before FFT
            #if DEBUG_PRINTS && !SKIP_VERIFICATION
            if (rank == 0 && global_y <= 2 && a == 0) {
                printf("[DEBUG Y=%d] Right before conjugate FFT - conj[0,0]=(%.3e,%.3e) prim[0,0]=(%.3e,%.3e)\n",
                       global_y, conj_array_start[0][0], conj_array_start[0][1],
                       prim_array_start[0][0], prim_array_start[0][1]);
            }
            #endif
            
            FFTW_EXECUTE_DFT(plan_2d, conj_array_start, conj_array_start);
            
            // DEBUG: Memory guard after conjugate FFT
            #if DEBUG_PRINTS
            if (rank < 4 && global_y <= 2) {
                // Verify buffer is still valid after FFT
                size_t conj_array_size = (size_t)N * N * sizeof(fftw_complex_t);
                volatile char checksum = 0;
                for (size_t i = 0; i < conj_array_size && i < 1024*1024; i += 4096) {
                    checksum ^= ((char*)conj_array_start)[i];
                }
                fprintf(stderr, "[Rank %d] After conjugate FFT (array %d, Y=%d): buffer OK (checksum=%d)\n",
                       rank, a, global_y, (int)checksum);
                fflush(stderr);
            }
            #endif
            
            // DEBUG: Check values right after FFT
            #if DEBUG_PRINTS && !SKIP_VERIFICATION
            if (rank == 0 && global_y <= 2 && a == 0) {
                printf("[DEBUG Y=%d] Right after conjugate FFT - conj[0,0]=(%.3e,%.3e) prim[0,0]=(%.3e,%.3e)\n",
                       global_y, conj_array_start[0][0], conj_array_start[0][1],
                       prim_array_start[0][0], prim_array_start[0][1]);
            }
            #endif
        }
    }
    
    // ========== DEBUG: Check imaginary parts AFTER 2D FFT (Real space - should be ~0) ==========
    #if DEBUG_PRINTS && !SKIP_VERIFICATION
    if (rank == 0 && global_y <= 2) {
        real_t debug_max_real_2d[4] = {0, 0, 0, 0};
        real_t debug_max_imag_2d[4] = {0, 0, 0, 0};
        for (int a = 0; a < narray; a++) {
            for (int x = 0; x < N; x++) {
                for (int z = 0; z < N; z++) {
                    double re = fabs_t(PRIM_SLICE(a, x, z)[0]);
                    double im = fabs_t(PRIM_SLICE(a, x, z)[1]);
                    debug_max_real_2d[a] = fmax_t(debug_max_real_2d[a], re);
                    debug_max_imag_2d[a] = fmax_t(debug_max_imag_2d[a], im);
                }
            }
        }
        printf("[DEBUG Y=%d] After 2D FFT (Real space) - PRIMARY slice - Max real/imag per array: ", global_y);
        for (int a = 0; a < narray; a++) {
            double ratio = debug_max_imag_2d[a] / (debug_max_real_2d[a] + 1e-16);
            printf("A%d: re=%.3e im=%.3e (ratio=%.3e) ", 
                   a, debug_max_real_2d[a], debug_max_imag_2d[a], ratio);
        }
        printf("(NOTE: Imag should be ~0 in real space!)\n");
        
        // Also check conjugate slice for conjugate pairs
        if (y_mirror != global_y) {
            real_t debug_max_real_conj[4] = {0, 0, 0, 0};
            real_t debug_max_imag_conj[4] = {0, 0, 0, 0};
            int num_different = 0;
            int num_identical = 0;
            for (int a = 0; a < narray; a++) {
                for (int x = 0; x < N; x++) {
                    for (int z = 0; z < N; z++) {
                        double re = fabs_t(CONJ_SLICE(a, x, z)[0]);
                        double im = fabs_t(CONJ_SLICE(a, x, z)[1]);
                        debug_max_real_conj[a] = fmax_t(debug_max_real_conj[a], re);
                        debug_max_imag_conj[a] = fmax_t(debug_max_imag_conj[a], im);
                        
                        // Check if this element is identical to primary
                        double prim_re = PRIM_SLICE(a, x, z)[0];
                        double prim_im = PRIM_SLICE(a, x, z)[1];
                        double conj_re = CONJ_SLICE(a, x, z)[0];
                        double conj_im = CONJ_SLICE(a, x, z)[1];
                        double diff_re = fabs_t(prim_re - conj_re);
                        double diff_im = fabs_t(prim_im - conj_im);
                        if (diff_re < 1e-10 && diff_im < 1e-10) {
                            num_identical++;
                        } else {
                            num_different++;
                        }
                    }
                }
            }
            printf("[DEBUG Y=%d] After 2D FFT (Real space) - CONJUGATE slice - Max real/imag per array: ", global_y);
            for (int a = 0; a < narray; a++) {
                double ratio = debug_max_imag_conj[a] / (debug_max_real_conj[a] + 1e-16);
                printf("A%d: re=%.3e im=%.3e (ratio=%.3e) ", 
                       a, debug_max_real_conj[a], debug_max_imag_conj[a], ratio);
            }
            printf("(NOTE: Imag should be ~0 in real space!)\n");
            printf("[DEBUG Y=%d] Array comparison: %d identical, %d different (total=%d)\n", 
                   global_y, num_identical, num_different, num_identical + num_different);
            
            // Show a few sample differences
            if (num_different > 0) {
                printf("[DEBUG Y=%d] Sample differences (first 3 non-identical): ", global_y);
                int shown = 0;
                for (int a = 0; a < narray && shown < 3; a++) {
                    for (int x = 0; x < N && shown < 3; x++) {
                        for (int z = 0; z < N && shown < 3; z++) {
                            double prim_re = PRIM_SLICE(a, x, z)[0];
                            double prim_im = PRIM_SLICE(a, x, z)[1];
                            double conj_re = CONJ_SLICE(a, x, z)[0];
                            double conj_im = CONJ_SLICE(a, x, z)[1];
                            double diff_re = fabs_t(prim_re - conj_re);
                            double diff_im = fabs_t(prim_im - conj_im);
                            if (diff_re > 1e-10 || diff_im > 1e-10) {
                                printf("A%d[%d,%d]: prim=(%.3e,%.3e) conj=(%.3e,%.3e) diff=(%.3e,%.3e) ", 
                                       a, x, z, prim_re, prim_im, conj_re, conj_im, diff_re, diff_im);
                                shown++;
                            }
                        }
                    }
                }
                printf("\n");
            }
        }
    }
    #endif
    
    // ========== Apply any remaining nskip at end of function ==========
    // This matches zeldovich.cpp: Pk.v2rng[y].advance(2 * nskip) at end
    // Ensures that exactly MAX_PPD * MAX_PPD complex numbers (2 * MAX_PPD * MAX_PPD real numbers)
    // were consumed/skipped for this y-row
    #if !PARALLELIZE_XZ_WITHIN_SLICE
    if (nskip > 0) {
        if (ps_handle != NULL && params_handle != NULL) {
            int64_t rng_index = global_y;
            #if DEBUG_RNG_SKIP
            if (global_y <= MAX_DEBUG_COORD || global_y == Nhalf - 1) {
                fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d: BEFORE advance_rng - nskip=%lld (0x%llx)\n",
                        N, global_y, (long long)nskip, (unsigned long long)nskip);
                fflush(stderr);
            }
            #endif
            zeldovich_ps_advance_rng(ps_handle, params_handle, rng_index, nskip);
            #if VERIFY_RNG_CALLS
            total_rng_skips += nskip;
            #endif
            #if DEBUG_RNG_SKIP
            if (global_y <= MAX_DEBUG_COORD || global_y == Nhalf - 1) {
                fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d: END OF FUNCTION - Applied remaining nskip=%lld (0x%llx)\n",
                        N, global_y, (long long)nskip, (unsigned long long)nskip);
                fflush(stderr);
            }
            #endif
        } else if (params_handle == NULL) {
            // Local PCG: advance by 2 * nskip (each complex number uses 2 random numbers)
            // Cast to uint64_t explicitly to avoid overflow in multiplication
            uint64_t advance_amount = (uint64_t)2 * (uint64_t)nskip;
            advance_pcg_global(global_y, advance_amount);
            #if VERIFY_RNG_CALLS
            total_rng_skips += nskip;
            #endif
            #if DEBUG_RNG_SKIP
            if (global_y <= MAX_DEBUG_COORD || global_y == Nhalf - 1) {
                fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d: END OF FUNCTION - Applied remaining nskip=%lld (PCG: 2*%lld)\n",
                        N, global_y, (long long)nskip, (long long)nskip);
                fflush(stderr);
            }
            #endif
        }
        nskip = 0;  // Reset after applying
    }
    #endif
    
    // ========== RNG VERIFICATION: Verify total RNG calls match expected ==========
    // Similar to zeldovich.cpp assertion: assert(Pk.v2rng[y] - checkpoint == 2 * MAX_PPD * MAX_PPD)
    // 
    // Why 2 * MAX_PPD * MAX_PPD?
    // - zeldovich.cpp processes ppd * ppd points per Y-slice (where ppd <= MAX_PPD)
    // - For each point: either call cgauss() (uses 2 random numbers) or skip (advance by 2 later)
    // - When ppd < MAX_PPD, missing (MAX_PPD-ppd) * MAX_PPD points are skipped
    // - Total = 2 * (actual points processed) + 2 * (missing points skipped)
    //         = 2 * ppd * ppd + 2 * (MAX_PPD * MAX_PPD - ppd * ppd)
    //         = 2 * MAX_PPD * MAX_PPD
    // This ensures RNG consistency: same Y-slice always consumes same number of random numbers
    // regardless of actual grid size, so overlapping regions match across different N values.
    #if VERIFY_RNG_CALLS
    {
        // Expected: 2 * MAX_PPD * MAX_PPD random numbers total (always, regardless of N)
        // - Each cgauss() call consumes 2 random numbers (tracked in total_rng_calls)
        // - Each skip advances by 2 * nskip random numbers (tracked in total_rng_skips, already accounts for the 2x)
        //   This includes: missing coordinates (x/z boundaries) and D=0 skips (DC/Nyquist modes)
        // Note: Cast to int64_t to avoid potential overflow issues
        int64_t max_ppd_val = (int64_t)MAX_PPD;
        int64_t expected_total = 2LL * max_ppd_val * max_ppd_val;
        int64_t actual_total = 2LL * total_rng_calls + 2LL * total_rng_skips;
        
        if (actual_total != expected_total) {
            fprintf(stderr, 
                    "[RNG-VERIFY ERROR] Rank %d, Y=%d: Expected %ld random numbers (MAX_PPD=%ld), got %ld "
                    "(calls=%ld * 2 = %ld, skips=%ld * 2 = %ld)\n",
                    rank, global_y, (long)expected_total, (long)max_ppd_val, (long)actual_total,
                    (long)total_rng_calls, (long)(2LL * total_rng_calls),
                    (long)total_rng_skips, (long)(2LL * total_rng_skips));
            fflush(stderr);
            // Don't abort in production, but warn
            #ifdef DEBUG
            assert(actual_total == expected_total);
            #endif
        } else if (rank == 0 && global_y <= 2) {
            // Print verification for first few slices in debug mode
            fprintf(stderr,
                    "[RNG-VERIFY OK] Rank %d, Y=%d: Total=%ld (calls=%ld, skips=%ld)\n",
                    rank, global_y, (long)actual_total,
                    (long)total_rng_calls, (long)total_rng_skips);
            fflush(stderr);
        }
    }
    #endif
    
    // Clean up local macros
    #undef PRIM_SLICE
    #undef CONJ_SLICE
}

