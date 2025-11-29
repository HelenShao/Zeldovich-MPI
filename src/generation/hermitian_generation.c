// ====================================================================================
// HERMITIAN 3D MATRIX MPI - HERMITIAN GENERATION MODULE
// ====================================================================================

#include "hermitian_generation.h"
#include "../utils/verification.h"
#include "../utils/zeldovich_wrapper.h"  // v15.2: For zeldovich-PLT integration
#include "../config.h"  // For DEBUG_PRINTS, SKIP_VERIFICATION, MAX_PPD
#include "../precision.h"  // For real_t, fabs_t, fmax_t
#include <stdio.h>
#include <stdlib.h>  // For malloc/free
#include <stdint.h>  // For int64_t
#include <math.h>    // For sqrt()
#include <omp.h>

// ====================================================================================
// FUNCTION IMPLEMENTATIONS
// ====================================================================================

void generate_hermitian_slice_pair_local(
    int N,
    int global_y,
    int y_mirror,
    fftw_complex_t *primary_slices,   // Now flat array (all narray arrays for primary slice)
    fftw_complex_t *conjugate_slices, // Now flat array (all narray arrays for conjugate slice)
    int narray,                      // NEW: number of arrays per slice
    fftw_plan_t plan_2d,
    int rank,
    const power_spectrum_params_t *ps_params,  // Legacy power spectrum parameters (NULL = use uniform RNG)
    PowerSpectrumHandle ps_handle,   // v15.2: zeldovich-PLT PowerSpectrum handle (NULL = use legacy or uniform RNG)
    ParametersHandle params_handle)  // v15.2: zeldovich-PLT Parameters handle (needed for fundamental wavenumber)
{
    // Debug: Log entry for segmentation fault investigation (enabled via DEBUG_PRINTS)
    #if DEBUG_PRINTS
    fprintf(stderr,
            "[Rank %d] ENTER generate_hermitian_slice_pair_local: "
            "Y_primary=%d, Y_mirror=%d, ps_handle=%p, params_handle=%p\n",
            rank, global_y, y_mirror, (void*)ps_handle, (void*)params_handle);
    fflush(stderr);
    #endif

    // Create local aliases for macro compatibility
    // primary_slices points to slice 0, conjugate_slices points to slice 1
    // We need to access them as if they're part of a larger buffer
    // For macro: Y_SLICE(0, ...) uses primary_slices, Y_SLICE(1, ...) uses conjugate_slices
    
    // Helper function to access primary slice (slice_idx=0)
    #define PRIM_SLICE(array_idx, x, z) \
        primary_slices[(int64_t)(x) + (N) * ((z) + (N) * (array_idx))]
    
    // Helper function to access conjugate slice (slice_idx=1)  
    #define CONJ_SLICE(array_idx, x, z) \
        conjugate_slices[(int64_t)(x) + (N) * ((z) + (N) * (array_idx))]
    
    // RNG consistency: Skip random numbers for missing grid points when N < MAX_PPD
    // This maintains consistency with what a full MAX_PPD × MAX_PPD grid would generate
    // Similar to zeldovich.cpp skip logic
    int64_t nskip = 0;
    if (N < MAX_PPD) {
        // Calculate skip for missing grid points
        // When crossing Nyquist boundary at z = N/2 + 1: skip (MAX_PPD - N) * MAX_PPD points (missing z-rows)
        // When crossing Nyquist boundary at x = N/2 + 1: skip (MAX_PPD - N) points (missing x-values in current z-row)
        // For parallel execution, calculate upfront; for sequential, track incrementally
        #if PARALLELIZE_XZ_WITHIN_SLICE
        // Parallel case: calculate total skip upfront (deterministic)
        // We cross z boundary once: (MAX_PPD - N) * MAX_PPD
        // We cross x boundary N times (once per z-row): N * (MAX_PPD - N)
        // Total: (MAX_PPD - N) * MAX_PPD + N * (MAX_PPD - N) = (MAX_PPD - N) * (MAX_PPD + N)
        // Actually, let's be more precise:
        // - When z == N/2 + 1: skip (MAX_PPD - N) * MAX_PPD (all missing z-rows)
        // - For each z from 0 to N-1, when x == N/2 + 1: skip (MAX_PPD - N)
        // Total: (MAX_PPD - N) * MAX_PPD + N * (MAX_PPD - N) = (MAX_PPD - N) * (MAX_PPD + N)
        nskip = (MAX_PPD - N) * MAX_PPD + N * (MAX_PPD - N);
        #else
        // Sequential case: track incrementally during iteration (like zeldovich.cpp)
        nskip = 0;  // Will be accumulated during loops
        #endif
    }
    
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
            // RNG consistency: Skip random numbers when crossing Nyquist boundary
            // When z == N/2 + 1, we enter the negative k region
            // If N < MAX_PPD, we need to skip missing z-rows (z = N to MAX_PPD-1)
            #if !PARALLELIZE_XZ_WITHIN_SLICE
            if (z == N/2 + 1 && N < MAX_PPD) {
                nskip += (MAX_PPD - N) * MAX_PPD;
            }
            #endif
            for (int x = 0; x < N; x++) {
                // RNG consistency: Skip random numbers when crossing Nyquist boundary
                // When x == N/2 + 1, we enter the negative k region
                // If N < MAX_PPD, we need to skip missing x-values (x = N to MAX_PPD-1) in current z-row
                #if !PARALLELIZE_XZ_WITHIN_SLICE
                if (x == N/2 + 1 && N < MAX_PPD) {
                    nskip += MAX_PPD - N;
                }
                #endif
                int x_mirror = (x == 0) ? 0 : N - x;
                int z_mirror = (z == 0) ? 0 : N - z;
                
                // ========== STEP 1: Calculate k-vector components ==========
                int kx = (x > N/2) ? x - N : x;
                int ky = (global_y > N/2) ? global_y - N : global_y;
                int kz = (z > N/2) ? z - N : z;
                double k2 = (double)(kx*kx + ky*ky + kz*kz);
                
                // ========== STEP 2: Generate D using RNG or cgauss() ==========
                fftw_complex D;
                if (k2 == 0.0) {
                    // DC mode: set to zero
                    D[0] = D[1] = 0.0;
                } else if (ps_handle != NULL && params_handle != NULL) {
                    // v15.2: Use zeldovich-PLT power spectrum-weighted Gaussian
                    // If ps_handle is available, we ALWAYS use power spectrum mode (cgauss)
                    
                    // TEMPORARY: Set Nyquist axis (Y = N/2) to zero
                    // TODO: Understand the "shifted by one location" issue in zeldovich-PLT
                    // (see zeldovich.cpp line 646-648). For now, matching zeldovich-PLT behavior
                    // which sets the Nyquist plane to zero after data reorganization.
                    // This should be reviewed and potentially changed after understanding the
                    // data shift operation in LoadBlock() and its relationship to Nyquist handling.
                    if (global_y == N/2) {
                        D[0] = D[1] = 0.0;
                    } else {
                        // Convert k indices to physical wavenumber: k_phys = k_index * fundamental
                        double fundamental = zeldovich_params_get_fundamental(params_handle);
                        double k2_phys = k2 * fundamental * fundamental;
                        double kmag = sqrt(k2_phys);
                        
                        // zeldovich_ps_cgauss returns double precision, convert to real_t
                        // zeldovich-PLT's v2rng array is sized to ppd/2, so valid indices are 0 to (N/2 - 1)
                        // Since we've already handled global_y == N/2 above, global_y is now < N/2
                        int64_t rng_index = global_y;
                        double D_real, D_imag;
                        #if PARALLELIZE_XZ_WITHIN_SLICE
                        // Lock protects generator access
                        #else
                        // Sequential access: no locks needed
                        #endif
                        zeldovich_ps_cgauss(ps_handle, kmag, rng_index, &D_real, &D_imag);
                        D[0] = (real_t)D_real;
                        D[1] = (real_t)D_imag;
                    }
                } else if (ps_params != NULL) {
                    // Legacy: Use standalone power spectrum-weighted Gaussian (cgauss)
                    // Convert k indices to physical wavenumber: k_phys = k_index * fundamental
                    // For now, fundamental = 1.0 (can be added to ps_params later)
                    double fundamental = 1.0;  // TODO: Add to power_spectrum_params_t
                    double k2_phys = k2 * fundamental * fundamental;
                    double kmag = sqrt(k2_phys);
                    
                    #if PARALLELIZE_XZ_WITHIN_SLICE
                    // Lock protects generator access
                    #else
                    // Sequential access: no locks needed
                    #endif
                    cgauss(ps_params, kmag, global_y, &D);
                } else {
                    // Fallback: use uniform random numbers (white noise mode, no power spectrum)
                    // This should only happen when ps_handle is NULL (no parameter file provided)
                    #if PARALLELIZE_XZ_WITHIN_SLICE
                    // Lock protects generator access: each call advances generator by +1
                    // Order is non-deterministic across threads, but generator advances correctly
                    #else
                    // Sequential access: no locks needed, generator advances deterministically
                    #endif
                    // Note: If ps_handle is available, we should have used cgauss() above
                    // This branch is only for backward compatibility (no power spectrum mode)
                    double D_re = random_real_pcg_global(global_y);
                    double D_im = random_real_pcg_global(global_y);
                    D[0] = D_re;
                    D[1] = D_im;
                }
                
                // ========== STEP 3: Compute F, G, H deterministically from D ==========
                fftw_complex F, G, H;
                if (k2 == 0.0) {
                    // DC mode: set everything to zero (D already set to 0 above)
                    F[0] = F[1] = G[0] = G[1] = H[0] = H[1] = 0.0;
                } else {
                    double ik2 = 1.0 / k2;
                    
                    // F = i × kx/k² × D = i × kx × ik2 × (D_re + i×D_im)
                    //   = i × kx × ik2 × D_re - kx × ik2 × D_im
                    //   = -kx × ik2 × D_im + i × kx × ik2 × D_re
                    F[0] = -kx * ik2 * D[1];  // Real part
                    F[1] =  kx * ik2 * D[0];  // Imaginary part
                    
                    G[0] = -ky * ik2 * D[1];
                    G[1] =  ky * ik2 * D[0];
                    
                    H[0] = -kz * ik2 * D[1];
                    H[1] =  kz * ik2 * D[0];
                }
                
                // ========== STEP 4: Store in arrays ==========
                // Array 0: D + i*F (density + X-displacement)
                PRIM_SLICE(0, x, z)[0] = D[0];  // Real = D_re
                PRIM_SLICE(0, x, z)[1] = F[0];  // Imag = F_re (X-displacement)
                
                // Array 1: G + i*H (Y-displacement + Z-displacement)
                PRIM_SLICE(1, x, z)[0] = G[0];  // Real = G_re (Y-displacement)
                PRIM_SLICE(1, x, z)[1] = H[0];  // Imag = H_re (Z-displacement)
                
                if (narray >= 4) {
                    // Array 2: 0 + i*F*f (X-velocity)
                    double f = 1.0;  // PLT growth rate (placeholder for now)
                    PRIM_SLICE(2, x, z)[0] = 0.0;
                    PRIM_SLICE(2, x, z)[1] = F[0] * f;
                    
                    // Array 3: G*f + i*H*f (Y-velocity + Z-velocity)
                    PRIM_SLICE(3, x, z)[0] = G[0] * f;
                    PRIM_SLICE(3, x, z)[1] = H[0] * f;
                }
                
                // ========== STEP 5: Store Hermitian conjugates ==========
                // For ALL arrays independently
                for (int a = 0; a < narray; a++) {
                    double re = PRIM_SLICE(a, x, z)[0];
                    double im = PRIM_SLICE(a, x, z)[1];
                    
                    // Store conjugate at mirror location
                    CONJ_SLICE(a, x_mirror, z_mirror)[0] = re;
                    CONJ_SLICE(a, x_mirror, z_mirror)[1] = -im;
                }
            }
        }
        
        // Advance RNG by 2 * nskip to account for missing grid points
        // Each complex number needs 2 random numbers (real and imaginary)
        if (nskip > 0 && N < MAX_PPD) {
            advance_pcg_global(global_y, 2 * nskip);
        }
    } else {
        // ========== SELF-CONJUGATE: Y=0 or Y=N/2 ==========
        // Same approach as conjugate pair: generate D, compute F,G,H
        // Reset nskip for self-conjugate case
        if (N < MAX_PPD) {
            #if PARALLELIZE_XZ_WITHIN_SLICE
            // Parallel case: calculate total skip upfront
            // For self-conjugate, we iterate z = 0 to N/2, x varies by z
            // We don't cross z boundary, but we do cross x boundary
            // For z = 0: x goes 0 to N/2+1, crosses at x = N/2+1: skip (MAX_PPD - N)
            // For z > 0: x goes 0 to N-1, crosses at x = N/2+1: skip (MAX_PPD - N) per z
            // Total: (MAX_PPD - N) + (N/2) * (MAX_PPD - N) = (MAX_PPD - N) * (1 + N/2)
            nskip = (MAX_PPD - N) * (1 + N/2);
            #else
            nskip = 0;  // Will be accumulated during loops
            #endif
        } else {
            nskip = 0;
        }
        #if USE_ZELDOVICH_METHOD
        // Zeldovich method: Fill half the plane, mirror the rest
        // PARALLELIZATION: Conditional based on PARALLELIZE_XZ_WITHIN_SLICE flag
        #if PARALLELIZE_XZ_WITHIN_SLICE
        // Parallel z-loop: requires locks for thread-safe RNG access
        #pragma omp parallel for
        #else
        // Sequential z-loop: no locks needed
        #endif
        for (int z = 0; z <= N/2; z++) {
            // RNG consistency: Skip random numbers when crossing Nyquist boundary
            // When z == N/2 + 1, we enter the negative k region
            // If N < MAX_PPD, we need to skip missing z-rows (z = N to MAX_PPD-1)
            // Note: For self-conjugate, we only iterate z = 0 to N/2, so we check at z = N/2 + 1
            // But since we stop at z = N/2, we need to account for missing points differently
            // Actually, for self-conjugate, we don't cross the Nyquist boundary in z (we stop at N/2)
            // But we still need to account for missing x-values when crossing x = N/2 + 1
            int x_max = (z == 0 ? N/2 + 1 : N);
            for (int x = 0; x < x_max; x++) {
                // RNG consistency: Skip random numbers when crossing Nyquist boundary
                // When x == N/2 + 1, we enter the negative k region
                // If N < MAX_PPD, we need to skip missing x-values (x = N to MAX_PPD-1) in current z-row
                #if !PARALLELIZE_XZ_WITHIN_SLICE
                if (x == N/2 + 1 && N < MAX_PPD) {
                    nskip += MAX_PPD - N;
                }
                #endif
                int x_mirror = (x == 0) ? 0 : N - x;
                int z_mirror = (z == 0) ? 0 : N - z;
                
                // Calculate k-vector components first (needed for both cgauss and uniform RNG)
                int kx = (x > N/2) ? x - N : x;
                int ky = (global_y > N/2) ? global_y - N : global_y;
                int kz = (z > N/2) ? z - N : z;
                double k2 = (double)(kx*kx + ky*ky + kz*kz);
                
                // Generate D using RNG or cgauss()
                // For self-conjugate slices, we still need to check for power spectrum mode
                #if PARALLELIZE_XZ_WITHIN_SLICE
                // Lock protects generator access
                #else
                // Sequential access: no locks needed
                #endif
                fftw_complex D;
                if (k2 == 0.0) {
                    // DC mode: set to zero
                    D[0] = D[1] = 0.0;
                } else if (ps_handle != NULL && params_handle != NULL) {
                    // v15.2: Use zeldovich-PLT power spectrum-weighted Gaussian
                    // If ps_handle is available, we ALWAYS use power spectrum mode (cgauss)
                    
                    // TEMPORARY: Set Nyquist axis (Y = N/2) to zero
                    // TODO: Understand the "shifted by one location" issue in zeldovich-PLT
                    // (see zeldovich.cpp line 646-648). For now, matching zeldovich-PLT behavior
                    // which sets the Nyquist plane to zero after data reorganization.
                    // This should be reviewed and potentially changed after understanding the
                    // data shift operation in LoadBlock() and its relationship to Nyquist handling.
                    if (global_y == N/2) {
                        D[0] = D[1] = 0.0;
                    } else {
                        // Convert k indices to physical wavenumber: k_phys = k_index * fundamental
                        double fundamental = zeldovich_params_get_fundamental(params_handle);
                        double k2_phys = k2 * fundamental * fundamental;
                        double kmag = sqrt(k2_phys);
                        
                        // zeldovich_ps_cgauss returns double precision, convert to real_t
                        // zeldovich-PLT's v2rng array is sized to ppd/2, so valid indices are 0 to (N/2 - 1)
                        // Since we've already handled global_y == N/2 above, global_y is now < N/2
                        int64_t rng_index = global_y;
                        double D_real, D_imag;
                        zeldovich_ps_cgauss(ps_handle, kmag, rng_index, &D_real, &D_imag);
                        D[0] = (real_t)D_real;
                        D[1] = (real_t)D_imag;
                    }
                } else if (ps_params != NULL) {
                    // Legacy: Use standalone power spectrum-weighted Gaussian (cgauss)
                    double fundamental = 1.0;  // TODO: Add to power_spectrum_params_t
                    double k2_phys = k2 * fundamental * fundamental;
                    double kmag = sqrt(k2_phys);
                    
                    cgauss(ps_params, kmag, global_y, &D);
                } else {
                    // Fallback: use uniform random numbers (white noise mode, no power spectrum)
                    // This should only happen when ps_handle is NULL (no parameter file provided)
                    // Note: If ps_handle is available, we should have used cgauss() above
                    double D_re = random_real_pcg_global(global_y);
                    double D_im = random_real_pcg_global(global_y);
                    D[0] = D_re;
                    D[1] = D_im;
                }
                
                // Compute F, G, H from D
                fftw_complex F, G, H;
                if (k2 == 0.0) {
                    F[0] = F[1] = G[0] = G[1] = H[0] = H[1] = 0.0;
                } else {
                    double ik2 = 1.0 / k2;
                    F[0] = -kx * ik2 * D[1];
                    F[1] =  kx * ik2 * D[0];
                    G[0] = -ky * ik2 * D[1];
                    G[1] =  ky * ik2 * D[0];
                    H[0] = -kz * ik2 * D[1];
                    H[1] =  kz * ik2 * D[0];
                }
                
                // Store in arrays (all arrays in primary_slices, self-conjugate uses same buffer)
                PRIM_SLICE(0, x, z)[0] = D[0];
                PRIM_SLICE(0, x, z)[1] = F[0];
                PRIM_SLICE(1, x, z)[0] = G[0];
                PRIM_SLICE(1, x, z)[1] = H[0];
                
                if (narray >= 4) {
                    double f = 1.0;  // PLT growth rate (placeholder)
                    PRIM_SLICE(2, x, z)[0] = 0.0;
                    PRIM_SLICE(2, x, z)[1] = F[0] * f;
                    PRIM_SLICE(3, x, z)[0] = G[0] * f;
                    PRIM_SLICE(3, x, z)[1] = H[0] * f;
                }
                
                // Mirror for self-conjugate
                if (x != x_mirror || z != z_mirror) {
                    for (int a = 0; a < narray; a++) {
                        double re = PRIM_SLICE(a, x, z)[0];
                        double im = PRIM_SLICE(a, x, z)[1];
                        PRIM_SLICE(a, x_mirror, z_mirror)[0] = re;
                        PRIM_SLICE(a, x_mirror, z_mirror)[1] = -im;
                    }
                }
            }
        }
        
        // Advance RNG by 2 * nskip to account for missing grid points
        // Each complex number needs 2 random numbers (real and imaginary)
        if (nskip > 0 && N < MAX_PPD) {
            advance_pcg_global(global_y, 2 * nskip);
        }
        
        // Enforce self-conjugate constraints (imaginary parts = 0 at special points)
        if (global_y == 0) {
            for (int a = 0; a < narray; a++) {
                PRIM_SLICE(a, 0, 0)[0] = 0.0;
                PRIM_SLICE(a, 0, 0)[1] = 0.0;
            }
        }
        // Set imaginary parts to 0 at special points
        for (int a = 0; a < narray; a++) {
            PRIM_SLICE(a, 0, 0)[1] = 0.0;
            PRIM_SLICE(a, N/2, 0)[1] = 0.0;
            PRIM_SLICE(a, 0, N/2)[1] = 0.0;
            PRIM_SLICE(a, N/2, N/2)[1] = 0.0;
        }
        #else
        // Fallback method (not recommended, kept for compatibility)
        // PARALLELIZATION: Conditional based on PARALLELIZE_XZ_WITHIN_SLICE flag
        #if PARALLELIZE_XZ_WITHIN_SLICE
        // Parallel (x,z) loops: requires locks for thread-safe RNG access
        #pragma omp parallel for collapse(2)
        #else
        // Sequential (x,z) loops: no locks needed
        #endif
        for (int x = 0; x < N; x++) {
            for (int z = 0; z < N; z++) {
                int x_mirror = (x == 0) ? 0 : N - x;
                int z_mirror = (z == 0) ? 0 : N - z;
                
                // Calculate k-vector components
                int kx = (x > N/2) ? x - N : x;
                int ky = (global_y > N/2) ? global_y - N : global_y;
                int kz = (z > N/2) ? z - N : z;
                double k2 = (double)(kx*kx + ky*ky + kz*kz);
                
                // Generate D using RNG or cgauss()
                fftw_complex D;
                if (k2 == 0.0) {
                    // DC mode: set to zero
                    D[0] = D[1] = 0.0;
                } else if (ps_handle != NULL && params_handle != NULL) {
                    // v15.2: Use zeldovich-PLT power spectrum-weighted Gaussian
                    // If ps_handle is available, we ALWAYS use power spectrum mode (cgauss)
                    // Convert k indices to physical wavenumber: k_phys = k_index * fundamental
                    // TEMPORARY: Set Nyquist axis (Y = N/2) to zero
                    // TODO: Understand the "shifted by one location" issue in zeldovich-PLT
                    // (see zeldovich.cpp line 646-648). For now, matching zeldovich-PLT behavior
                    // which sets the Nyquist plane to zero after data reorganization.
                    // This should be reviewed and potentially changed after understanding the
                    // data shift operation in LoadBlock() and its relationship to Nyquist handling.
                    if (global_y == N/2) {
                        D[0] = D[1] = 0.0;
                    } else {
                        double fundamental = zeldovich_params_get_fundamental(params_handle);
                        double k2_phys = k2 * fundamental * fundamental;
                        double kmag = sqrt(k2_phys);
                        
                        // zeldovich_ps_cgauss returns double precision, convert to real_t
                        // zeldovich-PLT's v2rng array is sized to ppd/2, so valid indices are 0 to (N/2 - 1)
                        // Since we've already handled global_y == N/2 above, global_y is now < N/2
                        int64_t rng_index = global_y;
                        double D_real, D_imag;
                        #if PARALLELIZE_XZ_WITHIN_SLICE
                        // Lock protects generator access
                        #else
                        // Sequential access: no locks needed
                        #endif
                        zeldovich_ps_cgauss(ps_handle, kmag, rng_index, &D_real, &D_imag);
                        D[0] = (real_t)D_real;
                        D[1] = (real_t)D_imag;
                    }
                } else if (ps_params != NULL) {
                    // Legacy: Use standalone power spectrum-weighted Gaussian (cgauss)
                    // Convert k indices to physical wavenumber: k_phys = k_index * fundamental
                    // For now, fundamental = 1.0 (can be added to ps_params later)
                    double fundamental = 1.0;  // TODO: Add to power_spectrum_params_t
                    double k2_phys = k2 * fundamental * fundamental;
                    double kmag = sqrt(k2_phys);
                    
                    #if PARALLELIZE_XZ_WITHIN_SLICE
                    // Lock protects generator access
                    #else
                    // Sequential access: no locks needed
                    #endif
                    cgauss(ps_params, kmag, global_y, &D);
                } else {
                    // Fallback: use uniform random numbers (white noise mode, no power spectrum)
                    // This should only happen when ps_handle is NULL (no parameter file provided)
                    #if PARALLELIZE_XZ_WITHIN_SLICE
                    // Lock protects generator access: each call advances generator by +1
                    // Order is non-deterministic across threads, but generator advances correctly
                    #else
                    // Sequential access: no locks needed, generator advances deterministically
                    #endif
                    // Note: If ps_handle is available, we should have used cgauss() above
                    // This branch is only for backward compatibility (no power spectrum mode)
                    double D_re = random_real_pcg_global(global_y);
                    double D_im = random_real_pcg_global(global_y);
                    D[0] = D_re;
                    D[1] = D_im;
                }
                
                fftw_complex F, G, H;
                if (k2 == 0.0) {
                    F[0] = F[1] = G[0] = G[1] = H[0] = H[1] = 0.0;
                } else {
                    double ik2 = 1.0 / k2;
                    F[0] = -kx * ik2 * D[1];
                    F[1] =  kx * ik2 * D[0];
                    G[0] = -ky * ik2 * D[1];
                    G[1] =  ky * ik2 * D[0];
                    H[0] = -kz * ik2 * D[1];
                    H[1] =  kz * ik2 * D[0];
                }
                
                // Handle special points
                if (global_y == 0 && x == 0 && z == 0) {
                    // DC mode: all zero
                    for (int a = 0; a < narray; a++) {
                        PRIM_SLICE(a, x, z)[0] = 0.0;
                        PRIM_SLICE(a, x, z)[1] = 0.0;
                    }
                } else if ((x == 0 || x == N/2) && (z == 0 || z == N/2)) {
                    // Self-symmetric points: imaginary = 0
                    PRIM_SLICE(0, x, z)[0] = D[0];
                    PRIM_SLICE(0, x, z)[1] = 0.0;
                    PRIM_SLICE(1, x, z)[0] = G[0];
                    PRIM_SLICE(1, x, z)[1] = 0.0;
                    if (narray >= 4) {
                        double f = 1.0;
                        PRIM_SLICE(2, x, z)[0] = 0.0;
                        PRIM_SLICE(2, x, z)[1] = 0.0;
                        PRIM_SLICE(3, x, z)[0] = G[0] * f;
                        PRIM_SLICE(3, x, z)[1] = 0.0;
                    }
                } else {
                    // Normal points
                    PRIM_SLICE(0, x, z)[0] = D[0];
                    PRIM_SLICE(0, x, z)[1] = F[0];
                    PRIM_SLICE(1, x, z)[0] = G[0];
                    PRIM_SLICE(1, x, z)[1] = H[0];
                    if (narray >= 4) {
                        double f = 1.0;
                        PRIM_SLICE(2, x, z)[0] = 0.0;
                        PRIM_SLICE(2, x, z)[1] = F[0] * f;
                        PRIM_SLICE(3, x, z)[0] = G[0] * f;
                        PRIM_SLICE(3, x, z)[1] = H[0] * f;
                    }
                }
                
                // Mirror
                if (x != x_mirror || z != z_mirror) {
                    for (int a = 0; a < narray; a++) {
                        double re = PRIM_SLICE(a, x, z)[0];
                        double im = PRIM_SLICE(a, x, z)[1];
                        PRIM_SLICE(a, x_mirror, z_mirror)[0] = re;
                        PRIM_SLICE(a, x_mirror, z_mirror)[1] = -im;
                    }
                }
            }
        }
        #endif
    }
    
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
    
    // Clean up local macros
    #undef PRIM_SLICE
    #undef CONJ_SLICE
}

