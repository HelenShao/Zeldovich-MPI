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

    // Create local aliases for macro compatibility
    // primary_slices points to slice 0, conjugate_slices points to slice 1
    // Need to access them as if they're part of a larger buffer

    // Helper function to access primary slice (slice_idx=0)
    #define PRIM_SLICE(array_idx, x, z) \
        primary_slices[(int64_t)(x) + (N) * ((z) + (N) * (array_idx))]
    
    // Helper function to access conjugate slice (slice_idx=1)  
    #define CONJ_SLICE(array_idx, x, z) \
        conjugate_slices[(int64_t)(x) + (N) * ((z) + (N) * (array_idx))]
    
    // RNG verification: Track total RNG calls and skips to verify consistency
    // Similar to zeldovich.cpp assertion: assert(Pk.v2rng[y] - checkpoint == 2 * MAX_PPD * MAX_PPD)
    #if VERIFY_RNG_CALLS
    int64_t total_rng_calls = 0;  // Count of cgauss() calls (each uses 2 random numbers)
    int64_t total_rng_skips = 0;  // Count of skipped random numbers (from nskip advances)
    int64_t total_d_zero_skips = 0;  // Count of D=0 skips (immediate advances)
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
        nskip = (MAX_PPD - N) * MAX_PPD + N * (MAX_PPD - N);
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
            advance_pcg_global(global_y, 2 * nskip);
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
                int64_t skip_amount = (MAX_PPD - N) * MAX_PPD;
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
                double k2 = (double)(kx*kx + ky*ky + kz*kz);
                
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
                if (k2 == 0.0) {
                    // DC mode: set to zero
                    D[0] = D[1] = 0.0;
                    // RNG consistency: When D=0, we skip the RNG call, so advance immediately
                    // This matches zeldovich.cpp line 361: advance when D is forced to zero
                    #if !PARALLELIZE_XZ_WITHIN_SLICE
                    #if DEBUG_RNG_SKIP
                    int log_d0 = (x <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD && z <= MAX_DEBUG_COORD) ||
                                  (x == Nhalf - 1 && global_y == Nhalf - 1 && z == Nhalf - 1);
                    if (log_d0) {
                        fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d (x,z)=(%d,%d): D=0 (DC mode), ADVANCE RNG by 1\n",
                                N, global_y, x, z);
                        fflush(stderr);
                    }
                    #endif
                    if (ps_handle != NULL && params_handle != NULL) {
                        int64_t rng_index = global_y;
                        zeldovich_ps_advance_rng(ps_handle, params_handle, rng_index, 1);
                        #if VERIFY_RNG_CALLS
                        total_d_zero_skips++;
                        #endif
                    } else if (ps_params != NULL || params_handle == NULL) {
                        // Local PCG: advance by 2 (each complex number uses 2 random numbers)
                        advance_pcg_global(global_y, 2);
                        #if VERIFY_RNG_CALLS
                        total_d_zero_skips++;
                        #endif
                    }
                    #endif
                } else if (is_nyquist) {
                    // Nyquist frequency: set to zero (self-conjugate, no separate partner)
                    // This ensures proper mirroring alignment between first and second halves
                    D[0] = D[1] = 0.0;
                    // RNG consistency: When D=0, we skip the RNG call, so advance immediately
                    // This matches zeldovich.cpp line 361: advance when D is forced to zero
                    #if !PARALLELIZE_XZ_WITHIN_SLICE
                    #if DEBUG_RNG_SKIP
                    int log_nyq = (x <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD && z <= MAX_DEBUG_COORD) ||
                                  (x == Nhalf - 1 && global_y == Nhalf - 1 && z == Nhalf - 1);
                    if (log_nyq) {
                        fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d (x,z)=(%d,%d): D=0 (Nyquist: kx=%d ky=%d kz=%d), ADVANCE RNG by 1\n",
                                N, global_y, x, z, kx, ky, kz);
                        fflush(stderr);
                    }
                    #endif
                    if (ps_handle != NULL && params_handle != NULL) {
                        int64_t rng_index = global_y;
                        zeldovich_ps_advance_rng(ps_handle, params_handle, rng_index, 1);
                        #if VERIFY_RNG_CALLS
                        total_d_zero_skips++;
                        #endif
                    } else if (ps_params != NULL || params_handle == NULL) {
                        // Local PCG: advance by 2 (each complex number uses 2 random numbers)
                        advance_pcg_global(global_y, 2);
                        #if VERIFY_RNG_CALLS
                        total_d_zero_skips++;
                        #endif
                    }
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
                }
                
                // ========== STEP 3: Compute F, G, H from D ==========
                fftw_complex F, G, H;
                double ik2 = 1.0 / k2;
                
                // F = i x kx/k^2 x D = i x kx x ik2 x (D_re + ixD_im)
                //   = i x kx x ik2 x D_re - kx x ik2 x D_im
                //   = -kx x ik2 x D_im + i x kx x ik2 x D_re
                F[0] = -kx * ik2 * D[1];  // Real part
                F[1] =  kx * ik2 * D[0];  // Imaginary part
                
                G[0] = -ky * ik2 * D[1];
                G[1] =  ky * ik2 * D[0];
                
                H[0] = -kz * ik2 * D[1];
                H[1] =  kz * ik2 * D[0];
                
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
                // For N=256 and N=512: overlapping region is x,y,z <= 128 (Nhalf for N=256)
                //   N=256 tests (127,127,127) which is in overlapping region
                //   N=512 should also test (127,127,127) to verify it matches N=256
                int test_boundary = (boundary_coord >= 0 && 
                                     x == boundary_coord && global_y == boundary_coord && z == boundary_coord);
                // Test coordinates up to max(MAX_DEBUG_COORD, min((N/2)-1, MAX_DEBUG_BOUNDARY_COORD))
                // This ensures overlapping region boundary is tested for large N comparisons
                // For N=256: (N/2)-1 = 127 ≤ 128, so tests coordinates ≤ 127 (includes (127,127,127))
                // For N=512: (N/2)-1 = 255 > 128, so tests coordinates ≤ 128 (includes (127,127,127))
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
                    // Also always print coordinates ≤ MAX_DEBUG_COORD (small coordinate cube)
                    if (x <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD && z <= MAX_DEBUG_COORD) {
                        should_print = 1;  // Always print small coordinate cube
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
            nskip = (MAX_PPD - N) * (1 + Nhalf);
            #else
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
        // PARALLELIZATION: Conditional based on PARALLELIZE_XZ_WITHIN_SLICE flag
        #if PARALLELIZE_XZ_WITHIN_SLICE
        // Parallel z-loop: requires locks for thread-safe RNG access
        #pragma omp parallel for
        #else
        // Sequential z-loop: no locks needed
        #endif
        for (int z = 0; z <= Nhalf; z++) {
            int x_max = (z == 0 ? Nhalf + 1 : N);
            for (int x = 0; x < x_max; x++) {
                // RNG consistency: When crossing Nyquist boundary (x == Nhalf + 1),
                // skip ALL missing x-values (x = N to MAX_PPD-1) in current z-row
                // This matches zeldovich.cpp: skip at Nyquist boundary before processing negative kx region
                // The missing frequencies are in the MIDDLE of the MAX_PPD array (high positive and negative k),
                // due to FFT ordering: [0, 1, ..., N/2, -N/2+1, ..., -1]
                // High frequencies (both positive and negative) are physically located in the middle,
                // so we skip them all at once when we first enter the mirrored region
                // Note: We skip MAX_PPD - N (all missing x-values), not MAX_PPD - x_max,
                // because the missing frequencies are determined by the full grid size N, not x_max
                #if !PARALLELIZE_XZ_WITHIN_SLICE
                if (x == Nhalf + 1 && N < MAX_PPD) {
                    // Skip ALL missing x-values: from x=N to x=MAX_PPD-1
                    // This accounts for the "gap" in frequency space (high frequencies in the middle of array)
                    int64_t skip_amount = MAX_PPD - N;
                    nskip += skip_amount; // skip missing x-values in this z-row
                    
                    #if DEBUG_RNG_SKIP
                    // Debug: Log skip accumulation for test coordinates
                    int log_skip = (z <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD) ||
                                   (z == Nhalf && global_y == Nhalf);
                    if (log_skip) {
                        fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d z=%d x=%d: ACCUMULATE skip at Nyquist boundary (self-conj): +%lld (missing x-values), total nskip=%lld\n",
                                N, global_y, z, x, (long long)skip_amount, (long long)nskip);
                        fflush(stderr);
                    }
                    #endif
                }
                #endif
                
                int x_mirror = (x == 0) ? 0 : N - x;
                int z_mirror = (z == 0) ? 0 : N - z;
                
                // Calculate k-vector components first (needed for both cgauss and uniform RNG)
                int kx = (x > Nhalf) ? x - N : x;
                int ky = (global_y > Nhalf) ? global_y - N : global_y;
                int kz = (z > Nhalf) ? z - N : z;
                double k2 = (double)(kx*kx + ky*ky + kz*kz);
                
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
                if (k2 == 0.0) {
                    // DC mode: set to zero
                    D[0] = D[1] = 0.0;
                    // RNG consistency: When D=0, we skip the RNG call, so advance immediately
                    // This matches zeldovich.cpp line 361: advance when D is forced to zero
                    #if !PARALLELIZE_XZ_WITHIN_SLICE
                    #if DEBUG_RNG_SKIP
                    int log_d0 = (x <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD && z <= MAX_DEBUG_COORD) ||
                                  (x == Nhalf - 1 && global_y == Nhalf - 1 && z == Nhalf - 1);
                    if (log_d0) {
                        fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d (x,z)=(%d,%d): D=0 (DC mode), ADVANCE RNG by 1\n",
                                N, global_y, x, z);
                        fflush(stderr);
                    }
                    #endif
                    if (ps_handle != NULL && params_handle != NULL) {
                        int64_t rng_index = global_y;
                        zeldovich_ps_advance_rng(ps_handle, params_handle, rng_index, 1);
                        #if VERIFY_RNG_CALLS
                        total_d_zero_skips++;
                        #endif
                    } else if (ps_params != NULL || params_handle == NULL) {
                        // Local PCG: advance by 2 (each complex number uses 2 random numbers)
                        advance_pcg_global(global_y, 2);
                        #if VERIFY_RNG_CALLS
                        total_d_zero_skips++;
                        #endif
                    }
                    #endif
                } else if (is_nyquist) {
                    // Nyquist frequency: set to zero (self-conjugate, no separate partner)
                    // This ensures proper mirroring alignment between first and second halves
                    D[0] = D[1] = 0.0;
                    // RNG consistency: When D=0, we skip the RNG call, so advance immediately
                    // This matches zeldovich.cpp line 361: advance when D is forced to zero
                    #if !PARALLELIZE_XZ_WITHIN_SLICE
                    #if DEBUG_RNG_SKIP
                    int log_nyq = (x <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD && z <= MAX_DEBUG_COORD) ||
                                  (x == Nhalf - 1 && global_y == Nhalf - 1 && z == Nhalf - 1);
                    if (log_nyq) {
                        fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d (x,z)=(%d,%d): D=0 (Nyquist: kx=%d ky=%d kz=%d), ADVANCE RNG by 1\n",
                                N, global_y, x, z, kx, ky, kz);
                        fflush(stderr);
                    }
                    #endif
                    if (ps_handle != NULL && params_handle != NULL) {
                        int64_t rng_index = global_y;
                        zeldovich_ps_advance_rng(ps_handle, params_handle, rng_index, 1);
                        #if VERIFY_RNG_CALLS
                        total_d_zero_skips++;
                        #endif
                    } else if (ps_params != NULL || params_handle == NULL) {
                        // Local PCG: advance by 2 (each complex number uses 2 random numbers)
                        advance_pcg_global(global_y, 2);
                        #if VERIFY_RNG_CALLS
                        total_d_zero_skips++;
                        #endif
                    }
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
                    }
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
                } else {
                    // Fallback: use uniform random numbers (white noise mode, no power spectrum)
                    // This should only happen when ps_handle is NULL (no parameter file provided)
                    // Note: If ps_handle is available, we should have used cgauss() above
                    
                    // RNG consistency: Advance our own RNG when crossing Nyquist boundaries
                    // This matches the other branches' behavior
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
                // For N=256 and N=512: overlapping region is x,y,z <= 128 (Nhalf for N=256)
                //   N=256 tests (127,127,127) which is in overlapping region
                //   N=512 should also test (127,127,127) to verify it matches N=256
                int test_boundary = (boundary_coord >= 0 && 
                                     x == boundary_coord && global_y == boundary_coord && z == boundary_coord);
                // Test coordinates up to max(MAX_DEBUG_COORD, min((N/2)-1, MAX_DEBUG_BOUNDARY_COORD))
                // This ensures overlapping region boundary is tested for large N comparisons
                // For N=256: (N/2)-1 = 127 ≤ 128, so tests coordinates ≤ 127 (includes (127,127,127))
                // For N=512: (N/2)-1 = 255 > 128, so tests coordinates ≤ 128 (includes (127,127,127))
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
                    // Also always print coordinates ≤ MAX_DEBUG_COORD (small coordinate cube)
                    if (x <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD && z <= MAX_DEBUG_COORD) {
                        should_print = 1;  // Always print small coordinate cube
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
            
            // RNG consistency: After processing last z-row (z == Nhalf),
            // skip missing z-rows (z = N to MAX_PPD-1, each containing MAX_PPD x-values)
            // This matches zeldovich.cpp: conceptually skip at Nyquist boundary (z == Nhalf + 1)
            // The missing frequencies are in the MIDDLE of the MAX_PPD array (high positive and negative k),
            // due to FFT ordering: [0, 1, ..., N/2, -N/2+1, ..., -1]
            // High frequencies (both positive and negative) are physically located in the middle,
            // so we skip them all at once when we would cross the boundary
            // Note: We process z = 0 to z = Nhalf, so we never reach z == Nhalf + 1 in the loop,
            // but we need to account for missing z-rows. The missing z-rows are determined by
            // the grid size N (z = N to MAX_PPD-1), same as conjugate branch, not by how many we process.
            #if !PARALLELIZE_XZ_WITHIN_SLICE
            if (z == Nhalf && N < MAX_PPD) {
                // Skip ALL missing z-rows: from z=N to z=MAX_PPD-1
                // Number of missing z-rows: (MAX_PPD - 1) - N + 1 = MAX_PPD - N
                // Each row has MAX_PPD x-values
                // This accounts for the "gap" in frequency space (high frequencies in the middle of array)
                // Same skip amount as conjugate branch, because missing z-rows are determined by N, not by processing range
                int64_t skip_amount = (MAX_PPD - N) * MAX_PPD;
                nskip += skip_amount; // skip missing z-rows
                
                // Advance immediately since we've finished all RNG calls for this z-row
                // nskip may include missing x-values from last x-loop + missing z-rows
                if (nskip > 0) {
                    if (ps_handle != NULL && params_handle != NULL) {
                        int64_t rng_index = global_y;
                        zeldovich_ps_advance_rng(ps_handle, params_handle, rng_index, nskip);
                        #if VERIFY_RNG_CALLS
                        total_rng_skips += nskip;
                        #endif
                    } else if (params_handle == NULL) {
                        // Local PCG: advance by 2 * nskip (each complex number uses 2 random numbers)
                        advance_pcg_global(global_y, 2 * nskip);
                        #if VERIFY_RNG_CALLS
                        total_rng_skips += nskip;
                        #endif
                    }
                    nskip = 0;  // Reset after advancing
                }
            }
            #endif
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
            PRIM_SLICE(a, Nhalf, 0)[1] = 0.0;
            PRIM_SLICE(a, 0, Nhalf)[1] = 0.0;
            PRIM_SLICE(a, Nhalf, Nhalf)[1] = 0.0;
        }
        #else
        // If not using zeldovich_method for filling in self-conjugate (buggy!)
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
                int kx = (x > Nhalf) ? x - N : x;
                int ky = (global_y > Nhalf) ? global_y - N : global_y;
                int kz = (z > Nhalf) ? z - N : z;
                double k2 = (double)(kx*kx + ky*ky + kz*kz);
                
                // Nyquist frequency zeroing: Force Nyquist elements to zero for all three axes
                // This matches zeldovich.cpp line 354: abs(kx)==kmax || abs(kz)==kmax || abs(ky)==kmax
                // The Nyquist frequency (k = N/2) is self-conjugate and doesn't have a separate partner.
                // Due to the Y-shift in the reflected shell, we need to zero these to align mirroring expectations.
                int abs_kx = (kx < 0) ? -kx : kx;
                int abs_ky = (ky < 0) ? -ky : ky;
                int abs_kz = (kz < 0) ? -kz : kz;
                int is_nyquist = (abs_kx == Nhalf || abs_ky == Nhalf || abs_kz == Nhalf);
                
                // Generate D using RNG or cgauss()
                fftw_complex D;
                if (k2 == 0.0) {
                    // DC mode: set to zero
                    D[0] = D[1] = 0.0;
                    // RNG consistency: When D=0, we skip the RNG call, so advance immediately
                    #if !PARALLELIZE_XZ_WITHIN_SLICE
                    #if DEBUG_RNG_SKIP
                    int log_d0 = (x <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD && z <= MAX_DEBUG_COORD) ||
                                  (x == Nhalf - 1 && global_y == Nhalf - 1 && z == Nhalf - 1);
                    if (log_d0) {
                        fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d (x,z)=(%d,%d): D=0 (DC mode, self-conj), ADVANCE RNG by 1\n",
                                N, global_y, x, z);
                        fflush(stderr);
                    }
                    #endif
                    if (ps_handle != NULL && params_handle != NULL) {
                        int64_t rng_index = global_y;
                        zeldovich_ps_advance_rng(ps_handle, params_handle, rng_index, 1);
                        #if VERIFY_RNG_CALLS
                        total_d_zero_skips++;
                        #endif
                    } else if (params_handle == NULL) {
                        // Local PCG: advance by 2 (each complex number uses 2 random numbers)
                        advance_pcg_global(global_y, 2);
                        #if VERIFY_RNG_CALLS
                        total_d_zero_skips++;
                        #endif
                    }
                    #endif
                } else if (is_nyquist) {
                    // Nyquist frequency: set to zero (self-conjugate, no separate partner)
                    // This ensures proper mirroring alignment between first and second halves
                    D[0] = D[1] = 0.0;
                    // RNG consistency: When D=0, we skip the RNG call, so advance immediately
                    #if !PARALLELIZE_XZ_WITHIN_SLICE
                    #if DEBUG_RNG_SKIP
                    int log_nyq = (x <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD && z <= MAX_DEBUG_COORD) ||
                                  (x == Nhalf - 1 && global_y == Nhalf - 1 && z == Nhalf - 1);
                    if (log_nyq) {
                        fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d (x,z)=(%d,%d): D=0 (Nyquist: kx=%d ky=%d kz=%d, self-conj), ADVANCE RNG by 1\n",
                                N, global_y, x, z, kx, ky, kz);
                        fflush(stderr);
                    }
                    #endif
                    if (ps_handle != NULL && params_handle != NULL) {
                        int64_t rng_index = global_y;
                        zeldovich_ps_advance_rng(ps_handle, params_handle, rng_index, 1);
                        #if VERIFY_RNG_CALLS
                        total_d_zero_skips++;
                        #endif
                    } else if (params_handle == NULL) {
                        // Local PCG: advance by 2 (each complex number uses 2 random numbers)
                        advance_pcg_global(global_y, 2);
                        #if VERIFY_RNG_CALLS
                        total_d_zero_skips++;
                        #endif
                    }
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
                        double D_real, D_imag;
                        #if DEBUG_RNG_SKIP
                        int log_cgauss_sc = (x <= MAX_DEBUG_COORD && global_y <= MAX_DEBUG_COORD && z <= MAX_DEBUG_COORD) ||
                                            (x == Nhalf - 1 && global_y == Nhalf - 1 && z == Nhalf - 1);
                        if (log_cgauss_sc) {
                            fprintf(stderr, "[SKIP-DEBUG] N=%d Y=%d (x,z)=(%d,%d): CALLING cgauss (self-conj, nskip=%lld, kmag=%.6e)\n",
                                    N, global_y, x, z, (long long)nskip, kmag);
                            fflush(stderr);
                        }
                        #endif
                        #if PARALLELIZE_XZ_WITHIN_SLICE
                        // Lock protects generator access
                        #else
                        // Sequential access: no locks needed
                        #endif
                        zeldovich_ps_cgauss(ps_handle, kmag, rng_index, &D_real, &D_imag);
                        #if VERIFY_RNG_CALLS
                        total_rng_calls++;  // Each cgauss() call uses 2 random numbers
                        #endif
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
                    #if VERIFY_RNG_CALLS
                    total_rng_calls++;  // Two random_real calls = 1 complex number = 2 random numbers
                    #endif
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
                } else if ((x == 0 || x == Nhalf) && (z == 0 || z == Nhalf)) {
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
        // - Each D=0 skip advances by 2 random numbers (tracked in total_d_zero_skips, already accounts for the 2x)
        // Note: Cast to int64_t to avoid potential overflow issues
        int64_t max_ppd_val = (int64_t)MAX_PPD;
        int64_t expected_total = 2LL * max_ppd_val * max_ppd_val;
        int64_t actual_total = 2LL * total_rng_calls + 2LL * total_rng_skips + 2LL * total_d_zero_skips;
        
        if (actual_total != expected_total) {
            fprintf(stderr, 
                    "[RNG-VERIFY ERROR] Rank %d, Y=%d: Expected %ld random numbers (MAX_PPD=%ld), got %ld "
                    "(calls=%ld * 2 = %ld, skips=%ld * 2 = %ld, d_zero=%ld * 2 = %ld)\n",
                    rank, global_y, (long)expected_total, (long)max_ppd_val, (long)actual_total,
                    (long)total_rng_calls, (long)(2LL * total_rng_calls),
                    (long)total_rng_skips, (long)(2LL * total_rng_skips),
                    (long)total_d_zero_skips, (long)(2LL * total_d_zero_skips));
            fflush(stderr);
            // Don't abort in production, but warn
            #ifdef DEBUG
            assert(actual_total == expected_total);
            #endif
        } else if (rank == 0 && global_y <= 2) {
            // Print verification for first few slices in debug mode
            fprintf(stderr,
                    "[RNG-VERIFY OK] Rank %d, Y=%d: Total=%ld (calls=%ld, skips=%ld, d_zero=%ld)\n",
                    rank, global_y, (long)actual_total,
                    (long)total_rng_calls, (long)total_rng_skips, (long)total_d_zero_skips);
            fflush(stderr);
        }
    }
    #endif
    
    // Clean up local macros
    #undef PRIM_SLICE
    #undef CONJ_SLICE
}

