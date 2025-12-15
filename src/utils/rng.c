// ====================================================================================
// PCG RANDOM NUMBER GENERATOR MODULE
// ====================================================================================
// PCG-based random number generation for Y-slices.
// Each rank has its own generator that has been advanced to correct state
// and advances with each call + for skips.
//
// Dependencies: pcg-rng/pcg_random.hpp
// ====================================================================================

#include "rng.h"
#include "../config.h"  // For PARALLELIZE_XZ_WITHIN_SLICE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>       // For ldexp()
#include "pcg-rng/pcg_random.hpp"

// Conditional OpenMP includes (only needed when using locks)
#if PARALLELIZE_XZ_WITHIN_SLICE
#include <omp.h>
#endif

// ====================================================================================
// GLOBAL PCG GENERATORS (shared among threads within a rank)
// ====================================================================================

// Global PCG random number generators array
static pcg64 *v2rng_global = NULL;
static int v2rng_global_initialized = 0;
static int v2rng_global_size = 0;

// OpenMP locks for thread-safe generator access (only when parallelizing (x,z) within slice)
#if PARALLELIZE_XZ_WITHIN_SLICE
static omp_lock_t *v2rng_locks = NULL;
#endif

// ====================================================================================
// FUNCTIONS
// ====================================================================================

// Initialize global PCG generators array (NxNxN)
// seed: Base seed for reproducibility
//
// One generator per Y-slice: Each Y-slice (Y = 0, 1, 2, ..., N/2) has its own RNG generator.
// Hermitian symmetry: Only Y values 0 to N/2 need to be generated; the rest are determined by symmetry.
// Deterministic sequences: Each generator produces an independent, reproducible sequence for its Y-slice.
// Array size: L/2 + 1 includes both Y=0 and Y=N/2 (both are self-conjugate).
void initialize_global_pcg(int L, int M, int N, uint64_t seed) {
    if (!v2rng_global_initialized) {
        // All threads share same v2rng array
        // Thread accesses v2rng[i] where i = assigned slice index
        
        // Allocate memory for v2rng array (one generator per i-slice)
        // Need L/2 + 1 to include index L/2 (self-conjugate case: Y=0 and Y=L/2)
        v2rng_global_size = L/2 + 1;  // Number of i-slices we need to process (0 to L/2)
        // Add padding to avoid buffer overflow in PCG state
        size_t alloc_size = v2rng_global_size * sizeof(pcg64) + 128;
        v2rng_global = (pcg64*)malloc(alloc_size);
        if (v2rng_global == NULL) {
            printf("ERROR: Failed to allocate v2rng_global array!\n");
            exit(1);
        }
        memset(v2rng_global, 0, alloc_size);
        
        // Initialize first generator
        v2rng_global[0] = pcg64(seed);
        
        // Initialize remaining generators with sequential copying
        // Each generator is advanced by 2*MAX_PPD*MAX_PPD to ensure independent sequences
        // This matches zeldovich.cpp approach: v2rng[i].advance(2 * MAX_PPD * MAX_PPD)
        // Using MAX_PPD (not actual grid size) maintains consistency across different N values
        // Cast to uint64_t to avoid integer overflow warning (MAX_PPD=65536 gives 8.5e9 > INT_MAX)
        for(int i = 1; i < v2rng_global_size; i++) {
            v2rng_global[i] = v2rng_global[i-1]; 
            v2rng_global[i].advance((uint64_t)2 * MAX_PPD * MAX_PPD); // Match zeldovich.cpp: each plane is MAX_PPD^2 complexes
        }
        
        // Initialize OpenMP locks for thread-safe access (only when parallelizing (x,z) within slice)
        #if PARALLELIZE_XZ_WITHIN_SLICE
        v2rng_locks = (omp_lock_t*)malloc(v2rng_global_size * sizeof(omp_lock_t));
        if (v2rng_locks == NULL) {
            printf("ERROR: Failed to allocate v2rng_locks array!\n");
            exit(1);
        }
        for (int i = 0; i < v2rng_global_size; i++) {
            omp_init_lock(&v2rng_locks[i]);
        }
        #endif
        
        v2rng_global_initialized = 1;
    }
}

// Clean up global PCG generators
void cleanup_global_pcg() {
    if (v2rng_global_initialized) {
        // Destroy OpenMP locks (only if they were allocated)
        #if PARALLELIZE_XZ_WITHIN_SLICE
        if (v2rng_locks != NULL) {
            for (int i = 0; i < v2rng_global_size; i++) {
                omp_destroy_lock(&v2rng_locks[i]);
            }
            free(v2rng_locks);
            v2rng_locks = NULL;
        }
        #endif
        
        free(v2rng_global);
        v2rng_global = NULL;
        v2rng_global_initialized = 0;
        v2rng_global_size = 0;
    }
}

// Generate random number using global v2rng array
// slice_index: Y-slice index (0 to N/2, inclusive)
// Returns: Random double in (0, 1] (matches zeldovich-PLT's one_rand<2>)
//
// Each call advances the generator for that Y-slice by +1
// Thread-safety: Uses OpenMP locks when PARALLELIZE_XZ_WITHIN_SLICE=1
// When PARALLELIZE_XZ_WITHIN_SLICE=0, sequential access doesn't need locks
//
// Implementation matches zeldovich-PLT's PowerSpectrum::one_rand<2>():
// - Shifts [0,1) to (0,1] to avoid log(0) in Box-Muller
// - Uses ldexp for precise floating-point conversion
double random_real_pcg_global(int slice_index) {
    if (!v2rng_global_initialized || v2rng_global == NULL) {
        fprintf(stderr, "ERROR: random_real_pcg_global called before initialize_global_pcg!\n");
        exit(1);
    }
    
    if (slice_index < 0 || slice_index >= v2rng_global_size) {
        fprintf(stderr, "ERROR: slice_index %d out of bounds [0, %d)!\n", 
                slice_index, v2rng_global_size);
        exit(1);
    }
    
    // Conditional thread-safety: use locks only when parallelizing (x,z) within slice
    #if PARALLELIZE_XZ_WITHIN_SLICE
    // Thread-safe access: acquire lock, generate number (advances generator), release lock
    // This ensures generator advances naturally with each call, just like zeldovich.cpp
    omp_set_lock(&v2rng_locks[slice_index]);
    uint64_t r = v2rng_global[slice_index]();
    omp_unset_lock(&v2rng_locks[slice_index]);
    #else
    // Sequential access: no locks needed, generator advances deterministically
    uint64_t r = v2rng_global[slice_index]();
    #endif
    
    // Match zeldovich-PLT's one_rand<2>() implementation:
    // Can't return 0! That will immediately break the log in Box-Muller
    // But 1.0 is a valid value, so we shift the domain from [0,1) to (0,1]
    // But first check if adding 1 would overflow
    if (r == UINT64_MAX) return 1.0;
    
    // Shift [0,1) to (0,1] by adding 1
    // This ensures we never return 0, which would cause log(0) in Box-Muller
    r += (uint64_t)1;
    
    // Convert using ldexp for precise floating-point conversion (matches zeldovich-PLT)
    // ldexp(r, -64) is equivalent to r * 2^-64, but more precise
    return ldexp((double)r, -64);
}

// Advance RNG generator for a specific Y-slice
// slice_index: Y-slice index (0 to N/2, inclusive)
// n: Number of steps to advance (each step = 1 random number)
//
// Used for skipping RNG calls when N < MAX_PPD to maintain consistency
// with what a full MAX_PPD × MAX_PPD grid would generate
void advance_pcg_global(int slice_index, uint64_t n) {
    if (!v2rng_global_initialized || v2rng_global == NULL) {
        fprintf(stderr, "ERROR: advance_pcg_global called before initialize_global_pcg!\n");
        exit(1);
    }
    
    if (slice_index < 0 || slice_index >= v2rng_global_size) {
        fprintf(stderr, "ERROR: slice_index %d out of bounds [0, %d)!\n", 
                slice_index, v2rng_global_size);
        exit(1);
    }
    
    if (n == 0) {
        return;  // Nothing to advance
    }
    
    // Conditional thread-safety: use locks only when parallelizing (x,z) within slice
    #if PARALLELIZE_XZ_WITHIN_SLICE
    // Thread-safe access: acquire lock, advance generator, release lock
    omp_set_lock(&v2rng_locks[slice_index]);
    v2rng_global[slice_index].advance(n);
    omp_unset_lock(&v2rng_locks[slice_index]);
    #else
    // Sequential access: no locks needed
    v2rng_global[slice_index].advance(n);
    #endif
}

