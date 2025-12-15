#ifndef RNG_H
#define RNG_H

// ====================================================================================
// PCG RANDOM NUMBER GENERATOR MODULE - HEADER
// ====================================================================================
// PCG-based random number generation for Hermitian matrix generation
//
// Dependencies: stdint.h (for uint64_t)
// ====================================================================================

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// ====================================================================================
// FUNCTION DECLARATIONS
// ====================================================================================

// Initialize global PCG generators array (NxNxN)
// seed: Base seed for reproducibility
//
// Allocates one generator per Y-slice index (0 to L/2, inclusive)
// Each generator is initialized to produce independent random sequences
void initialize_global_pcg(int L, int M, int N, uint64_t seed);

// Clean up global PCG generators
void cleanup_global_pcg();

// Generate random number using global v2rng array
// slice_index: Y-slice index (0 to N/2, inclusive)
// Returns: Random double in [0, 1)
//
// Each call advances the generator for that Y-slice by +1
double random_real_pcg_global(int slice_index);

// Advance RNG generator for a specific Y-slice
// slice_index: Y-slice index (0 to N/2, inclusive)
// n: Number of steps to advance (each step = 1 random number)
//
// Used for skipping RNG calls when N < MAX_PPD to maintain consistency
// with what a full MAX_PPD x MAX_PPD grid would generate
void advance_pcg_global(int slice_index, uint64_t n);

#ifdef __cplusplus
}
#endif

#endif

