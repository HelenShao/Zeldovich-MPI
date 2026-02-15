// ====================================================================================
// GRID DECOMPOSITION 
// ====================================================================================

#include "utils/decomposition.h"
#include <math.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

// ====================================================================================
// Calculate grid factors (grid_x * grid_z = num_ranks)
// ====================================================================================

void calculate_grid_factors(int num_ranks, int *grid_x_out, int *grid_z_out)
{
    int grid_x, grid_z;
    
    // AbacusAurora: 81^2 = 6561 nodes
    if (num_ranks == 6561) {
        grid_x = 81; grid_z = 81;
    } else if (num_ranks == 8) {
        grid_x = 2; grid_z = 4;
    } else if (num_ranks == 16) {
        grid_x = 4; grid_z = 4;
    } else if (num_ranks == 4) {
        grid_x = 2; grid_z = 2;
    } else if (num_ranks == 2) {
        grid_x = 1; grid_z = 2;
    } else {
        // Default: try to make square-ish
        grid_x = (int)sqrt_t((double)num_ranks);
        while (num_ranks % grid_x != 0) grid_x--;
        grid_z = num_ranks / grid_x;
    }
    
    *grid_x_out = grid_x;
    *grid_z_out = grid_z;
}

// Return 1 if n is prime, 0 otherwise. n must be >= 2.
static int is_prime(int n)
{
    if (n < 2) return 0;
    if (n == 2) return 1;
    if (n % 2 == 0) return 0;
    int d = 3;
    while (d * d <= n) {
        if (n % d == 0) return 0;
        d += 2;
    }
    return 1;
}

// When num_pencil_ranks is prime: use (num_pencil_ranks - 1) for grid (e.g. 13 -> 12 -> 3x4).
// Rank (num_pencil_ranks - 1) gets empty bounds (idle). Composite odd numbers (9, 15, etc.) use full num_pencil_ranks.
// Optional: when s_first_idle_rank >= 0, ranks [s_first_idle_rank, num_ranks) also get empty.
static int s_first_idle_rank = -1;

void decomposition_set_first_idle_rank(int first_idle_rank)
{
    s_first_idle_rank = first_idle_rank;
}

GridBounds get_grid_bounds(int dest, int N, int num_pencil_ranks)
{
    GridBounds bounds;

    if (s_first_idle_rank >= 0 && dest >= s_first_idle_rank) {
        bounds.x_start = bounds.x_end = 0;
        bounds.z_start = bounds.z_end = 0;
        return bounds;
    }
    // Last rank is idle for prime num_pencil_ranks
    if (is_prime(num_pencil_ranks) && dest == num_pencil_ranks - 1) {
        bounds.x_start = bounds.x_end = 0;
        bounds.z_start = bounds.z_end = 0;
        return bounds;
    }

    int effective_ranks = is_prime(num_pencil_ranks) ? num_pencil_ranks - 1 : num_pencil_ranks;
    int grid_x, grid_z;
    calculate_grid_factors(effective_ranks, &grid_x, &grid_z);

    int x_block = dest / grid_z;
    int z_block = dest % grid_z; 
    
    // X-range: Divide NxN slice into chunks with remainder handling
    int base_x = N / grid_x;
    int remainder_x = N % grid_x;
    if (x_block < remainder_x) {
        // First 'remainder_x' blocks get an extra element
        bounds.x_start = x_block * (base_x + 1);
        bounds.x_end = bounds.x_start + base_x + 1;
    } else {
        // Remaining blocks get base amount, need to offset by values in the first remainder_x blocks
        // which each contain (base_x + 1) values and the first (x_block - remainder_x) blocks get base_x values
        bounds.x_start = remainder_x * (base_x + 1) + (x_block - remainder_x) * base_x;
        bounds.x_end = bounds.x_start + base_x;
    }
    
    // Z-range: Divide NxN slice into chunks with remainder handling
    int base_z = N / grid_z;
    int remainder_z = N % grid_z;
    if (z_block < remainder_z) {
        // First 'remainder_z' blocks get an extra element
        bounds.z_start = z_block * (base_z + 1);
        bounds.z_end = bounds.z_start + base_z + 1;
    } else {
        // Remaining blocks get base amount
        bounds.z_start = remainder_z * (base_z + 1) + (z_block - remainder_z) * base_z;
        bounds.z_end = bounds.z_start + base_z;
    }
    
    return bounds;
}

ExtendedGridBounds get_extended_grid_bounds(int rank, int N, int num_ranks, int grid_x, int grid_z)
{
    ExtendedGridBounds ext_bounds;
    
    // Suppress unused parameter warnings (grid_x/grid_z kept for compatibility)
    (void)grid_x;
    (void)grid_z;
    
    // Step 1: Get core (non-overlapping) bounds using original logic
    GridBounds core = get_grid_bounds(rank, N, num_ranks);
    ext_bounds.core = core;
    
#if USE_X_PADDING
    // Step 3: Add X-padding (NO CLAMPING - V14 change!)
    // Allow negative and >N values for periodic boundary conditions
    ext_bounds.padded.x_start = core.x_start - X_PADDING;  // Can be negative!
    ext_bounds.padded.x_end = core.x_end + X_PADDING;      // Can exceed N!
    
    // V14: NO clamping! These logical coordinates will be mapped via PERIODIC_X() macro
    // - If x_start < 0: Data wraps from right side (x=-10 -> accesses x=N-10)
    // - If x_end > N: Data wraps to left side (x=N+5 -> accesses x=5)
#else
    // No padding: padded region equals core region
    ext_bounds.padded.x_start = core.x_start;
    ext_bounds.padded.x_end = core.x_end;
#endif
    
    // TODO: Add Z-periodic padding?
    ext_bounds.padded.z_start = core.z_start;
    ext_bounds.padded.z_end = core.z_end;
    
    // Step 4: Calculate number of pencils for both regions
    // Note: Padded region size is still (x_end - x_start) even if coordinates are out of [0,N)
    ext_bounds.num_pencils_core = (core.x_end - core.x_start) * (core.z_end - core.z_start);
    ext_bounds.num_pencils_padded = (ext_bounds.padded.x_end - ext_bounds.padded.x_start) * 
                                    (ext_bounds.padded.z_end - ext_bounds.padded.z_start);
    
    return ext_bounds;
}

// Unifed function to get padded or core bounds depending on USE_X_PADDING
GridBounds get_padded_bounds_simple(int dest, int N, int num_ranks)
{
    // Calculate grid factors using centralized function
    int grid_x, grid_z;
    calculate_grid_factors(num_ranks, &grid_x, &grid_z);
    
    // Get extended bounds and return chunk
    ExtendedGridBounds ext = get_extended_grid_bounds(dest, N, num_ranks, grid_x, grid_z);
#if USE_X_PADDING
    return ext.padded;  // Return padded chunk
#else
    return ext.core;    // Return core chunk 
#endif
}

#ifdef __cplusplus
}
#endif
