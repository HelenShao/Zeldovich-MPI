// ====================================================================================
// HERMITIAN 3D MATRIX MPI - GRID DECOMPOSITION UTILITIES
// ====================================================================================

#include "utils/decomposition.h"
#include <math.h>
#include <stdio.h>

// ====================================================================================
// FUNCTION IMPLEMENTATIONS
// ====================================================================================

// Calculate grid factors (grid_x × grid_z = num_ranks)
// Handles special cases including Abacus layout (81×81 = 6561 ranks)
void calculate_grid_factors(int num_ranks, int *grid_x_out, int *grid_z_out)
{
    int grid_x, grid_z;
    
    // Abacus-specific case: 81² = 6561 nodes
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

GridBounds get_grid_bounds(int dest, int N, int num_pencil_ranks)
{
    GridBounds bounds;
    
    // Determine 2D grid factors (grid_x × grid_z = num_pencil_ranks)
    int grid_x, grid_z;
    calculate_grid_factors(num_pencil_ranks, &grid_x, &grid_z);
    
    // Determine which block in the 2D grid this rank owns
    int x_block = dest / grid_z;  // Row in processor grid
    int z_block = dest % grid_z;  // Column in processor grid
    
    // X-dimension decomposition with remainder handling
    int base_x = N / grid_x;
    int remainder_x = N % grid_x;
    if (x_block < remainder_x) {
        // First 'remainder_x' blocks get an extra element
        bounds.x_start = x_block * (base_x + 1);
        bounds.x_end = bounds.x_start + base_x + 1;
    } else {
        // Remaining blocks get base amount
        bounds.x_start = remainder_x * (base_x + 1) + (x_block - remainder_x) * base_x;
        bounds.x_end = bounds.x_start + base_x;
    }
    
    // Z-dimension decomposition with remainder handling
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
    
    // Suppress unused parameter warnings (grid_x/grid_z kept for API consistency)
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
    // - If x_start < 0: Data wraps from right side (x=-10 → accesses x=N-10)
    // - If x_end > N: Data wraps to left side (x=N+5 → accesses x=5)
#else
    // No padding: padded region equals core region
    ext_bounds.padded.x_start = core.x_start;
    ext_bounds.padded.x_end = core.x_end;
#endif
    
    // Step 4: Z-direction stays unchanged (no padding in Z for now)
    // TODO: Add Z-direction periodic padding if needed
    ext_bounds.padded.z_start = core.z_start;
    ext_bounds.padded.z_end = core.z_end;
    
    // Step 5: Calculate number of pencils for both regions
    // Note: Padded region size is still (x_end - x_start) even if coordinates are out of [0,N)
    ext_bounds.num_pencils_core = (core.x_end - core.x_start) * (core.z_end - core.z_start);
    ext_bounds.num_pencils_padded = (ext_bounds.padded.x_end - ext_bounds.padded.x_start) * 
                                    (ext_bounds.padded.z_end - ext_bounds.padded.z_start);
    
    return ext_bounds;
}

GridBounds get_padded_bounds_simple(int dest, int N, int num_ranks)
{
    // Calculate grid factors using centralized function
    int grid_x, grid_z;
    calculate_grid_factors(num_ranks, &grid_x, &grid_z);
    
    // Get extended bounds and return appropriate region
    ExtendedGridBounds ext = get_extended_grid_bounds(dest, N, num_ranks, grid_x, grid_z);
#if USE_X_PADDING
    return ext.padded;  // Return padded region
#else
    return ext.core;    // Return core region (no padding)
#endif
}

// Validate domain decomposition for Abacus compatibility
int validate_abacus_compatibility(int N, int num_ranks, int grid_x, int grid_z)
{
    int valid = 1;  // 1 = valid, 0 = invalid
    
    // Check exact division (Abacus requirement)
    if (N % grid_x != 0) {
        fprintf(stderr, "[WARNING] N=%d not divisible by grid_x=%d (Abacus requires exact division)\n",
                N, grid_x);
        fprintf(stderr, "          Remainder: %d. This may cause uneven load distribution.\n",
                N % grid_x);
        valid = 0;
    }
    
    if (N % grid_z != 0) {
        fprintf(stderr, "[WARNING] N=%d not divisible by grid_z=%d (Abacus requires exact division)\n",
                N, grid_z);
        fprintf(stderr, "          Remainder: %d. This may cause uneven load distribution.\n",
                N % grid_z);
        valid = 0;
    }
    
    // Check if matches Abacus layout (81×81 nodes, N=6075)
    if (grid_x == 81 && grid_z == 81 && N == 6075) {
        int cells_per_node_x = N / grid_x;
        int cells_per_node_z = N / grid_z;
        if (cells_per_node_x == 75 && cells_per_node_z == 75) {
            printf("[INFO] Abacus-compatible decomposition detected:\n");
            printf("       Grid: %d×%d nodes (81×81)\n", grid_x, grid_z);
            printf("       Cells per node: %d (X) × %d (Y) × %d (Z)\n",
                   cells_per_node_x, N, cells_per_node_z);
            printf("       Total: %d×%d×%d cells/node\n",
                   cells_per_node_x, N, cells_per_node_z);
        }
    } else if (valid) {
        // Valid decomposition but not Abacus-specific
        int cells_per_node_x = N / grid_x;
        int cells_per_node_z = N / grid_z;
        if (N % grid_x == 0 && N % grid_z == 0) {
            printf("[INFO] Exact division decomposition: %d×%d nodes, %d×%d×%d cells/node\n",
                   grid_x, grid_z, cells_per_node_x, N, cells_per_node_z);
        }
    }
    
    return valid;
}

