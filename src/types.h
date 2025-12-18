#ifndef HERMITIAN_TYPES_H
#define HERMITIAN_TYPES_H

#include "precision.h"

// ====================================================================================
// GRID DECOMPOSITION
// ====================================================================================

// Grid bounds structure for 2D decomposition
// Rectangular region in (X,Z) 
// Ranges: [x_start, x_end) and [z_start, z_end) - half-open intervals
typedef struct {
    int x_start, x_end;  // X-direction range [x_start, x_end)
    int z_start, z_end;  // Z-direction range [z_start, z_end)
} GridBounds;

// Extended grid bounds with overlapping regions for Abacus compatibility
// Core region  : Non-overlapping, primary responsibility of this rank
// Padded region: Extended with X_PADDING on each side, may overlap with neighbors
//
// Example (N=1024, 3×3 grid, X_PADDING=10):
//   Rank 0 core: X=[0, 342), Z=[0, 342)
//   Rank 0 padded: X=[-10, 352), Z=[-10, 352)  // Wraps at boundaries
//
// The padded region uses periodic boundary conditions to handle coordinates
// outside [0, N). This ensures neighboring ranks have consistent ghost zones.
typedef struct {
    GridBounds core;              // Core region [x_start, x_end) - primary responsibility
    GridBounds padded;            // Padded region with overlap - actual storage
    int num_pencils_core;         // Number of pencils in core region
    int num_pencils_padded;       // Number of pencils in padded region (larger)
} ExtendedGridBounds;

// ====================================================================================
// INDEXING MACROS
// ====================================================================================
// Y-slice indexing: Access array 'array_idx' in slice 'slice_idx' at (x, z)
// Memory order: [Slice][Array][Z][X] (X is stride-1, fastest varying)
// Formula: x + N * (z + N * (array_idx + narray * slice_idx))
//
// Ex: N=8, narray=4, slice_idx=2, array_idx=1, x=3, z=5:
//   Index = 3 + 8*(5 + 8*(1 + 4*2)) = 3 + 8*(5 + 8*9) = 3 + 8*77 = 619
#define Y_SLICE(slice_idx, array_idx, x, z, N, narray) \
    local_y_slices[(int64_t)(x) + (N) * ((z) + (N) * ((array_idx) + (narray) * (slice_idx)))]

// Pencil indexing: Access array 'array_idx' in pencil 'pencil_idx' at Y position 'y'
// Memory order: [Pencil][Array][Y] (Y is stride-1 for FFT)
// Formula: y + N * (array_idx + narray * pencil_idx)
//
// Y-direction data is contiguous for 1D FFT along Y
#define PENCIL(pencil_idx, array_idx, y, N, narray) \
    local_pencils[(int64_t)(y) + (N) * ((array_idx) + (narray) * (pencil_idx))]

// Z-slab indexing macro: Access array 'array_idx' at (x_idx, y) for one Z-slab
// Memory order: [X][Array][Y] (Y is stride-1 for FFT)
// Formula: y + N * (array_idx + narray * x_idx)
#define ZSLAB(x_idx, array_idx, y, N, narray) \
    local_z_slab[(int64_t)(y) + (N) * ((array_idx) + (narray) * (x_idx))]

// ====================================================================================
// PERIODIC BOUNDARY CONDITION MACROS
// ====================================================================================
// Efficient periodic boundary condition mapping for coordinates.
// Maps logical coordinates (can be negative or >= N) to actual array indices [0, N).
// Optimized for the specific range: [-X_PADDING, N+X_PADDING).
//
// Performance: Fast-path for common case (x in [0,N)) - no computation needed!
// - If x in [0, N): return x directly (branch prediction makes this ~free)
// - If x >= N: return x - N (one subtraction, handles wrap to left)
// - If x < 0: use modulo arithmetic to handle all negative cases correctly
//
// Examples (N=1024, X_PADDING=10):
//   PERIODIC_X(50, 1024)   → 50    (common case, no-op)
//   PERIODIC_X(-5, 1024)   → 1019  (wraps from right: X=-5 → X=N-5)
//   PERIODIC_X(-10, 1024)  → 1014  (wraps from right: X=-10 → X=N-10)
//   PERIODIC_X(1030, 1024) → 6     (wraps to left: X=N+6 → X=6)

#if USE_X_PADDING
    #define PERIODIC_X(x, N) \
        ((x) >= 0 ? \
         ((x) < (N) ? (x) : (x) - (N)) : \
         (((x) % (N) + (N)) % (N)))
    
    #define PERIODIC_Z(z, N) \
        ((z) >= 0 ? \
         ((z) < (N) ? (z) : (z) - (N)) : \
         (((z) % (N) + (N)) % (N)))
#else
    // No padding: coordinates are always in [0, N), so no wrapping needed
    #define PERIODIC_X(x, N) (x)
    #define PERIODIC_Z(z, N) (z)
#endif

// ====================================================================================
// SIZE AND COUNT FUNCTIONS
// ====================================================================================
// Inline functions for computing sizes and counts based on grid configuration.

// Get number of elements in a GridBounds region
static inline int get_grid_bounds_size(const GridBounds *bounds) {
    int x_count = bounds->x_end - bounds->x_start;
    int z_count = bounds->z_end - bounds->z_start;
    return x_count * z_count;
}

// Get X-count from GridBounds
static inline int get_x_count(const GridBounds *bounds) {
    return bounds->x_end - bounds->x_start;
}

// Get Z-count from GridBounds
static inline int get_z_count(const GridBounds *bounds) {
    return bounds->z_end - bounds->z_start;
}

// Get total memory size for Y-slices
static inline size_t get_y_slice_memory_bytes(int num_slices, int N, int narray) {
    return (size_t)num_slices * narray * N * N * BYTES_PER_COMPLEX;
}

// Get total memory size for pencils
static inline size_t get_pencil_memory_bytes(int num_pencils, int N, int narray) {
    return (size_t)num_pencils * narray * N * BYTES_PER_COMPLEX;
}

// Get total memory size for Z-slab
static inline size_t get_z_slab_memory_bytes(int x_count, int N, int narray) {
    return (size_t)x_count * narray * N * BYTES_PER_COMPLEX;
}

// ====================================================================================
// GRID PRINTING (for debugging)
// ====================================================================================

// Print GridBounds structure
static inline void print_grid_bounds(const GridBounds *bounds, const char *label, int rank) {
    printf("[Rank %d] %s: X=[%d, %d), Z=[%d, %d), size=%d pencils\n",
           rank, label,
           bounds->x_start, bounds->x_end,
           bounds->z_start, bounds->z_end,
           get_grid_bounds_size(bounds));
}

// Print ExtendedGridBounds structure
static inline void print_extended_grid_bounds(const ExtendedGridBounds *ext_bounds, int rank) {
    printf("[Rank %d] Extended Grid Bounds:\n", rank);
    print_grid_bounds(&ext_bounds->core, "  Core", rank);
    print_grid_bounds(&ext_bounds->padded, "  Padded", rank);
    printf("[Rank %d]   Core pencils: %d, Padded pencils: %d\n",
           rank, ext_bounds->num_pencils_core, ext_bounds->num_pencils_padded);
}

#endif

