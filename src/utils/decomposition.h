#ifndef ZD_MPI_DECOMPOSITION_H
#define ZD_MPI_DECOMPOSITION_H

// ====================================================================================
// GRID DECOMPOSITION UTILITIES
// ====================================================================================
// Functions for calculating grid bounds & decompositions for
// the 2D (X,Z) pencil decomposition used for MPI comm. buffers
//
// Depends on: config.h, types.h
// ====================================================================================

#include "config.h"
#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

// Calculates cart decomposition dims from num_ranks (size_x * size_z = num_ranks)
// was: calculate_grid_factors(..., grid_x_out, grid_z_out)
void calculate_grid_factors(int num_ranks, int *size_x_out, int *size_z_out);

// CPD-aligned: size_x * size_z = num_ranks with size_x|cpd and size_z|cpd when possible.
void calculate_grid_factors_cpd_aligned(int num_ranks, int N, int cpd,
                                        int *size_x_out, int *size_z_out);

GridBounds get_grid_bounds(int dest, int N, int num_pencil_ranks);

// was: grid_x, grid_z (cart sizes; names were transposed vs Abacus)
GridBounds get_grid_bounds_CPD_aligned(int dest, int N, int num_ranks,
                                       int size_x, int size_z, int cpd);

// was: slab_to_rx(..., grid_x); divisor is NumZRanks = size_z
int slab_to_rx(int slab, int cpd, int size_z);

ExtendedGridBounds get_extended_grid_bounds_CPD_aligned(int rank, int N, int num_ranks,
                                                        int size_x, int size_z, int cpd);

ExtendedGridBounds get_extended_grid_bounds(int rank, int N, int num_ranks, int size_x, int size_z);

GridBounds get_padded_bounds_simple(int dest, int N, int num_ranks);

#ifdef __cplusplus
}
#endif

#endif
