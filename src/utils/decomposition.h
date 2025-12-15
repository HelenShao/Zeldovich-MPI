#ifndef HERMITIAN_DECOMPOSITION_H
#define HERMITIAN_DECOMPOSITION_H

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

// Calculate (X,Z) grid bounds for a given destination rank in pencil decomposition
// Returns GridBounds struct ("bounds") with x_start, x_end, z_start, z_end attributes
// Handles remainder distribution (first 'remainder' ranks get one extra element)
GridBounds get_grid_bounds(int dest, int N, int num_pencil_ranks);

// Calculate grid bounds with X-direction padding using PERIODIC boundary conditions
// Returns ExtendedGridBounds with both core (non-overlapping) and padded (overlapping) regions
//
// Core region: [x_start, x_end) - Primary responsibility, non-overlapping
// Padded region: [x_start - X_PADDING, x_end + X_PADDING) - Actual storage with periodic wrap
//
// If USE_X_PADDING=0: padded region equals core region (no padding)
ExtendedGridBounds get_extended_grid_bounds(int rank, int N, int num_ranks, int grid_x, int grid_z);

// Simple version: Get padded bounds for a destination rank
// Returns the padded GridBounds (or core if USE_X_PADDING=0)
GridBounds get_padded_bounds_simple(int dest, int N, int num_ranks);

// Validate domain decomposition for Abacus compatibility
// Checks that N is divisible by grid_x and grid_z (exact division required)
// Returns true if compatible, false otherwise
// Prints warnings/errors to stderr
int validate_abacus_compatibility(int N, int num_ranks, int grid_x, int grid_z);

// Calculate grid factors (grid_x x grid_z = num_ranks)
// Returns grid_x and grid_z via output parameters
void calculate_grid_factors(int num_ranks, int *grid_x_out, int *grid_z_out);

#endif 

