#ifndef OUTPUT_NEW_H
#define OUTPUT_NEW_H

// Header for output_new.cpp functions
// Provides WriteParticlesSlab_range (overloaded for both [array][x][y] and [y][x] layouts),
// WriteParticlesSlab_range_from_zslab (legacy, now equivalent to WriteParticlesSlab_range),
// and WriteParticlesSlab_new (Path 2: full-range reassembly)

#include <stdio.h>
// Undefine MAX_PPD macro from config.h to avoid conflict with zeldovich-PLT's const definition
#ifdef MAX_PPD
#undef MAX_PPD
#endif
// Include zeldovich-PLT headers to get Complx and Parameters definitions
#include <zeldovich.h>
#include <parameters.h>
#include <output.h>

// WriteParticlesSlab_range - C++ function overloads (same name, different signatures)
// The compiler selects the appropriate version based on the arguments provided.

// Overload 1: [array][x][y] layout (ZSLAB format) - NO TRANSPOSE NEEDED
// This is the preferred version for new code as it works directly with main.cpp's ZSLAB format
void WriteParticlesSlab_range(
   int rank,
   int i,                  // i = Z index (legacy z)
   int k_start_global,     // global X start for this rank's extent
   int k_extent,           // number of X values for this rank
   Complx *slab_data,      // Data in [array][x_local][y] layout (ZSLAB format)
   int N,                  // Grid size (ppd)
   int narray,             // Number of arrays (typically 4)
   Parameters &param
);

// Overload 2: [y][x] layout (JK format) - for backward compatibility
// This version accepts separate pointers for each array in [y][x] layout
void WriteParticlesSlab_range(
   int rank,
   int i,                  // i = Z index (legacy z)
   int k_start_global,     // global X start for this rank's extent
   int k_extent,           // number of X values for this rank
   Complx *slab1,          // Array 0 in [y][x] layout
   Complx *slab2,          // Array 1 in [y][x] layout
   Complx *slab3,          // Array 2 in [y][x] layout
   Complx *slab4,          // Array 3 in [y][x] layout
   Parameters &param
);

// WriteParticlesSlab_range_from_zslab - Legacy function (now equivalent to WriteParticlesSlab_range)
// This function is kept for backward compatibility but is functionally identical to
// WriteParticlesSlab_range overload 1. New code should use WriteParticlesSlab_range instead.
void WriteParticlesSlab_range_from_zslab(
   int rank,
   int i,                  // i = Z index (legacy z)
   int k_start_global,     // global X start for this rank's extent
   int k_extent,           // number of X values for this rank (x_count)
   Complx *slab_data,      // Data in [array][x_local][y] layout (ZSLAB format)
   int N,                  // Grid size (ppd)
   int narray,             // Number of arrays (typically 4)
   Parameters &param
);

// Write full-range particle ICs (all X, all Y) for one i-slab
// Used by Path 2 reassembly tool
void WriteParticlesSlab_new(
   FILE *output,
   int i,                  // i = Z index (legacy z)
   Complx *slab1,          // Array 0 in [y][x] layout
   Complx *slab2,          // Array 1 in [y][x] layout
   Complx *slab3,          // Array 2 in [y][x] layout
   Complx *slab4,          // Array 3 in [y][x] layout
   Parameters &param
);

// Setup and teardown functions
void SetupOutputDir(Parameters &param);
double InitOutputBuffers(Parameters &param);
void TeardownOutput();

// ====================================================================================
// MODE 3: Streaming-append CPD-ordered output (Bresenham division)
// ====================================================================================
//
// Assumes Bresenham division; layout must match Abacus subslab reader.
// Slab indices and x-ranges are computed on the fly (no preallocated buffer).

// Append one z-slab's particles to the rank's output file, grouped by CPD slab.
//
// For each slab s in [slab_x_start, slab_x_end), computes firstx/lastx on the fly,
// extracts the x-segment, converts complex -> RVZel particles in (y outer, x inner)
// order, and writes a contiguous segment.
//
// One z-block in the file = [slab s0 segment][slab s1 segment]...[slab sK segment].
void AppendZSlabParticles(
    FILE *fp,                 // Open particle file to append to (caller owns)
    FILE *fp_dens,            // Open density file, or NULL if !qdensity
    int slab_x_start,         // First CPD slab index this rank owns
    int slab_x_end,           // One past last (exclusive)
    int cpd,                  // Cells per dimension
    int z,                    // Global z index for this z-slab
    int k_start_global,       // This rank's x start (global): maps x_global -> x_local = x_global - k_start_global
    int k_extent,             // This rank's x count (slab_data uses local x in [0, k_extent))
    Complx *slab_data,        // local_z_slab in [array][x_local][y] layout
    int N,                    // Grid size (ppd)
    int narray,               // Number of arrays (4)
    Parameters &param         // For ICFormat, conversion factors, qdensity
);

#endif 

