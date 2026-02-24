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
// MODE 3: Streaming-append CPD-ordered output
// ====================================================================================

// Info for one CPD slab that this rank overlaps.
// Precomputed before the z-loop; used by AppendZSlabParticles() each z-iteration.
typedef struct {
    int slab_index;     // CPD slab index (0..cpd-1)
    int firstx;         // CPD slab x-range start (global)
    int lastx;          // CPD slab x-range end (global, exclusive)
    int ox_start;       // Overlap x start = max(rank_x_start, firstx)
    int ox_end;         // Overlap x end   = min(rank_x_end, lastx)
    int ox_count;       // = ox_end - ox_start (number of x-indices this rank contributes to this slab)
} CPDSlabInfo;

// Append one z-slab's particles to the rank's output file, grouped by CPD slab.
//
// For each CPD slab in slab_infos (ascending order), extracts the x-segment
// that overlaps this rank, converts complex -> RVZel particles in (y outer, x inner)
// order, and writes a contiguous segment to fp (and optionally fp_dens).
//
// One z-block in the file = [slab s0 segment][slab s1 segment]...[slab sK segment].
// Each slab segment = N * ox_count particles (y outer, x inner).
//
// Uses OpenMP parallel for over y within each slab segment (Option 2: direct to
// shared buffer, disjoint y-stripes, no contention).
void AppendZSlabParticles(
    FILE *fp,                 // Open particle file to append to (caller owns)
    FILE *fp_dens,            // Open density file, or NULL if !qdensity
    CPDSlabInfo *slab_infos,  // Array of CPD slabs this rank overlaps (sorted by slab_index)
    int num_slabs,            // Length of slab_infos array
    int z,                    // Global z index for this z-slab
    int k_start_global,       // This rank's x start (global)
    int k_extent,             // This rank's x count
    Complx *slab_data,        // local_z_slab in [array][x_local][y] layout
    int N,                    // Grid size (ppd)
    int narray,               // Number of arrays (4)
    Parameters &param         // For ICFormat, conversion factors, qdensity
);

#endif  // OUTPUT_NEW_H

