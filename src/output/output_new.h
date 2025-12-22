#ifndef OUTPUT_NEW_H
#define OUTPUT_NEW_H

// Header for output_new.cpp functions
// Provides WriteParticlesSlab_range (Option A) and WriteParticlesSlab_range_from_zslab (Option B)

#include <stdio.h>
// Undefine MAX_PPD macro from config.h to avoid conflict with zeldovich-PLT's const definition
#ifdef MAX_PPD
#undef MAX_PPD
#endif
// Include zeldovich-PLT headers to get Complx and Parameters definitions
#include <zeldovich.h>
#include <parameters.h>

// Option A: Transpose [array][x][y] → [y][x] and call WriteParticlesSlab_range
// Expects separate pointers for each array in [y][x] layout
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

// Option B: Direct access with [array][x][y] layout (no transpose)
// Works directly with ZSLAB format from main.cpp
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

// Setup and teardown functions
void SetupOutputDir(Parameters &param);
double InitOutputBuffers(Parameters &param);
void TeardownOutput();

#endif  // OUTPUT_NEW_H

