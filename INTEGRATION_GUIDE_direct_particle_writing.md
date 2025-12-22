# Direct Particle Writing from main.cpp

## Overview

Instead of writing complex FFT output and post-processing to particles, you can now write particle ICs directly from `main.cpp` using `WriteParticlesSlab_range_from_zslab()`.

## Integration into main.cpp

### Step 1: Include output_new header

Add to includes section (around line 85):
```cpp
#include "output/output_new.h"  // For WriteParticlesSlab_range_from_zslab
```

### Step 2: Initialize output system (before Z-loop)

Add after line 1406 (after allocating local_z_slab):
```cpp
// Initialize zeldovich-PLT parameter and output system
Parameters* zel_param = nullptr;
if (param_file != NULL) {
    zel_param = new Parameters(fs::path(param_file));
    SetupOutputDir(*zel_param);
    InitOutputBuffers(*zel_param);
} else {
    // Create minimal parameters for particle output
    std::string tmp_param_file = "/tmp/main_particle_output.par";
    FILE* tmp_fp = fopen(tmp_param_file.c_str(), "w");
    if (tmp_fp) {
        fprintf(tmp_fp, "BoxSize = 1000.0\n");
        fprintf(tmp_fp, "ZD_Version = 2\n");
        fprintf(tmp_fp, "ZD_PPD = %d\n", N);
        fprintf(tmp_fp, "ZD_Pk_scale = 1.0\n");
        fprintf(tmp_fp, "NP = %d\n", N * N * N);
        fprintf(tmp_fp, "ZD_NumBlock = 2\n");
        fprintf(tmp_fp, "CPD = %d\n", N);
        fprintf(tmp_fp, "ZD_Pk_norm = 1.0\n");
        fprintf(tmp_fp, "ZD_Pk_sigma = 1.0\n");
        fprintf(tmp_fp, "ZD_Pk_smooth = 0.0\n");
        fprintf(tmp_fp, "ZD_Pk_powerlaw_index = -2.0\n");
        fprintf(tmp_fp, "InitialConditionsDirectory = \"./output_particles\"\n");
        fprintf(tmp_fp, "InitialRedshift = 0.0\n");
        fprintf(tmp_fp, "ICFormat = \"RVZel\"\n");
        fprintf(tmp_fp, "ZD_Seed = 12345\n");
        fprintf(tmp_fp, "ZD_qPLT = 0\n");
        fprintf(tmp_fp, "ZD_qdensity = 1\n");
        fprintf(tmp_fp, "ZD_qascii = 0\n");
        fclose(tmp_fp);
        zel_param = new Parameters(fs::path(tmp_param_file));
        SetupOutputDir(*zel_param);
        InitOutputBuffers(*zel_param);
    }
}
```

### Step 3: Call particle writer after 1D FFT (in Z-loop)

Replace or add after the current output writing section (lines 1574-1620):

**Option A: Replace current .bin writing**
```cpp
// =======================================================================================
// PARTICLE OUTPUT WRITING (direct from Z-slab)
// =======================================================================================
// Write particle ICs directly instead of complex FFT output
// Uses WriteParticlesSlab_range_from_zslab with [Array][X][Y] layout
// =======================================================================================

#ifndef SKIP_FILE_WRITE
if (zel_param != NULL) {
    int i = z;  // i = Z coordinate (Zeldovich i)
    
    // Get this rank's X-range
    int k_start_global = my_extended_bounds.core.x_start;
    int k_extent = x_count;  // x_count already calculated above
    
    // Call particle writer directly with local_z_slab
    WriteParticlesSlab_range_from_zslab(
        rank,              // MPI rank
        i,                 // i index (Z slab)
        k_start_global,    // global X start
        k_extent,          // number of X values
        local_z_slab,      // [narray][x_count][N] data
        N,                 // Grid size
        narray,            // Number of arrays
        *zel_param         // Parameters
    );
    
    files_written++;
    if (DEBUG_PRINTS && rank == 0 && files_written <= 3) {
        printf("[DEBUG] Rank %d: Wrote particles for i=%d, X=[%d,%d)\n", 
               rank, i, k_start_global, k_start_global + k_extent);
    }
}
#endif
```

**Option B: Write both .bin and particles**
```cpp
// =======================================================================================
// OUTPUT WRITING (i,j,k notation)
// =======================================================================================
// Write BOTH complex FFT output (.bin) AND particle ICs
// =======================================================================================

#ifndef SKIP_FILE_WRITE
int i = z;  // i = Z coordinate (Zeldovich i)

// 1. Write complex FFT output (.bin files) - EXISTING CODE
char filename[256];
snprintf(filename, sizeof(filename), "rank_%d/i%d_slab_N%d.bin", 
        rank, i, N);

FILE *fp = fopen(filename, "wb");
if (fp) {
    // ... existing .bin writing code ...
    fclose(fp);
    files_written++;
}

// 2. Write particle ICs directly - NEW CODE
if (zel_param != NULL) {
    int k_start_global = my_extended_bounds.core.x_start;
    int k_extent = x_count;
    
    WriteParticlesSlab_range_from_zslab(
        rank, i, k_start_global, k_extent,
        local_z_slab, N, narray, *zel_param
    );
}
#endif
```

### Step 4: Cleanup (after Z-loop)

Add before `MPI_Finalize()`:
```cpp
// Cleanup particle output system
if (zel_param != NULL) {
    TeardownOutput();
    delete zel_param;
}
```

## Output Files

The function writes per-rank particle files:
- **Particles**: `ic_rank{rank}_i{i}_k{k_start}_{k_end}` (where i = Z-slab index, k = X extent)
- **Density** (if enabled): `dens_rank{rank}_i{i}_k{k_start}_{k_end}` (where i = Z-slab index, k = X extent)

These files contain particle data in the same format as Zeldovich-PLT output:
- `RVZel`: 6 floats per particle (3 displ + 3 vel)
- `RVdoubleZel`: 6 doubles per particle
- `Zeldovich`: 3 floats per particle (displ only)
- `ZelSimple`: 3 floats per particle (no grid indices)

## Advantages of Direct Writing

1. **Memory efficient**: No intermediate `.bin` storage
2. **Faster**: Single-pass processing (FFT → particles)
3. **Simpler**: No post-processing scripts needed
4. **Same output**: Identical particle format to WriteParticlesSlab_range

## When to Use Each Approach

| Approach | Use When |
|----------|----------|
| **Path 1** (write .bin, post-process) | Debugging FFT output, need FFT data for analysis |
| **Path 2** (reassemble, full ICs) | Need complete `ic_{i}` files (Abacus compatibility) |
| **Path 3** (direct writing) | Production runs, memory/storage limited |

## Example Workflow

### Compile with output_new
```bash
# Add to Makefile
SRCS += src/output/output_new.cpp
INCLUDES += -I$(ZELDOVICH_PLT_DIR)/src
LIBS += -L$(ZELDOVICH_PLT_DIR) -lzeldovich
```

### Run with parameter file
```bash
mpirun -np 4 ./main 256 param_N256.par
```

### Output
```
./output_particles/
├── ic_rank0_i0_k0_64      (rank 0, i=0, k=[0,64))
├── ic_rank1_i0_k64_128    (rank 1, i=0, k=[64,128))
├── ic_rank2_i0_k128_192   (rank 2, i=0, k=[128,192))
└── ic_rank3_i0_k192_256   (rank 3, i=0, k=[192,256))
```

## Notes

- The function uses the **same physics** as `WriteParticlesSlab_new` and `WriteParticlesSlab_range`
- **No transpose overhead**: Works directly with `[array][x][y]` layout
- **Global k indices**: Particles have global `(i,j,k)` coordinates
- **Thread-safe**: Uses output_mutex for statistics updates

