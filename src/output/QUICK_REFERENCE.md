# Output Module Quick Reference

## Axis Convention

```
axis0 = Z  (perpendicular to slabs)
axis1 = Y
axis2 = X

Legacy mapping: z→axis0, y→axis1, x→axis2
```

## Key Function

```cpp
int WriteParticlesSlabMPI(
    int axis0_index,              // Z value [0, N)
    const real_t *density_data,   // [axis1][axis2_local]
    const real_t *displacement_axis0,  // Z-displacement
    const real_t *displacement_axis1,  // Y-displacement
    const real_t *displacement_axis2,  // X-displacement
    const real_t *velocity_axis0,      // Z-velocity (or NULL)
    const real_t *velocity_axis1,      // Y-velocity (or NULL)
    const real_t *velocity_axis2,      // X-velocity (or NULL)
    int N,                        // Grid size
    int axis2_start,              // Starting X for this rank
    int axis2_count,              // Number of X values for this rank
    const char *output_dir,       // Output directory
    OutputFormat format,          // OUTPUT_RVZEL, etc.
    int rank,                     // MPI rank
    int write_density             // 1 = write density file
);
```

## Basic Usage

```cpp
#include "output/output_MPI.h"
#include "output/slab_extraction.h"

// Setup
CreateOutputDirForRank(output_dir, rank);
if (rank == 0) WriteMetadata(output_dir, N, num_ranks, OUTPUT_RVZEL, 1);

// Allocate slice buffers
real_t *slice_buffers[7];  // D, F, G, H, Fv, Gv, Hv
for (int i = 0; i < 7; i++) 
    slice_buffers[i] = malloc(N * axis2_count * sizeof(real_t));

// Write Z-slabs
for (int axis0 = 0; axis0 < N; axis0++) {
    // Extract Z-slice from 3D data
    ExtractZSlice(density_3d, slice_buffers[0], axis0, N, axis2_count);
    ExtractZSlice(displ_F_3d, slice_buffers[1], axis0, N, axis2_count);
    ExtractZSlice(displ_G_3d, slice_buffers[2], axis0, N, axis2_count);
    ExtractZSlice(displ_H_3d, slice_buffers[3], axis0, N, axis2_count);
    // ... extract velocities if needed ...
    
    // Write
    WriteParticlesSlabMPI(
        axis0, slice_buffers[0], slice_buffers[1], 
        slice_buffers[2], slice_buffers[3],
        slice_buffers[4], slice_buffers[5], slice_buffers[6],
        N, axis2_start, axis2_count, output_dir, 
        OUTPUT_RVZEL, rank, 1
    );
}

// Cleanup
for (int i = 0; i < 7; i++) free(slice_buffers[i]);
```

## Output Files

```
output_dir/rank_{rank}/z{z}_slab_N{N}.bin
```

Each file contains particles for:
- One axis0 (Z) value
- All axis1 (Y) values [0, N)
- Rank's axis2 (X) range [axis2_start, axis2_start+axis2_count)

## Particle Struct (RVZel)

```cpp
struct {
    int32_t i;      // axis0 (Z)
    int32_t j;      // axis1 (Y)
    int32_t k;      // axis2 (X)
    float displ[3]; // [0]=Z, [1]=Y, [2]=X
    float vel[3];   // [0]=Z, [1]=Y, [2]=X
}
```

## Data Layouts

**Input (after FFT):** `[axis0][axis1][axis2_local]`  
**Z-slice:** `[axis1][axis2_local]`  
**Iteration:** `for (axis1) for (axis2)`

## Formats

- `OUTPUT_RVDOUBLEZEL` - double precision, positions + displ + vel
- `OUTPUT_RVZEL` - float precision, positions + displ + vel  
- `OUTPUT_ZEL` - float, positions + displ
- `OUTPUT_ZEL_SIMPLE` - float, displ only

## Verification

```cpp
#include "output/output_verification.h"

VerifyParticleFile(filename, OUTPUT_RVZEL, N, 1);
VerifyAxisConvention(filename, OUTPUT_RVZEL, N, 1);
PrintOutputSummary(output_dir, rank, N, OUTPUT_RVZEL);
```

## Convert to Chunked

```bash
python scripts/convert_slabs_to_chunks.py output_dir \
    --chunk-size-x 256 --chunk-size-z 256
```

## Files

- **Header:** `src/output/output_MPI.h`
- **Implementation:** `src/output/output_MPI.cpp`
- **Extraction:** `src/output/slab_extraction.{h,cpp}`
- **Verification:** `src/output/output_verification.{h,cpp}`
- **Example:** `src/output/integration_example.cpp`
- **Test:** `src/output/test_output.c`
- **Docs:** `src/output/{README,USAGE_GUIDE,IMPLEMENTATION_SUMMARY}.md`
- **Script:** `scripts/convert_slabs_to_chunks.py`


