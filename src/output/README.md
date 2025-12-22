#  Output Module

Output functions for the MPI IC code, with backward compatibility to the Zeldovich-PLT output format.

## Axis Convention

### New Notation (axis0/axis1/axis2)

- **axis0** = Z coordinate (slab index, perpendicular to slabs)
- **axis1** = Y coordinate
- **axis2** = X coordinate

### Legacy Zeldovich Mapping

For backward compatibility with existing Abacus loadIC:

| Legacy Variable |  axis | Physical Meaning |
|----------------|----------------|------------------|
| `z` | axis0 | Z coordinate |
| `y` | axis1 | Y coordinate |
| `x` | axis2 | X coordinate |
| `out.displ[0]` | axis0 displacement | Z-displacement |
| `out.displ[1]` | axis1 displacement | Y-displacement |
| `out.displ[2]` | axis2 displacement | X-displacement |

## Data Organization

### Slab Structure

- **Output format**: axis0-distributed slabs (Z-slabs)
- **Each file contains**: All axis2 (X) × All axis1 (Y) for one axis0 (Z) value
- **File naming**: `rank_{rank}/z{z}_slab_N{N}.bin`

### MPI Distribution

After 3D FFT, data is X-distributed (axis2-distributed):
- Each rank owns: `axis2_range × all_axis1 × all_axis0`
- For each Z-slab: each rank writes its `axis2_range` for all `axis1`

### Iteration Order

For a single axis0 (Z) value:
```c
for (axis1 = 0; axis1 < N; axis1++)     // All Y values
    for (axis2 = 0; axis2 < N; axis2++)   // All X values
        Particle at (axis2, axis1, axis0)
```

Output particle struct:
```c
particle.i = axis0  // Z coordinate
particle.j = axis1  // Y coordinate  
particle.k = axis2  // X coordinate
particle.displ[0] = axis0_displacement  // Z-displacement
particle.displ[1] = axis1_displacement  // Y-displacement
particle.displ[2] = axis2_displacement  // X-displacement
```

## Usage

### 1. Initialize Output Directory

```c
#include "output/output_.h"

// Create output directory structure
CreateOutputDirForRank(output_dir, rank);

// Write metadata (rank 0 only)
if (rank == 0) {
    WriteMetadata(output_dir, N, num_ranks, OUTPUT_RVZEL, write_density);
}
```

### 2. Extract Z-Slices and Write

```c
#include "output/slab_extraction.h"

// Allocate buffers for 2D Z-slices [axis1][axis2_local]
int64_t slice_size = N * axis2_count;
real_t *density_slice = malloc(slice_size * sizeof(real_t));
real_t *displ_axis0_slice = malloc(slice_size * sizeof(real_t));
real_t *displ_axis1_slice = malloc(slice_size * sizeof(real_t));
real_t *displ_axis2_slice = malloc(slice_size * sizeof(real_t));
// ... allocate velocity slices if needed ...

// Loop over all Z values
for (int axis0 = 0; axis0 < N; axis0++) {
    // Extract Z-slice from 3D X-distributed data
    ExtractZSlice(density_3d, density_slice, axis0, N, axis2_count);
    ExtractZSlice(displ_axis0_3d, displ_axis0_slice, axis0, N, axis2_count);
    ExtractZSlice(displ_axis1_3d, displ_axis1_slice, axis0, N, axis2_count);
    ExtractZSlice(displ_axis2_3d, displ_axis2_slice, axis0, N, axis2_count);
    // ... extract velocity slices if needed ...
    
    // Write Z-slab for this rank's X-range
    WriteParticlesSlab(
        axis0,                  // Z index
        density_slice,
        displ_axis0_slice,      // Z-displacement
        displ_axis1_slice,      // Y-displacement
        displ_axis2_slice,      // X-displacement
        vel_axis0_slice,        // Z-velocity (or NULL)
        vel_axis1_slice,        // Y-velocity (or NULL)
        vel_axis2_slice,        // X-velocity (or NULL)
        N,
        axis2_start,            // Starting X index for this rank
        axis2_count,            // Number of X values this rank owns
        output_dir,
        OUTPUT_RVZEL,           // Particle format
        rank,
        write_density           // Write separate density file
    );
}

// Clean up
free(density_slice);
free(displ_axis0_slice);
// ... free other slices ...
```

## Supported Output Formats

| Format | Description | Precision | Contains |
|--------|-------------|-----------|----------|
| `OUTPUT_RVDOUBLEZEL` | Positions + Displacements + Velocities | double | i,j,k + 3 displ + 3 vel |
| `OUTPUT_RVZEL` | Positions + Displacements + Velocities | float | i,j,k + 3 displ + 3 vel |
| `OUTPUT_ZEL` | Positions + Displacements | float | i,j,k + 3 displ |
| `OUTPUT_ZEL_SIMPLE` | Displacements only | float | 3 displ |

## Output Files

### Per-Rank Output

Each MPI rank writes its own files:
```
{output_dir}/
├── rank_0/
│   ├── z0_slab_N1024.bin
│   ├── z1_slab_N1024.bin
│   ├── ...
│   └── z1023_slab_N1024.bin
├── rank_1/
│   ├── z0_slab_N1024.bin
│   └── ...
└── ic_metadata.json
```

### Metadata File

`ic_metadata.json` contains:
- Format version
- Axis convention documentation
- Grid size and number of slabs
- Particle format
- Number of MPI ranks
- File naming pattern

### Reassembly

To reassemble complete Z-slabs from per-rank files, use the existing reassembly scripts or tools. Each complete Z-slab will contain data from all ranks (all X values).

## Backward Compatibility

The output format is designed to be backward compatible with existing Zeldovich-PLT output:

1. **Same particle struct layout**: i,j,k fields and displ/vel arrays in same order
2. **Same iteration order**: axis1 then axis2 (Y then X) for each Z value
3. **Same file semantics**: Each file represents one Z-slab
4. **Compatible with loadIC.cpp**: After reassembly, files can be read by existing Abacus loadIC

The key difference is clearer naming (axis0/1/2) and per-rank file organization before final assembly.

## Future Enhancements

Planned enhancements (not yet implemented):

1. **Chunked format**: Alternative output with `ic_chunk_{x_block}_{z_block}.bin` files
2. **MPI-IO collective writes**: All ranks write to same file simultaneously
3. **HDF5 output**: Alternative hierarchical format for easier analysis
4. **Compression**: On-the-fly compression of output files

## Notes

- The current implementation writes per-rank files separately (no MPI collective I/O yet)
- Each rank writes its X-range contribution to each Z-slab
- Final reassembly required to create complete Z-slab files readable by loadIC
- Memory layout is optimized for extraction: Z-slices are contiguous in memory


