# Output Module Implementation - Complete Guide

## Summary

The output module for the MPI IC generation code has been successfully implemented with the following features:

### ✅ Completed Components

1. **Core Output Module** (`src/output/`)
   - `output_MPI.h/cpp` - Main output functions with axis0/1/2 notation
   - `output_verification.h/cpp` - Verification and testing utilities
   - `integration_example.cpp` - Example integration code
   - `README.md` - Comprehensive documentation

2. **Axis Convention**
   - **axis0** = Z coordinate (slab index, perpendicular to slabs)
   - **axis1** = Y coordinate
   - **axis2** = X coordinate
   - Maps to legacy Zeldovich: z→axis0, y→axis1, x→axis2
   - Maintains backward compatibility with existing loadIC

3. **Output Strategy**
   - Z-distributed slabs (axis0-distributed)
   - Each rank writes its X-range (axis2_range) for all Y values
   - File pattern: `rank_{rank}/z{z}_slab_N{N}.bin`
   - Supports all legacy formats: RVdoubleZel, RVZel, Zeldovich, ZelSimple

4. **Conversion Utility** (`scripts/`)
   - `convert_slabs_to_chunks.py` - Converts Z-slabs to chunked format
   - Supports custom chunk sizes
   - Includes metadata generation

## File Structure

```
InitialConditions/hermitian_3d_matrix_production/
├── src/
│   └── output/
│       ├── output_MPI.h           # Main header
│       ├── output_MPI.cpp         # Implementation
│       ├── output_verification.h        # Verification header
│       ├── output_verification.cpp      # Verification implementation
│       ├── integration_example.cpp      # Integration example
│       ├── README.md                    # Documentation
│       └── USAGE_GUIDE.md              # This file
└── scripts/
    └── convert_slabs_to_chunks.py      # Conversion utility
```

## Integration Steps

### Step 1: Add to Makefile

Add output module to your existing Makefile:

```makefile
# Output module source files
OUTPUT_SRCS = src/output/output_MPI.cpp \
              src/output/output_verification.cpp

OUTPUT_OBJS = $(OUTPUT_SRCS:.cpp=.o)

# Add to main compilation
hermitian_3d_matrix: $(OBJS) $(OUTPUT_OBJS)
	$(MPICXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Add output object compilation rules
src/output/%.o: src/output/%.cpp
	$(MPICXX) $(CXXFLAGS) -c $< -o $@
```

### Step 2: Include Headers in Main Code

In your main.cpp or wherever you handle FFT output:

```cpp
#include "output/output_MPI.h"
```

### Step 3: Call Output Function After FFT

After your 3D FFT is complete and data is in real space:

```cpp
// After FFT, data layout: [axis0][axis1][axis2_local]
// Where axis2_local is the X-distributed local range

// Prepare array of pointers to FFT output
real_t *fft_arrays[7] = {
    global_output_D,    // Density
    global_output_F,    // F → axis0 (Z) displacement
    global_output_G,    // G → axis1 (Y) displacement
    global_output_H,    // H → axis2 (X) displacement
    global_output_Fv,   // Fv → axis0 (Z) velocity (NULL if not PLT)
    global_output_Gv,   // Gv → axis1 (Y) velocity (NULL if not PLT)
    global_output_Hv    // Hv → axis2 (X) velocity (NULL if not PLT)
};

// Call output function
WriteOutputAfterFFT(
    fft_arrays,
    N,                  // Grid size
    x_start,            // Starting X index for this rank (axis2_start)
    local_x_count,      // Number of X values this rank owns (axis2_count)
    output_dir,         // Output directory path
    "RVZel",            // Particle format (from config)
    qdensity,           // Write density flag
    qPLT,               // Has velocities flag
    rank,               // MPI rank
    num_ranks,          // Total MPI ranks
    MPI_COMM_WORLD
);
```

### Step 4: Output Files

After running, each rank will have written:

```
output_dir/
├── rank_0/
│   ├── z0_slab_N1024.bin
│   ├── z1_slab_N1024.bin
│   ├── ...
│   └── z1023_slab_N1024.bin
├── rank_1/
│   └── ...
└── ic_metadata.json
```

### Step 5: (Optional) Convert to Chunked Format

For future Abacus loadIC with region loading:

```bash
python scripts/convert_slabs_to_chunks.py output_dir \
    --chunk-size-x 256 \
    --chunk-size-z 256 \
    --output output_dir_chunked
```

## Data Layout Reference

### Input to Output Module (from FFT)

After 3D FFT, data is X-distributed (axis2-distributed):

```
Each rank owns: [axis0][axis1][axis2_local]
- axis0 ∈ [0, N)              (all Z values)
- axis1 ∈ [0, N)              (all Y values)
- axis2_local ∈ [0, axis2_count)  (rank's X-range)
```

### Output Z-Slabs

Each Z-slab file contains:

```
File: rank_{rank}/z{z}_slab_N{N}.bin
Layout: [axis1][axis2_local] (row-major)
- axis1 ∈ [0, N)              (all Y values)
- axis2_local ∈ [0, axis2_count)  (rank's X-range)
```

### Particle Struct (e.g., RVZel)

```cpp
struct RVZelParticle {
    int32_t i;        // axis0 (Z)
    int32_t j;        // axis1 (Y)
    int32_t k;        // axis2 (X)
    float displ[3];   // [0]=Z, [1]=Y, [2]=X
    float vel[3];     // [0]=Z, [1]=Y, [2]=X
};
```

## Verification

### Check Output Integrity

```cpp
#include "output/output_verification.h"

// Verify a Z-slab file
VerifyParticleFile(
    "output_dir/rank_0/z0_slab_N1024.bin",
    OUTPUT_RVZEL,
    1024,
    1  // verbose
);

// Check axis convention
VerifyAxisConvention(
    "output_dir/rank_0/z0_slab_N1024.bin",
    OUTPUT_RVZEL,
    1024,
    1  // verbose
);

// Compute statistics
double stats[8];
ComputeParticleStats(
    "output_dir/rank_0/z0_slab_N1024.bin",
    OUTPUT_RVZEL,
    stats
);
printf("Displacement: min=%.3e, max=%.3e, mean=%.3e, std=%.3e\n",
       stats[0], stats[1], stats[2], stats[3]);
```

### Print Summary

```cpp
PrintOutputSummary(output_dir, rank, N, OUTPUT_RVZEL);
```

## Testing Strategy

### 1. Single Rank Test

Test with 1 MPI rank first:

```bash
mpirun -n 1 ./hermitian_3d_matrix 256 param_file.par
```

Verify:
- All Z-slabs written (256 files)
- File sizes correct
- No NaN/Inf values

### 2. Multi-Rank Test

Test with multiple ranks (e.g., 4):

```bash
mpirun -n 4 ./hermitian_3d_matrix 256 param_file.par
```

Verify:
- Each rank writes all Z-slabs
- X-ranges are non-overlapping
- Combined output covers all X values

### 3. Backward Compatibility Test

Compare against original Zeldovich output:

```bash
# Run Zeldovich (legacy)
./zeldovich seed=12345 N=256

# Run MPI code with same seed
./hermitian_3d_matrix 256 seed=12345

# Reassemble and compare
python scripts/reassemble_and_compare.py
```

## Future Enhancements (Not Yet Implemented)

1. **MPI-IO Collective Writes**
   - All ranks write to same file simultaneously
   - Better performance for large-scale runs

2. **HDF5 Output**
   - Alternative hierarchical format
   - Easier analysis and visualization

3. **Compression**
   - On-the-fly compression (blosc, zstd)
   - Reduce storage requirements

4. **Direct Chunked Output**
   - Write chunked format directly (skip conversion step)
   - Controlled by parameter

## Troubleshooting

### Issue: File not found errors

Check that output directory exists and rank subdirectories are created:

```cpp
CreateOutputDirForRank(output_dir, rank);
```

### Issue: Incorrect axis mapping

Verify data layout matches expected:
- Input: [axis0][axis1][axis2_local]
- axis0 = Z, axis1 = Y, axis2 = X

### Issue: NaN/Inf values in output

Check FFT output before writing:
- Verify inverse FFT completed
- Check for division by zero
- Verify power spectrum normalization

### Issue: File size mismatch

Expected size per Z-slab:
```
size = N × axis2_count × sizeof(ParticleStruct)
```

For RVZel with N=1024, axis2_count=256:
```
size = 1024 × 256 × 40 bytes = 10.5 MB
```

## Performance Notes

- Z-slice extraction is very fast (contiguous memory copy)
- Each rank writes independently (no communication)
- Bottleneck is typically disk I/O
- Consider using fast scratch filesystem for output
- Post-processing (reassembly) can be done offline

## Contact

For issues or questions about the output module, refer to:
- README.md in src/output/ directory
- integration_example.cpp for usage examples
- MPI IC generation documentation

---

**Implementation Complete!** All phases from the original plan have been implemented and are ready for integration.


