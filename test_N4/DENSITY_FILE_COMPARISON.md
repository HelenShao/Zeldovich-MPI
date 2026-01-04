# Density File Creation Comparison: Zeldovich-PLT vs Hermitian Code

## ⚠️ IMPORTANT: Parameter Name Issue

**The parameter name in zeldovich-PLT is `ZD_qdensity`, NOT `qdensity`!**

If you use `qdensity = 1` in your parameter file, it will **NOT** be read by zeldovich-PLT, and the value will default to `0`, meaning **no density files will be created**.

**Correct parameter name:** `ZD_qdensity = 1`

**Incorrect (won't work):** `qdensity = 1`

This is why you haven't seen density files even though you've been using `qdensity = 1` in your parameter files.

## Overview

Both codes create density files when `param.qdensity` is non-zero (typically `ZD_qdensity = 1`).
The density field is extracted from `real(slab1)` and written as `float` values.

## Zeldovich-PLT Code

### Initialization (`InitOutputBuffers` in `output.cpp`)

**When:** Called in `zeldovich.cpp` line 894 (before `ZeldovichXY`)

**Condition:** `if (param.qdensity)` (line 282)

**File Creation:**
```cpp
fs::path path = param.output_dir / fmt::format(fmt::runtime(param.density_filename.string()), param.ppd);
densfp = fopen(path.c_str(), "wb");
assert(densfp != NULL);
densoutput_tmp = new float[param.ppd * param.ppd];
```

**File Name Format:**
- Default: `"density{:d}"` (from `parameters.cpp` line 29)
- Final name: `density{ppd}` (e.g., `density4` for N=4)
- Location: `param.output_dir / density{ppd}`

**File Mode:** `"wb"` (write binary, creates new file)

### Writing (`WriteParticlesSlab` in `output.cpp`)

**When:** Called for each Z-slab in `ZeldovichXY` (line 677)

**Condition:** `if (param.qdensity)` (line 217)

**Process:**
1. Density values stored in buffer: `densoutput_tmp[i] = dens` (line 196)
2. After processing all particles in slab, append to file:
   ```cpp
   fwrite(densoutput_tmp, sizeof(*densoutput_tmp) * array.ppd * array.ppd, 1, densfp);
   ```
3. File is written in Z-order (one slab at a time, appended sequentially)

**Data Format:**
- Type: `float` (4 bytes per value)
- Layout: One 2D plane (ppd × ppd) per Z-slab
- Total size: `ppd × ppd × sizeof(float) × num_z_slabs`

### Cleanup (`TeardownOutput` in `output.cpp`)

**When:** Called at end of `main()` (line 1022)

**Process:**
```cpp
if (densoutput_tmp != NULL) {
    delete[] densoutput_tmp;
    fclose(densfp);
}
```

## Hermitian Code

### Initialization (`InitOutputBuffers` in `output_new.cpp`)

**When:** Called in `main.cpp` line 380 (before particle writing)

**Condition:** `if (param.qdensity)` (line 874)

**File Creation:**
```cpp
fs::path path = param.output_dir / fmt::format(fmt::runtime(param.density_filename.string()), param.ppd);
densfp = fopen(path.c_str(), "wb");
assert(densfp != NULL);
densoutput_tmp = new float[param.ppd * param.ppd];
```

**File Name Format:**
- Same as Zeldovich-PLT: `density{ppd}`
- Location: `param.output_dir / density{ppd}`

**File Mode:** `"wb"` (write binary, creates new file)

### Writing - Full Range Mode (`WriteParticlesSlab_new`)

**When:** Called for each Z-slab when writing full-range particle ICs

**Condition:** `if (param.qdensity && density_buffer != NULL)` (line 334)

**Process:**
1. Uses global `densoutput_tmp` buffer (line 132)
2. Density values stored: `density_buffer[count] = (float)dens` (line 287)
3. Appends to global file:
   ```cpp
   fwrite(density_buffer, sizeof(*density_buffer) * num_particles, 1, densfp);
   ```
4. Same behavior as Zeldovich-PLT (single global file, appended in Z-order)

### Writing - Local Range Mode (`WriteParticlesSlab_range` / `WriteParticlesSlab_range_from_zslab`)

**When:** Called for each Z-slab when writing per-rank particle ICs (MPI mode)

**Condition:** `if (param.qdensity && local_dens_tmp != NULL)` (line 773)

**Process:**
1. Allocates local buffer: `local_dens_tmp = new float[num_particles]` (line 588)
2. Density values stored: `local_dens_tmp[count] = (float)dens` (line 728)
3. Writes to **per-rank file** (NOT global file):
   ```cpp
   fs::path fn_dens = param.output_dir /
                      fmt::format("dens_rank{:d}_i{:d}_k{:d}_{:d}",
                                  rank, i, k_start_global, k_end_global);
   FILE *fp_d = fopen(fn_dens.c_str(), "wb");
   fwrite(local_dens_tmp, sizeof(float), (size_t)num_particles, fp_d);
   fclose(fp_d);
   ```

**File Name Format (Local Range):**
- Pattern: `dens_rank{rank}_i{i}_k{k_start}_{k_end}`
- Example: `dens_rank0_i0_k0_512` (rank 0, Z-slab 0, X-range 0-512)
- Location: `param.output_dir / dens_rank{rank}_i{i}_k{k_start}_{k_end}`

**Data Format:**
- Type: `float` (4 bytes per value)
- Layout: One 2D plane (ppd × k_extent) per file
- Total size per file: `ppd × k_extent × sizeof(float)`

### Cleanup (`TeardownOutput` in `output_new.cpp`)

**When:** Called at end of particle writing

**Process:**
```cpp
if (densoutput_tmp != NULL) {
    delete[] densoutput_tmp;
    fclose(densfp);
}
```

**Note:** Local buffers in `WriteParticlesSlab_range_from_zslab` are freed immediately after writing (line 821-823), not in `TeardownOutput`.

## Key Differences

| Aspect | Zeldovich-PLT | Hermitian (Full Range) | Hermitian (Local Range) |
|--------|---------------|------------------------|------------------------|
| **File Count** | 1 global file | 1 global file | Multiple per-rank files |
| **File Naming** | `density{ppd}` | `density{ppd}` | `dens_rank{rank}_i{i}_k{k_start}_{k_end}` |
| **File Mode** | `"wb"` (append) | `"wb"` (append) | `"wb"` (new file per slab) |
| **Buffer** | Global `densoutput_tmp` | Global `densoutput_tmp` | Local `local_dens_tmp` |
| **Writing** | Sequential append | Sequential append | Per-rank, per-slab files |
| **Use Case** | Serial execution | Serial/full-range reassembly | MPI parallel execution |

## Conditions for Density File Creation

**Both codes:**
- `param.qdensity != 0` (typically `ZD_qdensity = 1` in parameter file)
- `param.qdensity == 2` means "density only" mode (no displacements), but density files are still created if `ZD_qdensity = 1`

**⚠️ Parameter Name:**
- **zeldovich-PLT**: Must use `ZD_qdensity = 1` (NOT `qdensity = 1`)
- **Hermitian code**: Uses `qdensity = 1` (no `ZD_` prefix needed)

**Zeldovich-PLT:**
- Always creates single global file (serial execution)

**Hermitian:**
- Full range mode: Creates single global file (same as Zeldovich-PLT)
- Local range mode: Creates per-rank, per-slab files (MPI parallel execution)

## Summary

1. **Zeldovich-PLT**: Single global density file, written sequentially in Z-order
2. **Hermitian (Full Range)**: Same as Zeldovich-PLT - single global file
3. **Hermitian (Local Range)**: Multiple per-rank files, one per Z-slab per rank

The Hermitian code's local range mode is designed for MPI parallel execution where each rank writes its own density data to separate files, which can later be reassembled if needed.

