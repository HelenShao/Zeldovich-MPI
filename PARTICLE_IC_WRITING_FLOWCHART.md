# Particle Initial Conditions Writing Flowchart

This document illustrates all available paths for writing particle initial conditions from the MPI-decomposed z-slabs.

## Overview

After computing the 3D FFT and z-streaming, each MPI rank has:
- **Data**: `local_z_slab` in `[array][x_local][y]` layout
- **Content**: One z-slab (i-slab) with local X-range and all Y values
- **Format**: Complex FFT output (4 arrays: density, displacement/velocity components)

## All Available Options (Summary)

### Currently Implemented:

1. **Option A** (Direct, in `main.cpp`): Transpose `[array][x][y]` → `[y][x]`, then call `WriteParticlesSlab_range()`
   - Output: Per-rank IC files
   - Flag: `USE_PARTICLE_OUTPUT_OPTION_B = 0` (default)

2. **Option B** (Direct, in `main.cpp`): Use `local_z_slab` directly, call `WriteParticlesSlab_range_from_zslab()`
   - Output: Per-rank IC files
   - Flag: `USE_PARTICLE_OUTPUT_OPTION_B = 1`

3. **Fallback** (Direct, in `main.cpp`): Write complex .bin files to disk
   - Output: `rank_{rank}/i{i}_slab_N{N}.bin`
   - Condition: No parameter file provided

4. **Path 1** (Post-process, separate script): Read .bin files, call `WriteParticlesSlab_range()` per rank
   - Script: `conj-test/write_particles_from_ranks.cpp`
   - Output: Per-rank IC files

5. **Path 2** (Post-process, separate script): Read .bin files, reassemble, call `WriteParticlesSlab_new()`
   - Script: `conj-test/write_particles_from_reassembled_mpi.cpp`
   - Output: Full IC files (all X, all Y)

### Not Implemented (Potential Future):

6. **Option D** (Potential): Gather to rank 0 in `main.cpp`, write full IC files
   - Would integrate Path 2 logic directly into `main.cpp`

## Decision Tree

```
┌─────────────────────────────────────────────────────────────────────────┐
│                    START: After Z-Streaming Loop                        │
│              Each rank has local_z_slab[array][x_local][y]               │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
                    ┌─────────────────────────────────┐
                    │  param_file provided?           │
                    │  (params != NULL)                │
                    └─────────────────────────────────┘
                            │                    │
                    ┌───────┘                    └───────┐
                    │ YES                                │ NO
                    ▼                                    ▼
    ┌───────────────────────────────────┐   ┌──────────────────────────────┐
    │  DIRECT PARTICLE WRITING          │   │  FALLBACK: Write .bin files   │
    │  (Option A or B)                  │   │  rank_%d/i%d_slab_N%d.bin    │
    └───────────────────────────────────┘   └──────────────────────────────┘
                    │                                    │
                    ▼                                    │
    ┌───────────────────────────────────┐               │
    │  USE_PARTICLE_OUTPUT_OPTION_B?    │               │
    └───────────────────────────────────┘               │
            │                    │                       │
    ┌───────┘                    └───────┐              │
    │ = 0 (Option A)                     │ = 1 (Option B)│
    ▼                                     ▼              │
    ┌─────────────────────────┐  ┌──────────────────────┐│
    │ OPTION A: Transpose     │  │ OPTION B: Direct     ││
    │ [array][x][y] → [y][x]  │  │ [array][x][y]        ││
    │ Allocate 4 transposed   │  │ No transpose         ││
    │ slabs (2× memory peak)  │  │ No extra allocation  ││
    └─────────────────────────┘  └──────────────────────┘│
            │                             │              │
            ▼                             ▼              │
    ┌─────────────────────────┐  ┌──────────────────────┐│
    │ WriteParticlesSlab_range│  │WriteParticlesSlab_   ││
    │ (rank, i, k_start,      │  │range_from_zslab      ││
    │  k_extent, T_slab1..4)  │  │(rank, i, k_start,    ││
    └─────────────────────────┘  │ k_extent, local_z_   ││
            │                     │ slab, N, narray)     ││
            │                     └──────────────────────┘│
            │                             │              │
            └─────────────┬───────────────┘              │
                          │                              │
                          ▼                              │
            ┌─────────────────────────────┐             │
            │ OUTPUT: Per-Rank IC Files    │             │
            │ ic_rank{rank}_i{i}_k{k1}_{k2}│             │
            │ dens_rank{rank}_i{i}_k{k1}_{k2}│          │
            └─────────────────────────────┘             │
                                                         │
                                                         ▼
                                    ┌────────────────────────────────────┐
                                    │  SEPARATE POST-PROCESSING SCRIPTS │
                                    └────────────────────────────────────┘
                                                    │
                        ┌───────────────────────────┼───────────────────────────┐
                        │                           │                           │
                        ▼                           ▼                           ▼
        ┌──────────────────────────┐  ┌──────────────────────────┐  ┌──────────────┐
        │ PATH 1:                  │  │ PATH 2:                  │  │ (Not used)   │
        │ write_particles_from_    │  │ write_particles_from_     │  │              │
        │ ranks.cpp                │  │ reassembled_mpi.cpp      │  │              │
        └──────────────────────────┘  └──────────────────────────┘  └──────────────┘
                        │                           │
                        ▼                           ▼
        ┌──────────────────────────┐  ┌──────────────────────────┐
        │ Read per-rank .bin files │  │ Read per-rank .bin files │
        │ rank_%d/i%d_slab_N%d.bin │  │ rank_%d/i%d_slab_N%d.bin │
        └──────────────────────────┘  └──────────────────────────┘
                        │                           │
                        ▼                           ▼
        ┌──────────────────────────┐  ┌──────────────────────────┐
        │ Extract pointers to      │  │ ReassembleISlabFromRanks │
        │ 4 arrays (per rank)      │  │ Combine all ranks into   │
        └──────────────────────────┘  │ full [narray][Y][X] slab │
                        │               └──────────────────────────┘
                        ▼                           │
        ┌──────────────────────────┐                ▼
        │ WriteParticlesSlab_range│  ┌──────────────────────────┐
        │ (rank, i, k_start,       │  │ WriteParticlesSlab_new   │
        │  k_extent, slab1..4)    │  │ (i, slab1..4, all X, Y)  │
        └──────────────────────────┘  └──────────────────────────┘
                        │                           │
                        ▼                           ▼
        ┌──────────────────────────┐  ┌──────────────────────────┐
        │ OUTPUT: Per-Rank IC Files│  │ OUTPUT: Full IC Files     │
        │ ic_rank{rank}_i{i}_      │  │ ic_{ic_index}            │
        │ x{k1}_{k2}              │  │ (all X, all Y, per i)    │
        │ dens_rank{rank}_i{i}_    │  │ dens (appended)          │
        │ x{k1}_{k2}              │  │ (matches serial output)  │
        └──────────────────────────┘  └──────────────────────────┘
```

## Detailed Path Descriptions

### Path A: Direct Writing - Option A (Transpose)

**Location**: `main.cpp` lines 1616-1683  
**Flag**: `USE_PARTICLE_OUTPUT_OPTION_B = 0` (default)

**Steps**:
1. **Allocate transposed slabs**: 4 buffers of size `x_count × N` (layout: `[y][x]`)
2. **Transpose data**: Copy from `local_z_slab[array][x][y]` → `T_slab[y][x]`
3. **Call**: `WriteParticlesSlab_range(rank, i, k_start_global, k_extent, T_slab1..4, params)`
   - This internally calls `WriteParticlesSlab_unified()` with local-range parameters
   - Uses `YX_LOCAL` macro for indexing (stride = `k_extent`)
4. **Free transposed slabs**

**Memory**: 2× peak (original + transposed)  
**Time**: ~100-200 ms transpose overhead per z-slab  
**Output**: `ic_rank{rank}_i{i}_k{k_start}_{k_end}` per rank

**Best for**: N ≤ 8192 (when transpose overhead is acceptable)

---

### Path B: Direct Writing - Option B (No Transpose)

**Location**: `main.cpp` lines 1586-1614  
**Flag**: `USE_PARTICLE_OUTPUT_OPTION_B = 1`

**Steps**:
1. **No allocation**: Use `local_z_slab` directly
2. **No transpose**: Access data in `[array][x][y]` layout
3. **Call**: `WriteParticlesSlab_range_from_zslab(rank, i, k_start_global, k_extent, local_z_slab, N, narray, params)`

**Memory**: 1× (no extra allocation)  
**Time**: No transpose overhead  
**Output**: `ic_rank{rank}_i{i}_k{k_start}_{k_end}` per rank

**Best for**: N ≥ 8192 (saves ~0.85 GB + 60s per rank for N=32K)

---

### Path C: Fallback - Write .bin Files

**Location**: `main.cpp` lines 1684-1720  
**Condition**: `param_file == NULL` (no parameter file provided)

**Steps**:
1. **Write complex FFT output**: `rank_{rank}/i{i}_slab_N{N}.bin`
2. **Format**: Binary file with 4 arrays in `[Array][j][k]` layout (transposed during write)
3. **No particle conversion**: Raw complex data only

**Purpose**: Allows post-processing with separate scripts (Path 1 or Path 2)

---

### Path 1: Post-Process from .bin Files (Per-Rank)

**Location**: `conj-test/write_particles_from_ranks.cpp`

**Steps**:
1. **Read per-rank .bin files**: `rank_{rank}/i{i}_slab_N{N}.bin`
2. **Extract pointers**: Get pointers to 4 arrays from loaded data
   - Data layout: `[narray][Y][x_local]` in memory
   - Extract: `slab1 = &local_slab[0 * N * k_extent]`, `slab2 = &local_slab[1 * N * k_extent]`, etc.
3. **Determine X-range**: Use grid decomposition to find rank's `k_start_global` and `k_extent`
4. **Call**: `WriteParticlesSlab_range(rank, i, k_start_global, k_extent, slab1..4, params)`
   - This internally calls `WriteParticlesSlab_unified()` with local-range parameters

**Output**: `ic_rank{rank}_i{i}_k{k1}_{k2}` per rank (same as Path A/B)

**Use case**: When you want to convert .bin files to particle ICs without reassembly

---

### Path 2: Post-Process from .bin Files (Reassemble)

**Location**: `conj-test/write_particles_from_reassembled_mpi.cpp`

**Steps**:
1. **Read per-rank .bin files**: `rank_{rank}/i{i}_slab_N{N}.bin` from all ranks
2. **Reassemble in memory**: `ReassembleISlabFromRanks()` combines all ranks into full `[narray][Y][X]` slab
   - **Storage**: Reassembled data stored in `std::vector<Complx> full_slab` (in memory only)
   - **NOT saved to disk**: The reassembled data is temporary and only exists in memory
3. **Extract pointers**: Get pointers to 4 arrays from reassembled data
   - Data layout: `[narray][Y][X]` in memory
   - Extract: `slab1 = &full_slab[0 * N * N]`, `slab2 = &full_slab[1 * N * N]`, etc.
4. **Call**: `WriteParticlesSlab_new(i, slab1..4, params)` (writes all X, all Y)
   - This internally calls `WriteParticlesSlab_unified()` with full-range parameters

**Output Files** (written by `WriteParticlesSlab_new`):
- **Location**: `param->output_dir` (from parameter file's `InitialConditionsDirectory`)
  - Default: `./output_complete` if creating minimal parameters
  - Can be set in parameter file: `InitialConditionsDirectory = "./output"`
- **Files**:
  - `ic_{ic_index}` where `ic_index = i * cpd / ppd` (full particle IC files, matches serial Zeldovich output)
  - `dens` (appended density file, if `qdensity` is set)

**Important**: The reassembled `[narray][Y][X]` data itself is **NOT saved as a separate file**. It exists only temporarily in memory and is immediately passed to `WriteParticlesSlab_new()` which writes the final particle IC files.

**Use case**: When you need full IC files (all X, all Y) matching serial output format

---

## Implementation Details

### Unified Function: `WriteParticlesSlab_unified()`

**Refactoring**: Both `WriteParticlesSlab_new()` and `WriteParticlesSlab_range()` now call a unified internal function `WriteParticlesSlab_unified()` that handles both full-range and local-range cases.

**Key Differences**:
- **Full-range mode** (`WriteParticlesSlab_new`):
  - Uses `YX` macro with `ppd` stride: `YX(slab, j, k, param.ppd)`
  - Uses global buffers (`output_tmp`, `densoutput_tmp`)
  - Writes to `ic_{ic_index}` files (full slab)
  - Appends density to global `dens` file
  - Loop: `j ∈ [0, ppd)`, `k ∈ [0, ppd)`
  
- **Local-range mode** (`WriteParticlesSlab_range`):
  - Uses `YX_LOCAL` macro with `k_extent` stride: `YX_LOCAL(slab, j, k_local, k_extent)`
  - Allocates local buffers per call
  - Writes to per-rank files: `ic_rank{rank}_i{i}_k{k_start}_{k_end}` (where i = Z-slab index, k = X extent)
  - Writes per-rank density files
  - Loop: `j ∈ [0, ppd)`, `k_local ∈ [0, k_extent)`
  - Calculates `k_global = k_start_global + k_local`

**Benefits**:
- Single code path for physics calculations (pos/vel extraction)
- Consistent behavior between full and local range modes
- Easier maintenance (one function to update)
- Same physics logic, different indexing and I/O strategy

## Comparison Table

| Path | Location | Memory | Time | Output Format | Best For |
|------|----------|--------|------|---------------|----------|
| **A: Option A** | `main.cpp` (direct) | 2× peak | +100-200ms | Per-rank IC files | N ≤ 8192 |
| **B: Option B** | `main.cpp` (direct) | 1× | Fastest | Per-rank IC files | N ≥ 8192 |
| **C: Fallback** | `main.cpp` | 1× | Fast | .bin files | No param file |
| **1: Post-rank** | `write_particles_from_ranks.cpp` | 1× | Medium | Per-rank IC files | Convert .bin → IC |
| **2: Post-reassemble** | `write_particles_from_reassembled_mpi.cpp` | 2× | Slowest | Full IC files | Need serial format |

## Memory Layouts

### Input (all paths start here):
```
local_z_slab[array_idx][x_idx][y]
  = local_z_slab[array_idx * x_count * N + x_idx * N + y]
```

### Option A (after transpose):
```
T_slab[y][x]
  = T_slab[y * x_count + x_idx]
```

### Option B (direct access):
```
local_z_slab[array_idx][x_idx][y]  (same as input)
```

### Path 2 (after reassembly):
```
reassembled[array_idx][j][k]
  = reassembled[array_idx * N * N + j * N + k]
```

## File Naming Conventions

### Per-Rank IC Files (Paths A, B, 1):
```
ic_rank{rank}_i{i}_k{k_start}_{k_end}
dens_rank{rank}_i{i}_k{k_start}_{k_end}
```

### Full IC Files (Path 2):
```
{output_dir}/ic_{ic_index}  (where ic_index = i * cpd / ppd)
{output_dir}/dens  (appended density file)
```
**Note**: `output_dir` comes from parameter file's `InitialConditionsDirectory` (default: `./output_complete`)

### Intermediate .bin Files (Path C):
```
rank_{rank}/i{i}_slab_N{N}.bin
```

## Potential Future Option (Not Currently Implemented)

### Option D: Gather to Rank 0 and Write Full IC Files

**Concept**: Similar to Path 2, but integrated directly into `main.cpp`:
1. **Gather**: All ranks send their local X-ranges to rank 0
2. **Reassemble**: Rank 0 assembles full `[narray][Y][X]` slab
3. **Write**: Rank 0 calls `WriteParticlesSlab_new()` to write full IC files

**Pros**:
- Produces full IC files matching serial output
- Integrated into main computation (no separate script)
- All ranks can participate in writing

**Cons**:
- Requires MPI communication (gather overhead)
- Rank 0 needs memory for full slab (2× memory on one rank)
- All ranks wait for rank 0 to finish writing

**Status**: Not implemented. Use Path 2 (separate script) instead.

---

## Recommendations

1. **For production runs with parameter file**: Use **Option B** (set `USE_PARTICLE_OUTPUT_OPTION_B = 1`)
   - Most memory efficient
   - Fastest (no transpose)
   - Direct per-rank output

2. **For small N (≤ 8192)**: Use **Option A** (default)
   - Simpler code path
   - Transpose overhead is small
   - Same output format

3. **For debugging/comparison**: Use **Path 2** (reassemble)
   - Produces full IC files matching serial output
   - Useful for validation

4. **For flexibility**: Use **Path C + Path 1**
   - Write .bin files during computation
   - Convert to ICs later with different parameters
   - No need to re-run FFT computation

