# Reassemble Z-Slabs and Call WriteParticlesSlab

This program reads distributed rank files (`rank_*/z*_slab_N*.bin`) and reassembles them into full Z-slabs, then calls `WriteParticlesSlab` to write particle data in the format expected by zeldovich-PLT.

## Building

```bash
cd conj-test
make -f Makefile.reassemble
```

This will create the executable `reassemble_z_slabs`.

**Requirements:**
- zeldovich-PLT must be built (see main Makefile for details)
- C++17 compiler (g++ or compatible)
- zeldovich-PLT libraries must be available

## Usage

```bash
./reassemble_z_slabs <output_dir> <N> [param_file] [z_start] [z_end]
```

### Arguments

- `output_dir`: Directory containing `rank_*/z*_slab_N*.bin` files
- `N`: Grid size (e.g., 256, 512, 1024)
- `param_file`: (Optional) Parameter file for zeldovich-PLT. If not provided, minimal parameters will be created.
- `z_start`: (Optional) First Z-slab to process (default: 0)
- `z_end`: (Optional) Last Z-slab to process (default: N)

### Examples

```bash
# Process all Z-slabs from N=256 run
./reassemble_z_slabs ./output 256

# Process with parameter file
./reassemble_z_slabs ./output 256 param.par

# Process only Z-slabs 0-10
./reassemble_z_slabs ./output 256 param.par 0 10
```

## How It Works

1. **File Reading**: Reads all `rank_*/z*_slab_N*.bin` files for each Z-slab
2. **Reassembly**: Combines distributed X-extents from all ranks into full `N × N` slabs
3. **Memory Layout**: Creates contiguous 2D arrays `[Y][X]` for each of the 4 arrays
4. **WriteParticlesSlab**: Calls zeldovich-PLT's `WriteParticlesSlab` function to write particle data

## File Format

Input files (`rank_*/z*_slab_N*.bin`):
- Format: Binary `[Array][Y][X]` order
- Each file contains: `narray` arrays × `N` Y-values × `x_count` X-values
- Precision: Single (8 bytes/complex) or double (16 bytes/complex) - auto-detected

Output file (`particles_output.bin`):
- Format: Determined by `WriteParticlesSlab` (based on parameter file settings)
- Contains particle positions, velocities, and density

## Memory Layout

The program ensures the correct memory layout for `WriteParticlesSlab`:

- **File format**: `[Array][Y][X]` (row-major, X varies quickly)
- **Memory layout**: `full_slab[array_idx * N * N + y * N + x]`
- **WriteParticlesSlab expects**: `slab[x + ppd * y]` (same layout!)

The layouts match perfectly, so no data copying/reordering is needed beyond reassembly.

## Notes

- **N must be even**: BlockArray requires `numblock` to be even, and we use `numblock=2`
- **Precision**: Automatically detects single vs double precision from file size
- **Missing files**: Skips Z-slabs where files are not found (prints warning every 10th missing slab)
- **X-count verification**: Warns if total X-count from ranks doesn't match N (allows small rounding differences)

## Troubleshooting

**Error: "Cannot open file"**
- Check that `output_dir` contains `rank_*/z*_slab_N*.bin` files
- Verify file permissions

**Error: "Failed to create BlockArray"**
- Ensure N is even
- Check that zeldovich-PLT is built correctly

**Error: "Failed to load parameter file"**
- Parameter file is optional - program will create minimal parameters
- If provided, ensure it's a valid zeldovich-PLT parameter file

**Linking errors**
- Ensure zeldovich-PLT is built: `cd ../zeldovich-PLT && meson setup build && ninja -C build`
- Check that library paths in `Makefile.reassemble` are correct
