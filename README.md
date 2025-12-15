# Hermitian 3D Matrix MPI - Production Version

Generation of Hermitian 3D matrices in Fourier space with real-space output using hybrid MPI+OpenMP parallelization.

## Location

**Current working location**: `/scratch/gpfs/hshao/C_Bible/InitialConditions/hermitian_3d_matrix_production`

The previous location:
- `/scratch/gpfs/hshao/C_Bible/FFT_Hermitian/fftw_openmp/mpi_mock/hermitian_3d_matrix_production`

## Related Repositories

- **zeldovich-PLT**: Cloned from https://github.com/abacusorg/zeldovich-PLT
  - Located at: `/scratch/gpfs/hshao/C_Bible/InitialConditions/zeldovich-PLT`
  - Contains reference implementation and PowerSpectrum class with `cgauss()` method

## Overview

This code generates large-scale 3D Hermitian matrices in Fourier space and transforms them to real space using a distributed 3D FFT. The Hermitian symmetry ensures the result is purely real, as required for cosmological simulations (density fields, displacement fields).

### Key Features

- **Hermitian Symmetry**: Ensures real-space output (imaginary parts ≈ 0)
- **Hybrid MPI+OpenMP**: Scales to thousands of cores across multiple nodes
- **Memory Efficient**: Streaming Z-slab processing reduces memory by 80× compared to naive approaches
- **Flexible Precision**: Single (float) or double precision compile-time selection
- **Periodic Boundaries**: Optional X-direction padding with periodic wrapping for Abacus compatibility
- **Zeldovich Compatible**: Output format matches Zeldovich code (AZYX layout)
- **Production Ready**: Extensive verification, error handling, and performance monitoring
- **PCG RNG**: Version 15+ includes PCG random number generation with nskip tracking for consistency

### Performance Highlights

- **Scaling**: Tested up to N=32,768³ on 8,192 MPI ranks
- **Memory**: ~2 GB per rank for N=32K (single precision)
- **Speed**: ~50 seconds for N=1024³ on 16 nodes (16 ranks/node)

## Quick Start

### Prerequisites

- MPI (OpenMPI 4.x or MPICH 3.x)
- FFTW3 3.3.10 with OpenMP support
- GCC 7.x or later (C++11 support)
- 8 GB+ RAM per node

See [docs/DEPENDENCIES.md](docs/DEPENDENCIES.md) for detailed requirements.

### Build

```bash
# Navigate to the repository
cd /scratch/gpfs/hshao/C_Bible/InitialConditions/hermitian_3d_matrix_production

# Build with default settings (single precision, with padding)
make

# Build with double precision
make CFLAGS="-DUSE_DOUBLE_PRECISION"

# Build for production (minimal output)
make CFLAGS="-DPRODUCTION_MODE"

# Build without X-padding (saves ~10% memory)
make CFLAGS="-DUSE_X_PADDING=0"
```

For more configuration options, see [docs/CONFIGURATION.md](docs/CONFIGURATION.md).

### Run

```bash
# Small test (N=256, 14 ranks)
mpirun -np 14 ./hermitian_3d_matrix 256

# See test scripts in v15_test/ for examples
cd v15_test
./test_N256_13ranks.sh
```

### Output

The code produces Z-slab files for each MPI rank:

```
rank_0/
├── z0_slab_N256.bin
├── z1_slab_N256.bin
└── ...
rank_1/
├── z37_slab_N256.bin
├── z38_slab_N256.bin
└── ...
```

Each file contains one Z-slab in Zeldovich AZYX format:
- Layout: `[Array][Y][X]` where Array ∈ {0,1,2,3}
- Format: Binary `fftw_complex` (2×float or 2×double)
- Size: `x_count × narray × N × sizeof(fftw_complex)` bytes per file

## Architecture

### Version History

- **v15**: PCG RNG with nskip tracking and parallelize_within_slice option
  - MAX_PPD constant (65536) for RNG consistency across different grid sizes
  - nskip tracking for missing grid points (similar to zeldovich.cpp)
  - advance_pcg_global() function for RNG fast-forwarding
  - PARALLELIZE_XZ_WITHIN_SLICE option (default: 0, sequential like zeldovich)
- **v14**: Periodic boundary conditions for X-padding
- **v13**: Overlapping X-regions for Abacus compatibility
- **v12**: Z-slab streaming for Zeldovich compatibility
- **v11**: Persistent recv buffer with source-grouped layout
- **v10**: X-row streaming for memory efficiency

## Documentation

1. **Generation**: Each MPI rank generates Y-slice pairs with Hermitian symmetry
2. **2D FFT**: Apply 2D FFT (X-Z plane) to each Y-slice
3. **Redistribution**: MPI_Ialltoallv to convert Y-slices → (X,Z) pencils
4. **1D FFT**: Apply 1D FFT along Y-direction for each pencil
5. **Output**: Write Z-slabs in Zeldovich format

### Memory Layout

```
Generation (Y-slices):  [Slice][Array][Z][X]
After redistribution:   [Pencil][Array][Y]
Output (Z-slabs):       [Z][Array][Y][X]
```

## Configuration

### Compile-Time Options

Edit `src/config.h` or use `-D` flags:

| Option | Default | Description |
|--------|---------|-------------|
| `USE_DOUBLE_PRECISION` | 0 (float) | 1=double, 0=float |
| `USE_X_PADDING` | 1 | 1=periodic padding, 0=core grid only |
| `X_PADDING` | 10 | Padding size (if enabled) |
| `NARRAY` | 4 | Number of arrays per Y-slice |
| `DEBUG_PRINTS` | 0 | 1=verbose debug output |
| `VERIFY_HERMITIAN_SYMMETRY` | 0 | 1=enable verification checks |
| `USE_ZELDOVICH_METHOD` | 1 | 1=Zeldovich method for self-conjugate |

### Runtime Parameters

```bash
# Syntax: mpirun -np <ranks> ./executable <N>
mpirun -np 33 ./hermitian_3d_matrix_mpi 64

# Environment variables
export OMP_NUM_THREADS=16  # Threads per MPI rank
export FFTW_WISDOM_FILE=fftw_wisdom.dat  # FFTW optimization
```

## Project Structure

```
hermitian_3d_matrix_production/
├── README.md              # This file
├── Makefile               # Build system
├── .gitignore             # Git ignore patterns
├── src/                   # Source code
│   ├── main_original.c    # Original monolithic version (backup)
│   ├── config.h           # Configuration (Phase 2)
│   ├── types.h            # Data structures (Phase 3)
│   ├── precision.h        # Precision abstraction (Phase 3)
│   ├── generation/        # Y-slice generation (Phase 5)
│   ├── decomposition/     # Grid decomposition (Phase 5)
│   ├── communication/     # MPI communication (Phase 5)
│   ├── fft/               # FFT operations (Phase 5)
│   ├── io/                # File I/O (Phase 5)
│   └── utils/             # Utilities (Phase 4)
├── deps/                  # External dependencies
│   ├── pcg-rng/           # PCG random number generator
│   └── STimer/            # Timing utilities
├── tests/                 # Test suite
├── scripts/               # Utility scripts
│   ├── submit_job.slurm   # SLURM job submission
│   └── analyze_results.py # Output analysis
├── docs/                  # Documentation
│   ├── DEPENDENCIES.md    # Dependency documentation
│   ├── ARCHITECTURE.md    # Architecture details
│   ├── API.md             # API reference
│   └── USAGE.md           # Usage guide
└── examples/              # Example usage
```

## Testing

```bash
# Run small test suite
make test

# Run verification tests
make test_verification

# Run scaling tests
make test_scaling
```

## Performance Tuning

### Thread Configuration

```bash
# Balanced MPI+OpenMP (recommended)
export OMP_NUM_THREADS=16
mpirun -np 8 --bind-to socket ./executable 1024

# MPI-only (more ranks, fewer threads)
export OMP_NUM_THREADS=1
mpirun -np 128 ./executable 1024
```

### FFTW Wisdom

```bash
# Generate wisdom file (one-time setup)
./scripts/generate_fftw_wisdom.sh 1024

# Use wisdom in production
export FFTW_WISDOM_FILE=$PWD/fftw_wisdom.dat
```

### Memory Considerations

| N | Single Precision | Double Precision |
|---|------------------|------------------|
| 1024 | ~2 GB/rank | ~4 GB/rank |
| 4096 | ~32 GB/rank | ~64 GB/rank |
| 32768 | ~2 TB total | ~4 TB total |

## Troubleshooting

### Common Issues

1. **MPI rank mismatch**
   - Need exactly `(N/2 + 1)` ranks minimum
   - Use `N/2 + 1 + k` ranks for load balancing (k ≥ 0)

2. **Out of memory**
   - Reduce `OMP_NUM_THREADS` to increase MPI ranks
   - Use single precision instead of double
   - Check ulimit: `ulimit -s unlimited`

3. **Slow communication**
   - Use high-speed interconnect (InfiniBand)
   - Enable UCX transport: `--mca pml ucx`
   - Check MPI bindings: `--report-bindings`

4. **FFTW errors**
   - Ensure FFTW compiled with same precision
   - Check LD_LIBRARY_PATH includes FFTW lib
   - Verify OpenMP support: `ldd executable | grep fftw`

See [docs/TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md) for more details.

## Contributing

This is a production code under active development. Contributions welcome!

1. Fork the repository
2. Create a feature branch
3. Make changes with tests
4. Submit pull request

## Status

**Phase 1 Complete**:  Directory structure and documentation created

**Phase 2 Complete**:  Configuration extracted to `src/config.h` with comprehensive documentation

**Phase 3 Complete**:  Types and precision extracted to `src/types.h` and `src/precision.h`

**Phase 4 Complete**:  Utility functions extracted to `src/utils/` modules

**Current Version**: v14 (periodic boundary conditions) - Production refactoring in progress

**Code Reduction**: Original 3095 lines → Current ~2100 lines (995 lines removed, 32% reduction)

**Next Steps**:
- Phase 4: Extract utilities (utils/)
- Phase 5: Extract modules (generation/, decomposition/, fft/, io/, communication/)
- Phase 6: Final integration and testing

See [docs/PHASE1_COMPLETE.md](docs/PHASE1_COMPLETE.md), [docs/PHASE2_COMPLETE.md](docs/PHASE2_COMPLETE.md), and [docs/PHASE3_COMPLETE.md](docs/PHASE3_COMPLETE.md) for detailed completion reports.

See [../PRODUCTION_ORGANIZATION_PLAN.md](../PRODUCTION_ORGANIZATION_PLAN.md) for complete roadmap.

