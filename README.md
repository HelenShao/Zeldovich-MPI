# Hermitian 3D Matrix MPI - Production Version

Generation of Hermitian 3D matrices in Fourier space with real-space output using hybrid MPI+OpenMP parallelization.

## Location
/home/helenshao/InitialConditions/hermitian_3d_matrix_production

## Related Repositories

- **zeldovich-PLT**: Cloned from https://github.com/abacusorg/zeldovich-PLT
  - Located at: `/scratch/gpfs/hshao/C_Bible/InitialConditions/zeldovich-PLT`
  - Contains reference implementation and PowerSpectrum class with `cgauss()` method

## Overview

Generates large-scale 3D Hermitian matrices in Fourier space and transforms them to real space using a distributed 3D FFT. The Hermitian symmetry ensures the result is purely real, as required for N-body sims (density fields, displacement fields).

### Key Features

- **Hermitian Symmetry**: Ensures real-space output (imaginary parts ≈ 0)
- **Hybrid MPI+OpenMP**: Scales to thousands of cores across multiple nodes
- **Memory Efficient**: Streaming Z-slab processing reduces memory by 80× compared to naive approaches
- **Flexible Precision**: Single (float) or double precision compile-time selection
- **Periodic Boundaries**: Optional X-direction padding with periodic wrapping for Abacus compatibility
- **Zeldovich Compatible**: Output format matches Zeldovich code (AZYX layout)
- **PCG RNG**: Version 15+ includes PCG random number generation with nskip tracking for consistency

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
