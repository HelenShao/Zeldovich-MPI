# Paths are relative to this directory unless prefixed with ../
# Source list follows Makefile groups, sorted per Meson style guide.
# =============================================================================
# Zeldovich_MPI — C sources built as C++ (they include C++ headers)
# =============================================================================
ZD_MPI_generation_cpp = configure_file(
  input: 'generation/ZD_MPI_generation.c',
  output: 'ZD_MPI_generation.cxx',
  copy: true,
)
rng_cpp = configure_file(
  input: 'utils/rng.c',
  output: 'rng.cxx',
  copy: true,
)
# =============================================================================
# Zeldovich_MPI — translation units (files())
# =============================================================================
zd_mpi_sources = files(
  '../deps/ParseHeader/generated/phParser.tab.cc',
  '../deps/ParseHeader/generated/phScanner.cc',
  '../deps/ParseHeader/src/HeaderStream.cc',
  '../deps/ParseHeader/src/ParseHeader.cc',
  '../deps/ParseHeader/src/phDriver.cc',
  '../deps/ParseHeader/src/stringutil.cc',
  '../deps/zeldovich_core/src/parameters.cpp',
  '../deps/zeldovich_core/src/power_spectrum.cpp',
  '../deps/zeldovich_core/src/STimer.cc',
  'communication/mpi_exchange.c',
  'fft/fft_setup.c',
  'fft/fft_wisdom.c',
  'ic_stage_api.cpp',
  'mpi_topology.cpp',
  'output/output_new.cpp',
  'streaming/z_streaming.c',
  'utils/batch_helpers.c',
  'utils/decomposition.c',
  'utils/plt_eigenmodes.c',
  'utils/power_spectrum.c',
  'utils/printing.c',
  'utils/verification.c',
  'utils/zeldovich_wrapper.cpp',
  'zeldovich_mpi_driver.cpp',
)
zeldovich_ic_lib_sources = zd_mpi_sources + [ZD_MPI_generation_cpp, rng_cpp]
# =============================================================================
# Static library — IC logic + deps (reused by standalone exe and Abacus subproject)
# =============================================================================
zeldovich_ic_lib = static_library(
  'zeldovich_ic',
  zeldovich_ic_lib_sources,
  dependencies: zd_mpi_common_dep,
  install: false,
)
zeldovich_ic_dep = declare_dependency(
  link_with: zeldovich_ic_lib,
  include_directories: zd_mpi_inc,
  dependencies: [mpi_dep, omp_dep, fftw_dep, fmt_dep],
)
# =============================================================================
# Executable — Zeldovich_MPI (main only; links IC static lib)
# =============================================================================
executable(
  'Zeldovich_MPI',
  ['main.cpp'],
  link_with: zeldovich_ic_lib,
  dependencies: zd_mpi_common_dep,
  install: true,
)
# =============================================================================
# Executable — wisdom_rank0 (FFTW preflight; lighter link than main)
# =============================================================================
wisdom_inc = include_directories('.', 'fft')
executable(
  'wisdom_rank0',
  ['fft/wisdom_rank0.c', 'fft/wisdom_rank0_main.c'],
  include_directories: wisdom_inc,
  dependencies: [fftw_dep, omp_dep],
  link_language: 'cpp',
  install: true,
)
