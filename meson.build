# Root Meson build for Zeldovich_MPI (see plan ZD_MPI_meson)
# Style: https://mesonbuild.com/Style-guide.html
# =============================================================================
# Project
# =============================================================================
project(
  'Zeldovich_MPI',
  ['c', 'cpp'],
  version: '0.1.0',
  meson_version: '>=1.7.0',
  default_options: [
    'c_std=gnu11',
    'cpp_std=c++17',
    'warning_level=2',
    'buildtype=release',
    'optimization=3',
  ],
)
# =============================================================================
# Compiler (toolchains + global flags)
# =============================================================================
cc = meson.get_compiler('c')
cxx = meson.get_compiler('cpp')
# Intel icx often stays on strict C; match glibc + Make for M_PI / posix_memalign.
add_project_arguments('-D_GNU_SOURCE', language: 'c')
_supported_cpp_warns = cxx.get_supported_arguments(['-Wno-cast-function-type'])
add_project_arguments(_supported_cpp_warns, language: 'cpp')
# =============================================================================
# User options → preprocessor (see meson.options)
# =============================================================================
_double_prec = get_option('double_precision')
proj_defs = []
if _double_prec
  proj_defs += '-DUSE_DOUBLE_PRECISION'
endif
if get_option('production_mode')
  proj_defs += '-DPRODUCTION_MODE'
endif
proj_defs += '-DSKIP_FILE_WRITE=@0@'.format(get_option('skip_file_write'))
proj_defs += '-DPARALLELIZE_Z_LOOP=@0@'.format(get_option('parallelize_z_loop'))
if get_option('debug_prints')
  proj_defs += '-DDEBUG_PRINTS=1'
else
  proj_defs += '-DDEBUG_PRINTS=0'
endif
if get_option('debug_rng_consistency')
  proj_defs += '-DDEBUG_RNG_CONSISTENCY=1'
else
  proj_defs += '-DDEBUG_RNG_CONSISTENCY=0'
endif
if get_option('debug_rng_skip')
  proj_defs += '-DDEBUG_RNG_SKIP=1'
else
  proj_defs += '-DDEBUG_RNG_SKIP=0'
endif
foreach _a : get_option('extra_cpp_args')
  proj_defs += _a
endforeach
add_project_arguments(proj_defs, language: ['c', 'cpp'])
# Optional sanitizer (compile + link)
if get_option('enable_asan')
  _asan = ['-fsanitize=address', '-g', '-fno-omit-frame-pointer']
  add_project_arguments(_asan, language: ['c', 'cpp'])
  add_project_link_arguments(_asan, language: ['c', 'cpp'])
endif
# =============================================================================
# Dependencies (external packages)
# =============================================================================
fmt_dep = dependency(
  'fmt',
  fallback: ['fmt', 'fmt_dep'],
  required: true,
  default_options: ['default_library=static'],
)
mpi_dep = dependency('mpi', language: 'cpp')
omp_dep = dependency('openmp')
# =============================================================================
# Libraries — FFTW (pkg-config; e.g. module load fftw on Aurora), or oneMKL's FFTW3 interface
# Threaded plan API lives in libfftw3_omp / libfftw3f_omp. Aurora does not
# ship fftw3_omp.pc. Search the fftw3 libdir, and do not require that
# pkg-config variable: a cached dep still says "found" after the prefix is gone.
# =============================================================================
_fft = get_option('fft')
_mkl_gnu = cxx.get_id() == 'gcc'
use_mkl_fft = false
if _fft != 'fftw'
  # The 2D plan needs MKL's OpenMP threading layer (the sequential layer ignores
  # FFTW_PLAN_WITH_NTHREADS), built for the compiler's OpenMP runtime.
  _mkl_dep = dependency(_mkl_gnu ? 'mkl-dynamic-lp64-gomp' : 'mkl-dynamic-lp64-iomp', required: _fft == 'mkl')
  use_mkl_fft = _mkl_dep.found()
endif
if use_mkl_fft
  fftw_dep = declare_dependency(
    compile_args: ['-I' + (_mkl_dep.get_variable('includedir') / 'fftw')],
    dependencies: _mkl_dep.partial_dependency(link_args: true),
  )
else
  if _double_prec
    _fftw_base_dep = dependency('fftw3', required: true)
    _fftw_omp_name = 'fftw3_omp'
  else
    _fftw_base_dep = dependency('fftw3f', required: true)
    _fftw_omp_name = 'fftw3f_omp'
  endif
  _fftw_libdir = _fftw_base_dep.get_variable(pkgconfig: 'libdir', default_value: '')
  _fftw_omp_found = cxx.find_library(
    _fftw_omp_name,
    dirs: _fftw_libdir == '' ? [] : [_fftw_libdir],
    required: true,
  )
  fftw_dep = declare_dependency(
    dependencies: [_fftw_base_dep, _fftw_omp_found],
  )
endif
if get_option('use_fftw_wisdom')
  if use_mkl_fft
    message('use_fftw_wisdom is ignored: oneMKL\'s FFTW3 interface has no wisdom')
  else
    add_project_arguments('-DUSE_FFTW_WISDOM', language: ['c', 'cpp'])
  endif
endif
# =============================================================================
# Include directories (repo layout)
# =============================================================================
_zd_mpi_inc_base = [
  'deps',
  'deps/zeldovich_core/include',
  'deps/zeldovich_core/pcg',
  'src',
  'src/communication',
  'src/fft',
  'src/generation',
  'src/output',
  'src/streaming',
  'src/utils',
]
# ParseHeader via Meson subproject (dedupes with parent Abacus when abacus_embed=true).
ph_proj = subproject('ParseHeader', default_options: ['b_sanitize=none'])
parseheader_dep = ph_proj.get_variable('parseheader_dep')
zd_mpi_inc = include_directories(_zd_mpi_inc_base)
# =============================================================================
# Composite dependency (includes + link deps for Zeldovich_MPI)
# =============================================================================
zd_mpi_common_dep = declare_dependency(
  include_directories: zd_mpi_inc,
  dependencies: [mpi_dep, omp_dep, fftw_dep, fmt_dep],
)
# =============================================================================
# Tests
# =============================================================================
# No root-level tests; add e.g. subdir('tests') when defined.
# PLT eigenmodes staged next to the subproject build dir (same pattern as zeldovich-PLT).
fs = import('fs')
fs.copyfile('eigmodes128')

# =============================================================================
# Subdirectories (executables in src/meson.build)
# =============================================================================
subdir('src')