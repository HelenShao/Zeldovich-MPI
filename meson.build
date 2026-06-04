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
proj_defs += '-DPARTICLE_OUTPUT_MODE=@0@'.format(get_option('particle_output_mode'))
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
if get_option('use_fftw_wisdom')
  proj_defs += '-DUSE_FFTW_WISDOM'
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
# Libraries — FFTW (pkg-config or fftw_root + find_library)
# =============================================================================
_fftw_root = get_option('fftw_root')
if _double_prec
  _omp_lib_base = 'fftw3'
  _fftw_pc = dependency('fftw3', required: false)
else
  _omp_lib_base = 'fftw3f'
  _fftw_pc = dependency('fftw3f', required: false)
endif
if _fftw_pc.found()
  # find_library does not use pkg-config -L paths; pass fftw libdir (needed when
  # fftw3f / fftw3f_omp live under a non-default prefix, e.g. local float build).
  _fftw_libdir = _fftw_pc.get_variable(pkgconfig: 'libdir')
  _fftw_omp_lib = cxx.find_library(
    _omp_lib_base + '_omp',
    dirs: _fftw_libdir,
    required: true,
  )
  fftw_dep = declare_dependency(
    dependencies: [_fftw_pc, _fftw_omp_lib],
  )
elif _fftw_root != ''
  fftw_inc = include_directories(join_paths(_fftw_root, 'include'))
  fftw_libdir = join_paths(_fftw_root, 'lib')
  fftw_omp = cxx.find_library(_omp_lib_base + '_omp', dirs: fftw_libdir, required: true)
  fftw_base = cxx.find_library(_omp_lib_base, dirs: fftw_libdir, required: true)
  fftw_dep = declare_dependency(
    include_directories: fftw_inc,
    dependencies: [fftw_omp, fftw_base],
  )
else
  error(
    'FFTW not found via pkg-config and fftw_root is empty. ' +
    'Load an fftw module or pass -Dfftw_root=/path/to/fftw/prefix',
  )
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
_abacus_embed = get_option('abacus_embed')
if _abacus_embed
  # Parent Abacus must call meson.override_dependency('parseheader', ...) first.
  parseheader_dep = dependency('parseheader')
  zd_mpi_inc = include_directories(_zd_mpi_inc_base)
else
  parseheader_dep = disabler()
  zd_mpi_inc = include_directories(
    _zd_mpi_inc_base + [
      'deps/ParseHeader/generated',
      'deps/ParseHeader/include',
    ],
  )
endif
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
# =============================================================================
# Subdirectories (executables in src/meson.build)
# =============================================================================
subdir('src')