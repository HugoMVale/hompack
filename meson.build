project(
  'hompack',
  ['fortran', 'c'],
  version: '0.1.0',
  meson_version: '>=1.1',
  license: 'MIT',
  default_options: [
  'buildtype=release',
  'fortran_std=none'
  ],
)

cc = meson.get_compiler('c')
fc = meson.get_compiler('fortran')
os = host_machine.system()

# Compiler flags
if fc.get_id() == 'gcc'
  add_project_arguments('-Wno-compare-reals', '-Wno-unused-dummy-argument', language: ['fortran'])
elif fc.get_id() == 'llvm-flang'
  add_project_arguments('', language: ['fortran'])
else
  warning('Compiler is @0@. If the build fails, try using the GCC compiler.'.format(fc.get_id()))
endif

# BLAS
blas_option = get_option('blas_option')
blas_dep = disabler()
blas_source = files()
ref_blas_lapack_sources = files(
  'external/blas/daxpy.f',
  'external/blas/dcopy.f',
  'external/blas/ddot.f',
  'external/blas/dgemm.f',
  'external/blas/dgemv.f',
  'external/blas/dger.f',
  'external/blas/dnrm2.f90',
  'external/blas/dscal.f',
  'external/blas/dswap.f',
  'external/blas/dtpmv.f',
  'external/blas/dtpsv.f',
  'external/blas/dtrmm.f',
  'external/blas/dtrsv.f',
  'external/blas/idamax.f',
  'external/blas/lsame.f',
  'external/lapack/dgeqpf.f',
  'external/lapack/dgeqr2.f',
  'external/lapack/dgeqrf.f',
  'external/lapack/disnan.f',
  'external/lapack/dlacpy.f',
  'external/lapack/dlaic1.f',
  'external/lapack/dlaisnan.f',
  'external/lapack/dlamch.f',
  'external/lapack/dlapy2.f',
  'external/lapack/dlarf.f',
  'external/lapack/dlarf1f.f',
  'external/lapack/dlarfb.f',
  'external/lapack/dlarfg.f',
  'external/lapack/dlarft.f',
  'external/lapack/dorg2r.f',
  'external/lapack/dorgqr.f',
  'external/lapack/dorm2r.f',
  'external/lapack/dormqr.f',
  'external/lapack/ieeeck.f',
  'external/lapack/iladlc.f',
  'external/lapack/iladlr.f',
  'external/lapack/ilaenv.f',
  'external/lapack/iparmq.f',
  'external/lapack/xerbla.f',
)

if blas_option == 'refblas'
  message('`Reference BLAS` selected.')
  blas_source = ref_blas_lapack_sources

elif blas_option == 'openblas'
  message('`OpenBLAS` selected.')
  blas_dep = dependency('openblas', required: false)
  if not blas_dep.found()
    blas_dep = fc.find_library('openblas', required: false)
  endif
  if not blas_dep.found()
    error('OpenBLAS library not found. Install OpenBLAS or select another BLAS option.')
  endif

elif blas_option == 'auto'
  message('`Auto` BLAS selection enabled.')
  
  if os == 'darwin'
    blas_dep = dependency('appleframeworks', modules: ['Accelerate'], required: false)
    if blas_dep.found()
      message('Using Apple Accelerate framework (BLAS/LAPACK).')
    endif
  endif
  
  if not blas_dep.found()
    foreach lib : ['openblas']
      blas_dep = dependency(lib, required: false)
      if not blas_dep.found()
        blas_dep = fc.find_library(lib, required: false)
      endif
      if blas_dep.found()
        message('Detected BLAS library: @0@'.format(lib))
        break
      endif
    endforeach
  endif

  if not blas_dep.found()
    warning('No BLAS library detected. Falling back to reference BLAS.')
    blas_source = ref_blas_lapack_sources
  endif

else
  error('Invalid BLAS option: @0@'.format(blas_option))
endif

# Sources
hompack_source = files(
  'src/blas_interfaces.f90',
  'src/lapack_interfaces.f90',
  'src/hompack_kinds.f90',
  'src/hompack_global.f90',
  'src/hompack_interfaces.f90',
  'src/hompack.f',
  'src/hompack_core.f',
   'src/hompack_nf.f90',
)

all_sources = hompack_source + blas_source


# Library
alls_deps = []
if blas_dep.found()
  alls_deps += blas_dep
endif

if get_option('build_shared')
  hompack_lib = shared_library(
    meson.project_name(),
    sources: all_sources,
    dependencies: alls_deps,
    install: true
  )
else
  hompack_lib = static_library(
    meson.project_name(),
    sources: all_sources,
    dependencies: alls_deps
  )
endif

incdir = include_directories('c/include')
hompack_dep = declare_dependency(
  include_directories: [incdir],
  link_with: [hompack_lib]
)

# Examples
if get_option('build_examples')
  subdir('example')
endif

# Tests
if get_option('build_tests')
  subdir('test')
endif
