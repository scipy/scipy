odrpack = static_library('odrpack',
  [
    'odrpack/d_lpk.f',
    'odrpack/d_mprec.f',
    'odrpack/d_odr.f',
    'odrpack/dlunoc.f'
  ],
  fortran_args: '-Wno-conversion',  # silence "conversion from REAL(8) to INTEGER(4)"
)

# TODO: add support for multiple BLAS libraries!
# Note that this works as is with conda (on Linx at least), because
# openblas is detected properly via pkg-config. For multiple BLAS/LAPACK
# libraries, we'll need a custom dependency in Meson (like is done for
# scalapack).
# conda: pkgconfig.variables('openblas')
# https://github.com/mesonbuild/meson/issues/2835
# TODO: 64-bit BLAS
blas = dependency('openblas')

py3.extension_module('__odrpack',
  '__odrpack.c',
  link_with: odrpack,
  include_directories: inc_np,
  c_args: numpy_nodepr_api,
  dependencies: [py3_dep, blas],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/odr'
)


python_sources = [
  '__init__.py',
  '_add_newdocs.py',
  '_models.py',
  '_odrpack.py',
  'models.py',
  'odrpack.py'
]

py3.install_sources(
  python_sources,
  pure: false,
  subdir: 'scipy/odr'
)

subdir('tests')
