arpack_sources = [
  'ARPACK/_arpack.h',
  'ARPACK/_arpack_n_single.c',
  'ARPACK/_arpack_n_double.c',
  'ARPACK/_arpack_n_single_complex.c',
  'ARPACK/_arpack_n_double_complex.c',
  'ARPACK/_arpack_s_single.c',
  'ARPACK/_arpack_s_double.c',
  'ARPACK/_arpack_n_single.h',
  'ARPACK/_arpack_n_double.h',
  'ARPACK/_arpack_n_single_complex.h',
  'ARPACK/_arpack_n_double_complex.h',
  'ARPACK/_arpack_s_double.h',
  'ARPACK/_arpack_s_single.h'
]

arpack_lib = static_library('_arpack',
  arpack_sources,
  c_args: ['-Wall', '-Wextra', '-Wpedantic'],
  dependencies: [m_dep, lapack_dep, blas_dep],
)

py3.extension_module('_arpacklib',
  '_arpackmodule.c',
  link_with: arpack_lib,
  dependencies: [np_dep],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/sparse/linalg/_eigen/arpack',
)

py3.install_sources([
    '__init__.py',
    'arpack.py',
  ],
  subdir: 'scipy/sparse/linalg/_eigen/arpack'
)

subdir('tests')
