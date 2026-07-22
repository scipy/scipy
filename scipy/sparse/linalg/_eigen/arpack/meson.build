arnaud_sources = [
  'arnaud/src/arnaud_n_single.c',
  'arnaud/src/arnaud_n_double.c',
  'arnaud/src/arnaud_n_single_complex.c',
  'arnaud/src/arnaud_n_double_complex.c',
  'arnaud/src/arnaud_s_single.c',
  'arnaud/src/arnaud_s_double.c'
]

arnaud_lib = static_library('arnaud',
  arnaud_sources,
  include_directories: ['arnaud/include'],
  dependencies: [m_dep, lapack_dep, np_dep, py3_dep],
  c_args: ['-DARNAUD_HAS_BLAS_CONFIG'],
)

py3.extension_module('_arpacklib',
  '_arpackmodule.c',
  link_with: arnaud_lib,
  dependencies: [np_dep, lapack_dep],
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
