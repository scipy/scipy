py3.extension_module('_trlib',
  [
    lib_cython_gen.process('_trlib.pyx'),
    'trlib_krylov.c',
    'trlib_eigen_inverse.c',
    'trlib_leftmost.c',
    'trlib_quadratic_zero.c',
    'trlib_tri_factor.c'
  ],
  c_args: cython_c_args,
  include_directories: [
    '../../_lib',
    '../../_build_utils/src'
  ],
  dependencies: [lapack, blas, np_dep],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/optimize/_trlib'
)

py3.install_sources(
  ['__init__.py'],
  subdir: 'scipy/optimize/_trlib'
)
