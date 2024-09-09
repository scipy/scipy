python_sources = [
  '__init__.py',
  'test_basic.py',
  'test_blas.py',
  'test_cython_blas.py',
  'test_cython_lapack.py',
  'test_cythonized_array_utils.py',
  'test_decomp.py',
  'test_decomp_cholesky.py',
  'test_decomp_cossin.py',
  'test_decomp_ldl.py',
  'test_decomp_lu.py',
  'test_decomp_polar.py',
  'test_decomp_update.py',
  'test_extending.py',
  'test_fblas.py',
  'test_interpolative.py',
  'test_lapack.py',
  'test_matfuncs.py',
  'test_matmul_toeplitz.py',
  'test_procrustes.py',
  'test_sketches.py',
  'test_solve_toeplitz.py',
  'test_solvers.py',
  'test_special_matrices.py'
]


py3.install_sources(
  python_sources,
  subdir: 'scipy/linalg/tests',
  install_tag: 'tests'
)

py3.install_sources(
  ['_cython_examples/extending.pyx', '_cython_examples/meson.build'],
  subdir: 'scipy/linalg/tests/_cython_examples',
  install_tag: 'tests'
)

subdir('data')
