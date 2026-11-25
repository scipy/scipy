# Check if scipy's cython_blas exports blas_int (ILP64 support).
# This is the pattern that downstream packages like scikit-learn and
# statsmodels use to detect and adapt to ILP64 builds.
scipy_has_blas_int = cy.compiles(
  'from scipy.linalg.cython_blas cimport blas_int',
  name: 'scipy cython_blas blas_int'
)

_blas_int_conf = configuration_data()
if scipy_has_blas_int
  _blas_int_conf.set('BLAS_INT_DEF',
    'from scipy.linalg.cython_blas cimport blas_int')
else
  # Backwards compat for scipy<1.18.0
  _blas_int_conf.set('BLAS_INT_DEF', 'ctypedef int blas_int')
endif

_blas_int_pxi = configure_file(
  input: '_blas_int.pxi.in',
  output: '_blas_int.pxi',
  configuration: _blas_int_conf,
)

# Part 1: Internal .pxd + wrappers (sklearn-like pattern)
py.extension_module(
  '_blas_lapack_wrappers',
  '_blas_lapack_wrappers.pyx',
  _blas_int_pxi,
  install: true,
  subdir: 'ilp64_test_package',
  cython_args: cython_args,
  c_args: ['-DCYTHON_CCOMPLEX=0'],
)

# Part 1: Consumer that cimports from the internal .pxd
py.extension_module(
  '_consumer',
  '_consumer.pyx',
  install: true,
  subdir: 'ilp64_test_package',
  cython_args: cython_args,
  c_args: ['-DCYTHON_CCOMPLEX=0'],
)

# Part 2: Direct use of scipy's blas_int
py.extension_module(
  '_direct_blas',
  '_direct_blas.pyx',
  install: true,
  subdir: 'ilp64_test_package',
  cython_args: cython_args,
  c_args: ['-DCYTHON_CCOMPLEX=0'],
)

py.install_sources(
  '__init__.py',
  subdir: 'ilp64_test_package',
)
