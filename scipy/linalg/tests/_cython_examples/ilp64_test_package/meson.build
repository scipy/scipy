project('ilp64_test_package', 'c', 'cython')

fs = import('fs')
py = import('python').find_installation(pure: false)
cy = meson.get_compiler('cython')

if not cy.version().version_compare('>=3.2.0')
  error('Requires Cython >= 3.2.0')
endif

cython_args = ['-Xfreethreading_compatible=True']

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
  input: 'src/ilp64_test_package/_blas_int.pxi.in',
  output: '_blas_int.pxi',
  configuration: _blas_int_conf,
)

subdir('src/ilp64_test_package')
