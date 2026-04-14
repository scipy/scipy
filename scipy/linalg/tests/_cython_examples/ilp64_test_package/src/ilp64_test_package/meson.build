_wrappers_pxd = fs.copyfile('_blas_lapack_wrappers.pxd')

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
  _wrappers_pxd,
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
  '_blas_lapack_wrappers.pxd',
  subdir: 'ilp64_test_package',
)
