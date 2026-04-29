py3.install_sources([
    '__init__.py',
    '_rbfinterp_kernel_pythran.py',
    'test_bary_rational.py',
    'test_bsplines.py',
    'test_fitpack.py',
    'test_fitpack2.py',
    'test_gil.py',
    'test_interpnd.py',
    'test_interpolate.py',
    'test_ndgriddata.py',
    'test_pade.py',
    'test_polyint.py',
    'test_rbf.py',
    'test_rbfinterp.py',
    'test_rgi.py'
  ],
  subdir: 'scipy/interpolate/tests',
  install_tag: 'tests'
)

py3.install_sources([
    'data/bug-1310.npz',
    'data/estimate_gradients_hang.npy',
    'data/gcvspl.npz'
  ],
  subdir: 'scipy/interpolate/tests/data',
  install_tag: 'tests'
)

if use_pythran
  py3.extension_module('_rbfinterp_kernel_pythran',
    pythran_gen.process('_rbfinterp_kernel_pythran.py'),
    cpp_args: cpp_args_pythran,
    dependencies: [pythran_dep, np_dep],
    link_args: version_link_args,
    install: true,
    subdir: 'scipy/interpolate/tests'
  )
else
  py3.install_sources(
    ['_rbfinterp_kernel_pythran.py'],
    subdir: 'scipy/interpolate/tests'
  )
endif