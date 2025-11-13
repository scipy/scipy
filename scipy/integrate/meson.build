vode_src = [
  'odepack/vode.f',
  'odepack/zvode.f'
]

quadpack_test_src = [
  'tests/_test_multivariate.c'
]

odepacklib = static_library('odepacklib',
  'src/lsoda.c',
)

vode_lib = static_library('vode_lib',
  vode_src,
  fortran_args: fortran_ignore_warnings,
  override_options: ['b_lto=false'],
  gnu_symbol_visibility: 'hidden',
)

dop_lib = static_library('dop_lib',
  'src/dop.c',
)

py3.extension_module('_quadpack',
  ['__quadpack.h', '__quadpack.c'],
  dependencies: [np_dep, ccallback_dep],
  install: true,
  subdir: 'scipy/integrate'
)

py3.extension_module('_odepack',
  '_odepackmodule.c',
  link_with: [odepacklib],
  link_args: version_link_args,
  dependencies: [lapack_lp64_dep, np_dep, ccallback_dep],
  install: true,
  subdir: 'scipy/integrate'
)

py3.extension_module('_vode',
  f2py_gen.process('vode.pyf'),
  link_with: [vode_lib],
  c_args: [Wno_unused_variable],
  link_args: version_link_args,
  dependencies: [lapack_lp64_dep, fortranobject_dep],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/integrate'
)

py3.extension_module('_dop',
  ['_dopmodule.c'],
  link_with: [dop_lib],
  dependencies: [np_dep, ccallback_dep],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/integrate'
)

py3.extension_module('_test_multivariate',
  [quadpack_test_src],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/integrate',
  install_tag: 'tests'
)


subdir('_ivp')
subdir('_rules')
subdir('tests')

py3.install_sources([
    '__init__.py',
    '_bvp.py',
    '_lebedev.py',
    '_ode.py',
    '_quad_vec.py',
    '_quadrature.py',
    '_odepack_py.py',
    '_quadpack_py.py',
    '_tanhsinh.py',
    '_cubature.py',
    'dop.py',
    'lsoda.py',
    'odepack.py',
    'quadpack.py',
    'vode.py',
    'src/LICENSE_DOP'
  ],
  subdir: 'scipy/integrate'
)
