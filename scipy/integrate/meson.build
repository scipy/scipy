py3.extension_module('_quadpack',
  ['__quadpack.h', '__quadpack.c'],
  dependencies: [np_dep, ccallback_dep],
  install: true,
  subdir: 'scipy/integrate'
)

py3.extension_module('_odepack',
  ['_odepackmodule.c', 'src/lsoda.c'],
  link_args: version_link_args,
  dependencies: [lapack_lp64_dep, np_dep, ccallback_dep],
  install: true,
  subdir: 'scipy/integrate'
)

py3.extension_module('_vode',
  ['_dzvodemodule.c', 'src/vode.c', 'src/zvode.c'],
  link_args: version_link_args,
  dependencies: [lapack_lp64_dep, np_dep, ccallback_dep],
  install: true,
  subdir: 'scipy/integrate'
)

py3.extension_module('_dop',
  ['_dopmodule.c', 'src/dop.c'],
  dependencies: [np_dep, ccallback_dep],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/integrate'
)

py3.extension_module('_test_multivariate',
  ['tests/_test_multivariate.c'],
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
