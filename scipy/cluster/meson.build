py3.extension_module('_hierarchy',
  linalg_cython_gen.process('_hierarchy.pyx'),
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/cluster'
)

py3.extension_module('_optimal_leaf_ordering',
  linalg_cython_gen.process('_optimal_leaf_ordering.pyx'),
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/cluster'
)

py3.extension_module('_vq',
  linalg_cython_gen.process('_vq.pyx'),
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/cluster'
)


py3.install_sources([
    '__init__.py',
    'hierarchy.py',
    'vq.py'
  ],
  subdir: 'scipy/cluster'
)

subdir('tests')
