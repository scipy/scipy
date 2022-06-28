slib = py3.extension_module('_hierarchy',
  linalg_cython_gen.process('_hierarchy.pyx'),
  c_args: cython_c_args,
  include_directories: [inc_np],
  dependencies: [py3_dep],
  install: true,
  subdir: 'scipy/cluster'
)

slib = py3.extension_module('_optimal_leaf_ordering',
  linalg_cython_gen.process('_optimal_leaf_ordering.pyx'),
  c_args: cython_c_args,
  include_directories: [inc_np],
  dependencies: [py3_dep],
  install: true,
  subdir: 'scipy/cluster'
)

slib = py3.extension_module('_vq',
  linalg_cython_gen.process('_vq.pyx'),
  c_args: cython_c_args,
  include_directories: [inc_np],
  dependencies: [py3_dep],
  install: true,
  subdir: 'scipy/cluster'
)


py3.install_sources([
    '__init__.py',
    'hierarchy.py',
    'vq.py'
  ],
  pure: false,             # Will be installed next to binaries
  subdir: 'scipy/cluster'  # Folder relative to site-packages to install to
)

subdir('tests')
