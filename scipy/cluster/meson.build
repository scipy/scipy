slib = py3.extension_module('_hierarchy',
  linalg_cython_gen.process('_hierarchy.pyx'),
  c_args: [
   '-Wno-cpp',  # TODO: replace with numpy_nodepr_api
  ],
  include_directories: [inc_np],
  dependencies: [py3_dep],
  install: true,
  subdir: 'scipy/cluster'
)

slib = py3.extension_module('_optimal_leaf_ordering',
  linalg_cython_gen.process('_optimal_leaf_ordering.pyx',
  c_args: [
   '-Wno-cpp',  # TODO: replace with numpy_nodepr_api
  ],
  include_directories: [inc_np],
  dependencies: [py3_dep],
  install: true,
  subdir: 'scipy/cluster'
)

# _vq
_vq_c = custom_target('_vq',
  output: '_vq.c',
  input: [
    '_vq.pyx',
    _cython_tree,
    # _cython_blas_pxd, cython_lapack_pxd
  ],
  command: [
    cython, '-3', '--fast-fail',
    '@INPUT0@', '-o', '@OUTPUT@', '-I', '@BUILD_ROOT@'
  ]
)

slib = py3.extension_module('_vq',
  _vq_c,
  c_args: [
   '-Wno-cpp'  # TODO: replace with numpy_nodepr_api
  ],
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
