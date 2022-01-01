# TODO: these all depend on scipy.linalg.cython_blas/lapack!
#       trying to add cython_blas_pxd will result in:
#
# scipy/cluster/meson.build:12: WARNING: Custom target input '['_hierarchy.pyx', <CustomTarget b917806@@_cython_tree@cus: ['/usr/bin/cp', '@INPUT@', '-t', '@OUTDIR@']>, <InternalDependency null: True>]' can't be converted to File object(s).
# ...
# AssertionError: Unknown source type: <InternalDependency null: True>

# _hierarchy
_hierarchy_c = custom_target('_hierarchy_c',
  output: '_hierarchy.c',
  input: [
    '_hierarchy.pyx',
    _cython_tree,
    #cython_blas_pxd, cython_lapack_pxd
  ],
  command: [
    cython, '-3', '--fast-fail',
    '@INPUT0@', '-o', '@OUTPUT@', '-I', '@BUILD_ROOT@'
  ]
)

slib = py3.extension_module('_hierarchy',
  _hierarchy_c,
  c_args: [
   '-Wno-cpp',  # TODO: replace with numpy_nodepr_api
  ],
  include_directories: [inc_np],
  dependencies: [py3_dep],
  install: true,
  subdir: 'scipy/cluster'
)

# _optimal_leaf_ordering
_optimal_leaf_ordering_c = custom_target('_optimal_leaf_ordering',
  output: '_optimal_leaf_ordering.c',
  input: [
    '_optimal_leaf_ordering.pyx',
    _cython_tree,
    # _cython_blas_pxd, cython_lapack_pxd
  ],
  command: [
    cython, '-3', '--fast-fail',
    '@INPUT0@', '-o', '@OUTPUT@', '-I', '@BUILD_ROOT@'
  ]
)

slib = py3.extension_module('_optimal_leaf_ordering',
  _optimal_leaf_ordering_c,
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
