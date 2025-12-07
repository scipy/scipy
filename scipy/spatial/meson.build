_spatial_pxd = [
  fs.copyfile('_qhull.pxd'),
  fs.copyfile('setlist.pxd'),
]

# pyx -> c, pyx -> cpp generators, depending on copied pxd files.
# _qhull.pyx depends also on _lib/messagestream
spt_cython_gen = generator(cython,
  arguments : cython_args,
  output : '@BASENAME@.c',
  depends : [_cython_tree, _spatial_pxd, _lib_pxd, cython_lapack_pxd])

py3.extension_module('_qhull',
  [spt_cython_gen.process('_qhull.pyx'),
    'qhull_misc.h', 'qhull_misc.c'],
  c_args: cython_c_args,
  include_directories: [
    '../_lib',
    '../_build_utils/src',
  ],
  link_args: version_link_args,
  dependencies: [np_dep, qhull_r_dep],
  install: true,
  subdir: 'scipy/spatial'
)

ckdtree_src = [
  'ckdtree/src/build.cxx',
  'ckdtree/src/count_neighbors.cxx',
  'ckdtree/src/query.cxx',
  'ckdtree/src/query_ball_point.cxx',
  'ckdtree/src/query_ball_tree.cxx',
  'ckdtree/src/query_pairs.cxx',
  'ckdtree/src/sparse_distances.cxx'
]

py3.extension_module('_ckdtree',
  ckdtree_src + [cython_gen_cpp.process('_ckdtree.pyx')],
  cpp_args: cython_cpp_args,
  include_directories: [
    '../_lib',
    '../_build_utils/src',
    'ckdtree/src'
  ],
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/spatial'
)

py3.extension_module('_distance_wrap',
  'src/distance_wrap.c',
  include_directories: '../_lib',
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/spatial'
)

py3.extension_module('_distance_pybind',
  ['src/distance_pybind.cpp'],
  include_directories: ['src/'],
  dependencies: [np_dep, pybind11_dep],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/spatial'
)

py3.extension_module('_voronoi',
  [cython_gen.process('_voronoi.pyx')],
  c_args: [cython_c_args, Wno_maybe_uninitialized],
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/spatial'
)

py3.extension_module('_hausdorff',
  [cython_gen.process('_hausdorff.pyx')],
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/spatial'
)

py3.install_sources([
    '_qhull.pyi',
    '_voronoi.pyi',
    'distance.pyi'
  ],
  subdir: 'scipy/spatial'
)

py3.install_sources([
    '__init__.py',
    '_geometric_slerp.py',
    '_kdtree.py',
    '_plotutils.py',
    '_procrustes.py',
    '_spherical_voronoi.py',
    'ckdtree.py',
    'distance.py',
    'kdtree.py',
    'qhull.py'
  ],
  subdir: 'scipy/spatial'
)

subdir('tests')
subdir('transform')
