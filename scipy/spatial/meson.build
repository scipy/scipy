_spatial_pxd = [
  fs.copyfile('_qhull.pxd'),
  fs.copyfile('setlist.pxd'),
]

# pyx -> c, pyx -> cpp generators, depending on copied pxd files.
# _qhull.pyx depends also on _lib/messagestream
spt_cython_gen = generator(cython,
  arguments : cython_args,
  output : '@BASENAME@.c',
  depends : [_cython_tree, _spatial_pxd, _lib_pxd])

qhull_src = [
  'qhull_src/src/geom2_r.c',
  'qhull_src/src/geom_r.c',
  'qhull_src/src/global_r.c',
  'qhull_src/src/io_r.c',
  'qhull_src/src/libqhull_r.c',
  'qhull_src/src/mem_r.c',
  'qhull_src/src/merge_r.c',
  'qhull_src/src/poly2_r.c',
  'qhull_src/src/poly_r.c',
  'qhull_src/src/qset_r.c',
  'qhull_src/src/random_r.c',
  'qhull_src/src/rboxlib_r.c',
  'qhull_src/src/stat_r.c',
  'qhull_src/src/user_r.c',
  'qhull_src/src/usermem_r.c',
  'qhull_src/src/userprintf_r.c',
  'qhull_src/src/userprintf_rbox_r.c'
]

qhull = py3.extension_module('_qhull',
  [spt_cython_gen.process('_qhull.pyx'),
    'qhull_misc.h', 'qhull_misc.c'] + qhull_src,
  c_args: cython_c_args,
  include_directories: [
    '../_lib',
    '../_build_utils/src',
    'qhull_src/src'
  ],
  link_args: version_link_args,
  dependencies: [lapack, np_dep],
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

ckdtree = py3.extension_module('_ckdtree',
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

_distance_wrap = py3.extension_module('_distance_wrap',
  'src/distance_wrap.c',
  include_directories: '../_lib',
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/spatial'
)

_distance_pybind = py3.extension_module('_distance_pybind',
  ['src/distance_pybind.cpp'],
  include_directories: 'src/',
  dependencies: [np_dep, pybind11_dep],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/spatial'
)

_voronoi = py3.extension_module('_voronoi',
  [cython_gen.process('_voronoi.pyx')],
  c_args: [cython_c_args, Wno_maybe_uninitialized],
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/spatial'
)

_hausdorff = py3.extension_module('_hausdorff',
  [cython_gen.process('_hausdorff.pyx')],
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/spatial'
)

py3.install_sources([
    'qhull_src/COPYING.txt'
  ],
  subdir: 'scipy/spatial/qhull_src'
)

py3.install_sources([
    '_ckdtree.pyi',
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
