_spatial_pxd = custom_target('_spatial_pxd',
  output: ['_qhull.pxd', 'setlist.pxd'],
  input: ['_qhull.pxd', 'setlist.pxd'],
  command: [copier, '@INPUT@', '@OUTDIR@']
)

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
  ['_qhull.pyx', 'qhull_misc.h', 'qhull_misc.c'] + qhull_src,
  c_args: '-Wno-cpp',
  include_directories: [
    incdir_numpy,
    '../_lib',
    '../_build_utils/src',
    'qhull_src/src'
  ],
  dependencies: [py3_dep, lapack],
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

cpp_args = []
compiler = meson.get_compiler('cpp')

if compiler.get_id() == 'msvc'
  cpp_args = cpp_args + ['/EHsc']
elif compiler.has_argument('-fvisibility=hidden')
  cpp_args = cpp_args + ['-fvisibility=hidden']
endif

ckdtree = py3.extension_module('_ckdtree',
  ckdtree_src + ['_ckdtree.pyx'],
  # Cannot use -fnovisibility=hidden for Cython cpp, because it
  # conceals PyInit_ckdtree symbol, and prevents module load.
  cpp_args: ['-Wno-cpp'], # TODO: replace with numpy_nodepr_api
  include_directories: [
    incdir_numpy,
    '../_lib',
    '../_build_utils/src',
    'ckdtree/src'
  ],
  dependencies: [py3_dep],
  override_options: ['cython_language=cpp'],
  install: true,
  subdir: 'scipy/spatial'
)

_distance_wrap = py3.extension_module('_distance_wrap',
  'src/distance_wrap.c',
  c_args: numpy_nodepr_api,
  include_directories: [incdir_numpy, '../_lib'],
  dependencies: [py3_dep, npymath_lib],
  install: true,
  subdir: 'scipy/spatial'
)

_distance_pybind = py3.extension_module('_distance_pybind',
  ['src/distance_pybind.cpp'],
  cpp_args: cpp_args + [numpy_nodepr_api],
  include_directories: [incdir_pybind11, incdir_numpy, 'src/'],
  dependencies: [py3_dep],
  install: true,
  subdir: 'scipy/spatial'
)

_voronoi = py3.extension_module('_voronoi',
  ['_voronoi.pyx'],
  c_args: numpy_nodepr_api,
  include_directories: [incdir_numpy],
  dependencies: [py3_dep],
  install: true,
  subdir: 'scipy/spatial'
)

_hausdorff = py3.extension_module('_hausdorff',
  ['_hausdorff.pyx'],
  c_args: numpy_nodepr_api,
  include_directories: [incdir_numpy],
  dependencies: [py3_dep],
  install: true,
  subdir: 'scipy/spatial'
)

py3.install_sources([
    'qhull_src/COPYING.txt'
  ],
  pure: false,
  subdir: 'scipy/spatial/qhull_src'
)

py3.install_sources([
    '_ckdtree.pyi',
    '_qhull.pyi',
    '_voronoi.pyi',
    'distance.pyi'
  ],
  pure: false,
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
  pure: false,
  subdir: 'scipy/spatial'
)

subdir('tests')
subdir('transform')
