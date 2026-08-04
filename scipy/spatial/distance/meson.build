py3.extension_module('_distance_wrap',
  'src/distance_wrap.c',
  include_directories: '../../_lib',
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/spatial/distance'
)

py3.extension_module('_distance_pybind',
  ['src/distance_pybind.cpp'],
  include_directories: ['src/'],
  dependencies: [np_dep, pybind11_dep],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/spatial/distance'
)

py3.extension_module('_hausdorff',
  [cython_gen.process('_hausdorff.pyx')],
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/spatial/distance'
)

py3.install_sources([
    '__init__.py',
    '_distance.pyi',
    '_distance.py'
  ],
  subdir: 'scipy/spatial/distance'
)

subdir('tests')
