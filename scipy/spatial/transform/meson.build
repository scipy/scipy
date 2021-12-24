rotation = py3.extension_module('_rotation',
  ['_rotation.pyx'],
  c_args: '-Wno-cpp', # TODO: replace with numpy_nodepr_api
  include_directories: [incdir_numpy],
  dependencies: [py3_dep],
  install: true,
  subdir: 'scipy/spatial/transform'
)

py3.install_sources(
  '_rotation.pyi',
  pure: false,
  subdir: 'scipy/spatial/transform'
)

py3.install_sources([
    '__init__.py',
    '_rotation_groups.py',
    '_rotation_spline.py',
    'rotation.py'
  ],
  pure: false,
  subdir: 'scipy/spatial/transform'
)

subdir('tests')
