rotation = py3.extension_module('rotation',
  ['rotation.pyx'],
  c_args: '-Wno-cpp', # TODO: replace with numpy_nodepr_api
  include_directories: [incdir_numpy],
  dependencies: [py3_dep],
  install: true,
  subdir: 'scipy/spatial/transform'
)

py3.install_sources(
  'rotation.pyi',
  pure: false,
  subdir: 'scipy/spatial/transform'
)

py3.install_sources([
    '__init__.py',
    '_rotation_groups.py',
    '_rotation_spline.py'
  ],
  pure: false,
  subdir: 'scipy/spatial/transform'
)

subdir('tests')
