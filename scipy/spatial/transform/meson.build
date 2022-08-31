rotation = py3.extension_module('_rotation',
  [cython_gen.process('_rotation.pyx')],
  c_args: [cython_c_args, Wno_maybe_uninitialized],
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
