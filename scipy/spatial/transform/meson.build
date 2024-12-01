py3.extension_module('_rotation',
  [cython_gen.process('_rotation.pyx')],
  c_args: [cython_c_args, Wno_maybe_uninitialized],
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/spatial/transform'
)

py3.extension_module('_translation',
  [cython_gen.process('_translation.pyx')],
  c_args: [cython_c_args, Wno_maybe_uninitialized],
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/spatial/transform'
)

py3.install_sources(
  '_rotation.pyi',
  subdir: 'scipy/spatial/transform'
)

py3.install_sources(
  '_translation.pyi',
  subdir: 'scipy/spatial/transform'
)

py3.install_sources([
    '__init__.py',
    '_rotation_groups.py',
    '_rotation_spline.py',
    'rotation.py',
    'translation.py'
  ],
  subdir: 'scipy/spatial/transform'
)

subdir('tests')
