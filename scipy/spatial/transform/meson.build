py3.extension_module('_cython_backend',
  [cython_gen.process('_cython_backend.pyx')],
  c_args: [cython_c_args, Wno_maybe_uninitialized],
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/spatial/transform'
)

#py3.extension_module('_rigid_transform',
#  [cython_gen.process('_rigid_transform.pyx')],
#  c_args: [cython_c_args, Wno_maybe_uninitialized],
#  dependencies: np_dep,
#  link_args: version_link_args,
#  install: true,
#  subdir: 'scipy/spatial/transform'
#)

py3.install_sources([
    '__init__.py',
    '_rotation_groups.py',
    '_rotation_spline.py',
    '_array_api_backend.py',
    'rotation.py'
  ],
  subdir: 'scipy/spatial/transform'
)

subdir('tests')
