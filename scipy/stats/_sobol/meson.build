py3.extension_module('__init__',
  cython_gen.process('__init__.pyx'),
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/stats/_sobol'
)

py3.install_sources([
    '__init__.pyi',
    '_sobol_direction_numbers.npz',
  ],
  subdir: 'scipy/stats/_sobol'
)
