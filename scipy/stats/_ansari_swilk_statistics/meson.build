py3.extension_module('__init__',
  cython_gen.process('__init__.pyx'),
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/stats/_ansari_swilk_statistics'
)
