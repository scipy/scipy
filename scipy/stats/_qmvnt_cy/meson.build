py3.extension_module('__init__',
  cython_gen_cpp.process('__init__.pyx'),
  cpp_args: cython_cpp_args,
  dependencies: [np_dep],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/stats/_qmvnt_cy'
)
