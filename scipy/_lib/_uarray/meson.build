py3.extension_module('_uarray',
  ['_uarray_dispatch.cxx', 'vectorcall.cxx'],
  dependencies: py3_dep,
  cpp_args: [_cpp_Wno_terminate, _cpp_Wno_unused_function],
  install: true,
  subdir: 'scipy/_lib/_uarray'
)


python_sources = [
  '__init__.py',
  '_backend.py',
  'LICENSE'
]

py3.install_sources(
  python_sources,
  pure: false,
  subdir: 'scipy/_lib/_uarray'
)
