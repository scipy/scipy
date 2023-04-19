py3.extension_module('_uarray',
  ['_uarray_dispatch.cxx', 'vectorcall.cxx'],
  cpp_args: [_cpp_Wno_terminate, _cpp_Wno_unused_function],
  link_args: version_link_args,
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
  subdir: 'scipy/_lib/_uarray'
)
