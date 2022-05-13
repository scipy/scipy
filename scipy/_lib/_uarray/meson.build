# TODO: add_data_files('license')

# TODO: this used -std=c++14 if available, add c++11 otherwise
#       can we now rely on c++14 unconditionally?

compiler = meson.get_compiler('cpp')

py3.extension_module('_uarray',
  ['_uarray_dispatch.cxx', 'vectorcall.cxx'],
  cpp_args: ['-Wno-terminate', '-Wno-unused-function'],
  dependencies: py3_dep,
  install: true,
  subdir: 'scipy/_lib/_uarray'
)


python_sources = [
  '__init__.py',
  '_backend.py',
]

py3.install_sources(
  python_sources,
  pure: false,
  subdir: 'scipy/_lib/_uarray'
)
