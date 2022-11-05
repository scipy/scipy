py3.install_sources([
    '__init__.py',
  ],
  pure: false,
  subdir: 'scipy/stats/_rcont'
)

_rcont = static_library('_rcont',
  ['rcont.c'],
  include_directories: [inc_np],
  dependencies: [npyrandom_lib, npymath_lib],
)

rcont = py3.extension_module('rcont',
  cython_gen.process('rcont.pyx'),
  c_args: numpy_nodepr_api,
  include_directories: [inc_np],
  link_args: version_link_args,
  link_with: _rcont,
  install: true,
  subdir: 'scipy/stats/_rcont', 
  dependencies: [npyrandom_lib, npymath_lib],
)
