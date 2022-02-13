py3.install_sources([
    '__init__.py',
  ],
  pure: false,
  subdir: 'scipy/stats/_levy_stable'
)

_levyst = static_library('_levyst',
  ['c_src/levyst.c', 'c_src/levyst.h'],
)

levyst = py3.extension_module('levyst',
  cython_gen.process('levyst.pyx'),
  c_args: numpy_nodepr_api,
  include_directories: [inc_np],
  dependencies: [py3_dep],
  link_with: _levyst,
  install: true,
  subdir: 'scipy/stats/_levy_stable'
)
