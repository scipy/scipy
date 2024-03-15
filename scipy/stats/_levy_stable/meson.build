py3.install_sources([
    '__init__.py',
  ],
  subdir: 'scipy/stats/_levy_stable'
)

_levyst = static_library('_levyst',
  ['c_src/levyst.c', 'c_src/levyst.h'],
)

py3.extension_module('levyst',
  cython_gen.process('levyst.pyx'),
  dependencies: np_dep,
  link_args: version_link_args,
  link_with: _levyst,
  install: true,
  subdir: 'scipy/stats/_levy_stable'
)
