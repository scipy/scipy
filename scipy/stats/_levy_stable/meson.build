py3.install_sources([
    '__init__.py',
  ],
  subdir: 'scipy/stats/_levy_stable'
)

_levyst = static_library('_levyst',
  'src/levyst.c',
  gnu_symbol_visibility: 'hidden',
)

py3.extension_module('levyst',
  cython_gen.process('src/levyst.pyx'),
  include_directories: ['src'],
  dependencies: np_dep,
  link_args: version_link_args,
  link_with: _levyst,
  install: true,
  subdir: 'scipy/stats/_levy_stable'
)
