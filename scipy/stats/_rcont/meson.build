py3.install_sources([
    '__init__.py',
  ],
  pure: false,
  subdir: 'scipy/stats/_rcont'
)

py3.extension_module('rcont',
  [
    '_rcont.c',
    '../libnpyrandom/logfactorial.c',
    '../libnpyrandom/distributions.c',
  ],
  cython_gen.process('rcont.pyx'),
  include_directories: ['../libnpyrandom'],
  dependencies: [np_dep],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/stats/_rcont',
)
