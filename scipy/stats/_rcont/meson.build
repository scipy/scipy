py3.install_sources([
    '__init__.py',
  ],
  pure: false,
  subdir: 'scipy/stats/_rcont'
)

py3.extension_module('rcont',
  ['_rcont.c', 'logfactorial.c'],
  cython_gen.process('rcont.pyx'),
  dependencies: [np_dep, npyrandom_lib, npymath_lib],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/stats/_rcont',
)
