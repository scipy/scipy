py3.extension_module('__init__',
  [
    'src/_rcont.c',
    '../src/libnpyrandom/logfactorial.c',
    '../src/libnpyrandom/distributions.c',
  ],
  cython_gen.process('__init__.pyx'),
  include_directories: ['../src/libnpyrandom'],
  dependencies: [np_dep],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/stats/_rcont',
)
