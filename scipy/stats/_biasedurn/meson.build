py3.extension_module('__init__',
  [
    cython_gen_cpp.process('__init__.pyx'),
    'src/fnchyppr.cpp',
    'src/impls.cpp',
    'src/stoc1.cpp',
    'src/stoc3.cpp',
    'src/stocR.h',
    'src/wnchyppr.cpp',
    '../src/libnpyrandom/distributions.c',
  ],
  cpp_args: ['-DR_BUILD', cython_cpp_args],
  include_directories: ['../src/libnpyrandom'],
  link_args: version_link_args,
  dependencies: [np_dep],
  install: true,
  subdir: 'scipy/stats/_biasedurn'
)

py3.install_sources('__init__.pxd',
  subdir: 'scipy/stats/_biasedurn'
)
