_unuran_pxd = [
  fs.copyfile('__init__.py'),
  fs.copyfile('unuran.pxd'),
]

# Cython pyx -> c generator with _lib_pxd dependency, because the unuran
# code depends on _lib/_ccallback_c
unuran_cython_gen = generator(cython,
  arguments : cython_args,
  output : '@BASENAME@.c',
  depends : [_cython_tree, _lib_pxd, _stats_pxd, _unuran_pxd]
)

unuran_sources = unuran.get_variable('unuran_sources')
unuran_include_dirs = unuran.get_variable('unuran_include_dirs')
unuran_defines = unuran.get_variable('unuran_defines')

py3.extension_module('unuran_wrapper',
  [
    unuran_cython_gen.process('unuran_wrapper.pyx'),
    unuran_sources
  ],
  c_args: [unuran_defines, cython_c_args],
  include_directories: unuran_include_dirs,
  dependencies: [np_dep, ccallback_dep],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/stats/_unuran'
)

py3.install_sources([
    '__init__.py',
    'unuran_wrapper.pyi',
  ],
  subdir: 'scipy/stats/_unuran'
)
