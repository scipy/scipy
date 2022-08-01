# Needed to trick Cython, it won't do a relative import outside a package
# (see https://github.com/mesonbuild/meson/issues/8961)
_dummy_init_cyoptimize = custom_target('_dummy_init_cyoptimize',
  output: [
    '__init__.py',
    'c_zeros.pxd',
    '_zeros.pxd'
  ],
  input: [
    '__init__.py',
    'c_zeros.pxd',
    '_zeros.pxd'
  ],
  command: [copier, '@INPUT@', '@OUTDIR@']
)

_zeros_pyx = custom_target('_zeros_pyx',
  output: '_zeros.pyx',
  input: '_zeros.pyx.in',
  command: [
    py3, tempita, '@INPUT@', '-o', '@OUTDIR@',
    '--ignore', _dummy_init_cyoptimize
  ]
)

cy_opt_gen = generator(cython_cli,
  arguments : ['@INPUT@', '@OUTPUT@'],
  output : '@BASENAME@.c',
  depends : [_cython_tree,
    _dummy_init_optimize,
    _dummy_init_cyoptimize])

_zeros = py3.extension_module('_zeros',
  cy_opt_gen.process(_zeros_pyx),
  c_args: cython_c_args,
  include_directories: inc_np,
  dependencies: [py3_dep],
  install: true,
  subdir: 'scipy/optimize/cython_optimize'
)

py3.install_sources([
    '__init__.py',
    '_zeros.pxd',
    'c_zeros.pxd'
  ],
  pure: false,
  subdir: 'scipy/optimize/cython_optimize'
)
