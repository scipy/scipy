_dummy_init_cyoptimize = [
  fs.copyfile('__init__.py'),
  fs.copyfile('c_zeros.pxd'),
  fs.copyfile('_zeros.pxd'),
]

_zeros_pyx = custom_target('_zeros_pyx',
  output: '_zeros.pyx',
  input: '_zeros.pyx.in',
  command: [
    tempita, '@INPUT@', '-o', '@OUTDIR@',
    '--ignore', _dummy_init_cyoptimize[0]
  ]
)

cy_opt_gen = generator(cython,
  arguments : cython_args,
  output : '@BASENAME@.c',
  depends : [_cython_tree,
    cython_optimize_pxd,
    _dummy_init_optimize,
    _dummy_init_cyoptimize])

py3.extension_module('_zeros',
  cy_opt_gen.process(_zeros_pyx),
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/optimize/cython_optimize'
)

py3.install_sources([
    '__init__.py',
    '_zeros.pxd',
    'c_zeros.pxd'
  ],
  subdir: 'scipy/optimize/cython_optimize'
)
