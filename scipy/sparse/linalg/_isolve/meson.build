# Note: 3 source files in iterative/ are unused
methods = [
  'iterative/BiCGREVCOM.f.src',
  'iterative/BiCGSTABREVCOM.f.src',
  'iterative/CGREVCOM.f.src',
  'iterative/CGSREVCOM.f.src',
  'iterative/GMRESREVCOM.f.src',
  'iterative/QMRREVCOM.f.src'
]

Util = ['iterative/getbreak.f.src']

methods_f_0 = custom_target('methods_f_0',
  output: ['BiCGREVCOM.f'],
  input: [methods[0]],
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

methods_f_1 = custom_target('methods_f_1',
  output: ['BiCGSTABREVCOM.f'],
  input: [methods[1]],
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

methods_f_2 = custom_target('methods_f_2',
  output: ['CGREVCOM.f'],
  input: [methods[2]],
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

methods_f_3 = custom_target('methods_f_3',
  output: ['CGSREVCOM.f'],
  input: [methods[3]],
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

methods_f_4 = custom_target('methods_f_4',
  output: ['GMRESREVCOM.f'],
  input: [methods[4]],
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

methods_f_5 = custom_target('methods_f_5',
  output: ['QMRREVCOM.f'],
  input: [methods[5]],
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

Util_f_0 = custom_target('Util_f_0',
  output: ['getbreak.f'],
  input: [Util[0]],
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

_iterative_module = custom_target('_iterative_module',
  output: ['_iterative.pyf', '_iterativemodule.c'],
  input: 'iterative/_iterative.pyf.src',
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

_iterative = py3.extension_module('_iterative',
  [
    _iterative_module,
    fortranobject_c,
    g77_abi_wrappers,
    methods_f_0,
    methods_f_1,
    methods_f_2,
    methods_f_3,
    methods_f_4,
    methods_f_5,
    Util_f_0
  ],
  c_args: numpy_nodepr_api,
  fortran_args: fortran_ignore_warnings,
  include_directories: [inc_np, inc_f2py],
  dependencies: [py3_dep, lapack],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/sparse/linalg/_isolve'
)

py3.install_sources([
    '__init__.py',
    '_gcrotmk.py',
    'iterative.py',
    'lgmres.py',
    'lsmr.py',
    'lsqr.py',
    'minres.py',
    'tfqmr.py',
    'utils.py',
  ],
  pure: false,
  subdir: 'scipy/sparse/linalg/_isolve'
)

subdir('tests')
