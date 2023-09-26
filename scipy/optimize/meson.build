_direct = py3.extension_module('_direct',
  ['_direct/direct_wrap.c',
   '_direct/direct-internal.h',
   '_direct/DIRect.c',
   '_direct/DIRserial.c',
   '_direct/DIRsubrout.c',
   '_directmodule.h',
   '_directmodule.c'
  ],
  c_args: [numpy_nodepr_api],
  include_directories: '../_lib/src',
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/optimize'
)

minpack_lib = static_library('minpack',
  [
    'minpack/chkder.f',
    'minpack/dogleg.f',
    'minpack/dpmpar.f',
    'minpack/enorm.f',
    'minpack/fdjac1.f',
    'minpack/fdjac2.f',
    'minpack/hybrd.f',
    'minpack/hybrd1.f',
    'minpack/hybrj.f',
    'minpack/hybrj1.f',
    'minpack/lmder.f',
    'minpack/lmder1.f',
    'minpack/lmdif.f',
    'minpack/lmdif1.f',
    'minpack/lmpar.f',
    'minpack/lmstr.f',
    'minpack/lmstr1.f',
    'minpack/qform.f',
    'minpack/qrfac.f',
    'minpack/qrsolv.f',
    'minpack/r1mpyq.f',
    'minpack/r1updt.f',
    'minpack/rwupdt.f'
  ]
)

_minpack = py3.extension_module('_minpack',
  ['minpack.h', '__minpack.h', '_minpackmodule.c'],
  link_with: minpack_lib,
  c_args: numpy_nodepr_api,
  include_directories: '../_lib/src',
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/optimize'
)

rectangular_lsap = static_library('rectangular_lsap',
  [
    'rectangular_lsap/rectangular_lsap.h',
    'rectangular_lsap/rectangular_lsap.cpp'
  ],
)

py3.extension_module('_lsap',
  ['_lsap.c'],
  link_with: rectangular_lsap,
  c_args: numpy_nodepr_api,
  include_directories: '../_lib/src',
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/optimize'
)

rootfind_lib = static_library('rootfind',
  [
    'Zeros/bisect.c',
    'Zeros/brenth.c',
    'Zeros/brentq.c',
    'Zeros/ridder.c',
    'Zeros/zeros.h'
  ]
)

py3.extension_module('_zeros',
  'zeros.c',
  link_with: rootfind_lib,
  c_args: numpy_nodepr_api,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/optimize'
)

lbfgsb_module = custom_target('lbfgsb_module',
  output: ['_lbfgsb-f2pywrappers.f', '_lbfgsbmodule.c'],
  input: 'lbfgsb_src/lbfgsb.pyf',
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

_lbfgsb = py3.extension_module('_lbfgsb',
  [
    'lbfgsb_src/lbfgsb.f',
    'lbfgsb_src/linpack.f',
    'lbfgsb_src/timer.f',
    lbfgsb_module,
  ],
  c_args: numpy_nodepr_api,
  fortran_args: fortran_ignore_warnings,
  link_args: version_link_args,
  dependencies: [lapack, blas, fortranobject_dep],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/optimize'
)

moduleTNC = py3.extension_module('_moduleTNC',
  ['tnc/tnc.h',
    cython_gen.process('tnc/_moduleTNC.pyx'),
    'tnc/tnc.c'],
  c_args: cython_c_args,
  include_directories: 'tnc',
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/optimize'
)

cobyla_module = custom_target('cobyla_module',
  output: ['_cobylamodule.c'],
  input: 'cobyla/cobyla.pyf',
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

_cobyla = py3.extension_module('_cobyla',
  [cobyla_module, 'cobyla/cobyla2.f', 'cobyla/trstlp.f'],
  c_args: [numpy_nodepr_api, Wno_unused_variable],
  fortran_args: fortran_ignore_warnings,
  link_args: version_link_args,
  dependencies: [fortranobject_dep],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/optimize'
)

minpack2_module = custom_target('minpack2_module',
  output: ['_minpack2module.c'],
  input: 'minpack2/minpack2.pyf',
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

minpack2 = py3.extension_module('_minpack2',
  [minpack2_module, 'minpack2/dcsrch.f', 'minpack2/dcstep.f'],
  c_args: numpy_nodepr_api,
  fortran_args: fortran_ignore_warnings,
  link_args: version_link_args,
  dependencies: [fortranobject_dep],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/optimize'
)

slsqp_module = custom_target('slsqp_module',
  output: ['_slsqpmodule.c'],
  input: 'slsqp/slsqp.pyf',
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

_slsqp = py3.extension_module('_slsqp',
  [slsqp_module, 'slsqp/slsqp_optmz.f'],
  c_args: numpy_nodepr_api,
  fortran_args: fortran_ignore_warnings,
  link_args: version_link_args,
  dependencies: [fortranobject_dep],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/optimize'
)

_pava_pybind = py3.extension_module('_pava_pybind',
  ['_pava/pava_pybind.cpp'],
  cpp_args: [numpy_nodepr_api],
  include_directories: '_pava',
  dependencies: [np_dep, pybind11_dep],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/optimize'
)

if use_pythran
  _group_columns_cpp = custom_target('_group_columns',
    output: ['_group_columns.cpp'],
    input: '_group_columns.py',
    command: [pythran, '-E', '@INPUT@', '-o', '@OUTDIR@/_group_columns.cpp']
  )

  _group_columns = py3.extension_module('_group_columns',
    [_group_columns_cpp],
    cpp_args: cpp_args_pythran,
    dependencies: [pythran_dep, np_dep],
    link_args: version_link_args,
    install: true,
    subdir: 'scipy/optimize'
  )
else
  _group_columns = py3.extension_module('_group_columns',
    cython_gen.process('_group_columns.pyx'),
    c_args: cython_c_args,
    dependencies: np_dep,
    link_args: version_link_args,
    install: true,
    subdir: 'scipy/optimize'
  )
endif

_dummy_init_optimize = fs.copyfile('__init__.py')

# Copying this .pxd file is only needed because of a Cython bug, see
# discussion on SciPy PR gh-18810.
cython_optimize_pxd = [
  fs.copyfile('cython_optimize.pxd'),
]

opt_gen = generator(cython,
  arguments : cython_args,
  output : '@BASENAME@.c',
  depends : [_cython_tree, cython_blas_pxd, cython_optimize_pxd, _dummy_init_optimize])

_bglu_dense_c = opt_gen.process('_bglu_dense.pyx')

py3.extension_module('_bglu_dense',
  _bglu_dense_c,
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/optimize'
)

subdir('_lsq')
subdir('_trlib')
subdir('_trustregion_constr')
subdir('cython_optimize')
subdir('_shgo_lib')
subdir('_highs')
subdir('tests')

py3.install_sources(
  ['lbfgsb_src/README'],
  subdir: 'scipy/optimize'
)

py3.install_sources([
  '__init__.py',
  '_basinhopping.py',
  '_chandrupatla.py',
  '_cobyla_py.py',
  '_constraints.py',
  '_differentiable_functions.py',
  '_differentialevolution.py',
  '_direct_py.py',
  '_dcsrch.py',
  '_dual_annealing.py',
  '_hessian_update_strategy.py',
  '_isotonic.py',
  '_lbfgsb_py.py',
  '_linesearch.py',
  '_linprog.py',
  '_linprog_doc.py',
  '_linprog_highs.py',
  '_linprog_ip.py',
  '_linprog_rs.py',
  '_linprog_simplex.py',
  '_linprog_util.py',
  '_milp.py',
  '_minimize.py',
  '_minpack_py.py',
  '_nnls.py',
  '_nonlin.py',
  '_numdiff.py',
  '_optimize.py',
  '_qap.py',
  '_remove_redundancy.py',
  '_root.py',
  '_root_scalar.py',
  '_shgo.py',
  '_slsqp_py.py',
  '_spectral.py',
  '_tnc.py',
  '_trustregion.py',
  '_trustregion_dogleg.py',
  '_trustregion_exact.py',
  '_trustregion_krylov.py',
  '_trustregion_ncg.py',
  '_tstutils.py',
  '_zeros_py.py',
  'cobyla.py',
  'cython_optimize.pxd',
  'lbfgsb.py',
  'linesearch.py',
  'minpack.py',
  'minpack2.py',
  'moduleTNC.py',
  'nonlin.py',
  'optimize.py',
  'slsqp.py',
  'tnc.py',
  'zeros.py'
],
  subdir: 'scipy/optimize'
)
