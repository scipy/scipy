py3.extension_module('_direct',
  ['_direct/direct_wrap.c',
   '_direct/direct-internal.h',
   '_direct/DIRect.c',
   '_direct/DIRserial.c',
   '_direct/DIRsubrout.c',
   '_directmodule.h',
   '_directmodule.c'
  ],
  include_directories: '../_lib/src',
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/optimize'
)

py3.extension_module('_minpack',
  [
    '__minpack.h',
    '__minpack.c',
  ],
  dependencies: [np_dep, ccallback_dep],
  install: true,
  subdir: 'scipy/optimize'
)

rectangular_lsap = static_library('rectangular_lsap',
  [
    'rectangular_lsap/rectangular_lsap.h',
    'rectangular_lsap/rectangular_lsap.cpp'
  ],
  gnu_symbol_visibility: 'hidden',
)

py3.extension_module('_lsap',
  ['_lsap.c'],
  link_with: rectangular_lsap,
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
  ],
  gnu_symbol_visibility: 'hidden',
)

py3.extension_module('_zeros',
  'zeros.c',
  link_with: rootfind_lib,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/optimize'
)

py3.extension_module('_lbfgsb',
  [
    '__lbfgsb.h',
    '__lbfgsb.c',
  ],
  dependencies: [lapack_dep, blas_dep, np_dep],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/optimize'
)

py3.extension_module('_moduleTNC',
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

py3.extension_module('_cobyla',
  [f2py_gen.process('cobyla/cobyla.pyf'), 'cobyla/cobyla2.f', 'cobyla/trstlp.f'],
  c_args: [Wno_unused_variable],
  fortran_args: fortran_ignore_warnings,
  link_args: version_link_args,
  dependencies: [fortranobject_dep],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/optimize'
)

py3.extension_module('_slsqp',
  [f2py_gen.process('slsqp/slsqp.pyf'), 'slsqp/slsqp_optmz.f'],
  fortran_args: fortran_ignore_warnings,
  link_args: version_link_args,
  dependencies: [fortranobject_dep],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/optimize'
)

py3.extension_module('_pava_pybind',
  ['_pava/pava_pybind.cpp'],
  include_directories: '_pava',
  dependencies: [np_dep, pybind11_dep],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/optimize'
)

if use_pythran
  py3.extension_module('_group_columns',
    pythran_gen.process('_group_columns.py'),
    cpp_args: cpp_args_pythran,
    dependencies: [pythran_dep, np_dep],
    link_args: version_link_args,
    install: true,
    subdir: 'scipy/optimize'
  )
else
  py3.extension_module('_group_columns',
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

py3.extension_module('_bglu_dense',
  opt_gen.process('_bglu_dense.pyx'),
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/optimize'
)

py3.extension_module('_cython_nnls',
  opt_gen.process('_cython_nnls.pyx'),
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


py3.install_sources([
  '__init__.py',
  '_basinhopping.py',
  '_bracket.py',
  '_chandrupatla.py',
  '_cobyla_py.py',
  '_cobyqa_py.py',
  '_constraints.py',
  '_differentiable_functions.py',
  '_differentialevolution.py',
  '_direct_py.py',
  '_dcsrch.py',
  '_dual_annealing.py',
  '_elementwise.py',
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
  'elementwise.py',
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
