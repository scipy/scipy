_stats_pxd = custom_target('_stats_pxd',
  output: [
    '__init__.py',
    '_stats.pxd',
    '_biasedurn.pxd',
    'unuran.pxd',
  ],
  input: [
    '__init__.py',
    '_stats.pxd',
    '_biasedurn.pxd',
    '_unuran/unuran.pxd',
  ],
  command: [copier, '@INPUT@', '@OUTDIR@']
)

stats_special_cython_gen = generator(cython_cli,
  arguments : ['@INPUT@', '@OUTPUT@'],
  output : '@BASENAME@.c',
  depends : [
    _cython_tree,
    _ufuncs_pxi_pxd_sources,
    _stats_pxd,
    cython_special])

statlib_lib = static_library('statlib_lib',
  [
    'statlib/ansari.f',
    'statlib/spearman.f',
    'statlib/swilk.f'
  ],
  fortran_args: _fflag_Wno_unused_variable
)

statlib_module = custom_target('statlib_module',
  output: ['_statlib-f2pywrappers.f', '_statlibmodule.c'],
  input: 'statlib.pyf',
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

statlib = py3.extension_module('_statlib',
  [statlib_module, fortranobject_c],
  c_args: numpy_nodepr_api,
  fortran_args: fortran_ignore_warnings,
  link_with: statlib_lib,
  dependencies: [py3_dep],
  include_directories: [inc_np, inc_f2py],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/stats'
)

_stats = py3.extension_module('_stats',
  stats_special_cython_gen.process('_stats.pyx'),
  c_args: cython_c_args,
  dependencies: [py3_dep],
  include_directories: [inc_np],
  install: true,
  subdir: 'scipy/stats'
)

mvn_module = custom_target('mvn_module',
  output: ['_mvn-f2pywrappers.f', '_mvnmodule.c'],
  input: 'mvn.pyf',
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

mvn = py3.extension_module('_mvn',
  [mvn_module, 'mvndst.f', fortranobject_c],
  c_args: numpy_nodepr_api,
  fortran_args: fortran_ignore_warnings,
  dependencies: [py3_dep],
  include_directories: [inc_np, inc_f2py],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/stats'
)

_sobol = py3.extension_module('_sobol',
  cython_gen.process('_sobol.pyx'),
  c_args: cython_c_args,
  dependencies: [py3_dep],
  include_directories: [inc_np],
  install: true,
  subdir: 'scipy/stats'
)

py3.install_sources([
    '_sobol_direction_numbers.npz'
  ],
  pure: false,
  subdir: 'scipy/stats'
)

pthread = dependency('threads')

_qmc_cy = py3.extension_module('_qmc_cy',
  cython_gen_cpp.process('_qmc_cy.pyx'),
  cpp_args: cython_cpp_args,
  dependencies: [py3_dep, pthread],
  include_directories: [inc_np],
  install: true,
  subdir: 'scipy/stats'
)

_stats_gen_pyx = custom_target('_stats_gen_pyx',
  output: [
    '_biasedurn.pyx',  # 0 (used in this directory)
    'beta_ufunc.pyx',  # 1 (used in stats/_boost - S_B)
    'binom_ufunc.pyx',  # 2 (S_B)
    'func_defs.pxd',  # 3 (S_B)
    'hypergeom_ufunc.pyx',  # 4 (S_B)
    'nbinom_ufunc.pyx',  # 5 (S_B)
    'templated_pyufunc.pxd',  # 6 (S_B)
    'unuran_wrapper.pyx',  # 7 (used in stats/_unuran)
    'ncf_ufunc.pyx'  # 8 (S_B)
  ],
  input: '_generate_pyx.py',
  command: [py3, '@INPUT@', '-o', '@OUTDIR@'],
  depends: _stats_pxd
)

cython_stats_gen_cpp = generator(cython_cli,
  arguments : ['@INPUT@', '@OUTPUT@', '--cplus'],
  output : '@BASENAME@.cpp',
  depends : [_cython_tree, _stats_gen_pyx])

# Build recipes defined here to get correct output path when used from
# other subdirs.
beta_ufunc_pyx = cython_stats_gen_cpp.process(_stats_gen_pyx[1])
binom_ufunc_pyx = cython_stats_gen_cpp.process(_stats_gen_pyx[2])
hypergeom_ufunc_pyx = cython_stats_gen_cpp.process(_stats_gen_pyx[4])
nbinom_ufunc_pyx = cython_stats_gen_cpp.process(_stats_gen_pyx[5])
ncf_ufunc_pyx = cython_stats_gen_cpp.process(_stats_gen_pyx[8])
# Extra dependency from _lib
unuran_wrap_pyx = lib_cython_gen.process(_stats_gen_pyx[7])

biasedurn = py3.extension_module('_biasedurn',
  [
    cython_gen_cpp.process(_stats_gen_pyx[0]),  # _biasedurn.pyx
    'biasedurn/fnchyppr.cpp',
    'biasedurn/impls.cpp',
    'biasedurn/stoc1.cpp',
    'biasedurn/stoc3.cpp',
    'biasedurn/stocR.h',
    'biasedurn/wnchyppr.cpp',
  ],
  cpp_args: ['-DR_BUILD', cython_cpp_args],
  include_directories: [inc_np],
  # Note that we're not using npymath directly, it is needed to use npyrandom
  dependencies: [py3_dep, npyrandom_lib, npymath_lib],
  install: true,
  subdir: 'scipy/stats'
)

if use_pythran
  _hypotests_pythran = custom_target('_hypotests_pythran',
    output: ['_hypotests_pythran.cpp'],
    input: '_hypotests_pythran.py',
    command: [pythran, '-E', '@INPUT@', '-o', '@OUTDIR@/_hypotests_pythran.cpp']
  )

  _group_columns = py3.extension_module('_hypotests_pythran',
    _hypotests_pythran,
    cpp_args: cpp_args_pythran,
    include_directories: [incdir_pythran, incdir_numpy],
    dependencies: [py3_dep],
    install: true,
    subdir: 'scipy/stats'
  )
endif


py3.install_sources([
    '__init__.py',
    '_axis_nan_policy.py',
    '_binned_statistic.py',
    '_binomtest.py',
    '_common.py',
    '_constants.py',
    '_continuous_distns.py',
    '_crosstab.py',
    '_discrete_distns.py',
    '_distn_infrastructure.py',
    '_distr_params.py',
    '_entropy.py',
    '_fit.py',
    '_generate_pyx.py',
    '_hypotests.py',
    '_kde.py',
    '_ksstats.py',
    '_mannwhitneyu.py',
    '_morestats.py',
    '_mstats_basic.py',
    '_mstats_extras.py',
    '_multivariate.py',
    '_page_trend_test.py',
    '_qmc.py',
    '_relative_risk.py',
    '_resampling.py',
    '_result_classes.py',
    '_rvs_sampling.py',
    '_stats_mstats_common.py',
    '_stats_py.py',
    '_tukeylambda_stats.py',
    '_variation.py',
    '_warnings_errors.py',
    'biasedurn.py',
    'contingency.py',
    'distributions.py',
    'kde.py',
    'morestats.py',
    'mstats.py',
    'mstats_basic.py',
    'mstats_extras.py',
    'mvn.py',
    'qmc.py',
    'sampling.py',
    'statlib.py',
    'stats.py'
  ],
  pure: false,
  subdir: 'scipy/stats'
)

py3.install_sources([
    '_biasedurn.pxd',
    '_stats.pxd'
  ],
  pure: false,
  subdir: 'scipy/stats'
)

py3.install_sources([
    '_qmc_cy.pyi',
    '_sobol.pyi'
  ],
  pure: false,
  subdir: 'scipy/stats'
)

subdir('_boost')
subdir('_levy_stable')
subdir('_unuran')

subdir('tests')
