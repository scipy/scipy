fortran_ignore_warnings = [
  '-Wno-tabs', '-Wno-conversion',
  '-Wno-argument-mismatch', '-Wno-unused-dummy-argument',
  '-Wno-maybe-uninitialized', '-Wno-unused-label',
  '-Wno-unused-variable'
]

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

statlib_lib = static_library('statlib_lib',
  [
    'statlib/ansari.f',
    'statlib/spearman.f',
    'statlib/swilk.f'
  ],
  fortran_args: '-Wno-unused-variable'
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

# Relies on cython_special
_stats_c = custom_target('_stats_c',
  output: '_stats.c',
  input: [
    '_stats.pyx',
    _cython_tree,
    _stats_pxd[0],  # __init__.py
    cython_special[7],  # cython_special.pxd
  ],
  command: [
    cython, '-3', '--fast-fail',
    '@INPUT0@', '-o', '@OUTPUT@', '-I', '@BUILD_ROOT@'
  ]
)

_stats = py3.extension_module('_stats',
  _stats_c,
  c_args: '-Wno-cpp',  # Replace with numpy_nodepr_api after Cython 3.0 is out
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
  ['_sobol.pyx'],
  c_args: numpy_nodepr_api,
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
  '_qmc_cy.pyx',
  cpp_args: '-Wno-cpp',  # TODO: replace with numpy_nodepr_api
  dependencies: [py3_dep, pthread],
  include_directories: [inc_np],
  override_options: ['cython_language=cpp'],
  install: true,
  subdir: 'scipy/stats'
)

_stats_gen_pyx = custom_target('_stats_gen_pyx',
  output: [
    '_biasedurn.pyx',
    'beta_ufunc.pyx',
    'binom_ufunc.pyx',
    'func_defs.pxd',
    'hypergeom_ufunc.pyx',
    'nbinom_ufunc.pyx',
    'templated_pyufunc.pxd',
    'unuran_wrapper.pyx'
  ],
  input: '_generate_pyx.py',
  command: [py3, '@INPUT@', '-o', '@OUTDIR@']
)

biasedurn = py3.extension_module('_biasedurn',
  [
    _stats_gen_pyx[0],
    'biasedurn/fnchyppr.cpp',
    'biasedurn/impls.cpp',
    'biasedurn/stoc1.cpp',
    'biasedurn/stoc3.cpp',
    'biasedurn/stocR.h',
    'biasedurn/wnchyppr.cpp',
    _stats_pxd
  ],
  cpp_args: ['-DR_BUILD', numpy_nodepr_api],
  include_directories: [inc_np],
  dependencies: [py3_dep, npyrandom_lib, npymath_lib, cython_tree],
  override_options: ['cython_language=cpp'],
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
    cpp_args: [
      '-Wno-unused-function', '-Wno-unused-variable',
      '-Wno-deprecated-declarations', '-Wno-int-in-bool-context'
    ] + cpp_args_pythran,
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
    '_bootstrap.py',
    '_common.py',
    '_constants.py',
    '_continuous_distns.py',
    '_crosstab.py',
    '_discrete_distns.py',
    '_distn_infrastructure.py',
    '_distr_params.py',
    '_entropy.py',
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
    '_result_classes.py',
    '_rvs_sampling.py',
    '_stats_mstats_common.py',
    '_stats_py.py',
    '_tukeylambda_stats.py',
    '_variation.py',
    '_wilcoxon_data.py',
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
