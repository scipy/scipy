_stats_pxd = [
  fs.copyfile('__init__.py'),
  fs.copyfile('_stats.pxd'),
  fs.copyfile('_biasedurn.pxd'),
]

stats_special_cython_gen = generator(cython,
  arguments : cython_args,
  output : '@BASENAME@.c',
  depends : [
    _cython_tree,
    _ufuncs_pxi_pxd_sources,
    _stats_pxd,
    cython_special_pxd])

py3.extension_module('_stats',
  stats_special_cython_gen.process('_stats.pyx'),
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/stats'
)

py3.extension_module('_ansari_swilk_statistics',
  cython_gen.process('_ansari_swilk_statistics.pyx'),
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/stats'
)

py3.extension_module('_mvn',
  [f2py_gen.process('mvn.pyf'), 'mvndst.f'],
  # Wno-surprising is to suppress a pointless warning with GCC 10-12
  # (see GCC bug 98411: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=98411)
  fortran_args: [fortran_ignore_warnings, _fflag_Wno_surprising],
  dependencies: [fortranobject_dep],
  link_args: version_link_args,
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/stats'
)

py3.extension_module('_sobol',
  cython_gen.process('_sobol.pyx'),
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/stats'
)

py3.install_sources([
    '_sobol_direction_numbers.npz'
  ],
  subdir: 'scipy/stats'
)

py3.extension_module('_qmc_cy',
  cython_gen_cpp.process('_qmc_cy.pyx'),
  cpp_args: cython_cpp_args,
  dependencies: [np_dep, thread_dep],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/stats'
)

py3.extension_module('_biasedurn',
  [
    cython_gen_cpp.process('_biasedurn.pyx'),
    'biasedurn/fnchyppr.cpp',
    'biasedurn/impls.cpp',
    'biasedurn/stoc1.cpp',
    'biasedurn/stoc3.cpp',
    'biasedurn/stocR.h',
    'biasedurn/wnchyppr.cpp',
    'libnpyrandom/distributions.c',
  ],
  cpp_args: ['-DR_BUILD', cython_cpp_args],
  include_directories: ['libnpyrandom'],
  link_args: version_link_args,
  dependencies: [np_dep],
  install: true,
  subdir: 'scipy/stats'
)

if use_pythran
  py3.extension_module('_stats_pythran',
    pythran_gen.process('_stats_pythran.py'),
    cpp_args: cpp_args_pythran,
    dependencies: [pythran_dep, np_dep],
    link_args: version_link_args,
    install: true,
    subdir: 'scipy/stats'
  )
else
  py3.install_sources(
    ['_stats_pythran.py'],
    subdir: 'scipy/stats'
  )
endif


py3.install_sources([
    '__init__.py',
    '_axis_nan_policy.py',
    '_binned_statistic.py',
    '_binomtest.py',
    '_bws_test.py',
    '_censored_data.py',
    '_common.py',
    '_constants.py',
    '_continuous_distns.py',
    '_covariance.py',
    '_crosstab.py',
    '_discrete_distns.py',
    '_distn_infrastructure.py',
    '_distr_params.py',
    '_entropy.py',
    '_fit.py',
    '_hypotests.py',
    '_kde.py',
    '_ksstats.py',
    '_mannwhitneyu.py',
    '_mgc.py',
    '_morestats.py',
    '_mstats_basic.py',
    '_mstats_extras.py',
    '_multicomp.py',
    '_multivariate.py',
    '_odds_ratio.py',
    '_page_trend_test.py',
    '_qmc.py',
    '_qmvnt.py',
    '_relative_risk.py',
    '_resampling.py',
    '_result_classes.py',
    '_sampling.py',
    '_sensitivity_analysis.py',
    '_stats_mstats_common.py',
    '_stats_py.py',
    '_survival.py',
    '_tukeylambda_stats.py',
    '_variation.py',
    '_warnings_errors.py',
    '_wilcoxon.py',
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
    'stats.py'
  ],
  subdir: 'scipy/stats'
)

py3.install_sources([
    '_biasedurn.pxd',
    '_stats.pxd'
  ],
  subdir: 'scipy/stats'
)

py3.install_sources([
    '_qmc_cy.pyi',
    '_sobol.pyi'
  ],
  subdir: 'scipy/stats'
)

subdir('_levy_stable')
subdir('_unuran')
subdir('_rcont')

subdir('tests')
