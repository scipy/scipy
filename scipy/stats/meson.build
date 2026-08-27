# cython

_stats_pxd = [
  fs.copyfile('__init__.py'),
  fs.copyfile('_biasedurn/__init__.pxd', '_biasedurn.pxd'),
  fs.copyfile('_stats/__init__.pxd', '_stats.pxd'),
]

stats_special_cython_gen = generator(cython,
  arguments : cython_args,
  output : '@BASENAME@.c',
  depends : [
    _cython_tree,
    _ufuncs_pxi_pxd_sources,
    _stats_pxd,
    cython_special_pxd])


# pythran
    
if use_pythran
  py3.extension_module('_stats_pythran',
    pythran_gen.process('_stats_pythran.py'),
    # `_cpp_Wno_unused_local_typedefs` can be removed when we require pythran >= 0.19.0
    cpp_args: [
      cpp_args_pythran,
      _cpp_Wno_unused_local_typedefs,
      _cpp_Wno_maybe_uninitialized,
    ],
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


# python

py3.install_sources([
    '__init__.py',
    '_axis_nan_policy.py',
    '_binned_statistic.py',
    '_binomtest.py',
    '_bws_test.py',
    '_censored_data.py',
    '_circstats.py',
    '_common.py',
    '_constants.py',
    '_continued_fraction.py',
    '_continuous_distns.py',
    '_correlation.py',
    '_covariance.py',
    '_crosstab.py',
    '_discrete_distns.py',
    '_distn_infrastructure.py',
    '_distribution_infrastructure.py',
    '_distr_params.py',
    '_entropy.py',
    '_finite_differences.py',
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
    '_new_distributions.py',
    '_odds_ratio.py',
    '_page_trend_test.py',
    '_probability_distribution.py',
    '_qmc.py',
    '_qmvnt.py',
    '_quantile.py',
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
    'contingency.py',
    'distributions.py',
    'mstats.py',
    'qmc.py',
    'sampling.py',
  ],
  subdir: 'scipy/stats'
)


subdir('_ansari_swilk_statistics')
subdir('_biasedurn')
subdir('_levy_stable')
subdir('_qmc_cy')
subdir('_qmvnt_cy')
subdir('_rcont')
subdir('_sobol')
subdir('_stats')
subdir('_unuran')

subdir('tests')
