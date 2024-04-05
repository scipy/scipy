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

mvn_module = custom_target('mvn_module',
  output: ['_mvn-f2pywrappers.f', '_mvnmodule.c'],
  input: 'mvn.pyf',
  command: [generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('_mvn',
  [mvn_module, 'mvndst.f'],
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
  ],
  cpp_args: ['-DR_BUILD', cython_cpp_args],
  link_args: version_link_args,
  # Note that we're not using npymath directly, it is needed to use npyrandom
  dependencies: [np_dep, npyrandom_lib, npymath_lib],
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
    '_generate_pyx_.py',
    '_hypotests.py',
    '_kde.py',
    '_ksstats.py',
    '_mannwhitneyu.py',
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
    '_rvs_sampling.py',
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

_ufuncs_pxi_pxd_sources = [
  fs.copyfile('_cython_stats.pxi'),
  fs.copyfile('_ufuncs_extra_code.pxi'),
  fs.copyfile('_ufuncs_extra_code_common.pxi')
]

# TODO: this installs all generated files, while we want to install only:
#       _ufuncs.pyi, _cython_stats.pxd
#       There's no good way around this currently, because py3.install_sources
#       doesn't accept generated targets. See TODO near the end of
#       scipy/linalg/meson.build for more details.
cython_stats = custom_target('cython_stats',
  output: [
    '_ufuncs.pyx',
    '_ufuncs_defs.h',
    '_ufuncs_cxx.pyx',
    '_ufuncs_cxx.pxd',
    '_ufuncs_cxx_defs.h',
    '_ufuncs.pyi',
    'cython_stats.pyx',
    'cython_stats.pxd'
  ],
  input: ['_generate_pyx_.py', 'functions.json', '_add_newdocs.py'],
  command: [py3, '@INPUT0@', '-o', '@OUTDIR@'],
  install: true,
  install_dir: py3.get_install_dir() / 'scipy/stats',
  install_tag: 'devel',
)

# Only needed to establish a dependency; see comments in linalg/meson.build
# for cython_blas_pxd
cython_stats_pxd = custom_target('_dummy_cython_stats.pxd',
  output: '_dummy_cython_stats.pxd',
  input: cython_stats[7],
  command: [py3, '-c', '"@INPUT@"'],
)

# pyx -> c, pyx -> cpp generators, depending on copied pxi, pxd files.
uf_cython_gen = generator(cython,
  arguments : cython_args,
  output : '@BASENAME@.c',
  depends : [_cython_tree, _ufuncs_pxi_pxd_sources, cython_stats_pxd, _stats_pxd[0]])

uf_cython_gen_cpp = generator(cython,
  arguments : cython_cplus_args,
  output : '@BASENAME@.cpp',
  depends : [_cython_tree, _ufuncs_pxi_pxd_sources, cython_stats_pxd, _stats_pxd[0]])

py3.extension_module('_ufuncs',
  [
    uf_cython_gen.process(cython_stats[0]),  # _ufuncs.pyx
  ],
  c_args: [cython_c_args, Wno_maybe_uninitialized],
  include_directories: ['../_lib', '../_build_utils/src'],
  dependencies: [
    lapack,
    npymath_lib,
    np_dep,
  ],
  install: true,
  subdir: 'scipy/stats'
)

ufuncs_cxx_cpp_args = [
  cython_cpp_args,
  '-DBOOST_MATH_STANDALONE=1',
  '-DCYTHON_EXTERN_C=extern "C"',
]

py3.extension_module('_ufuncs_cxx',
  [
    uf_cython_gen_cpp.process(cython_stats[2]),  # _ufuncs_cxx.pyx
    ],
  cpp_args: ufuncs_cxx_cpp_args,
  include_directories: ['../_lib/boost_math/include', '../_lib',
                        '../_build_utils/src'],
  link_args: version_link_args,
  dependencies: [np_dep],
  install: true,
  subdir: 'scipy/stats'
)

subdir('_levy_stable')
subdir('_unuran')
subdir('_rcont')

subdir('tests')
