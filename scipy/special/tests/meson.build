python_sources = [
  '__init__.py',
  'test_basic.py',
  'test_bdtr.py',
  'test_boxcox.py',
  'test_cdflib.py',
  'test_cdft_asymptotic.py',
  'test_cephes_intp_cast.py',
  'test_cosine_distr.py',
  'test_cython_special.py',
  'test_data.py',
  'test_dd.py',
  'test_digamma.py',
  'test_ellip_harm.py',
  'test_erfinv.py',
  'test_exponential_integrals.py',
  'test_extending.py',
  'test_faddeeva.py',
  'test_gamma.py',
  'test_gammainc.py',
  'test_hyp2f1.py',
  'test_hypergeometric.py',
  'test_iv_ratio.py',
  'test_kolmogorov.py',
  'test_lambertw.py',
  'test_legendre.py',
  'test_log_softmax.py',
  'test_loggamma.py',
  'test_logit.py',
  'test_logsumexp.py',
  'test_mpmath.py',
  'test_nan_inputs.py',
  'test_ndtr.py',
  'test_ndtri_exp.py',
  'test_orthogonal.py',
  'test_orthogonal_eval.py',
  'test_owens_t.py',
  'test_pcf.py',
  'test_pdtr.py',
  'test_powm1.py',
  'test_precompute_expn_asy.py',
  'test_precompute_gammainc.py',
  'test_precompute_utils.py',
  'test_round.py',
  'test_sf_error.py',
  'test_sici.py',
  'test_specfun.py',
  'test_spence.py',
  'test_spfun_stats.py',
  'test_sph_harm.py',
  'test_spherical_bessel.py',
  'test_support_alternative_backends.py',
  'test_trig.py',
  'test_ufunc_signatures.py',
  'test_wright_bessel.py',
  'test_wrightomega.py',
  'test_zeta.py',
  'test_boost_ufuncs.py',
  'test_xsf_cuda.py',
]


py3.install_sources(
  python_sources,
  subdir: 'scipy/special/tests',
  install_tag: 'tests'
)

py3.install_sources(
  'data/__init__.py',
  subdir: 'scipy/special/tests/data',
  install_tag: 'tests'
)

py3.install_sources(
  ['_cython_examples/extending.pyx', '_cython_examples/meson.build'],
  subdir: 'scipy/special/tests/_cython_examples',
  install_tag: 'tests'
)
