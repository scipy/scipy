import os

import numpy as np
from numpy import arccosh, arcsinh, arctanh
from scipy.special import (
    erf, erfc, log1p, expm1,
    jn, jv, yn, yv, iv, kv, kn, gamma, gammaln, digamma, beta, cbrt,
    ellipe, ellipeinc, ellipk, ellipkm1, ellipj, erfinv, erfcinv, exp1, expi,
    expn, zeta, gammaincinv, lpmv
)

from scipy.special._testutils import FuncData

DATASETS = np.load(os.path.join(os.path.dirname(__file__),
                                "data", "boost.npz"))

def data(func, dataname, *a, **kw):
    kw.setdefault('dataname', dataname)
    return FuncData(func, DATASETS[dataname], *a, **kw)

def ellipk_(k):
    return ellipk(k*k)
def ellipe_(k):
    return ellipe(k*k)
def ellipeinc_(f, k):
    return ellipeinc(f, k*k)
def ellipj_(k):
    return ellipj(k*k)
def zeta_(x):
    return zeta(x, 1.)
def assoc_legendre_p_boost_(nu, mu, x):
    # the boost test data is for integer orders only
    return lpmv(mu, nu.astype(int), x)
def legendre_p_via_assoc_(nu, x):
    return lpmv(0, nu, x)

def test_boost():
    TESTS = [
        data(arccosh, 'acosh_data_ipp-acosh_data', 0, 1, rtol=5e-13),
        data(arccosh, 'acosh_data_ipp-acosh_data', 0j, 1, rtol=5e-14),

        data(arcsinh, 'asinh_data_ipp-asinh_data', 0, 1, rtol=1e-11),
        data(arcsinh, 'asinh_data_ipp-asinh_data', 0j, 1, rtol=1e-11),

        data(arctanh, 'atanh_data_ipp-atanh_data', 0, 1, rtol=1e-11),
        data(arctanh, 'atanh_data_ipp-atanh_data', 0j, 1, rtol=1e-11),

        data(assoc_legendre_p_boost_, 'assoc_legendre_p_ipp-assoc_legendre_p',
             (0,1,2), 3, rtol=1e-11),
        data(legendre_p_via_assoc_, 'legendre_p_ipp-legendre_p',
             (0,1), 2, rtol=1e-11),

        data(beta, 'beta_exp_data_ipp-beta_exp_data', (0,1), 2, rtol=1e-13),
        data(beta, 'beta_exp_data_ipp-beta_exp_data', (0,1), 2, rtol=1e-13),
        data(beta, 'beta_small_data_ipp-beta_small_data', (0,1), 2),

        data(cbrt, 'cbrt_data_ipp-cbrt_data', 1, 0),

        data(digamma, 'digamma_data_ipp-digamma_data', 0, 1),
        data(digamma, 'digamma_data_ipp-digamma_data', 0j, 1),
        data(digamma, 'digamma_neg_data_ipp-digamma_neg_data', 0, 1, rtol=1e-13),
        data(digamma, 'digamma_neg_data_ipp-digamma_neg_data', 0j, 1, rtol=1e-13),
        data(digamma, 'digamma_root_data_ipp-digamma_root_data', 0, 1, rtol=1e-11),
        data(digamma, 'digamma_root_data_ipp-digamma_root_data', 0j, 1, rtol=1e-11),
        data(digamma, 'digamma_small_data_ipp-digamma_small_data', 0, 1),
        data(digamma, 'digamma_small_data_ipp-digamma_small_data', 0j, 1),

        data(ellipk_, 'ellint_k_data_ipp-ellint_k_data', 0, 1),
        data(ellipkm1, '-ellipkm1', 0, 1),
        data(ellipe_, 'ellint_e_data_ipp-ellint_e_data', 0, 1),
        data(ellipeinc_, 'ellint_e2_data_ipp-ellint_e2_data', (0,1), 2, rtol=1e-14),

        data(erf, 'erf_data_ipp-erf_data', 0, 1),
        data(erf, 'erf_data_ipp-erf_data', 0j, 1, rtol=1e-14),
        data(erfc, 'erf_data_ipp-erf_data', 0, 2),
        data(erf, 'erf_large_data_ipp-erf_large_data', 0, 1),
        data(erf, 'erf_large_data_ipp-erf_large_data', 0j, 1),
        data(erfc, 'erf_large_data_ipp-erf_large_data', 0, 2),
        data(erf, 'erf_small_data_ipp-erf_small_data', 0, 1),
        data(erf, 'erf_small_data_ipp-erf_small_data', 0j, 1),
        data(erfc, 'erf_small_data_ipp-erf_small_data', 0, 2),

        data(erfinv, 'erf_inv_data_ipp-erf_inv_data', 0, 1),
        data(erfcinv, 'erfc_inv_data_ipp-erfc_inv_data', 0, 1),
        #data(erfcinv, 'erfc_inv_big_data_ipp-erfc_inv_big_data', 0, 1),

        data(exp1, 'expint_1_data_ipp-expint_1_data', 1, 2),
        data(exp1, 'expint_1_data_ipp-expint_1_data', 1j, 2, rtol=5e-9),
        data(expi, 'expinti_data_ipp-expinti_data', 0, 1, rtol=1e-13),
        data(expi, 'expinti_data_double_ipp-expinti_data_double', 0, 1),

        data(expn, 'expint_small_data_ipp-expint_small_data', (0,1), 2),
        data(expn, 'expint_data_ipp-expint_data', (0,1), 2, rtol=1e-14),

        data(gamma, 'test_gamma_data_ipp-near_0', 0, 1),
        data(gamma, 'test_gamma_data_ipp-near_1', 0, 1),
        data(gamma, 'test_gamma_data_ipp-near_2', 0, 1),
        data(gamma, 'test_gamma_data_ipp-near_m10', 0, 1),
        data(gamma, 'test_gamma_data_ipp-near_m55', 0, 1),
        data(gamma, 'test_gamma_data_ipp-near_0', 0j, 1, rtol=2e-9),
        data(gamma, 'test_gamma_data_ipp-near_1', 0j, 1, rtol=2e-9),
        data(gamma, 'test_gamma_data_ipp-near_2', 0j, 1, rtol=2e-9),
        data(gamma, 'test_gamma_data_ipp-near_m10', 0j, 1, rtol=2e-9),
        data(gamma, 'test_gamma_data_ipp-near_m55', 0j, 1, rtol=2e-9),
        data(gammaln, 'test_gamma_data_ipp-near_0', 0, 2, rtol=5e-11),
        data(gammaln, 'test_gamma_data_ipp-near_1', 0, 2, rtol=5e-11),
        data(gammaln, 'test_gamma_data_ipp-near_2', 0, 2, rtol=2e-10),
        data(gammaln, 'test_gamma_data_ipp-near_m10', 0, 2, rtol=5e-11),
        data(gammaln, 'test_gamma_data_ipp-near_m55', 0, 2, rtol=5e-11),

        data(log1p, 'log1p_expm1_data_ipp-log1p_expm1_data', 0, 1),
        data(expm1, 'log1p_expm1_data_ipp-log1p_expm1_data', 0, 2),

        data(iv, 'bessel_i_data_ipp-bessel_i_data', (0,1), 2, rtol=1e-12),
        data(iv, 'bessel_i_data_ipp-bessel_i_data', (0,1j), 2, rtol=2e-10, atol=1e-306),
        data(iv, 'bessel_i_int_data_ipp-bessel_i_int_data', (0,1), 2, rtol=1e-9),
        data(iv, 'bessel_i_int_data_ipp-bessel_i_int_data', (0,1j), 2, rtol=2e-10),

        data(jn, 'bessel_j_int_data_ipp-bessel_j_int_data', (0,1), 2, rtol=1e-12),
        data(jn, 'bessel_j_int_data_ipp-bessel_j_int_data', (0,1j), 2, rtol=1e-12),
        data(jn, 'bessel_j_large_data_ipp-bessel_j_large_data', (0,1), 2, rtol=6e-11),
        data(jn, 'bessel_j_large_data_ipp-bessel_j_large_data', (0,1j), 2, rtol=6e-11),

        data(jv, 'bessel_j_int_data_ipp-bessel_j_int_data', (0,1), 2, rtol=1e-12),
        data(jv, 'bessel_j_int_data_ipp-bessel_j_int_data', (0,1j), 2, rtol=1e-12),
        data(jv, 'bessel_j_data_ipp-bessel_j_data', (0,1), 2, rtol=1e-12),
        data(jv, 'bessel_j_data_ipp-bessel_j_data', (0,1j), 2, rtol=1e-12),

        data(kn, 'bessel_k_int_data_ipp-bessel_k_int_data', (0,1), 2, rtol=1e-12,
             knownfailure="Known bug in Cephes kn implementation"),

        data(kv, 'bessel_k_int_data_ipp-bessel_k_int_data', (0,1), 2, rtol=1e-12),
        data(kv, 'bessel_k_int_data_ipp-bessel_k_int_data', (0,1j), 2, rtol=1e-12),
        data(kv, 'bessel_k_data_ipp-bessel_k_data', (0,1), 2, rtol=1e-12),
        data(kv, 'bessel_k_data_ipp-bessel_k_data', (0,1j), 2, rtol=1e-12),

        data(yn, 'bessel_y01_data_ipp-bessel_y01_data', (0,1), 2, rtol=1e-12),
        data(yn, 'bessel_yn_data_ipp-bessel_yn_data', (0,1), 2, rtol=1e-12),

        data(yv, 'bessel_yn_data_ipp-bessel_yn_data', (0,1), 2, rtol=1e-12),
        data(yv, 'bessel_yn_data_ipp-bessel_yn_data', (0,1j), 2, rtol=1e-12),
        data(yv, 'bessel_yv_data_ipp-bessel_yv_data', (0,1), 2, rtol=1e-12,
             knownfailure="Known bug in Cephes yv implementation"),
        data(yv, 'bessel_yv_data_ipp-bessel_yv_data', (0,1j), 2, rtol=1e-10),

        data(zeta_, 'zeta_data_ipp-zeta_data', 0, 1, param_filter=(lambda s: s > 1)),
        data(zeta_, 'zeta_neg_data_ipp-zeta_neg_data', 0, 1, param_filter=(lambda s: s > 1)),
        data(zeta_, 'zeta_1_up_data_ipp-zeta_1_up_data', 0, 1, param_filter=(lambda s: s > 1)),
        data(zeta_, 'zeta_1_below_data_ipp-zeta_1_below_data', 0, 1, param_filter=(lambda s: s > 1)),

        data(gammaincinv, 'gamma_inv_data_ipp-gamma_inv_data', (0,1), 2,
             rtol=1e-12),
        data(gammaincinv, 'gamma_inv_big_data_ipp-gamma_inv_big_data',
             (0,1), 2, rtol=1e-11),

        # XXX: the data file needs reformatting...
        #data(gammaincinv, 'gamma_inv_small_data_ipp-gamma_inv_small_data',
        #     (0,1), 2),

        # -- not used yet:
        # assoc_legendre_p.txt
        # binomial_data.txt
        # binomial_large_data.txt
        # binomial_quantile_data.txt
        # ellint_f_data.txt
        # ellint_pi2_data.txt
        # ellint_pi3_data.txt
        # ellint_pi3_large_data.txt
        # ellint_rc_data.txt
        # ellint_rd_data.txt
        # ellint_rf_data.txt
        # ellint_rj_data.txt
        # expinti_data_long.txt
        # factorials.txt
        # gammap1m1_data.txt
        # hermite.txt
        # ibeta_data.txt
        # ibeta_int_data.txt
        # ibeta_inv_data.txt
        # ibeta_inva_data.txt
        # ibeta_large_data.txt
        # ibeta_small_data.txt
        # igamma_big_data.txt
        # igamma_int_data.txt
        # igamma_inva_data.txt
        # igamma_med_data.txt
        # igamma_small_data.txt
        # laguerre2.txt
        # laguerre3.txt
        # legendre_p.txt
        # legendre_p_large.txt
        # ncbeta.txt
        # ncbeta_big.txt
        # nccs.txt
        # near_0.txt
        # near_1.txt
        # near_2.txt
        # near_m10.txt
        # near_m55.txt
        # negative_binomial_quantile_data.txt
        # poisson_quantile_data.txt
        # sph_bessel_data.txt
        # sph_neumann_data.txt
        # spherical_harmonic.txt
        # tgamma_delta_ratio_data.txt
        # tgamma_delta_ratio_int.txt
        # tgamma_delta_ratio_int2.txt
        # tgamma_ratio_data.txt
    ]

    for test in TESTS:
        yield _test_factory, test

def _test_factory(test, dtype=np.double):
    """Boost test"""
    olderr = np.seterr(all='ignore')
    try:
        test.check(dtype=dtype)
    finally:
        np.seterr(**olderr)
