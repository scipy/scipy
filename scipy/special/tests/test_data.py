import os

import numpy as np
from numpy.testing import *
from scipy.special import (
    arccosh, arcsinh, arctanh, erf, erfc, log1p, expm1, 
    jn, jv, yn, yv, iv, kv, kn, gamma, gammaln, digamma, beta, cbrt,
    ellipe, ellipeinc, ellipk, ellipj, erfinv, erfcinv, exp1, expi, expn,
    zeta,
)
from numpy.testing.noseclasses import KnownFailureTest

DATASETS = [
    os.path.join(os.path.dirname(__file__), "data", "boost.npz")
]

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

def test_all():

    TESTS = [
        Data(arccosh, 'acosh_data_ipp-acosh_data', 0, 1),
        Data(arccosh, 'acosh_data_ipp-acosh_data', 0j, 1, rtol=5e-14),

        Data(arcsinh, 'asinh_data_ipp-asinh_data', 0, 1, rtol=1e-11),
        Data(arcsinh, 'asinh_data_ipp-asinh_data', 0j, 1, rtol=1e-11),

        Data(arctanh, 'atanh_data_ipp-atanh_data', 0, 1, rtol=1e-13),
        Data(arctanh, 'atanh_data_ipp-atanh_data', 0j, 1, rtol=1e-13),

        Data(beta, 'beta_exp_data_ipp-beta_exp_data', (0,1), 2, rtol=1e-13),
        Data(beta, 'beta_exp_data_ipp-beta_exp_data', (0,1), 2, rtol=1e-13),
        Data(beta, 'beta_small_data_ipp-beta_small_data', (0,1), 2),

        Data(cbrt, 'cbrt_data_ipp-cbrt_data', 1, 0),

        Data(digamma, 'digamma_data_ipp-digamma_data', 0, 1),
        Data(digamma, 'digamma_data_ipp-digamma_data', 0j, 1),
        Data(digamma, 'digamma_neg_data_ipp-digamma_neg_data', 0, 1, rtol=1e-13),
        Data(digamma, 'digamma_neg_data_ipp-digamma_neg_data', 0j, 1, rtol=1e-13),
        Data(digamma, 'digamma_root_data_ipp-digamma_root_data', 0, 1, rtol=1e-11),
        Data(digamma, 'digamma_root_data_ipp-digamma_root_data', 0j, 1, rtol=1e-11),
        Data(digamma, 'digamma_small_data_ipp-digamma_small_data', 0, 1),
        Data(digamma, 'digamma_small_data_ipp-digamma_small_data', 0j, 1),

        Data(ellipk_, 'ellint_k_data_ipp-ellint_k_data', 0, 1),
        Data(ellipe_, 'ellint_e_data_ipp-ellint_e_data', 0, 1),
        Data(ellipeinc_, 'ellint_e2_data_ipp-ellint_e2_data', (0,1), 2, rtol=1e-14),

        Data(erf, 'erf_data_ipp-erf_data', 0, 1),
        Data(erf, 'erf_data_ipp-erf_data', 0j, 1, rtol=1e-14),
        Data(erfc, 'erf_data_ipp-erf_data', 0, 2),
        Data(erf, 'erf_large_data_ipp-erf_large_data', 0, 1),
        Data(erf, 'erf_large_data_ipp-erf_large_data', 0j, 1),
        Data(erfc, 'erf_large_data_ipp-erf_large_data', 0, 2),
        Data(erf, 'erf_small_data_ipp-erf_small_data', 0, 1),
        Data(erf, 'erf_small_data_ipp-erf_small_data', 0j, 1),
        Data(erfc, 'erf_small_data_ipp-erf_small_data', 0, 2),

        Data(erfinv, 'erf_inv_data_ipp-erf_inv_data', 0, 1),
        Data(erfcinv, 'erfc_inv_data_ipp-erfc_inv_data', 0, 1),
        #Data(erfcinv, 'erfc_inv_big_data_ipp-erfc_inv_big_data', 0, 1),

        Data(exp1, 'expint_1_data_ipp-expint_1_data', 1, 2),
        Data(exp1, 'expint_1_data_ipp-expint_1_data', 1j, 2, rtol=5e-9),
        Data(expi, 'expinti_data_ipp-expinti_data', 0, 1, rtol=1e-13),
        Data(expi, 'expinti_data_double_ipp-expinti_data_double', 0, 1),

        Data(expn, 'expint_small_data_ipp-expint_small_data', (0,1), 2),
        Data(expn, 'expint_data_ipp-expint_data', (0,1), 2, rtol=1e-14),

        Data(gamma, 'test_gamma_data_ipp-near_0', 0, 1),
        Data(gamma, 'test_gamma_data_ipp-near_1', 0, 1),
        Data(gamma, 'test_gamma_data_ipp-near_2', 0, 1),
        Data(gamma, 'test_gamma_data_ipp-near_m10', 0, 1),
        Data(gamma, 'test_gamma_data_ipp-near_m55', 0, 1),
        Data(gamma, 'test_gamma_data_ipp-near_0', 0j, 1, rtol=2e-9),
        Data(gamma, 'test_gamma_data_ipp-near_1', 0j, 1, rtol=2e-9),
        Data(gamma, 'test_gamma_data_ipp-near_2', 0j, 1, rtol=2e-9),
        Data(gamma, 'test_gamma_data_ipp-near_m10', 0j, 1, rtol=2e-9),
        Data(gamma, 'test_gamma_data_ipp-near_m55', 0j, 1, rtol=2e-9),
        Data(gammaln, 'test_gamma_data_ipp-near_0', 0, 2, rtol=5e-11),
        Data(gammaln, 'test_gamma_data_ipp-near_1', 0, 2, rtol=5e-11),
        Data(gammaln, 'test_gamma_data_ipp-near_2', 0, 2, rtol=2e-10),
        Data(gammaln, 'test_gamma_data_ipp-near_m10', 0, 2, rtol=5e-11),
        Data(gammaln, 'test_gamma_data_ipp-near_m55', 0, 2, rtol=5e-11),

        Data(log1p, 'log1p_expm1_data_ipp-log1p_expm1_data', 0, 1),
        Data(expm1, 'log1p_expm1_data_ipp-log1p_expm1_data', 0, 2),

        Data(iv, 'bessel_i_data_ipp-bessel_i_data', (0,1), 2, rtol=1e-12),
        Data(iv, 'bessel_i_data_ipp-bessel_i_data', (0,1j), 2, rtol=2e-10, atol=1e-306),
        Data(iv, 'bessel_i_int_data_ipp-bessel_i_int_data', (0,1), 2, rtol=1e-9),
        Data(iv, 'bessel_i_int_data_ipp-bessel_i_int_data', (0,1j), 2, rtol=2e-10),

        Data(jn, 'bessel_j_int_data_ipp-bessel_j_int_data', (0,1), 2, rtol=1e-12),
        Data(jn, 'bessel_j_int_data_ipp-bessel_j_int_data', (0,1j), 2, rtol=1e-12),
        Data(jn, 'bessel_j_large_data_ipp-bessel_j_large_data', (0,1), 2, rtol=6e-11),
        Data(jn, 'bessel_j_large_data_ipp-bessel_j_large_data', (0,1j), 2, rtol=6e-11),

        Data(jv, 'bessel_j_int_data_ipp-bessel_j_int_data', (0,1), 2, rtol=1e-12),
        Data(jv, 'bessel_j_int_data_ipp-bessel_j_int_data', (0,1j), 2, rtol=1e-12),
        Data(jv, 'bessel_j_data_ipp-bessel_j_data', (0,1), 2, rtol=1e-12),
        Data(jv, 'bessel_j_data_ipp-bessel_j_data', (0,1j), 2, rtol=1e-12),

        Data(kn, 'bessel_k_int_data_ipp-bessel_k_int_data', (0,1), 2, rtol=1e-12,
             knownfailure="Known bug in Cephes kn implementation"),

        Data(kv, 'bessel_k_int_data_ipp-bessel_k_int_data', (0,1), 2, rtol=1e-12),
        Data(kv, 'bessel_k_int_data_ipp-bessel_k_int_data', (0,1j), 2, rtol=1e-12),
        Data(kv, 'bessel_k_data_ipp-bessel_k_data', (0,1), 2, rtol=1e-12),
        Data(kv, 'bessel_k_data_ipp-bessel_k_data', (0,1j), 2, rtol=1e-12),

        Data(yn, 'bessel_y01_data_ipp-bessel_y01_data', (0,1), 2, rtol=1e-12),
        Data(yn, 'bessel_yn_data_ipp-bessel_yn_data', (0,1), 2, rtol=1e-12),

        Data(yv, 'bessel_yn_data_ipp-bessel_yn_data', (0,1), 2, rtol=1e-12),
        Data(yv, 'bessel_yn_data_ipp-bessel_yn_data', (0,1j), 2, rtol=1e-12),
        Data(yv, 'bessel_yv_data_ipp-bessel_yv_data', (0,1), 2, rtol=1e-12,
             knownfailure="Known bug in Cephes yv implementation"),
        Data(yv, 'bessel_yv_data_ipp-bessel_yv_data', (0,1j), 2, rtol=1e-10),

        Data(zeta_, 'zeta_data_ipp-zeta_data', 0, 1, param_filter=(lambda s: s > 1)),
        Data(zeta_, 'zeta_neg_data_ipp-zeta_neg_data', 0, 1, param_filter=(lambda s: s > 1)),
        Data(zeta_, 'zeta_1_up_data_ipp-zeta_1_up_data', 0, 1, param_filter=(lambda s: s > 1)),
        Data(zeta_, 'zeta_1_below_data_ipp-zeta_1_below_data', 0, 1, param_filter=(lambda s: s > 1)),

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
        # gamma_inv_big_data.txt
        # gamma_inv_data.txt
        # gamma_inv_small_data.txt
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
    test.check(dtype=dtype)


#------------------------------------------------------------------------------

class Data(object):
    """
    Data set for checking a special function.

    Parameters
    ----------
    func : function
        Function to test
    filename : str
        Input file name
    param_columns : int or tuple of ints
        Columns indices in which the parameters to `func` lie.
        Can be imaginary integers to indicate that the parameter
        should be cast to complex.
    result_columns : int or tuple of ints
        Column indices for expected results from `func`.
    rtol : float, optional
        Required relative tolerance. Default is 5*eps.
    atol : float, optional
        Required absolute tolerance. Default is 5*tiny.
    param_filter : function, or tuple of functions/Nones, optional
        Filter functions to exclude some parameter ranges.
        If omitted, no filtering is done.
    knownfailure : str, optional
        Known failure error message to raise when the test is run.
        If omitted, no exception is raised.

    """

    def __init__(self, func, dataname, param_columns, result_columns,
                 rtol=None, atol=None, param_filter=None, knownfailure=None):
        self.func = func
        self.dataname = dataname
        if not hasattr(param_columns, '__len__'):
            param_columns = (param_columns,)
        if not hasattr(result_columns, '__len__'):
            result_columns = (result_columns,)
        self.param_columns = tuple(param_columns)
        self.result_columns = tuple(result_columns)
        self.rtol = rtol
        self.atol = atol
        if not hasattr(param_filter, '__len__'):
            param_filter = (param_filter,)
        self.param_filter = param_filter
        self.knownfailure = knownfailure

    def get_tolerances(self, dtype):
        info = np.finfo(dtype)
        rtol, atol = self.rtol, self.atol
        if rtol is None:
            rtol = 5*info.eps
        if atol is None:
            atol = 5*info.tiny
        return rtol, atol

    @classmethod
    def load_data(cls, dataname):
        if isinstance(dataname, np.ndarray):
            return dataname
        if not hasattr(cls, '_datasets'):
            cls._datasets = {}
            for ds in DATASETS:
                cls._datasets.update(np.load(ds))
        return cls._datasets[dataname].astype(dtype)

    def check(self, data=None, dtype=None):
        """Check the special function against the data."""

        if self.knownfailure:
            raise KnownFailureTest(self.knownfailure)

        if data is None:
            data = Data.load_data(self.dataname)
        if dtype is None:
            dtype = data.dtype
        else:
            data = data.astype(dtype)

        rtol, atol = self.get_tolerances(dtype)

        # Apply given filter functions
        if self.param_filter:
            param_mask = np.ones((data.shape[0],), np.bool_)
            for j, filter in zip(self.param_columns, self.param_filter):
                if filter:
                    param_mask &= filter(data[:,j])
            data = data[param_mask]

        # Pick parameters and results from the correct columns
        params = []
        for j in self.param_columns:
            if np.iscomplexobj(j):
                j = int(j.imag)
                params.append(data[:,j].astype(np.complex))
            else:
                params.append(data[:,j])
        wanted = tuple([data[:,j] for j in self.result_columns])

        # Evaluate
        got = self.func(*params)
        if not isinstance(got, tuple):
            got = (got,)

        # Check the validity of each output returned

        assert len(got) == len(wanted)

        for output_num, (x, y) in enumerate(zip(got, wanted)):
            pinf_x = np.isinf(x) & (x > 0)
            pinf_y = np.isinf(y) & (x > 0)
            minf_x = np.isinf(x) & (x < 0)
            minf_y = np.isinf(y) & (x < 0)
            nan_x = np.isnan(x)
            nan_y = np.isnan(y)

            abs_y = np.absolute(y)
            abs_y[~np.isfinite(abs_y)] = 0
            diff = np.absolute(x - y)
            diff[~np.isfinite(diff)] = 0

            rdiff = diff / np.absolute(y)
            rdiff[~np.isfinite(rdiff)] = 0

            tol_mask = (diff < atol + rtol*abs_y)
            pinf_mask = (pinf_x == pinf_y)
            minf_mask = (minf_x == minf_y)
            nan_mask = (nan_x == nan_y)

            bad_j = ~(tol_mask & pinf_mask & minf_mask & nan_mask)

            if np.any(bad_j):
                # Some bad results: inform what, where, and how bad
                msg = [""]
                msg.append("Max |adiff|: %g" % diff.max())
                msg.append("Max |rdiff|: %g" % rdiff.max())
                msg.append("Bad results for the following points (in output %d):"
                           % output_num)
                for j in np.where(bad_j)[0]:
                    j = int(j)
                    fmt = lambda x: "%30s" % np.array2string(x[j], precision=18)
                    a = "  ".join(map(fmt, params))
                    b = "  ".join(map(fmt, got))
                    c = "  ".join(map(fmt, wanted))
                    d = fmt(rdiff)
                    msg.append("%s => %s != %s  (rdiff %s)" % (a, b, c, d))
                assert False, "\n".join(msg)

    def __repr__(self):
        """Pretty-printing, esp. for Nose output"""
        if np.any(map(np.iscomplexobj, self.param_columns)):
            is_complex = " (complex)"
        else:
            is_complex = ""
        return "<Boost test for %s%s: %s>" % (self.func.__name__, is_complex,
                                              os.path.basename(self.dataname))
