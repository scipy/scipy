# this program corresponds to special.py

### Means test is not done yet
# E   Means test is giving error (E)
# F   Means test is failing (F)
# EF  Means test is giving error and Failing
#!   Means test is segfaulting
# 8   Means test runs forever

###  test_besselpoly
###  test_mathieu_a
###  test_mathieu_even_coef
###  test_mathieu_odd_coef
###  test_modfresnelp
###  test_modfresnelm
#    test_pbdv_seq
###  test_pbvv_seq
###  test_sph_harm
#    test_sph_in
#    test_sph_jn
#    test_sph_kn

from __future__ import division, print_function, absolute_import

import itertools
import warnings

import numpy as np
from numpy import array, isnan, r_, arange, finfo, pi, sin, cos, tan, exp, \
        log, zeros, sqrt, asarray, inf, nan_to_num, real, arctan, float_

from numpy.testing import assert_equal, assert_almost_equal, \
        assert_array_equal, assert_array_almost_equal, assert_approx_equal, \
        assert_, rand, dec, TestCase, run_module_suite, assert_allclose, \
        assert_raises, assert_array_almost_equal_nulp

from scipy import special
import scipy.special._ufuncs as cephes
from scipy.special import ellipk

from scipy.special._testutils import assert_tol_equal, with_special_errors, \
     assert_func_equal


class TestCephes(TestCase):
    def test_airy(self):
        cephes.airy(0)

    def test_airye(self):
        cephes.airye(0)

    def test_binom(self):
        n = np.array([0.264, 4, 5.2, 17])
        k = np.array([2, 0.4, 7, 3.3])
        nk = np.array(np.broadcast_arrays(n[:,None], k[None,:])
                      ).reshape(2, -1).T
        rknown = np.array([[-0.097152, 0.9263051596159367, 0.01858423645695389,
            -0.007581020651518199],[6, 2.0214389119675666, 0, 2.9827344527963846],
            [10.92, 2.22993515861399, -0.00585728, 10.468891352063146],
            [136, 3.5252179590758828, 19448, 1024.5526916174495]])
        assert_func_equal(cephes.binom, rknown.ravel(), nk, rtol=1e-13)

        # Test branches in implementation
        np.random.seed(1234)
        n = np.r_[np.arange(-7, 30), 1000*np.random.rand(30) - 500]
        k = np.arange(0, 102)
        nk = np.array(np.broadcast_arrays(n[:,None], k[None,:])
                      ).reshape(2, -1).T

        assert_func_equal(cephes.binom,
                          cephes.binom(nk[:,0], nk[:,1] * (1 + 1e-15)),
                          nk,
                          atol=1e-10, rtol=1e-10)

    def test_binom_2(self):
        # Test branches in implementation
        np.random.seed(1234)
        n = np.r_[np.logspace(1, 300, 20)]
        k = np.arange(0, 102)
        nk = np.array(np.broadcast_arrays(n[:,None], k[None,:])
                      ).reshape(2, -1).T

        assert_func_equal(cephes.binom,
                          cephes.binom(nk[:,0], nk[:,1] * (1 + 1e-15)),
                          nk,
                          atol=1e-10, rtol=1e-10)

    def test_binom_exact(self):
        @np.vectorize
        def binom_int(n, k):
            n = int(n)
            k = int(k)
            num = int(1)
            den = int(1)
            for i in range(1, k+1):
                num *= i + n - k
                den *= i
            return float(num/den)

        np.random.seed(1234)
        n = np.arange(1, 15)
        k = np.arange(0, 15)
        nk = np.array(np.broadcast_arrays(n[:,None], k[None,:])
                      ).reshape(2, -1).T
        nk = nk[nk[:,0] >= nk[:,1]]
        assert_func_equal(cephes.binom,
                          binom_int(nk[:,0], nk[:,1]),
                          nk,
                          atol=0, rtol=0)

    def test_bdtr(self):
        assert_equal(cephes.bdtr(1,1,0.5),1.0)

    def test_bdtri(self):
        assert_equal(cephes.bdtri(1,3,0.5),0.5)

    def test_bdtrc(self):
        assert_equal(cephes.bdtrc(1,3,0.5),0.5)

    def test_bdtrin(self):
        assert_equal(cephes.bdtrin(1,0,1),5.0)

    def test_bdtrik(self):
        cephes.bdtrik(1,3,0.5)

    def test_bei(self):
        assert_equal(cephes.bei(0),0.0)

    def test_beip(self):
        assert_equal(cephes.beip(0),0.0)

    def test_ber(self):
        assert_equal(cephes.ber(0),1.0)

    def test_berp(self):
        assert_equal(cephes.berp(0),0.0)

    def test_besselpoly(self):
        assert_equal(cephes.besselpoly(0,0,0),1.0)

    def test_beta(self):
        assert_equal(cephes.beta(1,1),1.0)
        assert_allclose(cephes.beta(-100.3, 1e-200), cephes.gamma(1e-200))
        assert_allclose(cephes.beta(0.0342, 171), 24.070498359873497, rtol=1e-14, atol=0)

    def test_betainc(self):
        assert_equal(cephes.betainc(1,1,1),1.0)
        assert_allclose(cephes.betainc(0.0342, 171, 1e-10), 0.55269916901806648)

    def test_betaln(self):
        assert_equal(cephes.betaln(1,1),0.0)
        assert_allclose(cephes.betaln(-100.3, 1e-200), cephes.gammaln(1e-200))
        assert_allclose(cephes.betaln(0.0342, 170), 3.1811881124242447, rtol=1e-14, atol=0)

    def test_betaincinv(self):
        assert_equal(cephes.betaincinv(1,1,1),1.0)
        assert_allclose(cephes.betaincinv(0.0342, 171, 0.25), 8.4231316935498957e-21, rtol=1e-12, atol=0)

    def test_beta_inf(self):
        assert_(np.isinf(special.beta(-1, 2)))

    def test_btdtr(self):
        assert_equal(cephes.btdtr(1,1,1),1.0)

    def test_btdtri(self):
        assert_equal(cephes.btdtri(1,1,1),1.0)

    def test_btdtria(self):
        assert_equal(cephes.btdtria(1,1,1),5.0)

    def test_btdtrib(self):
        assert_equal(cephes.btdtrib(1,1,1),5.0)

    def test_cbrt(self):
        assert_approx_equal(cephes.cbrt(1),1.0)

    def test_chdtr(self):
        assert_equal(cephes.chdtr(1,0),0.0)

    def test_chdtrc(self):
        assert_equal(cephes.chdtrc(1,0),1.0)

    def test_chdtri(self):
        assert_equal(cephes.chdtri(1,1),0.0)

    def test_chdtriv(self):
        assert_equal(cephes.chdtriv(0,0),5.0)

    def test_chndtr(self):
        assert_equal(cephes.chndtr(0,1,0),0.0)
        p = cephes.chndtr(np.linspace(20, 25, 5), 2, 1.07458615e+02)
        assert_allclose(p, [1.21805009e-09, 2.81979982e-09, 6.25652736e-09,
                            1.33520017e-08, 2.74909967e-08],
                        rtol=1e-6, atol=0)
        assert_almost_equal(cephes.chndtr(np.inf, np.inf, 0), 2.0)
        assert_almost_equal(cephes.chndtr(2, 1, np.inf), 0.0)
        assert_(np.isnan(cephes.chndtr(np.nan, 1, 2)))
        assert_(np.isnan(cephes.chndtr(5, np.nan, 2)))
        assert_(np.isnan(cephes.chndtr(5, 1, np.nan)))

    def test_chndtridf(self):
        assert_equal(cephes.chndtridf(0,0,1),5.0)

    def test_chndtrinc(self):
        assert_equal(cephes.chndtrinc(0,1,0),5.0)

    def test_chndtrix(self):
        assert_equal(cephes.chndtrix(0,1,0),0.0)

    def test_cosdg(self):
        assert_equal(cephes.cosdg(0),1.0)

    def test_cosm1(self):
        assert_equal(cephes.cosm1(0),0.0)

    def test_cotdg(self):
        assert_almost_equal(cephes.cotdg(45),1.0)

    def test_dawsn(self):
        assert_equal(cephes.dawsn(0),0.0)
        assert_allclose(cephes.dawsn(1.23), 0.50053727749081767)

    def test_ellipe(self):
        assert_equal(cephes.ellipe(1),1.0)

    def test_ellipeinc(self):
        assert_equal(cephes.ellipeinc(0,1),0.0)

    def test_ellipj(self):
        cephes.ellipj(0,1)

    def test_ellipk(self):
        assert_allclose(ellipk(0), pi/2)

    def test_ellipkinc(self):
        assert_equal(cephes.ellipkinc(0,0),0.0)

    def test_erf(self):
        assert_equal(cephes.erf(0),0.0)

    def test_erfc(self):
        assert_equal(cephes.erfc(0),1.0)

    def test_exp1(self):
        cephes.exp1(1)

    def test_expi(self):
        cephes.expi(1)

    def test_expn(self):
        cephes.expn(1,1)

    def test_exp1_reg(self):
        # Regression for #834
        a = cephes.exp1(-complex(19.9999990))
        b = cephes.exp1(-complex(19.9999991))
        assert_array_almost_equal(a.imag, b.imag)

    def test_exp10(self):
        assert_approx_equal(cephes.exp10(2),100.0)

    def test_exp2(self):
        assert_equal(cephes.exp2(2),4.0)

    def test_expm1(self):
        assert_equal(cephes.expm1(0),0.0)

    def test_fdtr(self):
        assert_equal(cephes.fdtr(1,1,0),0.0)

    def test_fdtrc(self):
        assert_equal(cephes.fdtrc(1,1,0),1.0)

    def test_fdtri(self):
        # cephes.fdtri(1,1,0.5)  #BUG: gives NaN, should be 1
        assert_allclose(cephes.fdtri(1, 1, [0.499, 0.501]),
                        array([0.9937365, 1.00630298]), rtol=1e-6)

    def test_fdtridfd(self):
        assert_equal(cephes.fdtridfd(1,0,0),5.0)

    def test_fresnel(self):
        assert_equal(cephes.fresnel(0),(0.0,0.0))

    def test_gamma(self):
        assert_equal(cephes.gamma(5),24.0)

    def test_gammainc(self):
        assert_equal(cephes.gammainc(5,0),0.0)

    def test_gammaincc(self):
        assert_equal(cephes.gammaincc(5,0),1.0)

    def test_gammainccinv(self):
        assert_equal(cephes.gammainccinv(5,1),0.0)

    def test_gammaln(self):
        cephes.gammaln(10)

    def test_gammasgn(self):
        vals = np.array([-4, -3.5, -2.3, 1, 4.2], np.float64)
        assert_array_equal(cephes.gammasgn(vals), np.sign(cephes.rgamma(vals)))

    def test_gdtr(self):
        assert_equal(cephes.gdtr(1,1,0),0.0)

    def test_gdtrc(self):
        assert_equal(cephes.gdtrc(1,1,0),1.0)

    def test_gdtria(self):
        assert_equal(cephes.gdtria(0,1,1),0.0)

    def test_gdtrib(self):
        cephes.gdtrib(1,0,1)
        # assert_equal(cephes.gdtrib(1,0,1),5.0)

    def test_gdtrix(self):
        cephes.gdtrix(1,1,.1)

    def test_hankel1(self):
        cephes.hankel1(1,1)

    def test_hankel1e(self):
        cephes.hankel1e(1,1)

    def test_hankel2(self):
        cephes.hankel2(1,1)

    def test_hankel2e(self):
        cephes.hankel2e(1,1)

    def test_hyp1f1(self):
        assert_approx_equal(cephes.hyp1f1(1,1,1), exp(1.0))
        assert_approx_equal(cephes.hyp1f1(3,4,-6), 0.026056422099537251095)
        cephes.hyp1f1(1,1,1)

    def test_hyp1f2(self):
        cephes.hyp1f2(1,1,1,1)

    def test_hyp2f0(self):
        cephes.hyp2f0(1,1,1,1)

    def test_hyp2f1(self):
        assert_equal(cephes.hyp2f1(1,1,1,0),1.0)

    def test_hyp3f0(self):
        assert_equal(cephes.hyp3f0(1,1,1,0),(1.0,0.0))

    def test_hyperu(self):
        assert_equal(cephes.hyperu(0,1,1),1.0)

    def test_i0(self):
        assert_equal(cephes.i0(0),1.0)

    def test_i0e(self):
        assert_equal(cephes.i0e(0),1.0)

    def test_i1(self):
        assert_equal(cephes.i1(0),0.0)

    def test_i1e(self):
        assert_equal(cephes.i1e(0),0.0)

    def test_it2i0k0(self):
        cephes.it2i0k0(1)

    def test_it2j0y0(self):
        cephes.it2j0y0(1)

    def test_it2struve0(self):
        cephes.it2struve0(1)

    def test_itairy(self):
        cephes.itairy(1)

    def test_iti0k0(self):
        assert_equal(cephes.iti0k0(0),(0.0,0.0))

    def test_itj0y0(self):
        assert_equal(cephes.itj0y0(0),(0.0,0.0))

    def test_itmodstruve0(self):
        assert_equal(cephes.itmodstruve0(0),0.0)

    def test_itstruve0(self):
        assert_equal(cephes.itstruve0(0),0.0)

    def test_iv(self):
        assert_equal(cephes.iv(1,0),0.0)

    def _check_ive(self):
        assert_equal(cephes.ive(1,0),0.0)

    def test_j0(self):
        assert_equal(cephes.j0(0),1.0)

    def test_j1(self):
        assert_equal(cephes.j1(0),0.0)

    def test_jn(self):
        assert_equal(cephes.jn(0,0),1.0)

    def test_jv(self):
        assert_equal(cephes.jv(0,0),1.0)

    def _check_jve(self):
        assert_equal(cephes.jve(0,0),1.0)

    def test_k0(self):
        cephes.k0(2)

    def test_k0e(self):
        cephes.k0e(2)

    def test_k1(self):
        cephes.k1(2)

    def test_k1e(self):
        cephes.k1e(2)

    def test_kei(self):
        cephes.kei(2)

    def test_keip(self):
        assert_equal(cephes.keip(0),0.0)

    def test_ker(self):
        cephes.ker(2)

    def test_kerp(self):
        cephes.kerp(2)

    def _check_kelvin(self):
        cephes.kelvin(2)

    def test_kn(self):
        cephes.kn(1,1)

    def test_kolmogi(self):
        assert_equal(cephes.kolmogi(1),0.0)
        assert_(np.isnan(cephes.kolmogi(np.nan)))

    def test_kolmogorov(self):
        assert_equal(cephes.kolmogorov(0),1.0)

    def _check_kv(self):
        cephes.kv(1,1)

    def _check_kve(self):
        cephes.kve(1,1)

    def test_log1p(self):
        assert_equal(cephes.log1p(0),0.0)

    def test_lpmv(self):
        assert_equal(cephes.lpmv(0,0,1),1.0)

    def test_mathieu_a(self):
        assert_equal(cephes.mathieu_a(1,0),1.0)

    def test_mathieu_b(self):
        assert_equal(cephes.mathieu_b(1,0),1.0)

    def test_mathieu_cem(self):
        assert_equal(cephes.mathieu_cem(1,0,0),(1.0,0.0))

        # Test AMS 20.2.27
        @np.vectorize
        def ce_smallq(m, q, z):
            z *= np.pi/180
            if m == 0:
                return 2**(-0.5) * (1 - .5*q*cos(2*z))  # + O(q^2)
            elif m == 1:
                return cos(z) - q/8 * cos(3*z)  # + O(q^2)
            elif m == 2:
                return cos(2*z) - q*(cos(4*z)/12 - 1/4)  # + O(q^2)
            else:
                return cos(m*z) - q*(cos((m+2)*z)/(4*(m+1)) - cos((m-2)*z)/(4*(m-1)))  # + O(q^2)
        m = np.arange(0, 100)
        q = np.r_[0, np.logspace(-30, -9, 10)]
        assert_allclose(cephes.mathieu_cem(m[:,None], q[None,:], 0.123)[0],
                        ce_smallq(m[:,None], q[None,:], 0.123),
                        rtol=1e-14, atol=0)

    def test_mathieu_sem(self):
        assert_equal(cephes.mathieu_sem(1,0,0),(0.0,1.0))

        # Test AMS 20.2.27
        @np.vectorize
        def se_smallq(m, q, z):
            z *= np.pi/180
            if m == 1:
                return sin(z) - q/8 * sin(3*z)  # + O(q^2)
            elif m == 2:
                return sin(2*z) - q*sin(4*z)/12  # + O(q^2)
            else:
                return sin(m*z) - q*(sin((m+2)*z)/(4*(m+1)) - sin((m-2)*z)/(4*(m-1)))  # + O(q^2)
        m = np.arange(1, 100)
        q = np.r_[0, np.logspace(-30, -9, 10)]
        assert_allclose(cephes.mathieu_sem(m[:,None], q[None,:], 0.123)[0],
                        se_smallq(m[:,None], q[None,:], 0.123),
                        rtol=1e-14, atol=0)

    def test_mathieu_modcem1(self):
        assert_equal(cephes.mathieu_modcem1(1,0,0),(0.0,0.0))

    def test_mathieu_modcem2(self):
        cephes.mathieu_modcem2(1,1,1)

        # Test reflection relation AMS 20.6.19
        m = np.arange(0, 4)[:,None,None]
        q = np.r_[np.logspace(-2, 2, 10)][None,:,None]
        z = np.linspace(0, 1, 7)[None,None,:]

        y1 = cephes.mathieu_modcem2(m, q, -z)[0]

        fr = -cephes.mathieu_modcem2(m, q, 0)[0] / cephes.mathieu_modcem1(m, q, 0)[0]
        y2 = -cephes.mathieu_modcem2(m, q, z)[0] - 2*fr*cephes.mathieu_modcem1(m, q, z)[0]

        assert_allclose(y1, y2, rtol=1e-10)

    def test_mathieu_modsem1(self):
        assert_equal(cephes.mathieu_modsem1(1,0,0),(0.0,0.0))

    def test_mathieu_modsem2(self):
        cephes.mathieu_modsem2(1,1,1)

        # Test reflection relation AMS 20.6.20
        m = np.arange(1, 4)[:,None,None]
        q = np.r_[np.logspace(-2, 2, 10)][None,:,None]
        z = np.linspace(0, 1, 7)[None,None,:]

        y1 = cephes.mathieu_modsem2(m, q, -z)[0]
        fr = cephes.mathieu_modsem2(m, q, 0)[1] / cephes.mathieu_modsem1(m, q, 0)[1]
        y2 = cephes.mathieu_modsem2(m, q, z)[0] - 2*fr*cephes.mathieu_modsem1(m, q, z)[0]
        assert_allclose(y1, y2, rtol=1e-10)

    def test_mathieu_overflow(self):
        # Check that these return NaNs instead of causing a SEGV
        assert_equal(cephes.mathieu_cem(10000, 0, 1.3), (np.nan, np.nan))
        assert_equal(cephes.mathieu_sem(10000, 0, 1.3), (np.nan, np.nan))
        assert_equal(cephes.mathieu_cem(10000, 1.5, 1.3), (np.nan, np.nan))
        assert_equal(cephes.mathieu_sem(10000, 1.5, 1.3), (np.nan, np.nan))
        assert_equal(cephes.mathieu_modcem1(10000, 1.5, 1.3), (np.nan, np.nan))
        assert_equal(cephes.mathieu_modsem1(10000, 1.5, 1.3), (np.nan, np.nan))
        assert_equal(cephes.mathieu_modcem2(10000, 1.5, 1.3), (np.nan, np.nan))
        assert_equal(cephes.mathieu_modsem2(10000, 1.5, 1.3), (np.nan, np.nan))

    def test_mathieu_ticket_1847(self):
        # Regression test --- this call had some out-of-bounds access
        # and could return nan occasionally
        for k in range(60):
            v = cephes.mathieu_modsem2(2, 100, -1)
            # Values from ACM TOMS 804 (derivate by numerical differentiation)
            assert_allclose(v[0], 0.1431742913063671074347, rtol=1e-10)
            assert_allclose(v[1], 0.9017807375832909144719, rtol=1e-4)

    def test_modfresnelm(self):
        cephes.modfresnelm(0)

    def test_modfresnelp(self):
        cephes.modfresnelp(0)

    def _check_modstruve(self):
        assert_equal(cephes.modstruve(1,0),0.0)

    def test_nbdtr(self):
        assert_equal(cephes.nbdtr(1,1,1),1.0)

    def test_nbdtrc(self):
        assert_equal(cephes.nbdtrc(1,1,1),0.0)

    def test_nbdtri(self):
        assert_equal(cephes.nbdtri(1,1,1),1.0)

    def __check_nbdtrik(self):
        cephes.nbdtrik(1,.4,.5)

    def test_nbdtrin(self):
        assert_equal(cephes.nbdtrin(1,0,0),5.0)

    def test_ncfdtr(self):
        assert_equal(cephes.ncfdtr(1,1,1,0),0.0)

    def test_ncfdtri(self):
        assert_equal(cephes.ncfdtri(1,1,1,0),0.0)

    def test_ncfdtridfd(self):
        cephes.ncfdtridfd(1,0.5,0,1)

    def __check_ncfdtridfn(self):
        cephes.ncfdtridfn(1,0.5,0,1)

    def __check_ncfdtrinc(self):
        cephes.ncfdtrinc(1,0.5,0,1)

    def test_nctdtr(self):
        assert_equal(cephes.nctdtr(1,0,0),0.5)
        assert_equal(cephes.nctdtr(9, 65536, 45), 0.0)

        assert_approx_equal(cephes.nctdtr(np.inf, 1., 1.), 0.5, 5)
        assert_(np.isnan(cephes.nctdtr(2., np.inf, 10.)))
        assert_approx_equal(cephes.nctdtr(2., 1., np.inf), 1.)

        assert_(np.isnan(cephes.nctdtr(np.nan, 1., 1.)))
        assert_(np.isnan(cephes.nctdtr(2., np.nan, 1.)))
        assert_(np.isnan(cephes.nctdtr(2., 1., np.nan)))

    def __check_nctdtridf(self):
        cephes.nctdtridf(1,0.5,0)

    def test_nctdtrinc(self):
        cephes.nctdtrinc(1,0,0)

    def test_nctdtrit(self):
        cephes.nctdtrit(.1,0.2,.5)

    def test_ndtr(self):
        assert_equal(cephes.ndtr(0), 0.5)
        assert_almost_equal(cephes.ndtr(1), 0.84134474606)

    def test_ndtri(self):
        assert_equal(cephes.ndtri(0.5),0.0)

    def test_nrdtrimn(self):
        assert_approx_equal(cephes.nrdtrimn(0.5,1,1),1.0)

    def test_nrdtrisd(self):
        assert_tol_equal(cephes.nrdtrisd(0.5,0.5,0.5), 0.0,
                         atol=0, rtol=0)

    def test_obl_ang1(self):
        cephes.obl_ang1(1,1,1,0)

    def test_obl_ang1_cv(self):
        result = cephes.obl_ang1_cv(1,1,1,1,0)
        assert_almost_equal(result[0],1.0)
        assert_almost_equal(result[1],0.0)

    def _check_obl_cv(self):
        assert_equal(cephes.obl_cv(1,1,0),2.0)

    def test_obl_rad1(self):
        cephes.obl_rad1(1,1,1,0)

    def test_obl_rad1_cv(self):
        cephes.obl_rad1_cv(1,1,1,1,0)

    def test_obl_rad2(self):
        cephes.obl_rad2(1,1,1,0)

    def test_obl_rad2_cv(self):
        cephes.obl_rad2_cv(1,1,1,1,0)

    def test_pbdv(self):
        assert_equal(cephes.pbdv(1,0),(0.0,1.0))

    def test_pbvv(self):
        cephes.pbvv(1,0)

    def test_pbwa(self):
        cephes.pbwa(1,0)

    def test_pdtr(self):
        val = cephes.pdtr(0, 1)
        assert_almost_equal(val, np.exp(-1))
        # Edge case: m = 0.
        val = cephes.pdtr([0, 1, 2], 0.0)
        assert_array_equal(val, [1, 1, 1])

    def test_pdtrc(self):
        val = cephes.pdtrc(0, 1)
        assert_almost_equal(val, 1 - np.exp(-1))
        # Edge case: m = 0.
        val = cephes.pdtrc([0, 1, 2], 0.0)
        assert_array_equal(val, [0, 0, 0])

    def test_pdtri(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            cephes.pdtri(0.5,0.5)

    def test_pdtrik(self):
        k = cephes.pdtrik(0.5, 1)
        assert_almost_equal(cephes.gammaincc(k + 1, 1), 0.5)
        # Edge case: m = 0 or very small.
        k = cephes.pdtrik([[0], [0.25], [0.95]], [0, 1e-20, 1e-6])
        assert_array_equal(k, np.zeros((3, 3)))

    def test_pro_ang1(self):
        cephes.pro_ang1(1,1,1,0)

    def test_pro_ang1_cv(self):
        assert_array_almost_equal(cephes.pro_ang1_cv(1,1,1,1,0),
                                  array((1.0,0.0)))

    def _check_pro_cv(self):
        assert_equal(cephes.pro_cv(1,1,0),2.0)

    def test_pro_rad1(self):
        cephes.pro_rad1(1,1,1,0.1)

    def test_pro_rad1_cv(self):
        cephes.pro_rad1_cv(1,1,1,1,0)

    def test_pro_rad2(self):
        cephes.pro_rad2(1,1,1,0)

    def test_pro_rad2_cv(self):
        cephes.pro_rad2_cv(1,1,1,1,0)

    def test_psi(self):
        cephes.psi(1)

    def test_radian(self):
        assert_equal(cephes.radian(0,0,0),0)

    def test_rgamma(self):
        assert_equal(cephes.rgamma(1),1.0)

    def test_round(self):
        assert_equal(cephes.round(3.4),3.0)
        assert_equal(cephes.round(-3.4),-3.0)
        assert_equal(cephes.round(3.6),4.0)
        assert_equal(cephes.round(-3.6),-4.0)
        assert_equal(cephes.round(3.5),4.0)
        assert_equal(cephes.round(-3.5),-4.0)

    def test_shichi(self):
        cephes.shichi(1)

    def test_sici(self):
        cephes.sici(1)

        s, c = cephes.sici(np.inf)
        assert_almost_equal(s, np.pi * 0.5)
        assert_almost_equal(c, 0)

        s, c = cephes.sici(-np.inf)
        assert_almost_equal(s, -np.pi * 0.5)
        assert_(np.isnan(c), "cosine integral(-inf) is not nan")

    def test_sindg(self):
        assert_equal(cephes.sindg(90),1.0)

    def test_smirnov(self):
        assert_equal(cephes.smirnov(1,.1),0.9)
        assert_(np.isnan(cephes.smirnov(1,np.nan)))

    def test_smirnovi(self):
        assert_almost_equal(cephes.smirnov(1,cephes.smirnovi(1,0.4)),0.4)
        assert_almost_equal(cephes.smirnov(1,cephes.smirnovi(1,0.6)),0.6)
        assert_(np.isnan(cephes.smirnovi(1,np.nan)))

    def test_spence(self):
        assert_equal(cephes.spence(1),0.0)

    def test_stdtr(self):
        assert_equal(cephes.stdtr(1,0),0.5)
        assert_almost_equal(cephes.stdtr(1,1), 0.75)
        assert_almost_equal(cephes.stdtr(1,2), 0.852416382349)

    def test_stdtridf(self):
        cephes.stdtridf(0.7,1)

    def test_stdtrit(self):
        cephes.stdtrit(1,0.7)

    def test_struve(self):
        assert_equal(cephes.struve(0,0),0.0)

    def test_tandg(self):
        assert_equal(cephes.tandg(45),1.0)

    def test_tklmbda(self):
        assert_almost_equal(cephes.tklmbda(1,1),1.0)

    def test_y0(self):
        cephes.y0(1)

    def test_y1(self):
        cephes.y1(1)

    def test_yn(self):
        cephes.yn(1,1)

    def test_yv(self):
        cephes.yv(1,1)

    def _check_yve(self):
        cephes.yve(1,1)

    def test_zeta(self):
        cephes.zeta(2,2)

    def test_zetac(self):
        assert_equal(cephes.zetac(0),-1.5)

    def test_wofz(self):
        z = [complex(624.2,-0.26123), complex(-0.4,3.), complex(0.6,2.),
             complex(-1.,1.), complex(-1.,-9.), complex(-1.,9.),
             complex(-0.0000000234545,1.1234), complex(-3.,5.1),
             complex(-53,30.1), complex(0.0,0.12345),
             complex(11,1), complex(-22,-2), complex(9,-28),
             complex(21,-33), complex(1e5,1e5), complex(1e14,1e14)
             ]
        w = [
            complex(-3.78270245518980507452677445620103199303131110e-7,
                    0.000903861276433172057331093754199933411710053155),
            complex(0.1764906227004816847297495349730234591778719532788,
                    -0.02146550539468457616788719893991501311573031095617),
            complex(0.2410250715772692146133539023007113781272362309451,
                    0.06087579663428089745895459735240964093522265589350),
            complex(0.30474420525691259245713884106959496013413834051768,
                    -0.20821893820283162728743734725471561394145872072738),
            complex(7.317131068972378096865595229600561710140617977e34,
                    8.321873499714402777186848353320412813066170427e34),
            complex(0.0615698507236323685519612934241429530190806818395,
                    -0.00676005783716575013073036218018565206070072304635),
            complex(0.3960793007699874918961319170187598400134746631,
                    -5.593152259116644920546186222529802777409274656e-9),
            complex(0.08217199226739447943295069917990417630675021771804,
                    -0.04701291087643609891018366143118110965272615832184),
            complex(0.00457246000350281640952328010227885008541748668738,
                    -0.00804900791411691821818731763401840373998654987934),
            complex(0.8746342859608052666092782112565360755791467973338452,
                    0.),
            complex(0.00468190164965444174367477874864366058339647648741,
                    0.0510735563901306197993676329845149741675029197050),
            complex(-0.0023193175200187620902125853834909543869428763219,
                    -0.025460054739731556004902057663500272721780776336),
            complex(9.11463368405637174660562096516414499772662584e304,
                    3.97101807145263333769664875189354358563218932e305),
            complex(-4.4927207857715598976165541011143706155432296e281,
                    -2.8019591213423077494444700357168707775769028e281),
            complex(2.820947917809305132678577516325951485807107151e-6,
                    2.820947917668257736791638444590253942253354058e-6),
            complex(2.82094791773878143474039725787438662716372268e-15,
                    2.82094791773878143474039725773333923127678361e-15)
        ]
        assert_func_equal(cephes.wofz, w, z, rtol=1e-13)


class TestAiry(TestCase):
    def test_airy(self):
        # This tests the airy function to ensure 8 place accuracy in computation

        x = special.airy(.99)
        assert_array_almost_equal(x,array([0.13689066,-0.16050153,1.19815925,0.92046818]),8)
        x = special.airy(.41)
        assert_array_almost_equal(x,array([0.25238916,-.23480512,0.80686202,0.51053919]),8)
        x = special.airy(-.36)
        assert_array_almost_equal(x,array([0.44508477,-0.23186773,0.44939534,0.48105354]),8)

    def test_airye(self):
        a = special.airye(0.01)
        b = special.airy(0.01)
        b1 = [None]*4
        for n in range(2):
            b1[n] = b[n]*exp(2.0/3.0*0.01*sqrt(0.01))
        for n in range(2,4):
            b1[n] = b[n]*exp(-abs(real(2.0/3.0*0.01*sqrt(0.01))))
        assert_array_almost_equal(a,b1,6)

    def test_bi_zeros(self):
        bi = special.bi_zeros(2)
        bia = (array([-1.17371322, -3.2710930]),
        array([-2.29443968, -4.07315509]),
        array([-0.45494438, 0.39652284]),
        array([0.60195789, -0.76031014]))
        assert_array_almost_equal(bi,bia,4)

    def test_ai_zeros(self):
        ai = special.ai_zeros(1)
        assert_array_almost_equal(ai,(array([-2.33810741]),
                                     array([-1.01879297]),
                                     array([0.5357]),
                                     array([0.7012])),4)


class TestAssocLaguerre(TestCase):
    def test_assoc_laguerre(self):
        a1 = special.genlaguerre(11,1)
        a2 = special.assoc_laguerre(.2,11,1)
        assert_array_almost_equal(a2,a1(.2),8)
        a2 = special.assoc_laguerre(1,11,1)
        assert_array_almost_equal(a2,a1(1),8)


class TestBesselpoly(TestCase):
    def test_besselpoly(self):
        pass


class TestKelvin(TestCase):
    def test_bei(self):
        mbei = special.bei(2)
        assert_almost_equal(mbei, 0.9722916273066613,5)  # this may not be exact

    def test_beip(self):
        mbeip = special.beip(2)
        assert_almost_equal(mbeip,0.91701361338403631,5)  # this may not be exact

    def test_ber(self):
        mber = special.ber(2)
        assert_almost_equal(mber,0.75173418271380821,5)  # this may not be exact

    def test_berp(self):
        mberp = special.berp(2)
        assert_almost_equal(mberp,-0.49306712470943909,5)  # this may not be exact

    def test_bei_zeros(self):
        bi = special.bi_zeros(5)
        assert_array_almost_equal(bi[0],array([-1.173713222709127,
                                               -3.271093302836352,
                                               -4.830737841662016,
                                               -6.169852128310251,
                                               -7.376762079367764]),11)

        assert_array_almost_equal(bi[1],array([-2.294439682614122,
                                               -4.073155089071828,
                                               -5.512395729663599,
                                               -6.781294445990305,
                                               -7.940178689168587]),10)

        assert_array_almost_equal(bi[2],array([-0.454944383639657,
                                               0.396522836094465,
                                               -0.367969161486959,
                                               0.349499116831805,
                                               -0.336026240133662]),11)

        assert_array_almost_equal(bi[3],array([0.601957887976239,
                                               -0.760310141492801,
                                               0.836991012619261,
                                               -0.88947990142654,
                                               0.929983638568022]),11)

    def test_beip_zeros(self):
        bip = special.beip_zeros(5)
        assert_array_almost_equal(bip,array([3.772673304934953,
                                               8.280987849760042,
                                               12.742147523633703,
                                               17.193431752512542,
                                               21.641143941167325]),4)

    def test_ber_zeros(self):
        ber = special.ber_zeros(5)
        assert_array_almost_equal(ber,array([2.84892,
                                             7.23883,
                                             11.67396,
                                             16.11356,
                                             20.55463]),4)

    def test_berp_zeros(self):
        brp = special.berp_zeros(5)
        assert_array_almost_equal(brp,array([6.03871,
                                             10.51364,
                                             14.96844,
                                             19.41758,
                                             23.86430]),4)

    def test_kelvin(self):
        mkelv = special.kelvin(2)
        assert_array_almost_equal(mkelv,(special.ber(2) + special.bei(2)*1j,
                                         special.ker(2) + special.kei(2)*1j,
                                         special.berp(2) + special.beip(2)*1j,
                                         special.kerp(2) + special.keip(2)*1j),8)

    def test_kei(self):
        mkei = special.kei(2)
        assert_almost_equal(mkei,-0.20240006776470432,5)

    def test_keip(self):
        mkeip = special.keip(2)
        assert_almost_equal(mkeip,0.21980790991960536,5)

    def test_ker(self):
        mker = special.ker(2)
        assert_almost_equal(mker,-0.041664513991509472,5)

    def test_kerp(self):
        mkerp = special.kerp(2)
        assert_almost_equal(mkerp,-0.10660096588105264,5)

    def test_kei_zeros(self):
        kei = special.kei_zeros(5)
        assert_array_almost_equal(kei,array([3.91467,
                                              8.34422,
                                              12.78256,
                                              17.22314,
                                              21.66464]),4)

    def test_keip_zeros(self):
        keip = special.keip_zeros(5)
        assert_array_almost_equal(keip,array([4.93181,
                                                9.40405,
                                                13.85827,
                                                18.30717,
                                                22.75379]),4)

    # numbers come from 9.9 of A&S pg. 381
    def test_kelvin_zeros(self):
        tmp = special.kelvin_zeros(5)
        berz,beiz,kerz,keiz,berpz,beipz,kerpz,keipz = tmp
        assert_array_almost_equal(berz,array([2.84892,
                                               7.23883,
                                               11.67396,
                                               16.11356,
                                               20.55463]),4)
        assert_array_almost_equal(beiz,array([5.02622,
                                               9.45541,
                                               13.89349,
                                               18.33398,
                                               22.77544]),4)
        assert_array_almost_equal(kerz,array([1.71854,
                                               6.12728,
                                               10.56294,
                                               15.00269,
                                               19.44382]),4)
        assert_array_almost_equal(keiz,array([3.91467,
                                               8.34422,
                                               12.78256,
                                               17.22314,
                                               21.66464]),4)
        assert_array_almost_equal(berpz,array([6.03871,
                                                10.51364,
                                                14.96844,
                                                19.41758,
                                                23.86430]),4)
        assert_array_almost_equal(beipz,array([3.77267,
                 # table from 1927 had 3.77320
                 #  but this is more accurate
                                                8.28099,
                                                12.74215,
                                                17.19343,
                                                21.64114]),4)
        assert_array_almost_equal(kerpz,array([2.66584,
                                                7.17212,
                                                11.63218,
                                                16.08312,
                                                20.53068]),4)
        assert_array_almost_equal(keipz,array([4.93181,
                                                9.40405,
                                                13.85827,
                                                18.30717,
                                                22.75379]),4)

    def test_ker_zeros(self):
        ker = special.ker_zeros(5)
        assert_array_almost_equal(ker,array([1.71854,
                                               6.12728,
                                               10.56294,
                                               15.00269,
                                               19.44381]),4)

    def test_kerp_zeros(self):
        kerp = special.kerp_zeros(5)
        assert_array_almost_equal(kerp,array([2.66584,
                                                7.17212,
                                                11.63218,
                                                16.08312,
                                                20.53068]),4)


class TestBernoulli(TestCase):
    def test_bernoulli(self):
        brn = special.bernoulli(5)
        assert_array_almost_equal(brn,array([1.0000,
                                             -0.5000,
                                             0.1667,
                                             0.0000,
                                             -0.0333,
                                             0.0000]),4)


class TestBeta(TestCase):
    def test_beta(self):
        bet = special.beta(2,4)
        betg = (special.gamma(2)*special.gamma(4))/special.gamma(6)
        assert_almost_equal(bet,betg,8)

    def test_betaln(self):
        betln = special.betaln(2,4)
        bet = log(abs(special.beta(2,4)))
        assert_almost_equal(betln,bet,8)

    def test_betainc(self):
        btinc = special.betainc(1,1,.2)
        assert_almost_equal(btinc,0.2,8)

    def test_betaincinv(self):
        y = special.betaincinv(2,4,.5)
        comp = special.betainc(2,4,y)
        assert_almost_equal(comp,.5,5)


class TestCombinatorics(TestCase):
    def test_comb(self):
        assert_array_almost_equal(special.comb([10, 10], [3, 4]), [120., 210.])
        assert_almost_equal(special.comb(10, 3), 120.)
        assert_equal(special.comb(10, 3, exact=True), 120)
        assert_equal(special.comb(10, 3, exact=True, repetition=True), 220)

    def test_comb_with_np_int64(self):
        n = 70
        k = 30
        np_n = np.int64(n)
        np_k = np.int64(k)
        assert_equal(special.comb(np_n, np_k, exact=True),
                     special.comb(n, k, exact=True))

    def test_comb_zeros(self):
        assert_equal(special.comb(2, 3, exact=True), 0)
        assert_equal(special.comb(-1, 3, exact=True), 0)
        assert_equal(special.comb(2, -1, exact=True), 0)
        assert_equal(special.comb(2, -1, exact=False), 0)
        assert_array_almost_equal(special.comb([2, -1, 2, 10], [3, 3, -1, 3]),
                [0., 0., 0., 120.])

    def test_perm(self):
        assert_array_almost_equal(special.perm([10, 10], [3, 4]), [720., 5040.])
        assert_almost_equal(special.perm(10, 3), 720.)
        assert_equal(special.perm(10, 3, exact=True), 720)

    def test_perm_zeros(self):
        assert_equal(special.perm(2, 3, exact=True), 0)
        assert_equal(special.perm(-1, 3, exact=True), 0)
        assert_equal(special.perm(2, -1, exact=True), 0)
        assert_equal(special.perm(2, -1, exact=False), 0)
        assert_array_almost_equal(special.perm([2, -1, 2, 10], [3, 3, -1, 3]),
                [0., 0., 0., 720.])


class TestTrigonometric(TestCase):
    def test_cbrt(self):
        cb = special.cbrt(27)
        cbrl = 27**(1.0/3.0)
        assert_approx_equal(cb,cbrl)

    def test_cbrtmore(self):
        cb1 = special.cbrt(27.9)
        cbrl1 = 27.9**(1.0/3.0)
        assert_almost_equal(cb1,cbrl1,8)

    def test_cosdg(self):
        cdg = special.cosdg(90)
        cdgrl = cos(pi/2.0)
        assert_almost_equal(cdg,cdgrl,8)

    def test_cosdgmore(self):
        cdgm = special.cosdg(30)
        cdgmrl = cos(pi/6.0)
        assert_almost_equal(cdgm,cdgmrl,8)

    def test_cosm1(self):
        cs = (special.cosm1(0),special.cosm1(.3),special.cosm1(pi/10))
        csrl = (cos(0)-1,cos(.3)-1,cos(pi/10)-1)
        assert_array_almost_equal(cs,csrl,8)

    def test_cotdg(self):
        ct = special.cotdg(30)
        ctrl = tan(pi/6.0)**(-1)
        assert_almost_equal(ct,ctrl,8)

    def test_cotdgmore(self):
        ct1 = special.cotdg(45)
        ctrl1 = tan(pi/4.0)**(-1)
        assert_almost_equal(ct1,ctrl1,8)

    def test_specialpoints(self):
        assert_almost_equal(special.cotdg(45), 1.0, 14)
        assert_almost_equal(special.cotdg(-45), -1.0, 14)
        assert_almost_equal(special.cotdg(90), 0.0, 14)
        assert_almost_equal(special.cotdg(-90), 0.0, 14)
        assert_almost_equal(special.cotdg(135), -1.0, 14)
        assert_almost_equal(special.cotdg(-135), 1.0, 14)
        assert_almost_equal(special.cotdg(225), 1.0, 14)
        assert_almost_equal(special.cotdg(-225), -1.0, 14)
        assert_almost_equal(special.cotdg(270), 0.0, 14)
        assert_almost_equal(special.cotdg(-270), 0.0, 14)
        assert_almost_equal(special.cotdg(315), -1.0, 14)
        assert_almost_equal(special.cotdg(-315), 1.0, 14)
        assert_almost_equal(special.cotdg(765), 1.0, 14)

    def test_sinc(self):
        # the sinc implementation and more extensive sinc tests are in numpy
        assert_array_equal(special.sinc([0]), 1)
        assert_equal(special.sinc(0.0), 1.0)

    def test_sindg(self):
        sn = special.sindg(90)
        assert_equal(sn,1.0)

    def test_sindgmore(self):
        snm = special.sindg(30)
        snmrl = sin(pi/6.0)
        assert_almost_equal(snm,snmrl,8)
        snm1 = special.sindg(45)
        snmrl1 = sin(pi/4.0)
        assert_almost_equal(snm1,snmrl1,8)


class TestTandg(TestCase):

    def test_tandg(self):
        tn = special.tandg(30)
        tnrl = tan(pi/6.0)
        assert_almost_equal(tn,tnrl,8)

    def test_tandgmore(self):
        tnm = special.tandg(45)
        tnmrl = tan(pi/4.0)
        assert_almost_equal(tnm,tnmrl,8)
        tnm1 = special.tandg(60)
        tnmrl1 = tan(pi/3.0)
        assert_almost_equal(tnm1,tnmrl1,8)

    def test_specialpoints(self):
        assert_almost_equal(special.tandg(0), 0.0, 14)
        assert_almost_equal(special.tandg(45), 1.0, 14)
        assert_almost_equal(special.tandg(-45), -1.0, 14)
        assert_almost_equal(special.tandg(135), -1.0, 14)
        assert_almost_equal(special.tandg(-135), 1.0, 14)
        assert_almost_equal(special.tandg(180), 0.0, 14)
        assert_almost_equal(special.tandg(-180), 0.0, 14)
        assert_almost_equal(special.tandg(225), 1.0, 14)
        assert_almost_equal(special.tandg(-225), -1.0, 14)
        assert_almost_equal(special.tandg(315), -1.0, 14)
        assert_almost_equal(special.tandg(-315), 1.0, 14)


class TestEllip(TestCase):
    def test_ellipj_nan(self):
        """Regression test for #912."""
        special.ellipj(0.5, np.nan)

    def test_ellipj(self):
        el = special.ellipj(0.2,0)
        rel = [sin(0.2),cos(0.2),1.0,0.20]
        assert_array_almost_equal(el,rel,13)

    def test_ellipk(self):
        elk = special.ellipk(.2)
        assert_almost_equal(elk,1.659623598610528,11)

    def test_ellipkinc(self):
        elkinc = special.ellipkinc(pi/2,.2)
        elk = special.ellipk(0.2)
        assert_almost_equal(elkinc,elk,15)
        alpha = 20*pi/180
        phi = 45*pi/180
        m = sin(alpha)**2
        elkinc = special.ellipkinc(phi,m)
        assert_almost_equal(elkinc,0.79398143,8)
        # From pg. 614 of A & S

    def test_ellipkinc_2(self):
        # Regression test for gh-3550
        # ellipkinc(phi, mbad) was NaN and mvals[2:6] were twice the correct value
        mbad = 0.68359375000000011
        phi = 0.9272952180016123
        m = np.nextafter(mbad, 0)
        mvals = []
        for j in range(10):
            mvals.append(m)
            m = np.nextafter(m, 1)
        f = special.ellipkinc(phi, mvals)
        assert_array_almost_equal_nulp(f, 1.0259330100195334 * np.ones_like(f), 1)
        # this bug also appears at phi + n * pi for at least small n
        f1 = special.ellipkinc(phi + pi, mvals)
        assert_array_almost_equal_nulp(f1, 5.1296650500976675 * np.ones_like(f1), 2)

    def test_ellipe(self):
        ele = special.ellipe(.2)
        assert_almost_equal(ele,1.4890350580958529,8)

    def test_ellipeinc(self):
        eleinc = special.ellipeinc(pi/2,.2)
        ele = special.ellipe(0.2)
        assert_almost_equal(eleinc,ele,14)
        # pg 617 of A & S
        alpha, phi = 52*pi/180,35*pi/180
        m = sin(alpha)**2
        eleinc = special.ellipeinc(phi,m)
        assert_almost_equal(eleinc, 0.58823065, 8)

    def test_ellipeinc_2(self):
        # Regression test for gh-3550
        # ellipeinc(phi, mbad) was NaN and mvals[2:6] were twice the correct value
        mbad = 0.68359375000000011
        phi = 0.9272952180016123
        m = np.nextafter(mbad, 0)
        mvals = []
        for j in range(10):
            mvals.append(m)
            m = np.nextafter(m, 1)
        f = special.ellipeinc(phi, mvals)
        assert_array_almost_equal_nulp(f, 0.84442884574781019 * np.ones_like(f), 2)
        # this bug also appears at phi + n * pi for at least small n
        f1 = special.ellipeinc(phi + pi, mvals)
        assert_array_almost_equal_nulp(f1, 3.3471442287390509 * np.ones_like(f1), 4)


class TestErf(TestCase):

    def test_erf(self):
        er = special.erf(.25)
        assert_almost_equal(er,0.2763263902,8)

    def test_erf_zeros(self):
        erz = special.erf_zeros(5)
        erzr = array([1.45061616+1.88094300j,
                     2.24465928+2.61657514j,
                     2.83974105+3.17562810j,
                     3.33546074+3.64617438j,
                     3.76900557+4.06069723j])
        assert_array_almost_equal(erz,erzr,4)

    def _check_variant_func(self, func, other_func, rtol, atol=0):
        np.random.seed(1234)
        n = 10000
        x = np.random.pareto(0.02, n) * (2*np.random.randint(0, 2, n) - 1)
        y = np.random.pareto(0.02, n) * (2*np.random.randint(0, 2, n) - 1)
        z = x + 1j*y

        old_errors = np.seterr(all='ignore')
        try:
            w = other_func(z)
            w_real = other_func(x).real

            mask = np.isfinite(w)
            w = w[mask]
            z = z[mask]

            mask = np.isfinite(w_real)
            w_real = w_real[mask]
            x = x[mask]

            # test both real and complex variants
            assert_func_equal(func, w, z, rtol=rtol, atol=atol)
            assert_func_equal(func, w_real, x, rtol=rtol, atol=atol)
        finally:
            np.seterr(**old_errors)

    def test_erfc_consistent(self):
        self._check_variant_func(
            cephes.erfc,
            lambda z: 1 - cephes.erf(z),
            rtol=1e-12,
            atol=1e-14  # <- the test function loses precision
            )

    def test_erfcx_consistent(self):
        self._check_variant_func(
            cephes.erfcx,
            lambda z: np.exp(z*z) * cephes.erfc(z),
            rtol=1e-12
            )

    def test_erfi_consistent(self):
        self._check_variant_func(
            cephes.erfi,
            lambda z: -1j * cephes.erf(1j*z),
            rtol=1e-12
            )

    def test_dawsn_consistent(self):
        self._check_variant_func(
            cephes.dawsn,
            lambda z: sqrt(pi)/2 * np.exp(-z*z) * cephes.erfi(z),
            rtol=1e-12
            )

    def test_erfcinv(self):
        i = special.erfcinv(1)
        # Use assert_array_equal instead of assert_equal, so the comparsion
        # of -0.0 and 0.0 doesn't fail.
        assert_array_equal(i, 0)

    def test_erfinv(self):
        i = special.erfinv(0)
        assert_equal(i,0)

    def test_errprint(self):
        a = special.errprint()
        b = 1-a  # a is the state 1-a inverts state
        c = special.errprint(b)  # returns last state 'a'
        assert_equal(a,c)
        d = special.errprint(a)  # returns to original state
        assert_equal(d,b)  # makes sure state was returned
        # assert_equal(d,1-a)


class TestEuler(TestCase):
    def test_euler(self):
        eu0 = special.euler(0)
        eu1 = special.euler(1)
        eu2 = special.euler(2)   # just checking segfaults
        assert_almost_equal(eu0[0],1,8)
        assert_almost_equal(eu2[2],-1,8)
        eu24 = special.euler(24)
        mathworld = [1,1,5,61,1385,50521,2702765,199360981,
                     19391512145,2404879675441,
                     370371188237525,69348874393137901,
                     15514534163557086905]
        correct = zeros((25,),'d')
        for k in range(0,13):
            if (k % 2):
                correct[2*k] = -float(mathworld[k])
            else:
                correct[2*k] = float(mathworld[k])
        olderr = np.seterr(all='ignore')
        try:
            err = nan_to_num((eu24-correct)/correct)
            errmax = max(err)
        finally:
            np.seterr(**olderr)
        assert_almost_equal(errmax, 0.0, 14)


class TestExp(TestCase):
    def test_exp2(self):
        ex = special.exp2(2)
        exrl = 2**2
        assert_equal(ex,exrl)

    def test_exp2more(self):
        exm = special.exp2(2.5)
        exmrl = 2**(2.5)
        assert_almost_equal(exm,exmrl,8)

    def test_exp10(self):
        ex = special.exp10(2)
        exrl = 10**2
        assert_approx_equal(ex,exrl)

    def test_exp10more(self):
        exm = special.exp10(2.5)
        exmrl = 10**(2.5)
        assert_almost_equal(exm,exmrl,8)

    def test_expm1(self):
        ex = (special.expm1(2),special.expm1(3),special.expm1(4))
        exrl = (exp(2)-1,exp(3)-1,exp(4)-1)
        assert_array_almost_equal(ex,exrl,8)

    def test_expm1more(self):
        ex1 = (special.expm1(2),special.expm1(2.1),special.expm1(2.2))
        exrl1 = (exp(2)-1,exp(2.1)-1,exp(2.2)-1)
        assert_array_almost_equal(ex1,exrl1,8)


class TestFactorialFunctions(TestCase):
    def test_factorial(self):
        assert_array_almost_equal([6., 24., 120.],
                special.factorial([3, 4, 5], exact=False))
        assert_equal(special.factorial(5, exact=True), 120)

    def test_factorial2(self):
        assert_array_almost_equal([105., 384., 945.],
                special.factorial2([7, 8, 9], exact=False))
        assert_equal(special.factorial2(7, exact=True), 105)

    def test_factorialk(self):
        assert_equal(special.factorialk(5, 1, exact=True), 120)
        assert_equal(special.factorialk(5, 3, exact=True), 10)


class TestFresnel(TestCase):
    def test_fresnel(self):
        frs = array(special.fresnel(.5))
        assert_array_almost_equal(frs,array([0.064732432859999287, 0.49234422587144644]),8)

    # values from pg 329  Table 7.11 of A & S
    #  slightly corrected in 4th decimal place
    def test_fresnel_zeros(self):
        szo, czo = special.fresnel_zeros(5)
        assert_array_almost_equal(szo,
                                  array([2.0093+0.2885j,
                                          2.8335+0.2443j,
                                          3.4675+0.2185j,
                                          4.0026+0.2009j,
                                          4.4742+0.1877j]),3)
        assert_array_almost_equal(czo,
                                  array([1.7437+0.3057j,
                                          2.6515+0.2529j,
                                          3.3204+0.2240j,
                                          3.8757+0.2047j,
                                          4.3611+0.1907j]),3)
        vals1 = special.fresnel(szo)[0]
        vals2 = special.fresnel(czo)[1]
        assert_array_almost_equal(vals1,0,14)
        assert_array_almost_equal(vals2,0,14)

    def test_fresnelc_zeros(self):
        szo, czo = special.fresnel_zeros(6)
        frc = special.fresnelc_zeros(6)
        assert_array_almost_equal(frc,czo,12)

    def test_fresnels_zeros(self):
        szo, czo = special.fresnel_zeros(5)
        frs = special.fresnels_zeros(5)
        assert_array_almost_equal(frs,szo,12)


class TestGamma(TestCase):
    def test_gamma(self):
        gam = special.gamma(5)
        assert_equal(gam,24.0)

    def test_gammaln(self):
        gamln = special.gammaln(3)
        lngam = log(special.gamma(3))
        assert_almost_equal(gamln,lngam,8)

    def test_gammainc(self):
        gama = special.gammainc(.5,.5)
        assert_almost_equal(gama,.7,1)

    def test_gammaincnan(self):
        gama = special.gammainc(-1,1)
        assert_(isnan(gama))

    def test_gammainczero(self):
        # bad arg but zero integration limit
        gama = special.gammainc(-1,0)
        assert_equal(gama,0.0)

    def test_gammaincc(self):
        gicc = special.gammaincc(.5,.5)
        greal = 1 - special.gammainc(.5,.5)
        assert_almost_equal(gicc,greal,8)

    def test_gammainccnan(self):
        gama = special.gammaincc(-1,1)
        assert_(isnan(gama))

    def test_gammainccinv(self):
        gccinv = special.gammainccinv(.5,.5)
        gcinv = special.gammaincinv(.5,.5)
        assert_almost_equal(gccinv,gcinv,8)

    @with_special_errors
    def test_gammaincinv(self):
        y = special.gammaincinv(.4,.4)
        x = special.gammainc(.4,y)
        assert_almost_equal(x,0.4,1)
        y = special.gammainc(10, 0.05)
        x = special.gammaincinv(10, 2.5715803516000736e-20)
        assert_almost_equal(0.05, x, decimal=10)
        assert_almost_equal(y, 2.5715803516000736e-20, decimal=10)
        x = special.gammaincinv(50, 8.20754777388471303050299243573393e-18)
        assert_almost_equal(11.0, x, decimal=10)

    @with_special_errors
    def test_975(self):
        # Regression test for ticket #975 -- switch point in algorithm
        # check that things work OK at the point, immediately next floats
        # around it, and a bit further away
        pts = [0.25,
               np.nextafter(0.25, 0), 0.25 - 1e-12,
               np.nextafter(0.25, 1), 0.25 + 1e-12]
        for xp in pts:
            y = special.gammaincinv(.4, xp)
            x = special.gammainc(0.4, y)
            assert_tol_equal(x, xp, rtol=1e-12)

    def test_rgamma(self):
        rgam = special.rgamma(8)
        rlgam = 1/special.gamma(8)
        assert_almost_equal(rgam,rlgam,8)

    def test_infinity(self):
        assert_(np.isinf(special.gamma(-1)))
        assert_equal(special.rgamma(-1), 0)


class TestHankel(TestCase):

    def test_negv1(self):
        assert_almost_equal(special.hankel1(-3,2), -special.hankel1(3,2), 14)

    def test_hankel1(self):
        hank1 = special.hankel1(1,.1)
        hankrl = (special.jv(1,.1) + special.yv(1,.1)*1j)
        assert_almost_equal(hank1,hankrl,8)

    def test_negv1e(self):
        assert_almost_equal(special.hankel1e(-3,2), -special.hankel1e(3,2), 14)

    def test_hankel1e(self):
        hank1e = special.hankel1e(1,.1)
        hankrle = special.hankel1(1,.1)*exp(-.1j)
        assert_almost_equal(hank1e,hankrle,8)

    def test_negv2(self):
        assert_almost_equal(special.hankel2(-3,2), -special.hankel2(3,2), 14)

    def test_hankel2(self):
        hank2 = special.hankel2(1,.1)
        hankrl2 = (special.jv(1,.1) - special.yv(1,.1)*1j)
        assert_almost_equal(hank2,hankrl2,8)

    def test_neg2e(self):
        assert_almost_equal(special.hankel2e(-3,2), -special.hankel2e(3,2), 14)

    def test_hankl2e(self):
        hank2e = special.hankel2e(1,.1)
        hankrl2e = special.hankel2e(1,.1)
        assert_almost_equal(hank2e,hankrl2e,8)


class TestHyper(TestCase):
    def test_h1vp(self):
        h1 = special.h1vp(1,.1)
        h1real = (special.jvp(1,.1) + special.yvp(1,.1)*1j)
        assert_almost_equal(h1,h1real,8)

    def test_h2vp(self):
        h2 = special.h2vp(1,.1)
        h2real = (special.jvp(1,.1) - special.yvp(1,.1)*1j)
        assert_almost_equal(h2,h2real,8)

    def test_hyp0f1(self):
        # scalar input
        assert_allclose(special.hyp0f1(2.5, 0.5), 1.21482702689997, rtol=1e-12)
        assert_allclose(special.hyp0f1(2.5, 0), 1.0, rtol=1e-15)

        # float input, expected values match mpmath
        x = special.hyp0f1(3.0, [-1.5, -1, 0, 1, 1.5])
        expected = np.array([0.58493659229143, 0.70566805723127, 1.0,
                             1.37789689539747, 1.60373685288480])
        assert_allclose(x, expected, rtol=1e-12)

        # complex input
        x = special.hyp0f1(3.0, np.array([-1.5, -1, 0, 1, 1.5]) + 0.j)
        assert_allclose(x, expected.astype(np.complex), rtol=1e-12)

        # test broadcasting
        x1 = [0.5, 1.5, 2.5]
        x2 = [0, 1, 0.5]
        x = special.hyp0f1(x1, x2)
        expected = [1.0, 1.8134302039235093, 1.21482702689997]
        assert_allclose(x, expected, rtol=1e-12)
        x = special.hyp0f1(np.row_stack([x1] * 2), x2)
        assert_allclose(x, np.row_stack([expected] * 2), rtol=1e-12)
        assert_raises(ValueError, special.hyp0f1,
                      np.row_stack([x1] * 3), [0, 1])

    def test_hyp1f1(self):
        hyp1 = special.hyp1f1(.1,.1,.3)
        assert_almost_equal(hyp1, 1.3498588075760032,7)

        # test contributed by Moritz Deger (2008-05-29)
        # http://projects.scipy.org/scipy/scipy/ticket/659

        # reference data obtained from mathematica [ a, b, x, m(a,b,x)]:
        # produced with test_hyp1f1.nb
        ref_data = array([[-8.38132975e+00, -1.28436461e+01, -2.91081397e+01, 1.04178330e+04],
                          [2.91076882e+00, -6.35234333e+00, -1.27083993e+01, 6.68132725e+00],
                          [-1.42938258e+01, 1.80869131e-01, 1.90038728e+01, 1.01385897e+05],
                          [5.84069088e+00, 1.33187908e+01, 2.91290106e+01, 1.59469411e+08],
                          [-2.70433202e+01, -1.16274873e+01, -2.89582384e+01, 1.39900152e+24],
                          [4.26344966e+00, -2.32701773e+01, 1.91635759e+01, 6.13816915e+21],
                          [1.20514340e+01, -3.40260240e+00, 7.26832235e+00, 1.17696112e+13],
                          [2.77372955e+01, -1.99424687e+00, 3.61332246e+00, 3.07419615e+13],
                          [1.50310939e+01, -2.91198675e+01, -1.53581080e+01, -3.79166033e+02],
                          [1.43995827e+01, 9.84311196e+00, 1.93204553e+01, 2.55836264e+10],
                          [-4.08759686e+00, 1.34437025e+01, -1.42072843e+01, 1.70778449e+01],
                          [8.05595738e+00, -1.31019838e+01, 1.52180721e+01, 3.06233294e+21],
                          [1.81815804e+01, -1.42908793e+01, 9.57868793e+00, -2.84771348e+20],
                          [-2.49671396e+01, 1.25082843e+01, -1.71562286e+01, 2.36290426e+07],
                          [2.67277673e+01, 1.70315414e+01, 6.12701450e+00, 7.77917232e+03],
                          [2.49565476e+01, 2.91694684e+01, 6.29622660e+00, 2.35300027e+02],
                          [6.11924542e+00, -1.59943768e+00, 9.57009289e+00, 1.32906326e+11],
                          [-1.47863653e+01, 2.41691301e+01, -1.89981821e+01, 2.73064953e+03],
                          [2.24070483e+01, -2.93647433e+00, 8.19281432e+00, -6.42000372e+17],
                          [8.04042600e-01, 1.82710085e+01, -1.97814534e+01, 5.48372441e-01],
                          [1.39590390e+01, 1.97318686e+01, 2.37606635e+00, 5.51923681e+00],
                          [-4.66640483e+00, -2.00237930e+01, 7.40365095e+00, 4.50310752e+00],
                          [2.76821999e+01, -6.36563968e+00, 1.11533984e+01, -9.28725179e+23],
                          [-2.56764457e+01, 1.24544906e+00, 1.06407572e+01, 1.25922076e+01],
                          [3.20447808e+00, 1.30874383e+01, 2.26098014e+01, 2.03202059e+04],
                          [-1.24809647e+01, 4.15137113e+00, -2.92265700e+01, 2.39621411e+08],
                          [2.14778108e+01, -2.35162960e+00, -1.13758664e+01, 4.46882152e-01],
                          [-9.85469168e+00, -3.28157680e+00, 1.67447548e+01, -1.07342390e+07],
                          [1.08122310e+01, -2.47353236e+01, -1.15622349e+01, -2.91733796e+03],
                          [-2.67933347e+01, -3.39100709e+00, 2.56006986e+01, -5.29275382e+09],
                          [-8.60066776e+00, -8.02200924e+00, 1.07231926e+01, 1.33548320e+06],
                          [-1.01724238e-01, -1.18479709e+01, -2.55407104e+01, 1.55436570e+00],
                          [-3.93356771e+00, 2.11106818e+01, -2.57598485e+01, 2.13467840e+01],
                          [3.74750503e+00, 1.55687633e+01, -2.92841720e+01, 1.43873509e-02],
                          [6.99726781e+00, 2.69855571e+01, -1.63707771e+01, 3.08098673e-02],
                          [-2.31996011e+01, 3.47631054e+00, 9.75119815e-01, 1.79971073e-02],
                          [2.38951044e+01, -2.91460190e+01, -2.50774708e+00, 9.56934814e+00],
                          [1.52730825e+01, 5.77062507e+00, 1.21922003e+01, 1.32345307e+09],
                          [1.74673917e+01, 1.89723426e+01, 4.94903250e+00, 9.90859484e+01],
                          [1.88971241e+01, 2.86255413e+01, 5.52360109e-01, 1.44165360e+00],
                          [1.02002319e+01, -1.66855152e+01, -2.55426235e+01, 6.56481554e+02],
                          [-1.79474153e+01, 1.22210200e+01, -1.84058212e+01, 8.24041812e+05],
                          [-1.36147103e+01, 1.32365492e+00, -7.22375200e+00, 9.92446491e+05],
                          [7.57407832e+00, 2.59738234e+01, -1.34139168e+01, 3.64037761e-02],
                          [2.21110169e+00, 1.28012666e+01, 1.62529102e+01, 1.33433085e+02],
                          [-2.64297569e+01, -1.63176658e+01, -1.11642006e+01, -2.44797251e+13],
                          [-2.46622944e+01, -3.02147372e+00, 8.29159315e+00, -3.21799070e+05],
                          [-1.37215095e+01, -1.96680183e+01, 2.91940118e+01, 3.21457520e+12],
                          [-5.45566105e+00, 2.81292086e+01, 1.72548215e-01, 9.66973000e-01],
                          [-1.55751298e+00, -8.65703373e+00, 2.68622026e+01, -3.17190834e+16],
                          [2.45393609e+01, -2.70571903e+01, 1.96815505e+01, 1.80708004e+37],
                          [5.77482829e+00, 1.53203143e+01, 2.50534322e+01, 1.14304242e+06],
                          [-1.02626819e+01, 2.36887658e+01, -2.32152102e+01, 7.28965646e+02],
                          [-1.30833446e+00, -1.28310210e+01, 1.87275544e+01, -9.33487904e+12],
                          [5.83024676e+00, -1.49279672e+01, 2.44957538e+01, -7.61083070e+27],
                          [-2.03130747e+01, 2.59641715e+01, -2.06174328e+01, 4.54744859e+04],
                          [1.97684551e+01, -2.21410519e+01, -2.26728740e+01, 3.53113026e+06],
                          [2.73673444e+01, 2.64491725e+01, 1.57599882e+01, 1.07385118e+07],
                          [5.73287971e+00, 1.21111904e+01, 1.33080171e+01, 2.63220467e+03],
                          [-2.82751072e+01, 2.08605881e+01, 9.09838900e+00, -6.60957033e-07],
                          [1.87270691e+01, -1.74437016e+01, 1.52413599e+01, 6.59572851e+27],
                          [6.60681457e+00, -2.69449855e+00, 9.78972047e+00, -2.38587870e+12],
                          [1.20895561e+01, -2.51355765e+01, 2.30096101e+01, 7.58739886e+32],
                          [-2.44682278e+01, 2.10673441e+01, -1.36705538e+01, 4.54213550e+04],
                          [-4.50665152e+00, 3.72292059e+00, -4.83403707e+00, 2.68938214e+01],
                          [-7.46540049e+00, -1.08422222e+01, -1.72203805e+01, -2.09402162e+02],
                          [-2.00307551e+01, -7.50604431e+00, -2.78640020e+01, 4.15985444e+19],
                          [1.99890876e+01, 2.20677419e+01, -2.51301778e+01, 1.23840297e-09],
                          [2.03183823e+01, -7.66942559e+00, 2.10340070e+01, 1.46285095e+31],
                          [-2.90315825e+00, -2.55785967e+01, -9.58779316e+00, 2.65714264e-01],
                          [2.73960829e+01, -1.80097203e+01, -2.03070131e+00, 2.52908999e+02],
                          [-2.11708058e+01, -2.70304032e+01, 2.48257944e+01, 3.09027527e+08],
                          [2.21959758e+01, 4.00258675e+00, -1.62853977e+01, -9.16280090e-09],
                          [1.61661840e+01, -2.26845150e+01, 2.17226940e+01, -8.24774394e+33],
                          [-3.35030306e+00, 1.32670581e+00, 9.39711214e+00, -1.47303163e+01],
                          [7.23720726e+00, -2.29763909e+01, 2.34709682e+01, -9.20711735e+29],
                          [2.71013568e+01, 1.61951087e+01, -7.11388906e-01, 2.98750911e-01],
                          [8.40057933e+00, -7.49665220e+00, 2.95587388e+01, 6.59465635e+29],
                          [-1.51603423e+01, 1.94032322e+01, -7.60044357e+00, 1.05186941e+02],
                          [-8.83788031e+00, -2.72018313e+01, 1.88269907e+00, 1.81687019e+00],
                          [-1.87283712e+01, 5.87479570e+00, -1.91210203e+01, 2.52235612e+08],
                          [-5.61338513e-01, 2.69490237e+01, 1.16660111e-01, 9.97567783e-01],
                          [-5.44354025e+00, -1.26721408e+01, -4.66831036e+00, 1.06660735e-01],
                          [-2.18846497e+00, 2.33299566e+01, 9.62564397e+00, 3.03842061e-01],
                          [6.65661299e+00, -2.39048713e+01, 1.04191807e+01, 4.73700451e+13],
                          [-2.57298921e+01, -2.60811296e+01, 2.74398110e+01, -5.32566307e+11],
                          [-1.11431826e+01, -1.59420160e+01, -1.84880553e+01, -1.01514747e+02],
                          [6.50301931e+00, 2.59859051e+01, -2.33270137e+01, 1.22760500e-02],
                          [-1.94987891e+01, -2.62123262e+01, 3.90323225e+00, 1.71658894e+01],
                          [7.26164601e+00, -1.41469402e+01, 2.81499763e+01, -2.50068329e+31],
                          [-1.52424040e+01, 2.99719005e+01, -2.85753678e+01, 1.31906693e+04],
                          [5.24149291e+00, -1.72807223e+01, 2.22129493e+01, 2.50748475e+25],
                          [3.63207230e-01, -9.54120862e-02, -2.83874044e+01, 9.43854939e-01],
                          [-2.11326457e+00, -1.25707023e+01, 1.17172130e+00, 1.20812698e+00],
                          [2.48513582e+00, 1.03652647e+01, -1.84625148e+01, 6.47910997e-02],
                          [2.65395942e+01, 2.74794672e+01, 1.29413428e+01, 2.89306132e+05],
                          [-9.49445460e+00, 1.59930921e+01, -1.49596331e+01, 3.27574841e+02],
                          [-5.89173945e+00, 9.96742426e+00, 2.60318889e+01, -3.15842908e-01],
                          [-1.15387239e+01, -2.21433107e+01, -2.17686413e+01, 1.56724718e-01],
                          [-5.30592244e+00, -2.42752190e+01, 1.29734035e+00, 1.31985534e+00]])

        for a,b,c,expected in ref_data:
            result = special.hyp1f1(a,b,c)
            assert_(abs(expected - result)/expected < 1e-4)

    def test_hyp1f1_gh2957(self):
        hyp1 = special.hyp1f1(0.5, 1.5, -709.7827128933)
        hyp2 = special.hyp1f1(0.5, 1.5, -709.7827128934)
        assert_almost_equal(hyp1, hyp2, 12)

    def test_hyp1f2(self):
        pass

    def test_hyp2f0(self):
        pass

    def test_hyp2f1(self):
        # a collection of special cases taken from AMS 55
        values = [[0.5, 1, 1.5, 0.2**2, 0.5/0.2*log((1+0.2)/(1-0.2))],
                  [0.5, 1, 1.5, -0.2**2, 1./0.2*arctan(0.2)],
                  [1, 1, 2, 0.2, -1/0.2*log(1-0.2)],
                  [3, 3.5, 1.5, 0.2**2,
                      0.5/0.2/(-5)*((1+0.2)**(-5)-(1-0.2)**(-5))],
                  [-3, 3, 0.5, sin(0.2)**2, cos(2*3*0.2)],
                  [3, 4, 8, 1, special.gamma(8)*special.gamma(8-4-3)/special.gamma(8-3)/special.gamma(8-4)],
                  [3, 2, 3-2+1, -1, 1./2**3*sqrt(pi) *
                      special.gamma(1+3-2)/special.gamma(1+0.5*3-2)/special.gamma(0.5+0.5*3)],
                  [5, 2, 5-2+1, -1, 1./2**5*sqrt(pi) *
                      special.gamma(1+5-2)/special.gamma(1+0.5*5-2)/special.gamma(0.5+0.5*5)],
                  [4, 0.5+4, 1.5-2*4, -1./3, (8./9)**(-2*4)*special.gamma(4./3) *
                      special.gamma(1.5-2*4)/special.gamma(3./2)/special.gamma(4./3-2*4)],
                  # and some others
                  # ticket #424
                  [1.5, -0.5, 1.0, -10.0, 4.1300097765277476484],
                  # negative integer a or b, with c-a-b integer and x > 0.9
                  [-2,3,1,0.95,0.715],
                  [2,-3,1,0.95,-0.007],
                  [-6,3,1,0.95,0.0000810625],
                  [2,-5,1,0.95,-0.000029375],
                  # huge negative integers
                  (10, -900, 10.5, 0.99, 1.91853705796607664803709475658e-24),
                  (10, -900, -10.5, 0.99, 3.54279200040355710199058559155e-18),
                  ]
        for i, (a, b, c, x, v) in enumerate(values):
            cv = special.hyp2f1(a, b, c, x)
            assert_almost_equal(cv, v, 8, err_msg='test #%d' % i)

    def test_hyp3f0(self):
        pass

    def test_hyperu(self):
        val1 = special.hyperu(1,0.1,100)
        assert_almost_equal(val1,0.0098153,7)
        a,b = [0.3,0.6,1.2,-2.7],[1.5,3.2,-0.4,-3.2]
        a,b = asarray(a), asarray(b)
        z = 0.5
        hypu = special.hyperu(a,b,z)
        hprl = (pi/sin(pi*b))*(special.hyp1f1(a,b,z) /
                               (special.gamma(1+a-b)*special.gamma(b)) -
                               z**(1-b)*special.hyp1f1(1+a-b,2-b,z)
                               / (special.gamma(a)*special.gamma(2-b)))
        assert_array_almost_equal(hypu,hprl,12)

    def test_hyperu_gh2287(self):
        assert_almost_equal(special.hyperu(1, 1.5, 20.2),
                            0.048360918656699191, 12)


class TestBessel(TestCase):
    def test_itj0y0(self):
        it0 = array(special.itj0y0(.2))
        assert_array_almost_equal(it0,array([0.19933433254006822, -0.34570883800412566]),8)

    def test_it2j0y0(self):
        it2 = array(special.it2j0y0(.2))
        assert_array_almost_equal(it2,array([0.0049937546274601858, -0.43423067011231614]),8)

    def test_negv_iv(self):
        assert_equal(special.iv(3,2), special.iv(-3,2))

    def test_j0(self):
        oz = special.j0(.1)
        ozr = special.jn(0,.1)
        assert_almost_equal(oz,ozr,8)

    def test_j1(self):
        o1 = special.j1(.1)
        o1r = special.jn(1,.1)
        assert_almost_equal(o1,o1r,8)

    def test_jn(self):
        jnnr = special.jn(1,.2)
        assert_almost_equal(jnnr,0.099500832639235995,8)

    def test_negv_jv(self):
        assert_almost_equal(special.jv(-3,2), -special.jv(3,2), 14)

    def test_jv(self):
        values = [[0, 0.1, 0.99750156206604002],
                  [2./3, 1e-8, 0.3239028506761532e-5],
                  [2./3, 1e-10, 0.1503423854873779e-6],
                  [3.1, 1e-10, 0.1711956265409013e-32],
                  [2./3, 4.0, -0.2325440850267039],
                  ]
        for i, (v, x, y) in enumerate(values):
            yc = special.jv(v, x)
            assert_almost_equal(yc, y, 8, err_msg='test #%d' % i)

    def test_negv_jve(self):
        assert_almost_equal(special.jve(-3,2), -special.jve(3,2), 14)

    def test_jve(self):
        jvexp = special.jve(1,.2)
        assert_almost_equal(jvexp,0.099500832639235995,8)
        jvexp1 = special.jve(1,.2+1j)
        z = .2+1j
        jvexpr = special.jv(1,z)*exp(-abs(z.imag))
        assert_almost_equal(jvexp1,jvexpr,8)

    def test_jn_zeros(self):
        jn0 = special.jn_zeros(0,5)
        jn1 = special.jn_zeros(1,5)
        assert_array_almost_equal(jn0,array([2.4048255577,
                                              5.5200781103,
                                              8.6537279129,
                                              11.7915344391,
                                              14.9309177086]),4)
        assert_array_almost_equal(jn1,array([3.83171,
                                              7.01559,
                                              10.17347,
                                              13.32369,
                                              16.47063]),4)

        jn102 = special.jn_zeros(102,5)
        assert_tol_equal(jn102, array([110.89174935992040343,
                                       117.83464175788308398,
                                       123.70194191713507279,
                                       129.02417238949092824,
                                       134.00114761868422559]), rtol=1e-13)

        jn301 = special.jn_zeros(301,5)
        assert_tol_equal(jn301, array([313.59097866698830153,
                                       323.21549776096288280,
                                       331.22338738656748796,
                                       338.39676338872084500,
                                       345.03284233056064157]), rtol=1e-13)

    def test_jn_zeros_slow(self):
        jn0 = special.jn_zeros(0, 300)
        assert_tol_equal(jn0[260-1], 816.02884495068867280, rtol=1e-13)
        assert_tol_equal(jn0[280-1], 878.86068707124422606, rtol=1e-13)
        assert_tol_equal(jn0[300-1], 941.69253065317954064, rtol=1e-13)

        jn10 = special.jn_zeros(10, 300)
        assert_tol_equal(jn10[260-1], 831.67668514305631151, rtol=1e-13)
        assert_tol_equal(jn10[280-1], 894.51275095371316931, rtol=1e-13)
        assert_tol_equal(jn10[300-1], 957.34826370866539775, rtol=1e-13)

        jn3010 = special.jn_zeros(3010,5)
        assert_tol_equal(jn3010, array([3036.86590780927,
                                        3057.06598526482,
                                        3073.66360690272,
                                        3088.37736494778,
                                        3101.86438139042]), rtol=1e-8)

    def test_jnjnp_zeros(self):
        jn = special.jn

        def jnp(n, x):
            return (jn(n-1,x) - jn(n+1,x))/2
        for nt in range(1, 30):
            z, n, m, t = special.jnjnp_zeros(nt)
            for zz, nn, tt in zip(z, n, t):
                if tt == 0:
                    assert_allclose(jn(nn, zz), 0, atol=1e-6)
                elif tt == 1:
                    assert_allclose(jnp(nn, zz), 0, atol=1e-6)
                else:
                    raise AssertionError("Invalid t return for nt=%d" % nt)

    def test_jnp_zeros(self):
        jnp = special.jnp_zeros(1,5)
        assert_array_almost_equal(jnp, array([1.84118,
                                                5.33144,
                                                8.53632,
                                                11.70600,
                                                14.86359]),4)
        jnp = special.jnp_zeros(443,5)
        assert_tol_equal(special.jvp(443, jnp), 0, atol=1e-15)

    def test_jnyn_zeros(self):
        jnz = special.jnyn_zeros(1,5)
        assert_array_almost_equal(jnz,(array([3.83171,
                                                7.01559,
                                                10.17347,
                                                13.32369,
                                                16.47063]),
                                       array([1.84118,
                                                5.33144,
                                                8.53632,
                                                11.70600,
                                                14.86359]),
                                       array([2.19714,
                                                5.42968,
                                                8.59601,
                                                11.74915,
                                                14.89744]),
                                       array([3.68302,
                                                6.94150,
                                                10.12340,
                                                13.28576,
                                                16.44006])),5)

    def test_jvp(self):
        jvprim = special.jvp(2,2)
        jv0 = (special.jv(1,2)-special.jv(3,2))/2
        assert_almost_equal(jvprim,jv0,10)

    def test_k0(self):
        ozk = special.k0(.1)
        ozkr = special.kv(0,.1)
        assert_almost_equal(ozk,ozkr,8)

    def test_k0e(self):
        ozke = special.k0e(.1)
        ozker = special.kve(0,.1)
        assert_almost_equal(ozke,ozker,8)

    def test_k1(self):
        o1k = special.k1(.1)
        o1kr = special.kv(1,.1)
        assert_almost_equal(o1k,o1kr,8)

    def test_k1e(self):
        o1ke = special.k1e(.1)
        o1ker = special.kve(1,.1)
        assert_almost_equal(o1ke,o1ker,8)

    def test_jacobi(self):
        a = 5*rand() - 1
        b = 5*rand() - 1
        P0 = special.jacobi(0,a,b)
        P1 = special.jacobi(1,a,b)
        P2 = special.jacobi(2,a,b)
        P3 = special.jacobi(3,a,b)

        assert_array_almost_equal(P0.c,[1],13)
        assert_array_almost_equal(P1.c,array([a+b+2,a-b])/2.0,13)
        cp = [(a+b+3)*(a+b+4), 4*(a+b+3)*(a+2), 4*(a+1)*(a+2)]
        p2c = [cp[0],cp[1]-2*cp[0],cp[2]-cp[1]+cp[0]]
        assert_array_almost_equal(P2.c,array(p2c)/8.0,13)
        cp = [(a+b+4)*(a+b+5)*(a+b+6),6*(a+b+4)*(a+b+5)*(a+3),
              12*(a+b+4)*(a+2)*(a+3),8*(a+1)*(a+2)*(a+3)]
        p3c = [cp[0],cp[1]-3*cp[0],cp[2]-2*cp[1]+3*cp[0],cp[3]-cp[2]+cp[1]-cp[0]]
        assert_array_almost_equal(P3.c,array(p3c)/48.0,13)

    def test_kn(self):
        kn1 = special.kn(0,.2)
        assert_almost_equal(kn1,1.7527038555281462,8)

    def test_negv_kv(self):
        assert_equal(special.kv(3.0, 2.2), special.kv(-3.0, 2.2))

    def test_kv0(self):
        kv0 = special.kv(0,.2)
        assert_almost_equal(kv0, 1.7527038555281462, 10)

    def test_kv1(self):
        kv1 = special.kv(1,0.2)
        assert_almost_equal(kv1, 4.775972543220472, 10)

    def test_kv2(self):
        kv2 = special.kv(2,0.2)
        assert_almost_equal(kv2, 49.51242928773287, 10)

    def test_kn_largeorder(self):
        assert_allclose(special.kn(32, 1), 1.7516596664574289e+43)

    def test_kv_largearg(self):
        assert_equal(special.kv(0, 1e19), 0)

    def test_negv_kve(self):
        assert_equal(special.kve(3.0, 2.2), special.kve(-3.0, 2.2))

    def test_kve(self):
        kve1 = special.kve(0,.2)
        kv1 = special.kv(0,.2)*exp(.2)
        assert_almost_equal(kve1,kv1,8)
        z = .2+1j
        kve2 = special.kve(0,z)
        kv2 = special.kv(0,z)*exp(z)
        assert_almost_equal(kve2,kv2,8)

    def test_kvp_v0n1(self):
        z = 2.2
        assert_almost_equal(-special.kv(1,z), special.kvp(0,z, n=1), 10)

    def test_kvp_n1(self):
        v = 3.
        z = 2.2
        xc = -special.kv(v+1,z) + v/z*special.kv(v,z)
        x = special.kvp(v,z, n=1)
        assert_almost_equal(xc, x, 10)   # this function (kvp) is broken

    def test_kvp_n2(self):
        v = 3.
        z = 2.2
        xc = (z**2+v**2-v)/z**2 * special.kv(v,z) + special.kv(v+1,z)/z
        x = special.kvp(v, z, n=2)
        assert_almost_equal(xc, x, 10)

    def test_y0(self):
        oz = special.y0(.1)
        ozr = special.yn(0,.1)
        assert_almost_equal(oz,ozr,8)

    def test_y1(self):
        o1 = special.y1(.1)
        o1r = special.yn(1,.1)
        assert_almost_equal(o1,o1r,8)

    def test_y0_zeros(self):
        yo,ypo = special.y0_zeros(2)
        zo,zpo = special.y0_zeros(2,complex=1)
        all = r_[yo,zo]
        allval = r_[ypo,zpo]
        assert_array_almost_equal(abs(special.yv(0.0,all)),0.0,11)
        assert_array_almost_equal(abs(special.yv(1,all)-allval),0.0,11)

    def test_y1_zeros(self):
        y1 = special.y1_zeros(1)
        assert_array_almost_equal(y1,(array([2.19714]),array([0.52079])),5)

    def test_y1p_zeros(self):
        y1p = special.y1p_zeros(1,complex=1)
        assert_array_almost_equal(y1p,(array([0.5768+0.904j]), array([-0.7635+0.5892j])),3)

    def test_yn_zeros(self):
        an = special.yn_zeros(4,2)
        assert_array_almost_equal(an,array([5.64515, 9.36162]),5)
        an = special.yn_zeros(443,5)
        assert_tol_equal(an, [450.13573091578090314, 463.05692376675001542,
                              472.80651546418663566, 481.27353184725625838,
                              488.98055964441374646], rtol=1e-15)

    def test_ynp_zeros(self):
        ao = special.ynp_zeros(0,2)
        assert_array_almost_equal(ao,array([2.19714133, 5.42968104]),6)
        ao = special.ynp_zeros(43,5)
        assert_tol_equal(special.yvp(43, ao), 0, atol=1e-15)
        ao = special.ynp_zeros(443,5)
        assert_tol_equal(special.yvp(443, ao), 0, atol=1e-9)

    def test_ynp_zeros_large_order(self):
        ao = special.ynp_zeros(443,5)
        assert_tol_equal(special.yvp(443, ao), 0, atol=1e-14)

    def test_yn(self):
        yn2n = special.yn(1,.2)
        assert_almost_equal(yn2n,-3.3238249881118471,8)

    def test_negv_yv(self):
        assert_almost_equal(special.yv(-3,2), -special.yv(3,2), 14)

    def test_yv(self):
        yv2 = special.yv(1,.2)
        assert_almost_equal(yv2,-3.3238249881118471,8)

    def test_negv_yve(self):
        assert_almost_equal(special.yve(-3,2), -special.yve(3,2), 14)

    def test_yve(self):
        yve2 = special.yve(1,.2)
        assert_almost_equal(yve2,-3.3238249881118471,8)
        yve2r = special.yv(1,.2+1j)*exp(-1)
        yve22 = special.yve(1,.2+1j)
        assert_almost_equal(yve22,yve2r,8)

    def test_yvp(self):
        yvpr = (special.yv(1,.2) - special.yv(3,.2))/2.0
        yvp1 = special.yvp(2,.2)
        assert_array_almost_equal(yvp1,yvpr,10)

    def _cephes_vs_amos_points(self):
        """Yield points at which to compare Cephes implementation to AMOS"""
        # check several points, including large-amplitude ones
        for v in [-120, -100.3, -20., -10., -1., -.5,
                  0., 1., 12.49, 120., 301]:
            for z in [-1300, -11, -10, -1, 1., 10., 200.5, 401., 600.5,
                      700.6, 1300, 10003]:
                yield v, z

        # check half-integers; these are problematic points at least
        # for cephes/iv
        for v in 0.5 + arange(-60, 60):
            yield v, 3.5

    def check_cephes_vs_amos(self, f1, f2, rtol=1e-11, atol=0, skip=None):
        for v, z in self._cephes_vs_amos_points():
            if skip is not None and skip(v, z):
                continue
            c1, c2, c3 = f1(v, z), f1(v,z+0j), f2(int(v), z)
            if np.isinf(c1):
                assert_(np.abs(c2) >= 1e300, (v, z))
            elif np.isnan(c1):
                assert_(c2.imag != 0, (v, z))
            else:
                assert_tol_equal(c1, c2, err_msg=(v, z), rtol=rtol, atol=atol)
                if v == int(v):
                    assert_tol_equal(c3, c2, err_msg=(v, z),
                                     rtol=rtol, atol=atol)

    def test_jv_cephes_vs_amos(self):
        self.check_cephes_vs_amos(special.jv, special.jn, rtol=1e-10, atol=1e-305)

    def test_yv_cephes_vs_amos(self):
        self.check_cephes_vs_amos(special.yv, special.yn, rtol=1e-11, atol=1e-305)

    def test_yv_cephes_vs_amos_only_small_orders(self):
        skipper = lambda v, z: (abs(v) > 50)
        self.check_cephes_vs_amos(special.yv, special.yn, rtol=1e-11, atol=1e-305, skip=skipper)

    def test_iv_cephes_vs_amos(self):
        olderr = np.seterr(all='ignore')
        try:
            self.check_cephes_vs_amos(special.iv, special.iv, rtol=5e-9, atol=1e-305)
        finally:
            np.seterr(**olderr)

    @dec.slow
    def test_iv_cephes_vs_amos_mass_test(self):
        N = 1000000
        np.random.seed(1)
        v = np.random.pareto(0.5, N) * (-1)**np.random.randint(2, size=N)
        x = np.random.pareto(0.2, N) * (-1)**np.random.randint(2, size=N)

        imsk = (np.random.randint(8, size=N) == 0)
        v[imsk] = v[imsk].astype(int)

        old_err = np.seterr(all='ignore')
        try:
            c1 = special.iv(v, x)
            c2 = special.iv(v, x+0j)

            # deal with differences in the inf and zero cutoffs
            c1[abs(c1) > 1e300] = np.inf
            c2[abs(c2) > 1e300] = np.inf
            c1[abs(c1) < 1e-300] = 0
            c2[abs(c2) < 1e-300] = 0

            dc = abs(c1/c2 - 1)
            dc[np.isnan(dc)] = 0
        finally:
            np.seterr(**old_err)

        k = np.argmax(dc)

        # Most error apparently comes from AMOS and not our implementation;
        # there are some problems near integer orders there
        assert_(dc[k] < 2e-7, (v[k], x[k], special.iv(v[k], x[k]), special.iv(v[k], x[k]+0j)))

    def test_kv_cephes_vs_amos(self):
        self.check_cephes_vs_amos(special.kv, special.kn, rtol=1e-9, atol=1e-305)
        self.check_cephes_vs_amos(special.kv, special.kv, rtol=1e-9, atol=1e-305)

    def test_ticket_623(self):
        assert_tol_equal(special.jv(3, 4), 0.43017147387562193)
        assert_tol_equal(special.jv(301, 1300), 0.0183487151115275)
        assert_tol_equal(special.jv(301, 1296.0682), -0.0224174325312048)

    def test_ticket_853(self):
        """Negative-order Bessels"""
        # cephes
        assert_tol_equal(special.jv(-1, 1), -0.4400505857449335)
        assert_tol_equal(special.jv(-2, 1), 0.1149034849319005)
        assert_tol_equal(special.yv(-1, 1), 0.7812128213002887)
        assert_tol_equal(special.yv(-2, 1), -1.650682606816255)
        assert_tol_equal(special.iv(-1, 1), 0.5651591039924851)
        assert_tol_equal(special.iv(-2, 1), 0.1357476697670383)
        assert_tol_equal(special.kv(-1, 1), 0.6019072301972347)
        assert_tol_equal(special.kv(-2, 1), 1.624838898635178)
        assert_tol_equal(special.jv(-0.5, 1), 0.43109886801837607952)
        assert_tol_equal(special.yv(-0.5, 1), 0.6713967071418031)
        assert_tol_equal(special.iv(-0.5, 1), 1.231200214592967)
        assert_tol_equal(special.kv(-0.5, 1), 0.4610685044478945)
        # amos
        assert_tol_equal(special.jv(-1, 1+0j), -0.4400505857449335)
        assert_tol_equal(special.jv(-2, 1+0j), 0.1149034849319005)
        assert_tol_equal(special.yv(-1, 1+0j), 0.7812128213002887)
        assert_tol_equal(special.yv(-2, 1+0j), -1.650682606816255)

        assert_tol_equal(special.iv(-1, 1+0j), 0.5651591039924851)
        assert_tol_equal(special.iv(-2, 1+0j), 0.1357476697670383)
        assert_tol_equal(special.kv(-1, 1+0j), 0.6019072301972347)
        assert_tol_equal(special.kv(-2, 1+0j), 1.624838898635178)

        assert_tol_equal(special.jv(-0.5, 1+0j), 0.43109886801837607952)
        assert_tol_equal(special.jv(-0.5, 1+1j), 0.2628946385649065-0.827050182040562j)
        assert_tol_equal(special.yv(-0.5, 1+0j), 0.6713967071418031)
        assert_tol_equal(special.yv(-0.5, 1+1j), 0.967901282890131+0.0602046062142816j)

        assert_tol_equal(special.iv(-0.5, 1+0j), 1.231200214592967)
        assert_tol_equal(special.iv(-0.5, 1+1j), 0.77070737376928+0.39891821043561j)
        assert_tol_equal(special.kv(-0.5, 1+0j), 0.4610685044478945)
        assert_tol_equal(special.kv(-0.5, 1+1j), 0.06868578341999-0.38157825981268j)

        assert_tol_equal(special.jve(-0.5,1+0.3j), special.jv(-0.5, 1+0.3j)*exp(-0.3))
        assert_tol_equal(special.yve(-0.5,1+0.3j), special.yv(-0.5, 1+0.3j)*exp(-0.3))
        assert_tol_equal(special.ive(-0.5,0.3+1j), special.iv(-0.5, 0.3+1j)*exp(-0.3))
        assert_tol_equal(special.kve(-0.5,0.3+1j), special.kv(-0.5, 0.3+1j)*exp(0.3+1j))

        assert_tol_equal(special.hankel1(-0.5, 1+1j), special.jv(-0.5, 1+1j) + 1j*special.yv(-0.5,1+1j))
        assert_tol_equal(special.hankel2(-0.5, 1+1j), special.jv(-0.5, 1+1j) - 1j*special.yv(-0.5,1+1j))

    def test_ticket_854(self):
        """Real-valued Bessel domains"""
        assert_(isnan(special.jv(0.5, -1)))
        assert_(isnan(special.iv(0.5, -1)))
        assert_(isnan(special.yv(0.5, -1)))
        assert_(isnan(special.yv(1, -1)))
        assert_(isnan(special.kv(0.5, -1)))
        assert_(isnan(special.kv(1, -1)))
        assert_(isnan(special.jve(0.5, -1)))
        assert_(isnan(special.ive(0.5, -1)))
        assert_(isnan(special.yve(0.5, -1)))
        assert_(isnan(special.yve(1, -1)))
        assert_(isnan(special.kve(0.5, -1)))
        assert_(isnan(special.kve(1, -1)))
        assert_(isnan(special.airye(-1)[0:2]).all(), special.airye(-1))
        assert_(not isnan(special.airye(-1)[2:4]).any(), special.airye(-1))

    def test_ticket_503(self):
        """Real-valued Bessel I overflow"""
        assert_tol_equal(special.iv(1, 700), 1.528500390233901e302)
        assert_tol_equal(special.iv(1000, 1120), 1.301564549405821e301)

    def test_iv_hyperg_poles(self):
        assert_tol_equal(special.iv(-0.5, 1), 1.231200214592967)

    def iv_series(self, v, z, n=200):
        k = arange(0, n).astype(float_)
        r = (v+2*k)*log(.5*z) - special.gammaln(k+1) - special.gammaln(v+k+1)
        r[isnan(r)] = inf
        r = exp(r)
        err = abs(r).max() * finfo(float_).eps * n + abs(r[-1])*10
        return r.sum(), err

    def test_i0_series(self):
        for z in [1., 10., 200.5]:
            value, err = self.iv_series(0, z)
            assert_tol_equal(special.i0(z), value, atol=err, err_msg=z)

    def test_i1_series(self):
        for z in [1., 10., 200.5]:
            value, err = self.iv_series(1, z)
            assert_tol_equal(special.i1(z), value, atol=err, err_msg=z)

    def test_iv_series(self):
        for v in [-20., -10., -1., 0., 1., 12.49, 120.]:
            for z in [1., 10., 200.5, -1+2j]:
                value, err = self.iv_series(v, z)
                assert_tol_equal(special.iv(v, z), value, atol=err, err_msg=(v, z))

    def test_i0(self):
        values = [[0.0, 1.0],
                  [1e-10, 1.0],
                  [0.1, 0.9071009258],
                  [0.5, 0.6450352706],
                  [1.0, 0.4657596077],
                  [2.5, 0.2700464416],
                  [5.0, 0.1835408126],
                  [20.0, 0.0897803119],
                  ]
        for i, (x, v) in enumerate(values):
            cv = special.i0(x) * exp(-x)
            assert_almost_equal(cv, v, 8, err_msg='test #%d' % i)

    def test_i0e(self):
        oize = special.i0e(.1)
        oizer = special.ive(0,.1)
        assert_almost_equal(oize,oizer,8)

    def test_i1(self):
        values = [[0.0, 0.0],
                  [1e-10, 0.4999999999500000e-10],
                  [0.1, 0.0452984468],
                  [0.5, 0.1564208032],
                  [1.0, 0.2079104154],
                  [5.0, 0.1639722669],
                  [20.0, 0.0875062222],
                  ]
        for i, (x, v) in enumerate(values):
            cv = special.i1(x) * exp(-x)
            assert_almost_equal(cv, v, 8, err_msg='test #%d' % i)

    def test_i1e(self):
        oi1e = special.i1e(.1)
        oi1er = special.ive(1,.1)
        assert_almost_equal(oi1e,oi1er,8)

    def test_iti0k0(self):
        iti0 = array(special.iti0k0(5))
        assert_array_almost_equal(iti0,array([31.848667776169801, 1.5673873907283657]),5)

    def test_it2i0k0(self):
        it2k = special.it2i0k0(.1)
        assert_array_almost_equal(it2k,array([0.0012503906973464409, 3.3309450354686687]),6)

    def test_iv(self):
        iv1 = special.iv(0,.1)*exp(-.1)
        assert_almost_equal(iv1,0.90710092578230106,10)

    def test_negv_ive(self):
        assert_equal(special.ive(3,2), special.ive(-3,2))

    def test_ive(self):
        ive1 = special.ive(0,.1)
        iv1 = special.iv(0,.1)*exp(-.1)
        assert_almost_equal(ive1,iv1,10)

    def test_ivp0(self):
        assert_almost_equal(special.iv(1,2), special.ivp(0,2), 10)

    def test_ivp(self):
        y = (special.iv(0,2) + special.iv(2,2))/2
        x = special.ivp(1,2)
        assert_almost_equal(x,y,10)


class TestLaguerre(TestCase):
    def test_laguerre(self):
        lag0 = special.laguerre(0)
        lag1 = special.laguerre(1)
        lag2 = special.laguerre(2)
        lag3 = special.laguerre(3)
        lag4 = special.laguerre(4)
        lag5 = special.laguerre(5)
        assert_array_almost_equal(lag0.c,[1],13)
        assert_array_almost_equal(lag1.c,[-1,1],13)
        assert_array_almost_equal(lag2.c,array([1,-4,2])/2.0,13)
        assert_array_almost_equal(lag3.c,array([-1,9,-18,6])/6.0,13)
        assert_array_almost_equal(lag4.c,array([1,-16,72,-96,24])/24.0,13)
        assert_array_almost_equal(lag5.c,array([-1,25,-200,600,-600,120])/120.0,13)

    def test_genlaguerre(self):
        k = 5*rand()-0.9
        lag0 = special.genlaguerre(0,k)
        lag1 = special.genlaguerre(1,k)
        lag2 = special.genlaguerre(2,k)
        lag3 = special.genlaguerre(3,k)
        assert_equal(lag0.c,[1])
        assert_equal(lag1.c,[-1,k+1])
        assert_almost_equal(lag2.c,array([1,-2*(k+2),(k+1.)*(k+2.)])/2.0)
        assert_almost_equal(lag3.c,array([-1,3*(k+3),-3*(k+2)*(k+3),(k+1)*(k+2)*(k+3)])/6.0)


# Base polynomials come from Abrahmowitz and Stegan
class TestLegendre(TestCase):
    def test_legendre(self):
        leg0 = special.legendre(0)
        leg1 = special.legendre(1)
        leg2 = special.legendre(2)
        leg3 = special.legendre(3)
        leg4 = special.legendre(4)
        leg5 = special.legendre(5)
        assert_equal(leg0.c,[1])
        assert_equal(leg1.c,[1,0])
        assert_equal(leg2.c,array([3,0,-1])/2.0)
        assert_almost_equal(leg3.c,array([5,0,-3,0])/2.0)
        assert_almost_equal(leg4.c,array([35,0,-30,0,3])/8.0)
        assert_almost_equal(leg5.c,array([63,0,-70,0,15,0])/8.0)


class TestLambda(TestCase):
    def test_lmbda(self):
        lam = special.lmbda(1,.1)
        lamr = (array([special.jn(0,.1), 2*special.jn(1,.1)/.1]),
                array([special.jvp(0,.1), -2*special.jv(1,.1)/.01 + 2*special.jvp(1,.1)/.1]))
        assert_array_almost_equal(lam,lamr,8)


class TestLog1p(TestCase):
    def test_log1p(self):
        l1p = (special.log1p(10), special.log1p(11), special.log1p(12))
        l1prl = (log(11), log(12), log(13))
        assert_array_almost_equal(l1p,l1prl,8)

    def test_log1pmore(self):
        l1pm = (special.log1p(1), special.log1p(1.1), special.log1p(1.2))
        l1pmrl = (log(2),log(2.1),log(2.2))
        assert_array_almost_equal(l1pm,l1pmrl,8)


class TestLegendreFunctions(TestCase):
    def test_clpmn(self):
        z = 0.5+0.3j
        clp = special.clpmn(2, 2, z, 3)
        assert_array_almost_equal(clp,
                   (array([[1.0000, z, 0.5*(3*z*z-1)],
                           [0.0000, sqrt(z*z-1), 3*z*sqrt(z*z-1)],
                           [0.0000, 0.0000, 3*(z*z-1)]]),
                    array([[0.0000, 1.0000, 3*z],
                           [0.0000, z/sqrt(z*z-1), 3*(2*z*z-1)/sqrt(z*z-1)],
                           [0.0000, 0.0000, 6*z]])),
                    7)

    def test_clpmn_close_to_real_2(self):
        eps = 1e-10
        m = 1
        n = 3
        x = 0.5
        clp_plus = special.clpmn(m, n, x+1j*eps, 2)[0][m, n]
        clp_minus = special.clpmn(m, n, x-1j*eps, 2)[0][m, n]
        assert_array_almost_equal(array([clp_plus, clp_minus]),
                                  array([special.lpmv(m, n, x),
                                         special.lpmv(m, n, x)]),
                                  7)

    def test_clpmn_close_to_real_3(self):
        eps = 1e-10
        m = 1
        n = 3
        x = 0.5
        clp_plus = special.clpmn(m, n, x+1j*eps, 3)[0][m, n]
        clp_minus = special.clpmn(m, n, x-1j*eps, 3)[0][m, n]
        assert_array_almost_equal(array([clp_plus, clp_minus]),
                                  array([special.lpmv(m, n, x)*np.exp(-0.5j*m*np.pi),
                                         special.lpmv(m, n, x)*np.exp(0.5j*m*np.pi)]),
                                  7)

    def test_clpmn_across_unit_circle(self):
        eps = 1e-7
        m = 1
        n = 1
        x = 1j
        for type in [2, 3]:
            assert_almost_equal(special.clpmn(m, n, x+1j*eps, type)[0][m, n],
                            special.clpmn(m, n, x-1j*eps, type)[0][m, n], 6)

    def test_inf(self):
        for z in (1, -1):
            for n in range(4):
                for m in range(1, n):
                    lp = special.clpmn(m, n, z)
                    assert_(np.isinf(lp[1][1,1:]).all())
                    lp = special.lpmn(m, n, z)
                    assert_(np.isinf(lp[1][1,1:]).all())

    def test_deriv_clpmn(self):
        # data inside and outside of the unit circle
        zvals = [0.5+0.5j, -0.5+0.5j, -0.5-0.5j, 0.5-0.5j,
                 1+1j, -1+1j, -1-1j, 1-1j]
        m = 2
        n = 3
        for type in [2, 3]:
            for z in zvals:
                for h in [1e-3, 1e-3j]:
                    approx_derivative = (special.clpmn(m, n, z+0.5*h, type)[0]
                                         - special.clpmn(m, n, z-0.5*h, type)[0])/h
                    assert_allclose(special.clpmn(m, n, z, type)[1],
                                    approx_derivative,
                                    rtol=1e-4)

    def test_lpmn(self):
        lp = special.lpmn(0,2,.5)
        assert_array_almost_equal(lp,(array([[1.00000,
                                                      0.50000,
                                                      -0.12500]]),
                                      array([[0.00000,
                                                      1.00000,
                                                      1.50000]])),4)

    def test_lpn(self):
        lpnf = special.lpn(2,.5)
        assert_array_almost_equal(lpnf,(array([1.00000,
                                                        0.50000,
                                                        -0.12500]),
                                      array([0.00000,
                                                      1.00000,
                                                      1.50000])),4)

    def test_lpmv(self):
        lp = special.lpmv(0,2,.5)
        assert_almost_equal(lp,-0.125,7)
        lp = special.lpmv(0,40,.001)
        assert_almost_equal(lp,0.1252678976534484,7)

        # XXX: this is outside the domain of the current implementation,
        #      so ensure it returns a NaN rather than a wrong answer.
        olderr = np.seterr(all='ignore')
        try:
            lp = special.lpmv(-1,-1,.001)
        finally:
            np.seterr(**olderr)
        assert_(lp != 0 or np.isnan(lp))

    def test_lqmn(self):
        lqmnf = special.lqmn(0,2,.5)
        lqmnf = special.lqmn(0,2,.5)
        lqf = special.lqn(2,.5)
        assert_array_almost_equal(lqmnf[0][0],lqf[0],4)
        assert_array_almost_equal(lqmnf[1][0],lqf[1],4)

    def test_lqmn_shape(self):
        a, b = special.lqmn(4, 4, 1.1)
        assert_equal(a.shape, (5, 5))
        assert_equal(b.shape, (5, 5))

        a, b = special.lqmn(4, 0, 1.1)
        assert_equal(a.shape, (5, 1))
        assert_equal(b.shape, (5, 1))

    def test_lqn(self):
        lqf = special.lqn(2,.5)
        assert_array_almost_equal(lqf,(array([0.5493, -0.7253, -0.8187]),
                                       array([1.3333, 1.216, -0.8427])),4)


class TestMathieu(TestCase):

    def test_mathieu_a(self):
        pass

    def test_mathieu_even_coef(self):
        mc = special.mathieu_even_coef(2,5)
        # Q not defined broken and cannot figure out proper reporting order

    def test_mathieu_odd_coef(self):
        # same problem as above
        pass


class TestFresnelIntegral(TestCase):

    def test_modfresnelp(self):
        pass

    def test_modfresnelm(self):
        pass


class TestOblCvSeq(TestCase):
    def test_obl_cv_seq(self):
        obl = special.obl_cv_seq(0,3,1)
        assert_array_almost_equal(obl,array([-0.348602,
                                              1.393206,
                                              5.486800,
                                              11.492120]),5)


class TestParabolicCylinder(TestCase):
    def test_pbdn_seq(self):
        pb = special.pbdn_seq(1,.1)
        assert_array_almost_equal(pb,(array([0.9975,
                                              0.0998]),
                                      array([-0.0499,
                                             0.9925])),4)

    def test_pbdv(self):
        pbv = special.pbdv(1,.2)
        derrl = 1/2*(.2)*special.pbdv(1,.2)[0] - special.pbdv(0,.2)[0]

    def test_pbdv_seq(self):
        pbn = special.pbdn_seq(1,.1)
        pbv = special.pbdv_seq(1,.1)
        assert_array_almost_equal(pbv,(real(pbn[0]),real(pbn[1])),4)

    def test_pbdv_points(self):
        # simple case
        eta = np.linspace(-10, 10, 5)
        z = 2**(eta/2)*np.sqrt(np.pi)/special.gamma(.5-.5*eta)
        assert_tol_equal(special.pbdv(eta, 0.)[0], z, rtol=1e-14, atol=1e-14)

        # some points
        assert_tol_equal(special.pbdv(10.34, 20.44)[0], 1.3731383034455e-32, rtol=1e-12)
        assert_tol_equal(special.pbdv(-9.53, 3.44)[0], 3.166735001119246e-8, rtol=1e-12)

    def test_pbdv_gradient(self):
        x = np.linspace(-4, 4, 8)[:,None]
        eta = np.linspace(-10, 10, 5)[None,:]

        p = special.pbdv(eta, x)
        eps = 1e-7 + 1e-7*abs(x)
        dp = (special.pbdv(eta, x + eps)[0] - special.pbdv(eta, x - eps)[0]) / eps / 2.
        assert_tol_equal(p[1], dp, rtol=1e-6, atol=1e-6)

    def test_pbvv_gradient(self):
        x = np.linspace(-4, 4, 8)[:,None]
        eta = np.linspace(-10, 10, 5)[None,:]

        p = special.pbvv(eta, x)
        eps = 1e-7 + 1e-7*abs(x)
        dp = (special.pbvv(eta, x + eps)[0] - special.pbvv(eta, x - eps)[0]) / eps / 2.
        assert_tol_equal(p[1], dp, rtol=1e-6, atol=1e-6)


class TestPolygamma(TestCase):
    # from Table 6.2 (pg. 271) of A&S
    def test_polygamma(self):
        poly2 = special.polygamma(2,1)
        poly3 = special.polygamma(3,1)
        assert_almost_equal(poly2,-2.4041138063,10)
        assert_almost_equal(poly3,6.4939394023,10)

        # Test polygamma(0, x) == psi(x)
        x = [2, 3, 1.1e14]
        assert_almost_equal(special.polygamma(0, x), special.psi(x))

        # Test broadcasting
        n = [0, 1, 2]
        x = [0.5, 1.5, 2.5]
        expected = [-1.9635100260214238, 0.93480220054467933,
                    -0.23620405164172739]
        assert_almost_equal(special.polygamma(n, x), expected)
        expected = np.row_stack([expected]*2)
        assert_almost_equal(special.polygamma(n, np.row_stack([x]*2)),
                            expected)
        assert_almost_equal(special.polygamma(np.row_stack([n]*2), x),
                            expected)


class TestProCvSeq(TestCase):
    def test_pro_cv_seq(self):
        prol = special.pro_cv_seq(0,3,1)
        assert_array_almost_equal(prol,array([0.319000,
                                               2.593084,
                                               6.533471,
                                               12.514462]),5)


class TestPsi(TestCase):
    def test_psi(self):
        ps = special.psi(1)
        assert_almost_equal(ps,-0.57721566490153287,8)


class TestRadian(TestCase):
    def test_radian(self):
        rad = special.radian(90,0,0)
        assert_almost_equal(rad,pi/2.0,5)

    def test_radianmore(self):
        rad1 = special.radian(90,1,60)
        assert_almost_equal(rad1,pi/2+0.0005816135199345904,5)


class TestRiccati(TestCase):
    def test_riccati_jn(self):
        jnrl = (special.sph_jn(1,.2)[0]*.2,special.sph_jn(1,.2)[0]+special.sph_jn(1,.2)[1]*.2)
        ricjn = special.riccati_jn(1,.2)
        assert_array_almost_equal(ricjn,jnrl,8)

    def test_riccati_yn(self):
        ynrl = (special.sph_yn(1,.2)[0]*.2,special.sph_yn(1,.2)[0]+special.sph_yn(1,.2)[1]*.2)
        ricyn = special.riccati_yn(1,.2)
        assert_array_almost_equal(ricyn,ynrl,8)


class TestRound(TestCase):
    def test_round(self):
        rnd = list(map(int,(special.round(10.1),special.round(10.4),special.round(10.5),special.round(10.6))))

        # Note: According to the documentation, scipy.special.round is
        # supposed to round to the nearest even number if the fractional
        # part is exactly 0.5. On some platforms, this does not appear
        # to work and thus this test may fail. However, this unit test is
        # correctly written.
        rndrl = (10,10,10,11)
        assert_array_equal(rnd,rndrl)


def test_sph_harm():
    # Tests derived from tables in
    # http://en.wikipedia.org/wiki/Table_of_spherical_harmonics
    sh = special.sph_harm
    pi = np.pi
    exp = np.exp
    sqrt = np.sqrt
    sin = np.sin
    cos = np.cos
    yield (assert_array_almost_equal, sh(0,0,0,0),
           0.5/sqrt(pi))
    yield (assert_array_almost_equal, sh(-2,2,0.,pi/4),
           0.25*sqrt(15./(2.*pi)) *
           (sin(pi/4))**2.)
    yield (assert_array_almost_equal, sh(-2,2,0.,pi/2),
           0.25*sqrt(15./(2.*pi)))
    yield (assert_array_almost_equal, sh(2,2,pi,pi/2),
           0.25*sqrt(15/(2.*pi)) *
           exp(0+2.*pi*1j)*sin(pi/2.)**2.)
    yield (assert_array_almost_equal, sh(2,4,pi/4.,pi/3.),
           (3./8.)*sqrt(5./(2.*pi)) *
           exp(0+2.*pi/4.*1j) *
           sin(pi/3.)**2. *
           (7.*cos(pi/3.)**2.-1))
    yield (assert_array_almost_equal, sh(4,4,pi/8.,pi/6.),
           (3./16.)*sqrt(35./(2.*pi)) *
           exp(0+4.*pi/8.*1j)*sin(pi/6.)**4.)


class TestSpherical(TestCase):
    def test_sph_harm(self):
        # see test_sph_harm function
        pass

    def test_sph_in(self):
        i1n = special.sph_in(1,.2)
        inp0 = (i1n[0][1])
        inp1 = (i1n[0][0] - 2.0/0.2 * i1n[0][1])
        assert_array_almost_equal(i1n[0],array([1.0066800127054699381,
                                                0.066933714568029540839]),12)
        assert_array_almost_equal(i1n[1],[inp0,inp1],12)

    def test_sph_inkn(self):
        spikn = r_[special.sph_in(1,.2) + special.sph_kn(1,.2)]
        inkn = r_[special.sph_inkn(1,.2)]
        assert_array_almost_equal(inkn,spikn,10)

    def test_sph_in_kn_order0(self):
        x = 1.
        sph_i0 = special.sph_in(0, x)
        sph_i0_expected = np.array([np.sinh(x)/x,
                                    np.cosh(x)/x-np.sinh(x)/x**2])
        assert_array_almost_equal(r_[sph_i0], sph_i0_expected)
        sph_k0 = special.sph_kn(0, x)
        sph_k0_expected = np.array([0.5*pi*exp(-x)/x,
                                    -0.5*pi*exp(-x)*(1/x+1/x**2)])
        assert_array_almost_equal(r_[sph_k0], sph_k0_expected)
        sph_i0k0 = special.sph_inkn(0, x)
        assert_array_almost_equal(r_[sph_i0+sph_k0],
                                  r_[sph_i0k0],
                                  10)

    def test_sph_jn(self):
        s1 = special.sph_jn(2,.2)
        s10 = -s1[0][1]
        s11 = s1[0][0]-2.0/0.2*s1[0][1]
        s12 = s1[0][1]-3.0/0.2*s1[0][2]
        assert_array_almost_equal(s1[0],[0.99334665397530607731,
                                      0.066400380670322230863,
                                      0.0026590560795273856680],12)
        assert_array_almost_equal(s1[1],[s10,s11,s12],12)

    def test_sph_jnyn(self):
        jnyn = r_[special.sph_jn(1,.2) + special.sph_yn(1,.2)]  # tuple addition
        jnyn1 = r_[special.sph_jnyn(1,.2)]
        assert_array_almost_equal(jnyn1,jnyn,9)

    def test_sph_kn(self):
        kn = special.sph_kn(2,.2)
        kn0 = -kn[0][1]
        kn1 = -kn[0][0]-2.0/0.2*kn[0][1]
        kn2 = -kn[0][1]-3.0/0.2*kn[0][2]
        assert_array_almost_equal(kn[0],[6.4302962978445670140,
                                         38.581777787067402086,
                                         585.15696310385559829],12)
        assert_array_almost_equal(kn[1],[kn0,kn1,kn2],9)

    def test_sph_yn(self):
        sy1 = special.sph_yn(2,.2)[0][2]
        sy2 = special.sph_yn(0,.2)[0][0]
        sphpy = (special.sph_yn(1,.2)[0][0]-2*special.sph_yn(2,.2)[0][2])/3  # correct derivative value
        assert_almost_equal(sy1,-377.52483,5)  # previous values in the system
        assert_almost_equal(sy2,-4.9003329,5)
        sy3 = special.sph_yn(1,.2)[1][1]
        assert_almost_equal(sy3,sphpy,4)  # compare correct derivative val. (correct =-system val).


class TestStruve(object):
    def _series(self, v, z, n=100):
        """Compute Struve function & error estimate from its power series."""
        k = arange(0, n)
        r = (-1)**k * (.5*z)**(2*k+v+1)/special.gamma(k+1.5)/special.gamma(k+v+1.5)
        err = abs(r).max() * finfo(float_).eps * n
        return r.sum(), err

    def test_vs_series(self):
        """Check Struve function versus its power series"""
        for v in [-20, -10, -7.99, -3.4, -1, 0, 1, 3.4, 12.49, 16]:
            for z in [1, 10, 19, 21, 30]:
                value, err = self._series(v, z)
                assert_tol_equal(special.struve(v, z), value, rtol=0, atol=err), (v, z)

    def test_some_values(self):
        assert_tol_equal(special.struve(-7.99, 21), 0.0467547614113, rtol=1e-7)
        assert_tol_equal(special.struve(-8.01, 21), 0.0398716951023, rtol=1e-8)
        assert_tol_equal(special.struve(-3.0, 200), 0.0142134427432, rtol=1e-12)
        assert_tol_equal(special.struve(-8.0, -41), 0.0192469727846, rtol=1e-11)
        assert_equal(special.struve(-12, -41), -special.struve(-12, 41))
        assert_equal(special.struve(+12, -41), -special.struve(+12, 41))
        assert_equal(special.struve(-11, -41), +special.struve(-11, 41))
        assert_equal(special.struve(+11, -41), +special.struve(+11, 41))

        assert_(isnan(special.struve(-7.1, -1)))
        assert_(isnan(special.struve(-10.1, -1)))

    def test_regression_679(self):
        """Regression test for #679"""
        assert_tol_equal(special.struve(-1.0, 20 - 1e-8), special.struve(-1.0, 20 + 1e-8))
        assert_tol_equal(special.struve(-2.0, 20 - 1e-8), special.struve(-2.0, 20 + 1e-8))
        assert_tol_equal(special.struve(-4.3, 20 - 1e-8), special.struve(-4.3, 20 + 1e-8))


def test_chi2_smalldf():
    assert_almost_equal(special.chdtr(0.6,3), 0.957890536704110)


def test_chi2c_smalldf():
    assert_almost_equal(special.chdtrc(0.6,3), 1-0.957890536704110)


def test_chi2_inv_smalldf():
    assert_almost_equal(special.chdtri(0.6,1-0.957890536704110), 3)


def test_agm_simple():
    assert_allclose(special.agm(24, 6), 13.4581714817)
    assert_allclose(special.agm(1e30, 1), 2.2292230559453832047768593e28)


def test_legacy():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)

        # Legacy behavior: truncating arguments to integers
        assert_equal(special.bdtrc(1, 2, 0.3), special.bdtrc(1.8, 2.8, 0.3))
        assert_equal(special.bdtr(1, 2, 0.3), special.bdtr(1.8, 2.8, 0.3))
        assert_equal(special.bdtri(1, 2, 0.3), special.bdtri(1.8, 2.8, 0.3))
        assert_equal(special.expn(1, 0.3), special.expn(1.8, 0.3))
        assert_equal(special.hyp2f0(1, 2, 0.3, 1), special.hyp2f0(1, 2, 0.3, 1.8))
        assert_equal(special.nbdtrc(1, 2, 0.3), special.nbdtrc(1.8, 2.8, 0.3))
        assert_equal(special.nbdtr(1, 2, 0.3), special.nbdtr(1.8, 2.8, 0.3))
        assert_equal(special.nbdtri(1, 2, 0.3), special.nbdtri(1.8, 2.8, 0.3))
        assert_equal(special.pdtrc(1, 0.3), special.pdtrc(1.8, 0.3))
        assert_equal(special.pdtr(1, 0.3), special.pdtr(1.8, 0.3))
        assert_equal(special.pdtri(1, 0.3), special.pdtri(1.8, 0.3))
        assert_equal(special.kn(1, 0.3), special.kn(1.8, 0.3))
        assert_equal(special.yn(1, 0.3), special.yn(1.8, 0.3))
        assert_equal(special.smirnov(1, 0.3), special.smirnov(1.8, 0.3))
        assert_equal(special.smirnovi(1, 0.3), special.smirnovi(1.8, 0.3))


@with_special_errors
def test_error_raising():
    assert_raises(special.SpecialFunctionWarning, special.iv, 1, 1e99j)


def test_xlogy():
    def xfunc(x, y):
        if x == 0 and not np.isnan(y):
            return x
        else:
            return x*np.log(y)

    z1 = np.asarray([(0,0), (0, np.nan), (0, np.inf), (1.0, 2.0)], dtype=float)
    z2 = np.r_[z1, [(0, 1j), (1, 1j)]]

    w1 = np.vectorize(xfunc)(z1[:,0], z1[:,1])
    assert_func_equal(special.xlogy, w1, z1, rtol=1e-13, atol=1e-13)
    w2 = np.vectorize(xfunc)(z2[:,0], z2[:,1])
    assert_func_equal(special.xlogy, w2, z2, rtol=1e-13, atol=1e-13)


def test_xlog1py():
    def xfunc(x, y):
        if x == 0 and not np.isnan(y):
            return x
        else:
            return x * np.log1p(y)

    z1 = np.asarray([(0,0), (0, np.nan), (0, np.inf), (1.0, 2.0),
                     (1, 1e-30)], dtype=float)
    w1 = np.vectorize(xfunc)(z1[:,0], z1[:,1])
    assert_func_equal(special.xlog1py, w1, z1, rtol=1e-13, atol=1e-13)


def test_entr():
    def xfunc(x):
        if x < 0:
            return -np.inf
        else:
            return -special.xlogy(x, x)
    values = (0, 0.5, 1.0, np.inf)
    signs = [1]
    arr = []
    for sgn, v in itertools.product(signs, values):
        arr.append(sgn * v)
    z = np.array(arr, dtype=float)
    w = np.vectorize(xfunc, otypes=[np.float64])(z)
    assert_func_equal(special.entr, w, z, rtol=1e-13, atol=1e-13)


def test_kl_div():
    def xfunc(x, y):
        if x < 0 or y < 0 or (y == 0 and x != 0):
            # extension of natural domain to preserve convexity
            return np.inf
        elif np.isposinf(x) or np.isposinf(y):
            # limits within the natural domain
            return np.inf
        elif x == 0:
            return y
        else:
            return special.xlogy(x, x) - special.xlogy(x, y) - x + y
    values = (0, 0.5, 1.0)
    signs = [1]
    arr = []
    for sgna, va, sgnb, vb in itertools.product(signs, values, signs, values):
        arr.append((sgna*va, sgnb*vb))
    z = np.array(arr, dtype=float)
    w = np.vectorize(xfunc, otypes=[np.float64])(z[:,0], z[:,1])
    assert_func_equal(special.kl_div, w, z, rtol=1e-13, atol=1e-13)


def test_rel_entr():
    def xfunc(x, y):
        return special.xlogy(x, x) - special.xlogy(x, y)
    values = (0, 0.5, 1.0)
    signs = [1]
    arr = []
    for sgna, va, sgnb, vb in itertools.product(signs, values, signs, values):
        arr.append((sgna*va, sgnb*vb))
    z = np.array(arr, dtype=float)
    w = np.vectorize(xfunc, otypes=[np.float64])(z[:,0], z[:,1])
    assert_func_equal(special.rel_entr, w, z, rtol=1e-13, atol=1e-13)


if __name__ == "__main__":
    run_module_suite()
