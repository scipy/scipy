from __future__ import division, print_function, absolute_import

import warnings

import numpy as np
from numpy.testing import (TestCase, assert_array_almost_equal,
                           assert_array_equal, assert_raises, assert_equal,
                           assert_, run_module_suite)

from scipy.signal import (tf2zpk, zpk2tf, BadCoefficients, freqz, normalize,
                          buttord, cheby1, cheby2, ellip, cheb1ord, cheb2ord,
                          ellipord, butter, bessel, buttap, besselap,
                          cheb1ap, cheb2ap, ellipap)


class TestTf2zpk(TestCase):

    def test_simple(self):
        z_r = np.array([0.5, -0.5])
        p_r = np.array([1.j / np.sqrt(2), -1.j / np.sqrt(2)])
        # Sort the zeros/poles so that we don't fail the test if the order
        # changes
        z_r.sort()
        p_r.sort()
        b = np.poly(z_r)
        a = np.poly(p_r)

        z, p, k = tf2zpk(b, a)
        z.sort()
        p.sort()
        assert_array_almost_equal(z, z_r)
        assert_array_almost_equal(p, p_r)

    def test_bad_filter(self):
        """Regression test for #651: better handling of badly conditioned
        filter coefficients."""
        warnings.simplefilter("error", BadCoefficients)
        try:
            assert_raises(BadCoefficients, tf2zpk, [1e-15], [1.0, 1.0])
        finally:
            warnings.simplefilter("always", BadCoefficients)


class TestZpk2Tf(TestCase):

    def test_identity(self):
        """Test the identity transfer function."""
        z = []
        p = []
        k = 1.
        b, a = zpk2tf(z, p, k)
        b_r = np.array([1.])  # desired result
        a_r = np.array([1.])  # desired result
        # The test for the *type* of the return values is a regression
        # test for ticket #1095.  In the case p=[], zpk2tf used to
        # return the scalar 1.0 instead of array([1.0]).
        assert_array_equal(b, b_r)
        assert_(isinstance(b, np.ndarray))
        assert_array_equal(a, a_r)
        assert_(isinstance(a, np.ndarray))


class TestFreqz(TestCase):

    def test_ticket1441(self):
        """Regression test for ticket 1441."""
        # Because freqz previously used arange instead of linspace,
        # when N was large, it would return one more point than
        # requested.
        N = 100000
        w, h = freqz([1.0], worN=N)
        assert_equal(w.shape, (N,))

    def test_basic(self):
        w, h = freqz([1.0], worN=8)
        assert_array_almost_equal(w, np.pi * np.arange(8.0) / 8)
        assert_array_almost_equal(h, np.ones(8))

    def test_basic_whole(self):
        w, h = freqz([1.0], worN=8, whole=True)
        assert_array_almost_equal(w, 2 * np.pi * np.arange(8.0) / 8)
        assert_array_almost_equal(h, np.ones(8))

    def test_plot(self):

        def plot(w, h):
            assert_array_almost_equal(w, np.pi * np.arange(8.0) / 8)
            assert_array_almost_equal(h, np.ones(8))

        assert_raises(ZeroDivisionError,
                      freqz, [1.0], worN=8, plot=lambda w, h: 1 / 0)
        freqz([1.0], worN=8, plot=plot)


class TestNormalize(TestCase):

    def test_allclose(self):
        """Test for false positive on allclose in normalize() in
        filter_design.py"""
        # Test to make sure the allclose call within signal.normalize does not
        # choose false positives. Then check against a known output from MATLAB
        # to make sure the fix doesn't break anything.

        # These are the coefficients returned from
        #   `[b,a] = cheby1(8, 0.5, 0.048)'
        # in MATLAB. There are at least 15 significant figures in each
        # coefficient, so it makes sense to test for errors on the order of
        # 1e-13 (this can always be relaxed if different platforms have
        # different rounding errors)
        b_matlab = np.array([2.150733144728282e-11, 1.720586515782626e-10,
                             6.022052805239190e-10, 1.204410561047838e-09,
                             1.505513201309798e-09, 1.204410561047838e-09,
                             6.022052805239190e-10, 1.720586515782626e-10,
                             2.150733144728282e-11])
        a_matlab = np.array([1.000000000000000e+00, -7.782402035027959e+00,
                             2.654354569747454e+01, -5.182182531666387e+01,
                             6.334127355102684e+01, -4.963358186631157e+01,
                             2.434862182949389e+01, -6.836925348604676e+00,
                             8.412934944449140e-01])

        # This is the input to signal.normalize after passing through the
        # equivalent steps in signal.iirfilter as was done for MATLAB
        b_norm_in = np.array([1.5543135865293012e-06, 1.2434508692234413e-05,
                              4.3520780422820447e-05, 8.7041560845640893e-05,
                              1.0880195105705122e-04, 8.7041560845640975e-05,
                              4.3520780422820447e-05, 1.2434508692234413e-05,
                              1.5543135865293012e-06])
        a_norm_in = np.array([7.2269025909127173e+04, -5.6242661430467968e+05,
                              1.9182761917308895e+06, -3.7451128364682454e+06,
                              4.5776121393762771e+06, -3.5869706138592605e+06,
                              1.7596511818472347e+06, -4.9409793515707983e+05,
                              6.0799461347219651e+04])

        b_output, a_output = normalize(b_norm_in, a_norm_in)

        # The test on b works for decimal=14 but the one for a does not. For
        # the sake of consistency, both of these are decimal=13. If something
        # breaks on another platform, it is probably fine to relax this lower.
        assert_array_almost_equal(b_matlab, b_output, decimal=13)
        assert_array_almost_equal(a_matlab, a_output, decimal=13)


class TestBesselap(TestCase):

    def test_output_type(self):
        # Should consistently output arrays, not lists
        # https://github.com/scipy/scipy/pull/441
        for N in range(7):
            z, p, k = besselap(N)
            assert_(isinstance(z, np.ndarray))
            assert_(isinstance(p, np.ndarray))


class TestButtap(TestCase):

    def test_output_type(self):
        # Should consistently output arrays, not lists
        # https://github.com/scipy/scipy/pull/441
        for N in range(7):
            z, p, k = buttap(N)
            assert_(isinstance(z, np.ndarray))
            assert_(isinstance(p, np.ndarray))


class TestCheb1ap(TestCase):

    def test_output_type(self):
        # Should consistently output arrays, not lists
        # https://github.com/scipy/scipy/pull/441
        for N in range(7):
            z, p, k = cheb1ap(N, 1)
            assert_(isinstance(z, np.ndarray))
            assert_(isinstance(p, np.ndarray))


class TestCheb2ap(TestCase):

    def test_output_type(self):
        # Should consistently output arrays, not lists
        # https://github.com/scipy/scipy/pull/441
        for N in range(7):
            z, p, k = cheb2ap(N, 20)
            assert_(isinstance(z, np.ndarray))
            assert_(isinstance(p, np.ndarray))


class TestEllipap(TestCase):

    def test_output_type(self):
        # Should consistently output arrays, not lists
        # https://github.com/scipy/scipy/pull/441
        for N in range(7):
            z, p, k = ellipap(N, 1, 20)
            assert_(isinstance(z, np.ndarray))
            assert_(isinstance(p, np.ndarray))


class TestButtord(TestCase):

    def test_basic(self):
        # From http://dsp.etfbl.net/filtri/aproksimacija
        assert_array_almost_equal(buttord(1, 550/450, 1, 26, analog=True),
                                  [19, 1.0441379169150726])

        assert_equal(buttord(1, 1.2, 1, 80, analog=True)[0], 55)

        # From http://www.mathworks.com/help/signal/ref/buttord.html
        n, Wn = buttord(40/500, 150/500, 3, 60)
        assert_equal(n, 5)
        assert_array_almost_equal(Wn, 0.0810, decimal=4)

        n, Wn = buttord([60/500, 200/500], [50/500, 250/500], 3, 40)
        assert_equal(n, 16)
        assert_array_almost_equal(Wn, [0.1198, 0.4005], decimal=4)


class TestCheb1ord(TestCase):

    def test_basic(self):
        # From http://dsp.etfbl.net/filtri/aproksimacija
        assert_equal(cheb1ord(1, 1.2, 1, 80, analog=True)[0], 17)

        # From http://www.mathworks.com/help/signal/ref/cheb1ord.html
        n, Wp = cheb1ord(40/500, 150/500, 3, 60)
        assert_equal(n, 4)
        assert_array_almost_equal(Wp, 0.0800, decimal=4)

        n, Wp = cheb1ord([60/500, 200/500], [50/500, 250/500], 3, 40)
        assert_equal(n, 7)
        assert_array_almost_equal(Wp, [0.1200, 0.4000], decimal=4)


class TestCheb2ord(TestCase):

    def test_basic(self):
        # From http://www.mathworks.com/help/signal/ref/cheb2ord.html
        n, Ws = cheb2ord(40/500, 150/500, 3, 60)
        assert_equal(n, 4)
        # TODO: This fails.  Matlab docs wrong?
        # assert_array_almost_equal(Ws, 0.3000, decimal=4)

        n, Ws = cheb2ord([60/500, 200/500], [50/500, 250/500], 3, 40)
        assert_equal(n, 7)
        # TODO: This fails.  Matlab docs wrong?
        # assert_array_almost_equal(Ws, [0.1000, 0.5000], decimal=4)

        # http://cens.ioc.ee/local/man/matlab/toolbox/signal/cheb2ord.html
        n, Wn = cheb2ord(100/500, 150/500, 3, 15)
        assert_equal(n, 3)
        assert_array_almost_equal(Wn, 0.2609, decimal=4)

        n, Wn = cheb2ord([100/500, 200/500], [50/500, 250/500], 3, 30)
        assert_equal(n, 4)
        assert_array_almost_equal(Wn, [0.1633, 0.4665], decimal=4)


class TestEllipord(TestCase):

    def test_basic(self):
        # From http://dsp.etfbl.net/filtri/aproksimacija
        assert_equal(ellipord(1, 1.2, 1, 80, analog=True)[0], 9)

        # From http://www.mathworks.com/help/signal/ref/ellipord.html
        n, Wp = ellipord(40/500, 150/500, 3, 60)
        assert_equal(n, 4)
        assert_array_almost_equal(Wp, 0.0800, decimal=4)

        n, Wp = ellipord([60/500, 200/500], [50/500, 250/500], 3, 40)
        assert_equal(n, 5)
        assert_array_almost_equal(Wp, [0.1200, 0.4000], decimal=4)


class TestBessel(TestCase):

    def test_degenerate(self):
        # 0-order filter is just a passthrough
        b, a = bessel(0, 1, analog=True)
        assert_array_equal(b, [1])
        assert_array_equal(a, [1])

        # 1-order filter is same for all types
        b, a = bessel(1, 1, analog=True)
        assert_array_almost_equal(b, [1])
        assert_array_almost_equal(a, [1, 1])


class TestButter(TestCase):

    def test_degenerate(self):
        # 0-order filter is just a passthrough
        b, a = butter(0, 1, analog=True)
        assert_array_equal(b, [1])
        assert_array_equal(a, [1])

        # 1-order filter is same for all types
        b, a = butter(1, 1, analog=True)
        assert_array_almost_equal(b, [1])
        assert_array_almost_equal(a, [1, 1])

    def test_basic(self):
        # Requires https://github.com/scipy/scipy/pull/3085 to pass
        for N in range(25):
            wn = 0.01
            z, p, k = butter(N, wn, 'low', analog=True, output='zpk')
            assert_array_almost_equal([], z)
            assert_(len(p) == N)
            # All poles should be at distance wn from origin
            assert_array_almost_equal(wn, abs(p))
            assert_(all(np.real(p) <= 0)) # No poles in right half of S-plane
            assert_array_almost_equal(wn**N, k)

        for N in range(25):
            wn = 0.01
            z, p, k = butter(N, wn, 'high', analog=False, output='zpk')
            assert_array_equal(np.ones(N), z) # All zeros exactly at DC
            assert_(all(np.abs(p) <= 1)) # No poles outside unit circle

        # From http://dsp.etfbl.net/filtri/aproksimacija
        b1, a1 = butter(2, 1, analog=True)
        assert_array_almost_equal(b1, [1])
        assert_array_almost_equal(a1, [1, np.sqrt(2), 1])

        b2, a2 = butter(5, 1, analog=True)
        assert_array_almost_equal(b2, [1])
        assert_array_almost_equal(a2, [1, 3.2361, 5.2361,
                                       5.2361, 3.2361, 1], decimal=4)

        b3, a3 = butter(10, 1, analog=True)
        assert_array_almost_equal(b3, [1])
        assert_array_almost_equal(a3, [1, 6.3925, 20.4317, 42.8021, 64.8824,
                                       74.2334, 64.8824, 42.8021, 20.4317,
                                       6.3925, 1], decimal=4)

        b2, a2 = butter(19, 1.0441379169150726, analog=True)
        assert_array_almost_equal(b2, [2.2720], decimal=4)
        assert_array_almost_equal(a2, 1.0e+004 * np.array([
                        0.0001, 0.0013, 0.0080, 0.0335, 0.1045, 0.2570,
                        0.5164, 0.8669, 1.2338, 1.5010, 1.5672, 1.4044,
                        1.0759, 0.6986, 0.3791, 0.1681, 0.0588, 0.0153,
                        0.0026, 0.0002]), decimal=0)

        # From https://sites.google.com/site/nahums84/research/digital-filters---matlab
        b, a = butter(5, 0.4)
        assert_array_almost_equal(b, [0.0219, 0.1097, 0.2194,
                                      0.2194, 0.1097, 0.0219], decimal=4)
        assert_array_almost_equal(a, [1.0000, -0.9853, 0.9738,
                                     -0.3864, 0.1112, -0.0113], decimal=4)


class TestCheby1(TestCase):

    def test_degenerate(self):
        # 0-order filter is just a passthrough
        # Even-order filters have DC gain of -rp dB
        b, a = cheby1(0, 10*np.log10(2), 1, analog=True)
        assert_array_almost_equal(b, [1/np.sqrt(2)])
        assert_array_equal(a, [1])

        # 1-order filter is same for all types
        b, a = cheby1(1, 10*np.log10(2), 1, analog=True)
        assert_array_almost_equal(b, [1])
        assert_array_almost_equal(a, [1, 1])

    def test_basic(self):
        # Requires https://github.com/scipy/scipy/pull/3085 to pass
        for N in range(25):
            wn = 0.01
            z, p, k = cheby1(N, 1, wn, 'low', analog=True, output='zpk')
            assert_array_almost_equal([], z)
            assert_(len(p) == N)
            assert_(all(np.real(p) <= 0)) # No poles in right half of S-plane

        for N in range(25):
            wn = 0.01
            z, p, k = cheby1(N, 1, wn, 'high', analog=False, output='zpk')
            assert_array_equal(np.ones(N), z) # All zeros exactly at DC
            assert_(all(np.abs(p) <= 1)) # No poles outside unit circle

        # From TestNormalize
        b, a = cheby1(8, 0.5, 0.048)
        assert_array_almost_equal(b, np.array([
                             2.150733144728282e-11, 1.720586515782626e-10,
                             6.022052805239190e-10, 1.204410561047838e-09,
                             1.505513201309798e-09, 1.204410561047838e-09,
                             6.022052805239190e-10, 1.720586515782626e-10,
                             2.150733144728282e-11]), decimal=14)
        assert_array_almost_equal(a, np.array([
                             1.000000000000000e+00, -7.782402035027959e+00,
                             2.654354569747454e+01, -5.182182531666387e+01,
                             6.334127355102684e+01, -4.963358186631157e+01,
                             2.434862182949389e+01, -6.836925348604676e+00,
                             8.412934944449140e-01]), decimal=14)

        # From https://sites.google.com/site/vandankeuth/lab7-matlab
        b, a = cheby1(4, 1, [0.4, 0.7], btype='band')
        assert_array_almost_equal(b, [0.0084, 0, -0.0335, 0, 0.0502, 0,
                                      -0.0335, 0, 0.0084], decimal=4)
        assert_array_almost_equal(a, [1.0, 1.1191, 2.862, 2.2986, 3.4137,
                                      1.8653, 1.8982, 0.5676, 0.4103],
                                      decimal=4)

        # From http://dsp.etfbl.net/filtri/aproksimacija
        b2, a2 = cheby1(5, 3, 1, analog=True)
        assert_array_almost_equal(b2, [0.0626], decimal=4)
        assert_array_almost_equal(a2, [1, 0.5745, 1.4150, 0.5489, 0.4080,
                                       0.0626], decimal=4)

        # From http://searchcode.com/codesearch/raw/27092032
        b, a = cheby1(8, 0.5, 0.1)
        assert_array_almost_equal(b, 1.0e-006 * np.array([
             0.00703924326028, 0.05631394608227, 0.19709881128793,
             0.39419762257586, 0.49274702821983, 0.39419762257586,
             0.19709881128793, 0.05631394608227, 0.00703924326028]),
             decimal=13)
        assert_array_almost_equal(a, [
             1.00000000000000, -7.44912258934158,  24.46749067762108,
           -46.27560200466141, 55.11160187999928, -42.31640010161038,
            20.45543300484147, -5.69110270561444,   0.69770374759022],
            decimal=13)

        b, a = cheby1(8, 0.5, 0.25)
        assert_array_almost_equal(b, 1.0e-003 * np.array([
            0.00895261138923, 0.07162089111382, 0.25067311889837,
            0.50134623779673, 0.62668279724591, 0.50134623779673,
            0.25067311889837, 0.07162089111382, 0.00895261138923]),
            decimal=13)
        assert_array_almost_equal(a, [
            1.00000000000000, -5.97529229188545, 16.58122329202101,
            -27.71423273542923, 30.39509758355313, -22.34729670426879,
            10.74509800434910, -3.08924633697497, 0.40707685889802],
            decimal=13)


class TestCheby2(TestCase):

    def test_degenerate(self):
        # 0-order filter is just a passthrough
        # Stopband ripple factor doesn't matter
        b, a = cheby2(0, 123.456, 1, analog=True)
        assert_array_equal(b, [1])
        assert_array_equal(a, [1])

        # 1-order filter is same for all types
        b, a = cheby2(1, 10*np.log10(2), 1, analog=True)
        assert_array_almost_equal(b, [1])
        assert_array_almost_equal(a, [1, 1])

    def test_basic(self):
        # Requires https://github.com/scipy/scipy/pull/3085 to pass
        for N in range(25):
            wn = 0.01
            z, p, k = cheby2(N, 40, wn, 'low', analog=True, output='zpk')
            assert_(len(p) == N)
            assert_(all(np.real(p) <= 0)) # No poles in right half of S-plane

        for N in range(25):
            wn = 0.01
            z, p, k = cheby2(N, 40, wn, 'high', analog=False, output='zpk')
            assert_(all(np.abs(p) <= 1)) # No poles outside unit circle

        # From http://www.dsprelated.com/showmessage/20207/1.php
        B, A = cheby2(18, 100, 0.5)
        assert_array_almost_equal(B, [
            0.00167583914216, 0.01249479541868, 0.05282702120282,
            0.15939804265706, 0.37690207631117, 0.73227013789108,
            1.20191856962356, 1.69522872823393, 2.07598674519837,
            2.21972389625291, 2.07598674519838, 1.69522872823395,
            1.20191856962359, 0.73227013789110, 0.37690207631118,
            0.15939804265707, 0.05282702120282, 0.01249479541868,
            0.00167583914216], decimal=13)
        assert_array_almost_equal(A, [
            1.00000000000000, -0.27631970006174, 3.19751214254060,
            -0.15685969461355, 4.13926117356269, 0.60689917820044,
            2.95082770636540, 0.89016501910416, 1.32135245849798,
            0.51502467236824, 0.38906643866660, 0.15367372690642,
            0.07255803834919, 0.02422454070134, 0.00756108751837,
            0.00179848550988, 0.00033713574499, 0.00004258794833,
            0.00000281030149], decimal=13)


class TestEllip(TestCase):

    def test_degenerate(self):
        # 0-order filter is just a passthrough
        # Even-order filters have DC gain of -rp dB
        # Stopband ripple factor doesn't matter
        b, a = ellip(0, 10*np.log10(2), 123.456, 1, analog=True)
        assert_array_almost_equal(b, [1/np.sqrt(2)])
        assert_array_equal(a, [1])

        # 1-order filter is same for all types
        b, a = ellip(1, 10*np.log10(2), 1, 1, analog=True)
        assert_array_almost_equal(b, [1])
        assert_array_almost_equal(a, [1, 1])

    def test_basic(self):
        # Requires https://github.com/scipy/scipy/pull/3085 to pass
        for N in range(25):
            wn = 0.01
            z, p, k = ellip(N, 1, 40, wn, 'low', analog=True, output='zpk')
            assert_(len(p) == N)
            assert_(all(np.real(p) <= 0)) # No poles in right half of S-plane

        for N in range(25):
            wn = 0.01
            z, p, k = ellip(N, 1, 40, wn, 'high', analog=False, output='zpk')
            assert_(all(np.abs(p) <= 1)) # No poles outside unit circle

        # From http://dsp.etfbl.net/filtri/aproksimacija
        b3, a3 = ellip(5, 3, 26, 1, analog=True)
        assert_array_almost_equal(b3, [0.1420, 0, 0.3764, 0,
                                       0.2409], decimal=4)
        assert_array_almost_equal(a3, [1, 0.5686, 1.8061, 0.8017, 0.8012,
                                       0.2409], decimal=4)

        # From https://sites.google.com/site/michaeltsessa/lab-7
        b, a = ellip(3, 1, 60, [0.4, 0.7], 'stop')
        assert_array_almost_equal(b, [0.3310, 0.3469, 1.1042, 0.7044, 1.1042,
                                      0.3469, 0.3310], decimal=4)
        assert_array_almost_equal(a, [1.0000, 0.6973, 1.1441, 0.5878, 0.7323,
                                      0.1131, -0.0060], decimal=4)


if __name__ == "__main__":
    run_module_suite()
