from __future__ import division, print_function, absolute_import

import warnings

import numpy as np
from numpy.testing import (TestCase, assert_array_almost_equal,
                           assert_array_equal, assert_array_less,
                           assert_raises, assert_equal, assert_,
                           run_module_suite, assert_allclose, assert_warns)
from numpy import array, spacing, sin, pi, sort

from scipy.signal import (tf2zpk, zpk2tf, tf2sos, sos2tf, sos2zpk, zpk2sos,
                          BadCoefficients, freqz, normalize,
                          buttord, cheby1, cheby2, ellip, cheb1ord, cheb2ord,
                          ellipord, butter, bessel, buttap, besselap,
                          cheb1ap, cheb2ap, ellipap, iirfilter, freqs,
                          lp2lp, lp2hp, lp2bp, lp2bs, bilinear, group_delay,
                          firwin)
from scipy.signal.filter_design import _cplxreal, _cplxpair


class TestCplxPair(TestCase):

    def test_trivial_input(self):
        assert_equal(_cplxpair([]).size, 0)
        assert_equal(_cplxpair(1), 1)

    def test_output_order(self):
        assert_allclose(_cplxpair([1+1j, 1-1j]), [1-1j, 1+1j])

        a = [1+1j, 1+1j, 1, 1-1j, 1-1j, 2]
        b = [1-1j, 1+1j, 1-1j, 1+1j, 1, 2]
        assert_allclose(_cplxpair(a), b)

        # points spaced around the unit circle
        z = np.exp(2j*pi*array([4, 3, 5, 2, 6, 1, 0])/7)
        z1 = np.copy(z)
        np.random.shuffle(z)
        assert_allclose(_cplxpair(z), z1)
        np.random.shuffle(z)
        assert_allclose(_cplxpair(z), z1)
        np.random.shuffle(z)
        assert_allclose(_cplxpair(z), z1)

        # Should be able to pair up all the conjugates
        x = np.random.rand(10000) + 1j * np.random.rand(10000)
        y = x.conj()
        z = np.random.rand(10000)
        x = np.concatenate((x, y, z))
        np.random.shuffle(x)
        c = _cplxpair(x)

        # Every other element of head should be conjugates:
        assert_allclose(c[0:20000:2], np.conj(c[1:20000:2]))
        # Real parts of head should be in sorted order:
        assert_allclose(c[0:20000:2].real, np.sort(c[0:20000:2].real))
        # Tail should be sorted real numbers:
        assert_allclose(c[20000:], np.sort(c[20000:]))

    def test_real_integer_input(self):
        assert_array_equal(_cplxpair([2, 0, 1]), [0, 1, 2])

    def test_tolerances(self):
        eps = spacing(1)
        assert_allclose(_cplxpair([1j, -1j, 1+1j*eps], tol=2*eps),
                        [-1j, 1j, 1+1j*eps])

        # sorting close to 0
        assert_allclose(_cplxpair([-eps+1j, +eps-1j]), [-1j, +1j])
        assert_allclose(_cplxpair([+eps+1j, -eps-1j]), [-1j, +1j])
        assert_allclose(_cplxpair([+1j, -1j]), [-1j, +1j])

    def test_unmatched_conjugates(self):
        # 1+2j is unmatched
        assert_raises(ValueError, _cplxpair, [1+3j, 1-3j, 1+2j])

        # 1+2j and 1-3j are unmatched
        assert_raises(ValueError, _cplxpair, [1+3j, 1-3j, 1+2j, 1-3j])

        # 1+3j is unmatched
        assert_raises(ValueError, _cplxpair, [1+3j, 1-3j, 1+3j])

        # Not conjugates
        assert_raises(ValueError, _cplxpair, [4+5j, 4+5j])
        assert_raises(ValueError, _cplxpair, [1-7j, 1-7j])

        # No pairs
        assert_raises(ValueError, _cplxpair, [1+3j])
        assert_raises(ValueError, _cplxpair, [1-3j])


class TestCplxReal(TestCase):

    def test_trivial_input(self):
        assert_equal(_cplxreal([]), ([], []))
        assert_equal(_cplxreal(1), ([], [1]))

    def test_output_order(self):
        zc, zr = _cplxreal(np.roots(array([1, 0, 0, 1])))
        assert_allclose(np.append(zc, zr), [1/2 + 1j*sin(pi/3), -1])

        eps = spacing(1)

        a = [0+1j, 0-1j, eps + 1j, eps - 1j, -eps + 1j, -eps - 1j,
             1, 4, 2, 3, 0, 0,
             2+3j, 2-3j,
             1-eps + 1j, 1+2j, 1-2j, 1+eps - 1j,  # sorts out of order
             3+1j, 3+1j, 3+1j, 3-1j, 3-1j, 3-1j,
             2-3j, 2+3j]
        zc, zr = _cplxreal(a)
        assert_allclose(zc, [1j, 1j, 1j, 1+1j, 1+2j, 2+3j, 2+3j, 3+1j, 3+1j,
                             3+1j])
        assert_allclose(zr, [0, 0, 1, 2, 3, 4])

        z = array([1-eps + 1j, 1+2j, 1-2j, 1+eps - 1j, 1+eps+3j, 1-2*eps-3j,
                   0+1j, 0-1j, 2+4j, 2-4j, 2+3j, 2-3j, 3+7j, 3-7j, 4-eps+1j,
                   4+eps-2j, 4-1j, 4-eps+2j])

        zc, zr = _cplxreal(z)
        assert_allclose(zc, [1j, 1+1j, 1+2j, 1+3j, 2+3j, 2+4j, 3+7j, 4+1j,
                             4+2j])
        assert_equal(zr, [])

    def test_unmatched_conjugates(self):
        # 1+2j is unmatched
        assert_raises(ValueError, _cplxreal, [1+3j, 1-3j, 1+2j])

        # 1+2j and 1-3j are unmatched
        assert_raises(ValueError, _cplxreal, [1+3j, 1-3j, 1+2j, 1-3j])

        # 1+3j is unmatched
        assert_raises(ValueError, _cplxreal, [1+3j, 1-3j, 1+3j])

        # No pairs
        assert_raises(ValueError, _cplxreal, [1+3j])
        assert_raises(ValueError, _cplxreal, [1-3j])

    def test_real_integer_input(self):
        zc, zr = _cplxreal([2, 0, 1, 4])
        assert_array_equal(zc, [])
        assert_array_equal(zr, [0, 1, 2, 4])


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


class TestSos2Zpk(TestCase):

    def test_basic(self):
        sos = [[1, 0, 1, 1, 0, -0.81],
               [1, 0, 0, 1, 0, +0.49]]
        z, p, k = sos2zpk(sos)
        z2 = [1j, -1j, 0, 0]
        p2 = [0.9, -0.9, 0.7j, -0.7j]
        k2 = 1
        assert_array_almost_equal(sort(z), sort(z2), decimal=4)
        assert_array_almost_equal(sort(p), sort(p2), decimal=4)
        assert_array_almost_equal(k, k2)

        sos = [[1.00000, +0.61803, 1.0000, 1.00000, +0.60515, 0.95873],
               [1.00000, -1.61803, 1.0000, 1.00000, -1.58430, 0.95873],
               [1.00000, +1.00000, 0.0000, 1.00000, +0.97915, 0.00000]]
        z, p, k = sos2zpk(sos)
        z2 = [-0.3090 + 0.9511j, -0.3090 - 0.9511j, 0.8090 + 0.5878j,
              0.8090 - 0.5878j, -1.0000 + 0.0000j, 0]
        p2 = [-0.3026 + 0.9312j, -0.3026 - 0.9312j, 0.7922 + 0.5755j,
              0.7922 - 0.5755j, -0.9791 + 0.0000j, 0]
        k2 = 1
        assert_array_almost_equal(sort(z), sort(z2), decimal=4)
        assert_array_almost_equal(sort(p), sort(p2), decimal=4)

        sos = array([[1, 2, 3, 1, 0.2, 0.3],
                     [4, 5, 6, 1, 0.4, 0.5]])
        z = array([-1 - 1.41421356237310j, -1 + 1.41421356237310j,
                  -0.625 - 1.05326872164704j, -0.625 + 1.05326872164704j])
        p = array([-0.2 - 0.678232998312527j, -0.2 + 0.678232998312527j,
                  -0.1 - 0.538516480713450j, -0.1 + 0.538516480713450j])
        k = 4
        z2, p2, k2 = sos2zpk(sos)
        assert_allclose(_cplxpair(z2), z)
        assert_allclose(_cplxpair(p2), p)
        assert_allclose(k2, k)


class TestSos2Tf(TestCase):

    def test_basic(self):
        sos = [[1, 1, 1, 1, 0, -1],
               [-2, 3, 1, 1, 10, 1]]
        b, a = sos2tf(sos)
        assert_array_almost_equal(b, [-2, 1, 2, 4, 1])
        assert_array_almost_equal(a, [1, 10, 0, -10, -1])


class TestTf2Sos(TestCase):

    def test_basic(self):
        num = [2, 16, 44, 56, 32]
        den = [3, 3, -15, 18, -12]
        sos = tf2sos(num, den)
        sos2 = [[0.6667, 4.0000, 5.3333, 1.0000, +2.0000, -4.0000],
                [1.0000, 2.0000, 2.0000, 1.0000, -1.0000, +1.0000]]
        assert_array_almost_equal(sos, sos2, decimal=4)

        b = [1, -3, 11, -27, 18]
        a = [16, 12, 2, -4, -1]
        sos = tf2sos(b, a)
        sos2 = [[0.0625, -0.1875, 0.1250, 1.0000, -0.2500, -0.1250],
                [1.0000, +0.0000, 9.0000, 1.0000, +1.0000, +0.5000]]
        # assert_array_almost_equal(sos, sos2, decimal=4)


class TestZpk2Sos(TestCase):

    def test_basic(self):
        for pairing in ('nearest', 'keep_odd'):
            #
            # Cases that match octave
            #

            z = [-1, -1]
            p = [0.57149 + 0.29360j, 0.57149 - 0.29360j]
            k = 1
            sos = zpk2sos(z, p, k, pairing=pairing)
            sos2 = [[1, 2, 1, 1, -1.14298, 0.41280]]  # octave & MATLAB
            assert_array_almost_equal(sos, sos2, decimal=4)

            z = [1j, -1j]
            p = [0.9, -0.9, 0.7j, -0.7j]
            k = 1
            sos = zpk2sos(z, p, k, pairing=pairing)
            sos2 = [[1, 0, 1, 1, 0, +0.49],
                    [1, 0, 0, 1, 0, -0.81]]  # octave
            # sos2 = [[0, 0, 1, 1, -0.9, 0],
            #         [1, 0, 1, 1, 0.9, 0]]  # MATLAB
            assert_array_almost_equal(sos, sos2, decimal=4)

            z = []
            p = [0.8, -0.5+0.25j, -0.5-0.25j]
            k = 1.
            sos = zpk2sos(z, p, k, pairing=pairing)
            sos2 = [[1., 0., 0., 1., 1., 0.3125],
                    [1., 0., 0., 1., -0.8, 0.]]  # octave, MATLAB fails
            assert_array_almost_equal(sos, sos2, decimal=4)

            z = [1., 1., 0.9j, -0.9j]
            p = [0.99+0.01j, 0.99-0.01j, 0.1+0.9j, 0.1-0.9j]
            k = 1
            sos = zpk2sos(z, p, k, pairing=pairing)
            sos2 = [[1, 0, 0.81, 1, -0.2, 0.82],
                    [1, -2, 1, 1, -1.98, 0.9802]]  # octave
            # sos2 = [[1, -2, 1, 1,  -0.2, 0.82],
            #         [1, 0, 0.81, 1, -1.98, 0.9802]]  # MATLAB
            assert_array_almost_equal(sos, sos2, decimal=4)

            z = [0.9+0.1j, 0.9-0.1j, -0.9]
            p = [0.75+0.25j, 0.75-0.25j, 0.9]
            k = 1
            sos = zpk2sos(z, p, k, pairing=pairing)
            if pairing == 'keep_odd':
                sos2 = [[1, -1.8, 0.82, 1, -1.5, 0.625],
                        [1, 0.9, 0, 1, -0.9, 0]]  # octave; MATLAB fails
                assert_array_almost_equal(sos, sos2, decimal=4)
            else:  # pairing == 'nearest'
                sos2 = [[1, 0.9, 0, 1, -1.5, 0.625],
                        [1, -1.8, 0.82, 1, -0.9, 0]]  # our algorithm
                assert_array_almost_equal(sos, sos2, decimal=4)

            #
            # Cases that differ from octave:
            #

            z = [-0.3090 + 0.9511j, -0.3090 - 0.9511j, 0.8090 + 0.5878j,
                 +0.8090 - 0.5878j, -1.0000 + 0.0000j]
            p = [-0.3026 + 0.9312j, -0.3026 - 0.9312j, 0.7922 + 0.5755j,
                 +0.7922 - 0.5755j, -0.9791 + 0.0000j]
            k = 1
            sos = zpk2sos(z, p, k, pairing=pairing)
            # sos2 = [[1, 0.618, 1, 1, 0.6052, 0.95870],
            #         [1, -1.618, 1, 1, -1.5844, 0.95878],
            #         [1, 1, 0, 1, 0.9791, 0]]  # octave, MATLAB fails
            sos2 = [[1, 1, 0, 1, +0.97915, 0],
                    [1, 0.61803, 1, 1, +0.60515, 0.95873],
                    [1, -1.61803, 1, 1, -1.58430, 0.95873]]
            assert_array_almost_equal(sos, sos2, decimal=4)

            z = [-1 - 1.4142j, -1 + 1.4142j,
                 -0.625 - 1.0533j, -0.625 + 1.0533j]
            p = [-0.2 - 0.6782j, -0.2 + 0.6782j,
                 -0.1 - 0.5385j, -0.1 + 0.5385j]
            k = 4
            sos = zpk2sos(z, p, k, pairing=pairing)
            sos2 = [[4, 8, 12, 1, 0.2, 0.3],
                    [1, 1.25, 1.5, 1, 0.4, 0.5]]  # MATLAB
            # sos2 = [[4, 8, 12, 1, 0.4, 0.5],
            #         [1, 1.25, 1.5, 1, 0.2, 0.3]]  # octave
            assert_allclose(sos, sos2, rtol=1e-4, atol=1e-4)

            z = []
            p = [0.2, -0.5+0.25j, -0.5-0.25j]
            k = 1.
            sos = zpk2sos(z, p, k, pairing=pairing)
            sos2 = [[1., 0., 0., 1., -0.2, 0.],
                    [1., 0., 0., 1., 1., 0.3125]]
            # sos2 = [[1., 0., 0., 1., 1., 0.3125],
            #         [1., 0., 0., 1., -0.2, 0]]  # octave, MATLAB fails
            assert_array_almost_equal(sos, sos2, decimal=4)

            # The next two examples are adapted from Leland B. Jackson,
            # "Digital Filters and Signal Processing (1995) p.400:
            # http://books.google.com/books?id=VZ8uabI1pNMC&lpg=PA400&ots=gRD9pi8Jua&dq=Pole%2Fzero%20pairing%20for%20minimum%20roundoff%20noise%20in%20BSF.&pg=PA400#v=onepage&q=Pole%2Fzero%20pairing%20for%20minimum%20roundoff%20noise%20in%20BSF.&f=false

            deg2rad = np.pi / 180.
            k = 1.

            # first example
            thetas = [22.5, 45, 77.5]
            mags = [0.8, 0.6, 0.9]
            z = np.array([np.exp(theta * deg2rad * 1j) for theta in thetas])
            z = np.concatenate((z, np.conj(z)))
            p = np.array([mag * np.exp(theta * deg2rad * 1j)
                          for theta, mag in zip(thetas, mags)])
            p = np.concatenate((p, np.conj(p)))
            sos = zpk2sos(z, p, k)
            # sos2 = [[1, -0.43288, 1, 1, -0.38959, 0.81],  # octave,
            #         [1, -1.41421, 1, 1, -0.84853, 0.36],  # MATLAB fails
            #         [1, -1.84776, 1, 1, -1.47821, 0.64]]
            # Note that pole-zero pairing matches, but ordering is different
            sos2 = [[1, -1.41421, 1, 1, -0.84853, 0.36],
                    [1, -1.84776, 1, 1, -1.47821, 0.64],
                    [1, -0.43288, 1, 1, -0.38959, 0.81]]
            assert_array_almost_equal(sos, sos2, decimal=4)

            # second example
            z = np.array([np.exp(theta * deg2rad * 1j)
                          for theta in (85., 10.)])
            z = np.concatenate((z, np.conj(z), [1, -1]))
            sos = zpk2sos(z, p, k)

            # sos2 = [[1, -0.17431, 1, 1, -0.38959, 0.81],  # octave "wrong",
            #         [1, -1.96962, 1, 1, -0.84853, 0.36],  # MATLAB fails
            #         [1, 0, -1, 1, -1.47821, 0.64000]]
            # Our pole-zero pairing matches the text, Octave does not
            sos2 = [[1, 0, -1, 1, -0.84853, 0.36],
                    [1, -1.96962, 1, 1, -1.47821, 0.64],
                    [1, -0.17431, 1, 1, -0.38959, 0.81]]
            assert_array_almost_equal(sos, sos2, decimal=4)


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


class TestLp2lp(TestCase):

    def test_basic(self):
        b = [1]
        a = [1, np.sqrt(2), 1]
        b_lp, a_lp = lp2lp(b, a, 0.38574256627112119)
        assert_array_almost_equal(b_lp, [0.1488], decimal=4)
        assert_array_almost_equal(a_lp, [1, 0.5455, 0.1488], decimal=4)


class TestLp2hp(TestCase):

    def test_basic(self):
        b = [0.25059432325190018]
        a = [1, 0.59724041654134863, 0.92834805757524175, 0.25059432325190018]
        b_hp, a_hp = lp2hp(b, a, 2*np.pi*5000)
        assert_allclose(b_hp, [1, 0, 0, 0])
        assert_allclose(a_hp, [1, 1.1638e5, 2.3522e9, 1.2373e14], rtol=1e-4)


class TestLp2bp(TestCase):

    def test_basic(self):
        b = [1]
        a = [1, 2, 2, 1]
        b_bp, a_bp = lp2bp(b, a, 2*np.pi*4000, 2*np.pi*2000)
        assert_allclose(b_bp, [1.9844e12, 0, 0, 0], rtol=1e-6)
        assert_allclose(a_bp, [1, 2.5133e4, 2.2108e9, 3.3735e13,
                               1.3965e18, 1.0028e22, 2.5202e26], rtol=1e-4)


class TestLp2bs(TestCase):

    def test_basic(self):
        b = [1]
        a = [1, 1]
        b_bs, a_bs = lp2bs(b, a, 0.41722257286366754, 0.18460575326152251)
        assert_array_almost_equal(b_bs, [1, 0, 0.17407], decimal=5)
        assert_array_almost_equal(a_bs, [1, 0.18461, 0.17407], decimal=5)


class TestBilinear(TestCase):

    def test_basic(self):
        b = [0.14879732743343033]
        a = [1, 0.54552236880522209, 0.14879732743343033]
        b_z, a_z = bilinear(b, a, 0.5)
        assert_array_almost_equal(b_z, [0.087821, 0.17564, 0.087821],
                                  decimal=5)
        assert_array_almost_equal(a_z, [1, -1.0048, 0.35606], decimal=4)

        b = [1, 0, 0.17407467530697837]
        a = [1, 0.18460575326152251, 0.17407467530697837]
        b_z, a_z = bilinear(b, a, 0.5)
        assert_array_almost_equal(b_z, [0.86413, -1.2158, 0.86413],
                                  decimal=4)
        assert_array_almost_equal(a_z, [1, -1.2158, 0.72826],
                                  decimal=4)


class TestPrototypeType(TestCase):

    def test_output_type(self):
        # Prototypes should consistently output arrays, not lists
        # https://github.com/scipy/scipy/pull/441
        for func in (buttap,
                     besselap,
                     lambda N: cheb1ap(N, 1),
                     lambda N: cheb2ap(N, 20),
                     lambda N: ellipap(N, 1, 20)):
            for N in range(7):
                z, p, k = func(N)
                assert_(isinstance(z, np.ndarray))
                assert_(isinstance(p, np.ndarray))


def dB(x):
    # Return magnitude in decibels
    return 20 * np.log10(abs(x))


class TestButtord(TestCase):

    def test_lowpass(self):
        wp = 0.2
        ws = 0.3
        rp = 3
        rs = 60
        N, Wn = buttord(wp, ws, rp, rs, False)
        b, a = butter(N, Wn, 'lowpass', False)
        w, h = freqz(b, a)
        w /= np.pi
        assert_array_less(-rp, dB(h[w <= wp]))
        assert_array_less(dB(h[ws <= w]), -rs)

        assert_equal(N, 16)
        assert_allclose(Wn, 2.0002776782743284e-01, rtol=1e-15)

    def test_highpass(self):
        wp = 0.3
        ws = 0.2
        rp = 3
        rs = 70
        N, Wn = buttord(wp, ws, rp, rs, False)
        b, a = butter(N, Wn, 'highpass', False)
        w, h = freqz(b, a)
        w /= np.pi
        assert_array_less(-rp, dB(h[wp <= w]))
        assert_array_less(dB(h[w <= ws]), -rs)

        assert_equal(N, 18)
        assert_allclose(Wn, 2.9996603079132672e-01, rtol=1e-15)

    def test_bandpass(self):
        wp = [0.2, 0.5]
        ws = [0.1, 0.6]
        rp = 3
        rs = 80
        N, Wn = buttord(wp, ws, rp, rs, False)
        b, a = butter(N, Wn, 'bandpass', False)
        w, h = freqz(b, a)
        w /= np.pi
        assert_array_less(-rp - 0.1,
                          dB(h[np.logical_and(wp[0] <= w, w <= wp[1])]))
        assert_array_less(dB(h[np.logical_or(w <= ws[0], ws[1] <= w)]),
                          -rs + 0.1)

        assert_equal(N, 18)
        assert_allclose(Wn, [1.9998742411409134e-01, 5.0002139595676276e-01],
                        rtol=1e-15)

    def test_bandstop(self):
        wp = [0.1, 0.6]
        ws = [0.2, 0.5]
        rp = 3
        rs = 90
        N, Wn = buttord(wp, ws, rp, rs, False)
        b, a = butter(N, Wn, 'bandstop', False)
        w, h = freqz(b, a)
        w /= np.pi
        assert_array_less(-rp,
                          dB(h[np.logical_or(w <= wp[0], wp[1] <= w)]))
        assert_array_less(dB(h[np.logical_and(ws[0] <= w, w <= ws[1])]),
                          -rs)

        assert_equal(N, 20)
        assert_allclose(Wn, [1.4759432329294042e-01, 5.9997365985276407e-01],
                        rtol=1e-6)

    def test_analog(self):
        wp = 200
        ws = 600
        rp = 3
        rs = 60
        N, Wn = buttord(wp, ws, rp, rs, True)
        b, a = butter(N, Wn, 'lowpass', True)
        w, h = freqs(b, a)
        assert_array_less(-rp, dB(h[w <= wp]))
        assert_array_less(dB(h[ws <= w]), -rs)

        assert_equal(N, 7)
        assert_allclose(Wn, 2.0006785355671877e+02, rtol=1e-15)

        n, Wn = buttord(1, 550/450, 1, 26, analog=True)
        assert_equal(n, 19)
        assert_allclose(Wn, 1.0361980524629517, rtol=1e-15)

        assert_equal(buttord(1, 1.2, 1, 80, analog=True)[0], 55)


class TestCheb1ord(TestCase):

    def test_lowpass(self):
        wp = 0.2
        ws = 0.3
        rp = 3
        rs = 60
        N, Wn = cheb1ord(wp, ws, rp, rs, False)
        b, a = cheby1(N, rp, Wn, 'low', False)
        w, h = freqz(b, a)
        w /= np.pi
        assert_array_less(-rp - 0.1, dB(h[w <= wp]))
        assert_array_less(dB(h[ws <= w]), -rs + 0.1)

        assert_equal(N, 8)
        assert_allclose(Wn, 0.2, rtol=1e-15)

    def test_highpass(self):
        wp = 0.3
        ws = 0.2
        rp = 3
        rs = 70
        N, Wn = cheb1ord(wp, ws, rp, rs, False)
        b, a = cheby1(N, rp, Wn, 'high', False)
        w, h = freqz(b, a)
        w /= np.pi
        assert_array_less(-rp - 0.1, dB(h[wp <= w]))
        assert_array_less(dB(h[w <= ws]), -rs + 0.1)

        assert_equal(N, 9)
        assert_allclose(Wn, 0.3, rtol=1e-15)

    def test_bandpass(self):
        wp = [0.2, 0.5]
        ws = [0.1, 0.6]
        rp = 3
        rs = 80
        N, Wn = cheb1ord(wp, ws, rp, rs, False)
        b, a = cheby1(N, rp, Wn, 'band', False)
        w, h = freqz(b, a)
        w /= np.pi
        assert_array_less(-rp - 0.1,
                          dB(h[np.logical_and(wp[0] <= w, w <= wp[1])]))
        assert_array_less(dB(h[np.logical_or(w <= ws[0], ws[1] <= w)]),
                          -rs + 0.1)

        assert_equal(N, 9)
        assert_allclose(Wn, [0.2, 0.5], rtol=1e-15)

    def test_bandstop(self):
        wp = [0.1, 0.6]
        ws = [0.2, 0.5]
        rp = 3
        rs = 90
        N, Wn = cheb1ord(wp, ws, rp, rs, False)
        b, a = cheby1(N, rp, Wn, 'stop', False)
        w, h = freqz(b, a)
        w /= np.pi
        assert_array_less(-rp - 0.1,
                          dB(h[np.logical_or(w <= wp[0], wp[1] <= w)]))
        assert_array_less(dB(h[np.logical_and(ws[0] <= w, w <= ws[1])]),
                          -rs + 0.1)

        assert_equal(N, 10)
        assert_allclose(Wn, [0.14758232569947785, 0.6], rtol=1e-5)

    def test_analog(self):
        wp = 700
        ws = 100
        rp = 3
        rs = 70
        N, Wn = cheb1ord(wp, ws, rp, rs, True)
        b, a = cheby1(N, rp, Wn, 'high', True)
        w, h = freqs(b, a)
        assert_array_less(-rp - 0.1, dB(h[wp <= w]))
        assert_array_less(dB(h[w <= ws]), -rs + 0.1)

        assert_equal(N, 4)
        assert_allclose(Wn, 700, rtol=1e-15)

        assert_equal(cheb1ord(1, 1.2, 1, 80, analog=True)[0], 17)


class TestCheb2ord(TestCase):

    def test_lowpass(self):
        wp = 0.2
        ws = 0.3
        rp = 3
        rs = 60
        N, Wn = cheb2ord(wp, ws, rp, rs, False)
        b, a = cheby2(N, rs, Wn, 'lp', False)
        w, h = freqz(b, a)
        w /= np.pi
        assert_array_less(-rp - 0.1, dB(h[w <= wp]))
        assert_array_less(dB(h[ws <= w]), -rs + 0.1)

        assert_equal(N, 8)
        assert_allclose(Wn, 0.28647639976553163, rtol=1e-15)

    def test_highpass(self):
        wp = 0.3
        ws = 0.2
        rp = 3
        rs = 70
        N, Wn = cheb2ord(wp, ws, rp, rs, False)
        b, a = cheby2(N, rs, Wn, 'hp', False)
        w, h = freqz(b, a)
        w /= np.pi
        assert_array_less(-rp - 0.1, dB(h[wp <= w]))
        assert_array_less(dB(h[w <= ws]), -rs + 0.1)

        assert_equal(N, 9)
        assert_allclose(Wn, 0.20697492182903282, rtol=1e-15)

    def test_bandpass(self):
        wp = [0.2, 0.5]
        ws = [0.1, 0.6]
        rp = 3
        rs = 80
        N, Wn = cheb2ord(wp, ws, rp, rs, False)
        b, a = cheby2(N, rs, Wn, 'bp', False)
        w, h = freqz(b, a)
        w /= np.pi
        assert_array_less(-rp - 0.1,
                          dB(h[np.logical_and(wp[0] <= w, w <= wp[1])]))
        assert_array_less(dB(h[np.logical_or(w <= ws[0], ws[1] <= w)]),
                          -rs + 0.1)

        assert_equal(N, 9)
        assert_allclose(Wn, [0.14876937565923479, 0.59748447842351482],
                        rtol=1e-15)

    def test_bandstop(self):
        wp = [0.1, 0.6]
        ws = [0.2, 0.5]
        rp = 3
        rs = 90
        N, Wn = cheb2ord(wp, ws, rp, rs, False)
        b, a = cheby2(N, rs, Wn, 'bs', False)
        w, h = freqz(b, a)
        w /= np.pi
        assert_array_less(-rp - 0.1,
                          dB(h[np.logical_or(w <= wp[0], wp[1] <= w)]))
        assert_array_less(dB(h[np.logical_and(ws[0] <= w, w <= ws[1])]),
                          -rs + 0.1)

        assert_equal(N, 10)
        assert_allclose(Wn, [0.19926249974781743, 0.50125246585567362],
                        rtol=1e-6)

    def test_analog(self):
        wp = [20, 50]
        ws = [10, 60]
        rp = 3
        rs = 80
        N, Wn = cheb2ord(wp, ws, rp, rs, True)
        b, a = cheby2(N, rs, Wn, 'bp', True)
        w, h = freqs(b, a)
        assert_array_less(-rp - 0.1,
                          dB(h[np.logical_and(wp[0] <= w, w <= wp[1])]))
        assert_array_less(dB(h[np.logical_or(w <= ws[0], ws[1] <= w)]),
                          -rs + 0.1)

        assert_equal(N, 11)
        assert_allclose(Wn, [1.673740595370124e+01, 5.974641487254268e+01],
                        rtol=1e-15)


class TestEllipord(TestCase):

    def test_lowpass(self):
        wp = 0.2
        ws = 0.3
        rp = 3
        rs = 60
        N, Wn = ellipord(wp, ws, rp, rs, False)
        b, a = ellip(N, rp, rs, Wn, 'lp', False)
        w, h = freqz(b, a)
        w /= np.pi
        assert_array_less(-rp - 0.1, dB(h[w <= wp]))
        assert_array_less(dB(h[ws <= w]), -rs + 0.1)

        assert_equal(N, 5)
        assert_allclose(Wn, 0.2, rtol=1e-15)

    def test_highpass(self):
        wp = 0.3
        ws = 0.2
        rp = 3
        rs = 70
        N, Wn = ellipord(wp, ws, rp, rs, False)
        b, a = ellip(N, rp, rs, Wn, 'hp', False)
        w, h = freqz(b, a)
        w /= np.pi
        assert_array_less(-rp - 0.1, dB(h[wp <= w]))
        assert_array_less(dB(h[w <= ws]), -rs + 0.1)

        assert_equal(N, 6)
        assert_allclose(Wn, 0.3, rtol=1e-15)

    def test_bandpass(self):
        wp = [0.2, 0.5]
        ws = [0.1, 0.6]
        rp = 3
        rs = 80
        N, Wn = ellipord(wp, ws, rp, rs, False)
        b, a = ellip(N, rp, rs, Wn, 'bp', False)
        w, h = freqz(b, a)
        w /= np.pi
        assert_array_less(-rp - 0.1,
                          dB(h[np.logical_and(wp[0] <= w, w <= wp[1])]))
        assert_array_less(dB(h[np.logical_or(w <= ws[0], ws[1] <= w)]),
                          -rs + 0.1)

        assert_equal(N, 6)
        assert_allclose(Wn, [0.2, 0.5], rtol=1e-15)

    def test_bandstop(self):
        wp = [0.1, 0.6]
        ws = [0.2, 0.5]
        rp = 3
        rs = 90
        N, Wn = ellipord(wp, ws, rp, rs, False)
        b, a = ellip(N, rp, rs, Wn, 'bs', False)
        w, h = freqz(b, a)
        w /= np.pi
        assert_array_less(-rp - 0.1,
                          dB(h[np.logical_or(w <= wp[0], wp[1] <= w)]))
        assert_array_less(dB(h[np.logical_and(ws[0] <= w, w <= ws[1])]),
                          -rs + 0.1)

        assert_equal(N, 7)
        assert_allclose(Wn, [0.14758232794342988, 0.6], rtol=1e-5)

    def test_analog(self):
        wp = [1000, 6000]
        ws = [2000, 5000]
        rp = 3
        rs = 90
        N, Wn = ellipord(wp, ws, rp, rs, True)
        b, a = ellip(N, rp, rs, Wn, 'bs', True)
        w, h = freqs(b, a)
        assert_array_less(-rp - 0.1,
                          dB(h[np.logical_or(w <= wp[0], wp[1] <= w)]))
        assert_array_less(dB(h[np.logical_and(ws[0] <= w, w <= ws[1])]),
                          -rs + 0.1)

        assert_equal(N, 8)
        assert_allclose(Wn, [1666.6666, 6000])

        assert_equal(ellipord(1, 1.2, 1, 80, analog=True)[0], 9)


class TestBessel(TestCase):

    def test_degenerate(self):
        # 0-order filter is just a passthrough
        b, a = bessel(0, 1, analog=True)
        assert_array_equal(b, [1])
        assert_array_equal(a, [1])

        # 1-order filter is same for all types
        b, a = bessel(1, 1, analog=True)
        assert_array_equal(b, [1])
        assert_array_equal(a, [1, 1])

        z, p, k = bessel(1, 0.3, analog=True, output='zpk')
        assert_array_equal(z, [])
        assert_allclose(p, [-0.3], rtol=1e-14)
        assert_allclose(k, 0.3, rtol=1e-14)

    def test_high_order(self):
        # high even order
        z, p, k = bessel(24, 100, analog=True, output='zpk')
        z2 = []
        p2 = [
             -9.055312334014323e+01 + 4.844005815403969e+00j,
             -9.055312334014323e+01 - 4.844005815403969e+00j,
             -8.983105162681878e+01 + 1.454056170018573e+01j,
             -8.983105162681878e+01 - 1.454056170018573e+01j,
             -8.837357994162065e+01 + 2.426335240122282e+01j,
             -8.837357994162065e+01 - 2.426335240122282e+01j,
             -8.615278316179575e+01 + 3.403202098404543e+01j,
             -8.615278316179575e+01 - 3.403202098404543e+01j,
             -8.312326467067703e+01 + 4.386985940217900e+01j,
             -8.312326467067703e+01 - 4.386985940217900e+01j,
             -7.921695461084202e+01 + 5.380628489700191e+01j,
             -7.921695461084202e+01 - 5.380628489700191e+01j,
             -7.433392285433246e+01 + 6.388084216250878e+01j,
             -7.433392285433246e+01 - 6.388084216250878e+01j,
             -6.832565803501586e+01 + 7.415032695116071e+01j,
             -6.832565803501586e+01 - 7.415032695116071e+01j,
             -6.096221567378025e+01 + 8.470292433074425e+01j,
             -6.096221567378025e+01 - 8.470292433074425e+01j,
             -5.185914574820616e+01 + 9.569048385258847e+01j,
             -5.185914574820616e+01 - 9.569048385258847e+01j,
             -4.027853855197555e+01 + 1.074195196518679e+02j,
             -4.027853855197555e+01 - 1.074195196518679e+02j,
             -2.433481337524861e+01 + 1.207298683731973e+02j,
             -2.433481337524861e+01 - 1.207298683731973e+02j,
             ]
        k2 = 9.999999999999989e+47
        assert_array_equal(z, z2)
        assert_allclose(sorted(p, key=np.imag), sorted(p2, key=np.imag))
        assert_allclose(k, k2, rtol=1e-14)

        # high odd order
        z, p, k = bessel(23, 1000, analog=True, output='zpk')
        z2 = []
        p2 = [
             -2.497697202208956e+02 + 1.202813187870698e+03j,
             -2.497697202208956e+02 - 1.202813187870698e+03j,
             -4.126986617510172e+02 + 1.065328794475509e+03j,
             -4.126986617510172e+02 - 1.065328794475509e+03j,
             -5.304922463809596e+02 + 9.439760364018479e+02j,
             -5.304922463809596e+02 - 9.439760364018479e+02j,
             -9.027564978975828e+02 + 1.010534334242318e+02j,
             -9.027564978975828e+02 - 1.010534334242318e+02j,
             -8.909283244406079e+02 + 2.023024699647598e+02j,
             -8.909283244406079e+02 - 2.023024699647598e+02j,
             -8.709469394347836e+02 + 3.039581994804637e+02j,
             -8.709469394347836e+02 - 3.039581994804637e+02j,
             -8.423805948131370e+02 + 4.062657947488952e+02j,
             -8.423805948131370e+02 - 4.062657947488952e+02j,
             -8.045561642249877e+02 + 5.095305912401127e+02j,
             -8.045561642249877e+02 - 5.095305912401127e+02j,
             -7.564660146766259e+02 + 6.141594859516342e+02j,
             -7.564660146766259e+02 - 6.141594859516342e+02j,
             -6.965966033906477e+02 + 7.207341374730186e+02j,
             -6.965966033906477e+02 - 7.207341374730186e+02j,
             -6.225903228776276e+02 + 8.301558302815096e+02j,
             -6.225903228776276e+02 - 8.301558302815096e+02j,
             -9.066732476324988e+02]
        k2 = 9.999999999999983e+68
        assert_array_equal(z, z2)
        assert_allclose(sorted(p, key=np.imag), sorted(p2, key=np.imag))
        assert_allclose(k, k2, rtol=1e-14)


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

        z, p, k = butter(1, 0.3, output='zpk')
        assert_array_equal(z, [-1])
        assert_allclose(p, [3.249196962329063e-01], rtol=1e-14)
        assert_allclose(k, 3.375401518835469e-01, rtol=1e-14)

    def test_basic(self):
        # analog s-plane
        for N in range(25):
            wn = 0.01
            z, p, k = butter(N, wn, 'low', analog=True, output='zpk')
            assert_array_almost_equal([], z)
            assert_(len(p) == N)
            # All poles should be at distance wn from origin
            assert_array_almost_equal(wn, abs(p))
            assert_(all(np.real(p) <= 0))  # No poles in right half of S-plane
            assert_array_almost_equal(wn**N, k)

        # digital z-plane
        for N in range(25):
            wn = 0.01
            z, p, k = butter(N, wn, 'high', analog=False, output='zpk')
            assert_array_equal(np.ones(N), z)  # All zeros exactly at DC
            assert_(all(np.abs(p) <= 1))  # No poles outside unit circle

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

        b, a = butter(5, 0.4)
        assert_array_almost_equal(b, [0.0219, 0.1097, 0.2194,
                                      0.2194, 0.1097, 0.0219], decimal=4)
        assert_array_almost_equal(a, [1.0000, -0.9853, 0.9738,
                                     -0.3864, 0.1112, -0.0113], decimal=4)

    def test_highpass(self):
        # highpass, high even order
        z, p, k = butter(28, 0.43, 'high', output='zpk')
        z2 = np.ones(28)
        p2 = [
            2.068257195514592e-01 + 9.238294351481734e-01j,
            2.068257195514592e-01 - 9.238294351481734e-01j,
            1.874933103892023e-01 + 8.269455076775277e-01j,
            1.874933103892023e-01 - 8.269455076775277e-01j,
            1.717435567330153e-01 + 7.383078571194629e-01j,
            1.717435567330153e-01 - 7.383078571194629e-01j,
            1.588266870755982e-01 + 6.564623730651094e-01j,
            1.588266870755982e-01 - 6.564623730651094e-01j,
            1.481881532502603e-01 + 5.802343458081779e-01j,
            1.481881532502603e-01 - 5.802343458081779e-01j,
            1.394122576319697e-01 + 5.086609000582009e-01j,
            1.394122576319697e-01 - 5.086609000582009e-01j,
            1.321840881809715e-01 + 4.409411734716436e-01j,
            1.321840881809715e-01 - 4.409411734716436e-01j,
            1.262633413354405e-01 + 3.763990035551881e-01j,
            1.262633413354405e-01 - 3.763990035551881e-01j,
            1.214660449478046e-01 + 3.144545234797277e-01j,
            1.214660449478046e-01 - 3.144545234797277e-01j,
            1.104868766650320e-01 + 2.771505404367791e-02j,
            1.104868766650320e-01 - 2.771505404367791e-02j,
            1.111768629525075e-01 + 8.331369153155753e-02j,
            1.111768629525075e-01 - 8.331369153155753e-02j,
            1.125740630842972e-01 + 1.394219509611784e-01j,
            1.125740630842972e-01 - 1.394219509611784e-01j,
            1.147138487992747e-01 + 1.963932363793666e-01j,
            1.147138487992747e-01 - 1.963932363793666e-01j,
            1.176516491045901e-01 + 2.546021573417188e-01j,
            1.176516491045901e-01 - 2.546021573417188e-01j,
            ]
        k2 = 1.446671081817286e-06
        assert_array_equal(z, z2)
        assert_allclose(sorted(p, key=np.imag),
                        sorted(p2, key=np.imag), rtol=1e-7)
        assert_allclose(k, k2, rtol=1e-10)

        # highpass, high odd order
        z, p, k = butter(27, 0.56, 'high', output='zpk')
        z2 = np.ones(27)
        p2 = [
            -1.772572785680147e-01 + 9.276431102995948e-01j,
            -1.772572785680147e-01 - 9.276431102995948e-01j,
            -1.600766565322114e-01 + 8.264026279893268e-01j,
            -1.600766565322114e-01 - 8.264026279893268e-01j,
            -1.461948419016121e-01 + 7.341841939120078e-01j,
            -1.461948419016121e-01 - 7.341841939120078e-01j,
            -1.348975284762046e-01 + 6.493235066053785e-01j,
            -1.348975284762046e-01 - 6.493235066053785e-01j,
            -1.256628210712206e-01 + 5.704921366889227e-01j,
            -1.256628210712206e-01 - 5.704921366889227e-01j,
            -1.181038235962314e-01 + 4.966120551231630e-01j,
            -1.181038235962314e-01 - 4.966120551231630e-01j,
            -1.119304913239356e-01 + 4.267938916403775e-01j,
            -1.119304913239356e-01 - 4.267938916403775e-01j,
            -1.069237739782691e-01 + 3.602914879527338e-01j,
            -1.069237739782691e-01 - 3.602914879527338e-01j,
            -1.029178030691416e-01 + 2.964677964142126e-01j,
            -1.029178030691416e-01 - 2.964677964142126e-01j,
            -9.978747500816100e-02 + 2.347687643085738e-01j,
            -9.978747500816100e-02 - 2.347687643085738e-01j,
            -9.743974496324025e-02 + 1.747028739092479e-01j,
            -9.743974496324025e-02 - 1.747028739092479e-01j,
            -9.580754551625957e-02 + 1.158246860771989e-01j,
            -9.580754551625957e-02 - 1.158246860771989e-01j,
            -9.484562207782568e-02 + 5.772118357151691e-02j,
            -9.484562207782568e-02 - 5.772118357151691e-02j,
            -9.452783117928215e-02
            ]
        k2 = 9.585686688851069e-09
        assert_array_equal(z, z2)
        assert_allclose(sorted(p, key=np.imag),
                        sorted(p2, key=np.imag), rtol=1e-8)
        assert_allclose(k, k2)

    def test_bandpass(self):
        z, p, k = butter(8, [0.25, 0.33], 'band', output='zpk')
        z2 = [1, 1, 1, 1, 1, 1, 1, 1,
             -1, -1, -1, -1, -1, -1, -1, -1]
        p2 = [
            4.979909925436156e-01 + 8.367609424799387e-01j,
            4.979909925436156e-01 - 8.367609424799387e-01j,
            4.913338722555539e-01 + 7.866774509868817e-01j,
            4.913338722555539e-01 - 7.866774509868817e-01j,
            5.035229361778706e-01 + 7.401147376726750e-01j,
            5.035229361778706e-01 - 7.401147376726750e-01j,
            5.307617160406101e-01 + 7.029184459442954e-01j,
            5.307617160406101e-01 - 7.029184459442954e-01j,
            5.680556159453138e-01 + 6.788228792952775e-01j,
            5.680556159453138e-01 - 6.788228792952775e-01j,
            6.100962560818854e-01 + 6.693849403338664e-01j,
            6.100962560818854e-01 - 6.693849403338664e-01j,
            6.904694312740631e-01 + 6.930501690145245e-01j,
            6.904694312740631e-01 - 6.930501690145245e-01j,
            6.521767004237027e-01 + 6.744414640183752e-01j,
            6.521767004237027e-01 - 6.744414640183752e-01j,
            ]
        k2 = 3.398854055800844e-08
        assert_array_equal(z, z2)
        assert_allclose(sorted(p, key=np.imag),
                        sorted(p2, key=np.imag), rtol=1e-13)
        assert_allclose(k, k2, rtol=1e-13)

        # bandpass analog
        z, p, k = butter(4, [90.5, 110.5], 'bp', analog=True, output='zpk')
        z2 = np.zeros(4)
        p2 = [
            -4.179137760733086e+00 + 1.095935899082837e+02j,
            -4.179137760733086e+00 - 1.095935899082837e+02j,
            -9.593598668443835e+00 + 1.034745398029734e+02j,
            -9.593598668443835e+00 - 1.034745398029734e+02j,
            -8.883991981781929e+00 + 9.582087115567160e+01j,
            -8.883991981781929e+00 - 9.582087115567160e+01j,
            -3.474530886568715e+00 + 9.111599925805801e+01j,
            -3.474530886568715e+00 - 9.111599925805801e+01j,
            ]
        k2 = 1.600000000000001e+05
        assert_array_equal(z, z2)
        assert_allclose(sorted(p, key=np.imag), sorted(p2, key=np.imag))
        assert_allclose(k, k2, rtol=1e-15)

    def test_bandstop(self):
        z, p, k = butter(7, [0.45, 0.56], 'stop', output='zpk')
        z2 = [-1.594474531383421e-02 + 9.998728744679880e-01j,
              -1.594474531383421e-02 - 9.998728744679880e-01j,
              -1.594474531383421e-02 + 9.998728744679880e-01j,
              -1.594474531383421e-02 - 9.998728744679880e-01j,
              -1.594474531383421e-02 + 9.998728744679880e-01j,
              -1.594474531383421e-02 - 9.998728744679880e-01j,
              -1.594474531383421e-02 + 9.998728744679880e-01j,
              -1.594474531383421e-02 - 9.998728744679880e-01j,
              -1.594474531383421e-02 + 9.998728744679880e-01j,
              -1.594474531383421e-02 - 9.998728744679880e-01j,
              -1.594474531383421e-02 + 9.998728744679880e-01j,
              -1.594474531383421e-02 - 9.998728744679880e-01j,
              -1.594474531383421e-02 + 9.998728744679880e-01j,
              -1.594474531383421e-02 - 9.998728744679880e-01j]
        p2 = [-1.766850742887729e-01 + 9.466951258673900e-01j,
              -1.766850742887729e-01 - 9.466951258673900e-01j,
               1.467897662432886e-01 + 9.515917126462422e-01j,
               1.467897662432886e-01 - 9.515917126462422e-01j,
              -1.370083529426906e-01 + 8.880376681273993e-01j,
              -1.370083529426906e-01 - 8.880376681273993e-01j,
               1.086774544701390e-01 + 8.915240810704319e-01j,
               1.086774544701390e-01 - 8.915240810704319e-01j,
              -7.982704457700891e-02 + 8.506056315273435e-01j,
              -7.982704457700891e-02 - 8.506056315273435e-01j,
               5.238812787110331e-02 + 8.524011102699969e-01j,
               5.238812787110331e-02 - 8.524011102699969e-01j,
              -1.357545000491310e-02 + 8.382287744986582e-01j,
              -1.357545000491310e-02 - 8.382287744986582e-01j]
        k2 = 4.577122512960063e-01
        assert_allclose(sorted(z, key=np.imag), sorted(z2, key=np.imag))
        assert_allclose(sorted(p, key=np.imag), sorted(p2, key=np.imag))
        assert_allclose(k, k2, rtol=1e-14)

    def test_ba_output(self):
        b, a = butter(4, [100, 300], 'bandpass', analog=True)
        b2 = [1.6e+09, 0, 0, 0, 0]
        a2 = [1.000000000000000e+00, 5.226251859505511e+02,
              2.565685424949238e+05, 6.794127417357160e+07,
              1.519411254969542e+10, 2.038238225207147e+12,
              2.309116882454312e+14, 1.411088002066486e+16,
              8.099999999999991e+17]
        assert_allclose(b, b2, rtol=1e-14)
        assert_allclose(a, a2, rtol=1e-14)


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

        z, p, k = cheby1(1, 0.1, 0.3, output='zpk')
        assert_array_equal(z, [-1])
        assert_allclose(p, [-5.390126972799615e-01], rtol=1e-14)
        assert_allclose(k, 7.695063486399808e-01, rtol=1e-14)

    def test_basic(self):
        for N in range(25):
            wn = 0.01
            z, p, k = cheby1(N, 1, wn, 'low', analog=True, output='zpk')
            assert_array_almost_equal([], z)
            assert_(len(p) == N)
            assert_(all(np.real(p) <= 0))  # No poles in right half of S-plane

        for N in range(25):
            wn = 0.01
            z, p, k = cheby1(N, 1, wn, 'high', analog=False, output='zpk')
            assert_array_equal(np.ones(N), z)  # All zeros exactly at DC
            assert_(all(np.abs(p) <= 1))  # No poles outside unit circle

        # Same test as TestNormalize
        b, a = cheby1(8, 0.5, 0.048)
        assert_array_almost_equal(b, [
                             2.150733144728282e-11, 1.720586515782626e-10,
                             6.022052805239190e-10, 1.204410561047838e-09,
                             1.505513201309798e-09, 1.204410561047838e-09,
                             6.022052805239190e-10, 1.720586515782626e-10,
                             2.150733144728282e-11], decimal=14)
        assert_array_almost_equal(a, [
                             1.000000000000000e+00, -7.782402035027959e+00,
                             2.654354569747454e+01, -5.182182531666387e+01,
                             6.334127355102684e+01, -4.963358186631157e+01,
                             2.434862182949389e+01, -6.836925348604676e+00,
                             8.412934944449140e-01], decimal=14)

        b, a = cheby1(4, 1, [0.4, 0.7], btype='band')
        assert_array_almost_equal(b, [0.0084, 0, -0.0335, 0, 0.0502, 0,
                                     -0.0335, 0, 0.0084], decimal=4)
        assert_array_almost_equal(a, [1.0, 1.1191, 2.862, 2.2986, 3.4137,
                                      1.8653, 1.8982, 0.5676, 0.4103],
                                      decimal=4)

        b2, a2 = cheby1(5, 3, 1, analog=True)
        assert_array_almost_equal(b2, [0.0626], decimal=4)
        assert_array_almost_equal(a2, [1, 0.5745, 1.4150, 0.5489, 0.4080,
                                       0.0626], decimal=4)

        b, a = cheby1(8, 0.5, 0.1)
        assert_array_almost_equal(b, 1.0e-006 * np.array([
            0.00703924326028, 0.05631394608227, 0.19709881128793,
            0.39419762257586, 0.49274702821983, 0.39419762257586,
            0.19709881128793, 0.05631394608227, 0.00703924326028]),
            decimal=13)
        assert_array_almost_equal(a, [
              1.00000000000000, -7.44912258934158, 24.46749067762108,
              -46.27560200466141, 55.11160187999928, -42.31640010161038,
              20.45543300484147, -5.69110270561444, 0.69770374759022],
            decimal=13)

        b, a = cheby1(8, 0.5, 0.25)
        assert_array_almost_equal(b, 1.0e-003 * np.array([
            0.00895261138923, 0.07162089111382, 0.25067311889837,
            0.50134623779673, 0.62668279724591, 0.50134623779673,
            0.25067311889837, 0.07162089111382, 0.00895261138923]),
            decimal=13)
        assert_array_almost_equal(a, [1.00000000000000, -5.97529229188545,
                                      16.58122329202101, -27.71423273542923,
                                      30.39509758355313, -22.34729670426879,
                                      10.74509800434910, -3.08924633697497,
                                      0.40707685889802], decimal=13)

    def test_highpass(self):
        # high even order
        z, p, k = cheby1(24, 0.7, 0.2, 'high', output='zpk')
        z2 = np.ones(24)
        p2 = [-6.136558509657073e-01 + 2.700091504942893e-01j,
              -6.136558509657073e-01 - 2.700091504942893e-01j,
              -3.303348340927516e-01 + 6.659400861114254e-01j,
              -3.303348340927516e-01 - 6.659400861114254e-01j,
              8.779713780557169e-03 + 8.223108447483040e-01j,
              8.779713780557169e-03 - 8.223108447483040e-01j,
              2.742361123006911e-01 + 8.356666951611864e-01j,
              2.742361123006911e-01 - 8.356666951611864e-01j,
              4.562984557158206e-01 + 7.954276912303594e-01j,
              4.562984557158206e-01 - 7.954276912303594e-01j,
              5.777335494123628e-01 + 7.435821817961783e-01j,
              5.777335494123628e-01 - 7.435821817961783e-01j,
              6.593260977749194e-01 + 6.955390907990932e-01j,
              6.593260977749194e-01 - 6.955390907990932e-01j,
              7.149590948466562e-01 + 6.559437858502012e-01j,
              7.149590948466562e-01 - 6.559437858502012e-01j,
              7.532432388188739e-01 + 6.256158042292060e-01j,
              7.532432388188739e-01 - 6.256158042292060e-01j,
              7.794365244268271e-01 + 6.042099234813333e-01j,
              7.794365244268271e-01 - 6.042099234813333e-01j,
              7.967253874772997e-01 + 5.911966597313203e-01j,
              7.967253874772997e-01 - 5.911966597313203e-01j,
              8.069756417293870e-01 + 5.862214589217275e-01j,
              8.069756417293870e-01 - 5.862214589217275e-01j]
        k2 = 6.190427617192018e-04
        assert_array_equal(z, z2)
        assert_allclose(sorted(p, key=np.imag),
                        sorted(p2, key=np.imag), rtol=1e-10)
        assert_allclose(k, k2, rtol=1e-10)

        # high odd order
        z, p, k = cheby1(23, 0.8, 0.3, 'high', output='zpk')
        z2 = np.ones(23)
        p2 = [-7.676400532011010e-01,
              -6.754621070166477e-01 + 3.970502605619561e-01j,
              -6.754621070166477e-01 - 3.970502605619561e-01j,
              -4.528880018446727e-01 + 6.844061483786332e-01j,
              -4.528880018446727e-01 - 6.844061483786332e-01j,
              -1.986009130216447e-01 + 8.382285942941594e-01j,
              -1.986009130216447e-01 - 8.382285942941594e-01j,
              2.504673931532608e-02 + 8.958137635794080e-01j,
              2.504673931532608e-02 - 8.958137635794080e-01j,
              2.001089429976469e-01 + 9.010678290791480e-01j,
              2.001089429976469e-01 - 9.010678290791480e-01j,
              3.302410157191755e-01 + 8.835444665962544e-01j,
              3.302410157191755e-01 - 8.835444665962544e-01j,
              4.246662537333661e-01 + 8.594054226449009e-01j,
              4.246662537333661e-01 - 8.594054226449009e-01j,
              4.919620928120296e-01 + 8.366772762965786e-01j,
              4.919620928120296e-01 - 8.366772762965786e-01j,
              5.385746917494749e-01 + 8.191616180796720e-01j,
              5.385746917494749e-01 - 8.191616180796720e-01j,
              5.855636993537203e-01 + 8.060680937701062e-01j,
              5.855636993537203e-01 - 8.060680937701062e-01j,
              5.688812849391721e-01 + 8.086497795114683e-01j,
              5.688812849391721e-01 - 8.086497795114683e-01j]
        k2 = 1.941697029206324e-05
        assert_array_equal(z, z2)
        assert_allclose(sorted(p, key=np.imag),
                        sorted(p2, key=np.imag), rtol=1e-10)
        assert_allclose(k, k2, rtol=1e-10)

        z, p, k = cheby1(10, 1, 1000, 'high', analog=True, output='zpk')
        z2 = np.zeros(10)
        p2 = [-3.144743169501551e+03 + 3.511680029092744e+03j,
              -3.144743169501551e+03 - 3.511680029092744e+03j,
              -5.633065604514602e+02 + 2.023615191183945e+03j,
              -5.633065604514602e+02 - 2.023615191183945e+03j,
              -1.946412183352025e+02 + 1.372309454274755e+03j,
              -1.946412183352025e+02 - 1.372309454274755e+03j,
              -7.987162953085479e+01 + 1.105207708045358e+03j,
              -7.987162953085479e+01 - 1.105207708045358e+03j,
              -2.250315039031946e+01 + 1.001723931471477e+03j,
              -2.250315039031946e+01 - 1.001723931471477e+03j]
        k2 = 8.912509381337453e-01
        assert_array_equal(z, z2)
        assert_allclose(sorted(p, key=np.imag),
                        sorted(p2, key=np.imag), rtol=1e-13)
        assert_allclose(k, k2, rtol=1e-15)

    def test_bandpass(self):
        z, p, k = cheby1(8, 1, [0.3, 0.4], 'bp', output='zpk')
        z2 = [1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1]
        p2 = [3.077784854851463e-01 + 9.453307017592942e-01j,
              3.077784854851463e-01 - 9.453307017592942e-01j,
              3.280567400654425e-01 + 9.272377218689016e-01j,
              3.280567400654425e-01 - 9.272377218689016e-01j,
              3.677912763284301e-01 + 9.038008865279966e-01j,
              3.677912763284301e-01 - 9.038008865279966e-01j,
              4.194425632520948e-01 + 8.769407159656157e-01j,
              4.194425632520948e-01 - 8.769407159656157e-01j,
              4.740921994669189e-01 + 8.496508528630974e-01j,
              4.740921994669189e-01 - 8.496508528630974e-01j,
              5.234866481897429e-01 + 8.259608422808477e-01j,
              5.234866481897429e-01 - 8.259608422808477e-01j,
              5.844717632289875e-01 + 8.052901363500210e-01j,
              5.844717632289875e-01 - 8.052901363500210e-01j,
              5.615189063336070e-01 + 8.100667803850766e-01j,
              5.615189063336070e-01 - 8.100667803850766e-01j]
        k2 = 5.007028718074307e-09
        assert_array_equal(z, z2)
        assert_allclose(sorted(p, key=np.imag),
                        sorted(p2, key=np.imag), rtol=1e-13)
        assert_allclose(k, k2, rtol=1e-13)

    def test_bandstop(self):
        z, p, k = cheby1(7, 1, [0.5, 0.6], 'stop', output='zpk')
        z2 = [-1.583844403245361e-01 + 9.873775210440450e-01j,
              -1.583844403245361e-01 - 9.873775210440450e-01j,
              -1.583844403245361e-01 + 9.873775210440450e-01j,
              -1.583844403245361e-01 - 9.873775210440450e-01j,
              -1.583844403245361e-01 + 9.873775210440450e-01j,
              -1.583844403245361e-01 - 9.873775210440450e-01j,
              -1.583844403245361e-01 + 9.873775210440450e-01j,
              -1.583844403245361e-01 - 9.873775210440450e-01j,
              -1.583844403245361e-01 + 9.873775210440450e-01j,
              -1.583844403245361e-01 - 9.873775210440450e-01j,
              -1.583844403245361e-01 + 9.873775210440450e-01j,
              -1.583844403245361e-01 - 9.873775210440450e-01j,
              -1.583844403245361e-01 + 9.873775210440450e-01j,
              -1.583844403245361e-01 - 9.873775210440450e-01j]
        p2 = [-8.942974551472813e-02 + 3.482480481185926e-01j,
              -8.942974551472813e-02 - 3.482480481185926e-01j,
               1.293775154041798e-01 + 8.753499858081858e-01j,
               1.293775154041798e-01 - 8.753499858081858e-01j,
               3.399741945062013e-02 + 9.690316022705607e-01j,
               3.399741945062013e-02 - 9.690316022705607e-01j,
               4.167225522796539e-04 + 9.927338161087488e-01j,
               4.167225522796539e-04 - 9.927338161087488e-01j,
              -3.912966549550960e-01 + 8.046122859255742e-01j,
              -3.912966549550960e-01 - 8.046122859255742e-01j,
              -3.307805547127368e-01 + 9.133455018206508e-01j,
              -3.307805547127368e-01 - 9.133455018206508e-01j,
              -3.072658345097743e-01 + 9.443589759799366e-01j,
              -3.072658345097743e-01 - 9.443589759799366e-01j]
        k2 = 3.619438310405028e-01
        assert_allclose(sorted(z, key=np.imag),
                        sorted(z2, key=np.imag), rtol=1e-13)
        assert_allclose(sorted(p, key=np.imag),
                        sorted(p2, key=np.imag), rtol=1e-13)
        assert_allclose(k, k2, rtol=1e-15)

    def test_ba_output(self):
        # with transfer function conversion,  without digital conversion
        b, a = cheby1(5, 0.9, [210, 310], 'stop', analog=True)
        b2 = [1.000000000000006e+00, 0,
              3.255000000000020e+05, 0,
              4.238010000000026e+10, 0,
              2.758944510000017e+15, 0,
              8.980364380050052e+19, 0,
              1.169243442282517e+24
              ]
        a2 = [1.000000000000000e+00, 4.630555945694342e+02,
              4.039266454794788e+05, 1.338060988610237e+08,
              5.844333551294591e+10, 1.357346371637638e+13,
              3.804661141892782e+15, 5.670715850340080e+17,
              1.114411200988328e+20, 8.316815934908471e+21,
              1.169243442282517e+24
              ]
        assert_allclose(b, b2, rtol=1e-14)
        assert_allclose(a, a2, rtol=1e-14)


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

        z, p, k = cheby2(1, 50, 0.3, output='zpk')
        assert_array_equal(z, [-1])
        assert_allclose(p, [9.967826460175649e-01], rtol=1e-14)
        assert_allclose(k, 1.608676991217512e-03, rtol=1e-14)

    def test_basic(self):
        for N in range(25):
            wn = 0.01
            z, p, k = cheby2(N, 40, wn, 'low', analog=True, output='zpk')
            assert_(len(p) == N)
            assert_(all(np.real(p) <= 0))  # No poles in right half of S-plane

        for N in range(25):
            wn = 0.01
            z, p, k = cheby2(N, 40, wn, 'high', analog=False, output='zpk')
            assert_(all(np.abs(p) <= 1))  # No poles outside unit circle

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

    def test_highpass(self):
        # high even order
        z, p, k = cheby2(26, 60, 0.3, 'high', output='zpk')
        z2 = [9.981088955489852e-01 + 6.147058341984388e-02j,
              9.981088955489852e-01 - 6.147058341984388e-02j,
              9.832702870387426e-01 + 1.821525257215483e-01j,
              9.832702870387426e-01 - 1.821525257215483e-01j,
              9.550760158089112e-01 + 2.963609353922882e-01j,
              9.550760158089112e-01 - 2.963609353922882e-01j,
              9.162054748821922e-01 + 4.007087817803773e-01j,
              9.162054748821922e-01 - 4.007087817803773e-01j,
              8.700619897368064e-01 + 4.929423232136168e-01j,
              8.700619897368064e-01 - 4.929423232136168e-01j,
              5.889791753434985e-01 + 8.081482110427953e-01j,
              5.889791753434985e-01 - 8.081482110427953e-01j,
              5.984900456570295e-01 + 8.011302423760501e-01j,
              5.984900456570295e-01 - 8.011302423760501e-01j,
              6.172880888914629e-01 + 7.867371958365343e-01j,
              6.172880888914629e-01 - 7.867371958365343e-01j,
              6.448899971038180e-01 + 7.642754030030161e-01j,
              6.448899971038180e-01 - 7.642754030030161e-01j,
              6.804845629637927e-01 + 7.327624168637228e-01j,
              6.804845629637927e-01 - 7.327624168637228e-01j,
              8.202619107108660e-01 + 5.719881098737678e-01j,
              8.202619107108660e-01 - 5.719881098737678e-01j,
              7.228410452536148e-01 + 6.910143437705678e-01j,
              7.228410452536148e-01 - 6.910143437705678e-01j,
              7.702121399578629e-01 + 6.377877856007792e-01j,
              7.702121399578629e-01 - 6.377877856007792e-01j]
        p2 = [7.365546198286450e-01 + 4.842085129329526e-02j,
              7.365546198286450e-01 - 4.842085129329526e-02j,
              7.292038510962885e-01 + 1.442201672097581e-01j,
              7.292038510962885e-01 - 1.442201672097581e-01j,
              7.151293788040354e-01 + 2.369925800458584e-01j,
              7.151293788040354e-01 - 2.369925800458584e-01j,
              6.955051820787286e-01 + 3.250341363856910e-01j,
              6.955051820787286e-01 - 3.250341363856910e-01j,
              6.719122956045220e-01 + 4.070475750638047e-01j,
              6.719122956045220e-01 - 4.070475750638047e-01j,
              6.461722130611300e-01 + 4.821965916689270e-01j,
              6.461722130611300e-01 - 4.821965916689270e-01j,
              5.528045062872224e-01 + 8.162920513838372e-01j,
              5.528045062872224e-01 - 8.162920513838372e-01j,
              5.464847782492791e-01 + 7.869899955967304e-01j,
              5.464847782492791e-01 - 7.869899955967304e-01j,
              5.488033111260949e-01 + 7.520442354055579e-01j,
              5.488033111260949e-01 - 7.520442354055579e-01j,
              6.201874719022955e-01 + 5.500894392527353e-01j,
              6.201874719022955e-01 - 5.500894392527353e-01j,
              5.586478152536709e-01 + 7.112676877332921e-01j,
              5.586478152536709e-01 - 7.112676877332921e-01j,
              5.958145844148228e-01 + 6.107074340842115e-01j,
              5.958145844148228e-01 - 6.107074340842115e-01j,
              5.747812938519067e-01 + 6.643001536914696e-01j,
              5.747812938519067e-01 - 6.643001536914696e-01j]
        k2 = 9.932997786497189e-02
        assert_allclose(sorted(z, key=np.angle),
                        sorted(z2, key=np.angle), rtol=1e-13)
        assert_allclose(sorted(p, key=np.angle),
                        sorted(p2, key=np.angle), rtol=1e-12)
        assert_allclose(k, k2, rtol=1e-11)

        # high odd order
        z, p, k = cheby2(25, 80, 0.5, 'high', output='zpk')
        z2 = [9.690690376586687e-01 + 2.467897896011971e-01j,
              9.690690376586687e-01 - 2.467897896011971e-01j,
              9.999999999999492e-01,
              8.835111277191199e-01 + 4.684101698261429e-01j,
              8.835111277191199e-01 - 4.684101698261429e-01j,
              7.613142857900539e-01 + 6.483830335935022e-01j,
              7.613142857900539e-01 - 6.483830335935022e-01j,
              6.232625173626231e-01 + 7.820126817709752e-01j,
              6.232625173626231e-01 - 7.820126817709752e-01j,
              4.864456563413621e-01 + 8.737108351316745e-01j,
              4.864456563413621e-01 - 8.737108351316745e-01j,
              3.618368136816749e-01 + 9.322414495530347e-01j,
              3.618368136816749e-01 - 9.322414495530347e-01j,
              2.549486883466794e-01 + 9.669545833752675e-01j,
              2.549486883466794e-01 - 9.669545833752675e-01j,
              1.676175432109457e-01 + 9.858520980390212e-01j,
              1.676175432109457e-01 - 9.858520980390212e-01j,
              1.975218468277521e-03 + 9.999980492540941e-01j,
              1.975218468277521e-03 - 9.999980492540941e-01j,
              1.786959496651858e-02 + 9.998403260399917e-01j,
              1.786959496651858e-02 - 9.998403260399917e-01j,
              9.967933660557139e-02 + 9.950196127985684e-01j,
              9.967933660557139e-02 - 9.950196127985684e-01j,
              5.013970951219547e-02 + 9.987422137518890e-01j,
              5.013970951219547e-02 - 9.987422137518890e-01j]
        p2 = [4.218866331906864e-01,
              4.120110200127552e-01 + 1.361290593621978e-01j,
              4.120110200127552e-01 - 1.361290593621978e-01j,
              3.835890113632530e-01 + 2.664910809911026e-01j,
              3.835890113632530e-01 - 2.664910809911026e-01j,
              3.399195570456499e-01 + 3.863983538639875e-01j,
              3.399195570456499e-01 - 3.863983538639875e-01j,
              2.855977834508353e-01 + 4.929444399540688e-01j,
              2.855977834508353e-01 - 4.929444399540688e-01j,
              2.255765441339322e-01 + 5.851631870205766e-01j,
              2.255765441339322e-01 - 5.851631870205766e-01j,
              1.644087535815792e-01 + 6.637356937277153e-01j,
              1.644087535815792e-01 - 6.637356937277153e-01j,
              -7.293633845273095e-02 + 9.739218252516307e-01j,
              -7.293633845273095e-02 - 9.739218252516307e-01j,
              1.058259206358626e-01 + 7.304739464862978e-01j,
              1.058259206358626e-01 - 7.304739464862978e-01j,
              -5.703971947785402e-02 + 9.291057542169088e-01j,
              -5.703971947785402e-02 - 9.291057542169088e-01j,
              5.263875132656864e-02 + 7.877974334424453e-01j,
              5.263875132656864e-02 - 7.877974334424453e-01j,
              -3.007943405982616e-02 + 8.846331716180016e-01j,
              -3.007943405982616e-02 - 8.846331716180016e-01j,
              6.857277464483946e-03 + 8.383275456264492e-01j,
              6.857277464483946e-03 - 8.383275456264492e-01j]
        k2 = 6.507068761705037e-03
        assert_allclose(sorted(z, key=np.angle),
                        sorted(z2, key=np.angle), rtol=1e-13)
        assert_allclose(sorted(p, key=np.angle),
                        sorted(p2, key=np.angle), rtol=1e-12)
        assert_allclose(k, k2, rtol=1e-11)

    def test_bandpass(self):
        z, p, k = cheby2(9, 40, [0.07, 0.2], 'pass', output='zpk')
        z2 = [-9.999999999999999e-01,
               3.676588029658514e-01 + 9.299607543341383e-01j,
               3.676588029658514e-01 - 9.299607543341383e-01j,
               7.009689684982283e-01 + 7.131917730894889e-01j,
               7.009689684982283e-01 - 7.131917730894889e-01j,
               7.815697973765858e-01 + 6.238178033919218e-01j,
               7.815697973765858e-01 - 6.238178033919218e-01j,
               8.063793628819866e-01 + 5.913986160941200e-01j,
               8.063793628819866e-01 - 5.913986160941200e-01j,
               1.000000000000001e+00,
               9.944493019920448e-01 + 1.052168511576739e-01j,
               9.944493019920448e-01 - 1.052168511576739e-01j,
               9.854674703367308e-01 + 1.698642543566085e-01j,
               9.854674703367308e-01 - 1.698642543566085e-01j,
               9.762751735919308e-01 + 2.165335665157851e-01j,
               9.762751735919308e-01 - 2.165335665157851e-01j,
               9.792277171575134e-01 + 2.027636011479496e-01j,
               9.792277171575134e-01 - 2.027636011479496e-01j]
        p2 = [8.143803410489621e-01 + 5.411056063397541e-01j,
              8.143803410489621e-01 - 5.411056063397541e-01j,
              7.650769827887418e-01 + 5.195412242095543e-01j,
              7.650769827887418e-01 - 5.195412242095543e-01j,
              6.096241204063443e-01 + 3.568440484659796e-01j,
              6.096241204063443e-01 - 3.568440484659796e-01j,
              6.918192770246239e-01 + 4.770463577106911e-01j,
              6.918192770246239e-01 - 4.770463577106911e-01j,
              6.986241085779207e-01 + 1.146512226180060e-01j,
              6.986241085779207e-01 - 1.146512226180060e-01j,
              8.654645923909734e-01 + 1.604208797063147e-01j,
              8.654645923909734e-01 - 1.604208797063147e-01j,
              9.164831670444591e-01 + 1.969181049384918e-01j,
              9.164831670444591e-01 - 1.969181049384918e-01j,
              9.630425777594550e-01 + 2.317513360702271e-01j,
              9.630425777594550e-01 - 2.317513360702271e-01j,
              9.438104703725529e-01 + 2.193509900269860e-01j,
              9.438104703725529e-01 - 2.193509900269860e-01j]
        k2 = 9.345352824659604e-03
        assert_allclose(sorted(z, key=np.angle),
                        sorted(z2, key=np.angle), rtol=1e-13)
        assert_allclose(sorted(p, key=np.angle),
                        sorted(p2, key=np.angle), rtol=1e-13)
        assert_allclose(k, k2, rtol=1e-11)

    def test_bandstop(self):
        z, p, k = cheby2(6, 55, [0.1, 0.9], 'stop', output='zpk')
        z2 = [6.230544895101009e-01 + 7.821784343111114e-01j,
              6.230544895101009e-01 - 7.821784343111114e-01j,
              9.086608545660115e-01 + 4.175349702471991e-01j,
              9.086608545660115e-01 - 4.175349702471991e-01j,
              9.478129721465802e-01 + 3.188268649763867e-01j,
              9.478129721465802e-01 - 3.188268649763867e-01j,
              -6.230544895100982e-01 + 7.821784343111109e-01j,
              -6.230544895100982e-01 - 7.821784343111109e-01j,
              -9.086608545660116e-01 + 4.175349702472088e-01j,
              -9.086608545660116e-01 - 4.175349702472088e-01j,
              -9.478129721465784e-01 + 3.188268649763897e-01j,
              -9.478129721465784e-01 - 3.188268649763897e-01j]
        p2 = [-9.464094036167638e-01 + 1.720048695084344e-01j,
              -9.464094036167638e-01 - 1.720048695084344e-01j,
              -8.715844103386737e-01 + 1.370665039509297e-01j,
              -8.715844103386737e-01 - 1.370665039509297e-01j,
              -8.078751204586425e-01 + 5.729329866682983e-02j,
              -8.078751204586425e-01 - 5.729329866682983e-02j,
               9.464094036167665e-01 + 1.720048695084332e-01j,
               9.464094036167665e-01 - 1.720048695084332e-01j,
               8.078751204586447e-01 + 5.729329866683007e-02j,
               8.078751204586447e-01 - 5.729329866683007e-02j,
               8.715844103386721e-01 + 1.370665039509331e-01j,
               8.715844103386721e-01 - 1.370665039509331e-01j]
        k2 = 2.917823332763358e-03
        assert_allclose(sorted(z, key=np.angle),
                        sorted(z2, key=np.angle), rtol=1e-13)
        assert_allclose(sorted(p, key=np.angle),
                        sorted(p2, key=np.angle), rtol=1e-13)
        assert_allclose(k, k2, rtol=1e-11)

    def test_ba_output(self):
        # with transfer function conversion, without digital conversion
        b, a = cheby2(5, 20, [2010, 2100], 'stop', True)
        b2 = [1.000000000000000e+00, 0,  # Matlab: 6.683253076978249e-12,
              2.111512500000000e+07, 0,  # Matlab: 1.134325604589552e-04,
              1.782966433781250e+14, 0,  # Matlab: 7.216787944356781e+02,
              7.525901316990656e+20, 0,  # Matlab: 2.039829265789886e+09,
              1.587960565565748e+27, 0,  # Matlab: 2.161236218626134e+15,
              1.339913493808585e+33]
        a2 = [1.000000000000000e+00, 1.849550755473371e+02,
              2.113222918998538e+07, 3.125114149732283e+09,
              1.785133457155609e+14, 1.979158697776348e+16,
              7.535048322653831e+20, 5.567966191263037e+22,
              1.589246884221346e+27, 5.871210648525566e+28,
              1.339913493808590e+33]
        assert_allclose(b, b2, rtol=1e-14)
        assert_allclose(a, a2, rtol=1e-14)


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

        z, p, k = ellip(1, 1, 55, 0.3, output='zpk')
        assert_allclose(z, [-9.999999999999998e-01], rtol=1e-14)
        assert_allclose(p, [-6.660721153525525e-04], rtol=1e-10)
        assert_allclose(k, 5.003330360576763e-01, rtol=1e-14)

    def test_basic(self):
        for N in range(25):
            wn = 0.01
            z, p, k = ellip(N, 1, 40, wn, 'low', analog=True, output='zpk')
            assert_(len(p) == N)
            assert_(all(np.real(p) <= 0))  # No poles in right half of S-plane

        for N in range(25):
            wn = 0.01
            z, p, k = ellip(N, 1, 40, wn, 'high', analog=False, output='zpk')
            assert_(all(np.abs(p) <= 1))  # No poles outside unit circle

        b3, a3 = ellip(5, 3, 26, 1, analog=True)
        assert_array_almost_equal(b3, [0.1420, 0, 0.3764, 0,
                                       0.2409], decimal=4)
        assert_array_almost_equal(a3, [1, 0.5686, 1.8061, 0.8017, 0.8012,
                                       0.2409], decimal=4)

        b, a = ellip(3, 1, 60, [0.4, 0.7], 'stop')
        assert_array_almost_equal(b, [0.3310, 0.3469, 1.1042, 0.7044, 1.1042,
                                      0.3469, 0.3310], decimal=4)
        assert_array_almost_equal(a, [1.0000, 0.6973, 1.1441, 0.5878, 0.7323,
                                      0.1131, -0.0060], decimal=4)

    def test_highpass(self):
        # high even order
        z, p, k = ellip(24, 1, 80, 0.3, 'high', output='zpk')
        z2 = [9.761875332501075e-01 + 2.169283290099910e-01j,
              9.761875332501075e-01 - 2.169283290099910e-01j,
              8.413503353963494e-01 + 5.404901600661900e-01j,
              8.413503353963494e-01 - 5.404901600661900e-01j,
              7.160082576305009e-01 + 6.980918098681732e-01j,
              7.160082576305009e-01 - 6.980918098681732e-01j,
              6.456533638965329e-01 + 7.636306264739803e-01j,
              6.456533638965329e-01 - 7.636306264739803e-01j,
              6.127321820971366e-01 + 7.902906256703928e-01j,
              6.127321820971366e-01 - 7.902906256703928e-01j,
              5.983607817490196e-01 + 8.012267936512676e-01j,
              5.983607817490196e-01 - 8.012267936512676e-01j,
              5.922577552594799e-01 + 8.057485658286990e-01j,
              5.922577552594799e-01 - 8.057485658286990e-01j,
              5.896952092563588e-01 + 8.076258788449631e-01j,
              5.896952092563588e-01 - 8.076258788449631e-01j,
              5.886248765538837e-01 + 8.084063054565607e-01j,
              5.886248765538837e-01 - 8.084063054565607e-01j,
              5.881802711123132e-01 + 8.087298490066037e-01j,
              5.881802711123132e-01 - 8.087298490066037e-01j,
              5.879995719101164e-01 + 8.088612386766461e-01j,
              5.879995719101164e-01 - 8.088612386766461e-01j,
              5.879354086709576e-01 + 8.089078780868164e-01j,
              5.879354086709576e-01 - 8.089078780868164e-01j]
        p2 = [-3.184805259081650e-01 + 4.206951906775851e-01j,
              -3.184805259081650e-01 - 4.206951906775851e-01j,
               1.417279173459985e-01 + 7.903955262836452e-01j,
               1.417279173459985e-01 - 7.903955262836452e-01j,
               4.042881216964651e-01 + 8.309042239116594e-01j,
               4.042881216964651e-01 - 8.309042239116594e-01j,
               5.128964442789670e-01 + 8.229563236799665e-01j,
               5.128964442789670e-01 - 8.229563236799665e-01j,
               5.569614712822724e-01 + 8.155957702908510e-01j,
               5.569614712822724e-01 - 8.155957702908510e-01j,
               5.750478870161392e-01 + 8.118633973883931e-01j,
               5.750478870161392e-01 - 8.118633973883931e-01j,
               5.825314018170804e-01 + 8.101960910679270e-01j,
               5.825314018170804e-01 - 8.101960910679270e-01j,
               5.856397379751872e-01 + 8.094825218722543e-01j,
               5.856397379751872e-01 - 8.094825218722543e-01j,
               5.869326035251949e-01 + 8.091827531557583e-01j,
               5.869326035251949e-01 - 8.091827531557583e-01j,
               5.874697218855733e-01 + 8.090593298213502e-01j,
               5.874697218855733e-01 - 8.090593298213502e-01j,
               5.876904783532237e-01 + 8.090127161018823e-01j,
               5.876904783532237e-01 - 8.090127161018823e-01j,
               5.877753105317594e-01 + 8.090050577978136e-01j,
               5.877753105317594e-01 - 8.090050577978136e-01j]
        k2 = 4.918081266957108e-02
        assert_allclose(sorted(z, key=np.angle),
                        sorted(z2, key=np.angle), rtol=1e-4)
        assert_allclose(sorted(p, key=np.angle),
                        sorted(p2, key=np.angle), rtol=1e-4)
        assert_allclose(k, k2, rtol=1e-3)

        # high odd order
        z, p, k = ellip(23, 1, 70, 0.5, 'high', output='zpk')
        z2 = [9.999999999998661e-01,
              6.603717261750994e-01 + 7.509388678638675e-01j,
              6.603717261750994e-01 - 7.509388678638675e-01j,
              2.788635267510325e-01 + 9.603307416968041e-01j,
              2.788635267510325e-01 - 9.603307416968041e-01j,
              1.070215532544218e-01 + 9.942567008268131e-01j,
              1.070215532544218e-01 - 9.942567008268131e-01j,
              4.049427369978163e-02 + 9.991797705105507e-01j,
              4.049427369978163e-02 - 9.991797705105507e-01j,
              1.531059368627931e-02 + 9.998827859909265e-01j,
              1.531059368627931e-02 - 9.998827859909265e-01j,
              5.808061438534933e-03 + 9.999831330689181e-01j,
              5.808061438534933e-03 - 9.999831330689181e-01j,
              2.224277847754599e-03 + 9.999975262909676e-01j,
              2.224277847754599e-03 - 9.999975262909676e-01j,
              8.731857107534554e-04 + 9.999996187732845e-01j,
              8.731857107534554e-04 - 9.999996187732845e-01j,
              3.649057346914968e-04 + 9.999999334218996e-01j,
              3.649057346914968e-04 - 9.999999334218996e-01j,
              1.765538109802615e-04 + 9.999999844143768e-01j,
              1.765538109802615e-04 - 9.999999844143768e-01j,
              1.143655290967426e-04 + 9.999999934602630e-01j,
              1.143655290967426e-04 - 9.999999934602630e-01j]
        p2 = [-6.322017026545028e-01,
              -4.648423756662754e-01 + 5.852407464440732e-01j,
              -4.648423756662754e-01 - 5.852407464440732e-01j,
              -2.249233374627773e-01 + 8.577853017985717e-01j,
              -2.249233374627773e-01 - 8.577853017985717e-01j,
              -9.234137570557621e-02 + 9.506548198678851e-01j,
              -9.234137570557621e-02 - 9.506548198678851e-01j,
              -3.585663561241373e-02 + 9.821494736043981e-01j,
              -3.585663561241373e-02 - 9.821494736043981e-01j,
              -1.363917242312723e-02 + 9.933844128330656e-01j,
              -1.363917242312723e-02 - 9.933844128330656e-01j,
              -5.131505238923029e-03 + 9.975221173308673e-01j,
              -5.131505238923029e-03 - 9.975221173308673e-01j,
              -1.904937999259502e-03 + 9.990680819857982e-01j,
              -1.904937999259502e-03 - 9.990680819857982e-01j,
              -6.859439885466834e-04 + 9.996492201426826e-01j,
              -6.859439885466834e-04 - 9.996492201426826e-01j,
              -2.269936267937089e-04 + 9.998686250679161e-01j,
              -2.269936267937089e-04 - 9.998686250679161e-01j,
              -5.687071588789117e-05 + 9.999527573294513e-01j,
              -5.687071588789117e-05 - 9.999527573294513e-01j,
              -6.948417068525226e-07 + 9.999882737700173e-01j,
              -6.948417068525226e-07 - 9.999882737700173e-01j]
        k2 = 1.220910020289434e-02
        assert_allclose(sorted(z, key=np.angle),
                        sorted(z2, key=np.angle), rtol=1e-4)
        assert_allclose(sorted(p, key=np.angle),
                        sorted(p2, key=np.angle), rtol=1e-4)
        assert_allclose(k, k2, rtol=1e-3)

    def test_bandpass(self):
        z, p, k = ellip(7, 1, 40, [0.07, 0.2], 'pass', output='zpk')
        z2 = [-9.999999999999991e-01,
               6.856610961780020e-01 + 7.279209168501619e-01j,
               6.856610961780020e-01 - 7.279209168501619e-01j,
               7.850346167691289e-01 + 6.194518952058737e-01j,
               7.850346167691289e-01 - 6.194518952058737e-01j,
               7.999038743173071e-01 + 6.001281461922627e-01j,
               7.999038743173071e-01 - 6.001281461922627e-01j,
               9.999999999999999e-01,
               9.862938983554124e-01 + 1.649980183725925e-01j,
               9.862938983554124e-01 - 1.649980183725925e-01j,
               9.788558330548762e-01 + 2.045513580850601e-01j,
               9.788558330548762e-01 - 2.045513580850601e-01j,
               9.771155231720003e-01 + 2.127093189691258e-01j,
               9.771155231720003e-01 - 2.127093189691258e-01j]
        p2 = [8.063992755498643e-01 + 5.858071374778874e-01j,
              8.063992755498643e-01 - 5.858071374778874e-01j,
              8.050395347071724e-01 + 5.639097428109795e-01j,
              8.050395347071724e-01 - 5.639097428109795e-01j,
              8.113124936559144e-01 + 4.855241143973142e-01j,
              8.113124936559144e-01 - 4.855241143973142e-01j,
              8.665595314082394e-01 + 3.334049560919331e-01j,
              8.665595314082394e-01 - 3.334049560919331e-01j,
              9.412369011968871e-01 + 2.457616651325908e-01j,
              9.412369011968871e-01 - 2.457616651325908e-01j,
              9.679465190411238e-01 + 2.228772501848216e-01j,
              9.679465190411238e-01 - 2.228772501848216e-01j,
              9.747235066273385e-01 + 2.178937926146544e-01j,
              9.747235066273385e-01 - 2.178937926146544e-01j]
        k2 = 8.354782670263239e-03
        assert_allclose(sorted(z, key=np.angle),
                        sorted(z2, key=np.angle), rtol=1e-4)
        assert_allclose(sorted(p, key=np.angle),
                        sorted(p2, key=np.angle), rtol=1e-4)
        assert_allclose(k, k2, rtol=1e-3)

        z, p, k = ellip(5, 1, 75, [90.5, 110.5], 'pass', True, 'zpk')
        z2 = [-5.583607317695175e-14 + 1.433755965989225e+02j,
              -5.583607317695175e-14 - 1.433755965989225e+02j,
               5.740106416459296e-14 + 1.261678754570291e+02j,
               5.740106416459296e-14 - 1.261678754570291e+02j,
              -2.199676239638652e-14 + 6.974861996895196e+01j,
              -2.199676239638652e-14 - 6.974861996895196e+01j,
              -3.372595657044283e-14 + 7.926145989044531e+01j,
              -3.372595657044283e-14 - 7.926145989044531e+01j,
              0]
        p2 = [-8.814960004852743e-01 + 1.104124501436066e+02j,
              -8.814960004852743e-01 - 1.104124501436066e+02j,
              -2.477372459140184e+00 + 1.065638954516534e+02j,
              -2.477372459140184e+00 - 1.065638954516534e+02j,
              -3.072156842945799e+00 + 9.995404870405324e+01j,
              -3.072156842945799e+00 - 9.995404870405324e+01j,
              -2.180456023925693e+00 + 9.379206865455268e+01j,
              -2.180456023925693e+00 - 9.379206865455268e+01j,
              -7.230484977485752e-01 + 9.056598800801140e+01j,
              -7.230484977485752e-01 - 9.056598800801140e+01j]
        k2 = 3.774571622827070e-02
        assert_allclose(sorted(z, key=np.imag),
                        sorted(z2, key=np.imag), rtol=1e-4)
        assert_allclose(sorted(p, key=np.imag),
                        sorted(p2, key=np.imag), rtol=1e-6)
        assert_allclose(k, k2, rtol=1e-3)

    def test_bandstop(self):
        z, p, k = ellip(8, 1, 65, [0.2, 0.4], 'stop', output='zpk')
        z2 = [3.528578094286510e-01 + 9.356769561794296e-01j,
              3.528578094286510e-01 - 9.356769561794296e-01j,
              3.769716042264783e-01 + 9.262248159096587e-01j,
              3.769716042264783e-01 - 9.262248159096587e-01j,
              4.406101783111199e-01 + 8.976985411420985e-01j,
              4.406101783111199e-01 - 8.976985411420985e-01j,
              5.539386470258847e-01 + 8.325574907062760e-01j,
              5.539386470258847e-01 - 8.325574907062760e-01j,
              6.748464963023645e-01 + 7.379581332490555e-01j,
              6.748464963023645e-01 - 7.379581332490555e-01j,
              7.489887970285254e-01 + 6.625826604475596e-01j,
              7.489887970285254e-01 - 6.625826604475596e-01j,
              7.913118471618432e-01 + 6.114127579150699e-01j,
              7.913118471618432e-01 - 6.114127579150699e-01j,
              7.806804740916381e-01 + 6.249303940216475e-01j,
              7.806804740916381e-01 - 6.249303940216475e-01j]

        p2 = [-1.025299146693730e-01 + 5.662682444754943e-01j,
              -1.025299146693730e-01 - 5.662682444754943e-01j,
               1.698463595163031e-01 + 8.926678667070186e-01j,
               1.698463595163031e-01 - 8.926678667070186e-01j,
               2.750532687820631e-01 + 9.351020170094005e-01j,
               2.750532687820631e-01 - 9.351020170094005e-01j,
               3.070095178909486e-01 + 9.457373499553291e-01j,
               3.070095178909486e-01 - 9.457373499553291e-01j,
               7.695332312152288e-01 + 2.792567212705257e-01j,
               7.695332312152288e-01 - 2.792567212705257e-01j,
               8.083818999225620e-01 + 4.990723496863960e-01j,
               8.083818999225620e-01 - 4.990723496863960e-01j,
               8.066158014414928e-01 + 5.649811440393374e-01j,
               8.066158014414928e-01 - 5.649811440393374e-01j,
               8.062787978834571e-01 + 5.855780880424964e-01j,
               8.062787978834571e-01 - 5.855780880424964e-01j]
        k2 = 2.068622545291259e-01
        assert_allclose(sorted(z, key=np.angle),
                        sorted(z2, key=np.angle), rtol=1e-6)
        assert_allclose(sorted(p, key=np.angle),
                        sorted(p2, key=np.angle), rtol=1e-5)
        assert_allclose(k, k2, rtol=1e-5)

    def test_ba_output(self):
        # with transfer function conversion,  without digital conversion
        b, a = ellip(5, 1, 40, [201, 240], 'stop', True)
        b2 = [
             1.000000000000000e+00, 0,  # Matlab: 1.743506051190569e-13,
             2.426561778314366e+05, 0,  # Matlab: 3.459426536825722e-08,
             2.348218683400168e+10, 0,  # Matlab: 2.559179747299313e-03,
             1.132780692872241e+15, 0,  # Matlab: 8.363229375535731e+01,
             2.724038554089566e+19, 0,  # Matlab: 1.018700994113120e+06,
             2.612380874940186e+23
             ]
        a2 = [
             1.000000000000000e+00, 1.337266601804649e+02,
             2.486725353510667e+05, 2.628059713728125e+07,
             2.436169536928770e+10, 1.913554568577315e+12,
             1.175208184614438e+15, 6.115751452473410e+16,
             2.791577695211466e+19, 7.241811142725384e+20,
             2.612380874940182e+23
             ]
        assert_allclose(b, b2, rtol=1e-6)
        assert_allclose(a, a2, rtol=1e-4)


def test_sos_consistency():
    # Consistency checks of output='sos' for the specialized IIR filter
    # design functions.
    design_funcs = [(bessel, (0.1,)),
                    (butter, (0.1,)),
                    (cheby1, (45.0, 0.1)),
                    (cheby2, (0.087, 0.1)),
                    (ellip, (0.087, 45, 0.1))]
    for func, args in design_funcs:
        name = func.__name__

        b, a = func(2, *args, output='ba')
        sos = func(2, *args, output='sos')
        assert_allclose(sos, [np.hstack((b, a))], err_msg="%s(2,...)" % name)

        zpk = func(3, *args, output='zpk')
        sos = func(3, *args, output='sos')
        assert_allclose(sos, zpk2sos(*zpk), err_msg="%s(3,...)" % name)

        zpk = func(4, *args, output='zpk')
        sos = func(4, *args, output='sos')
        assert_allclose(sos, zpk2sos(*zpk), err_msg="%s(4,...)" % name)


class TestIIRFilter(TestCase):

    def test_symmetry(self):
        # All built-in IIR filters are real, so should have perfectly
        # symmetrical poles and zeros. Then ba representation (using
        # numpy.poly) will be purely real instead of having negligible
        # imaginary parts.
        for N in np.arange(1, 26):
            for ftype in ('butter', 'bessel', 'cheby1', 'cheby2', 'ellip'):
                z, p, k = iirfilter(N, 1.1, 1, 20, 'low', analog=True,
                                    ftype=ftype, output='zpk')
                assert_array_equal(sorted(z), sorted(z.conj()))
                assert_array_equal(sorted(p), sorted(p.conj()))
                assert_equal(k, np.real(k))

                b, a = iirfilter(N, 1.1, 1, 20, 'low', analog=True,
                                 ftype=ftype, output='ba')
                assert_(issubclass(b.dtype.type, np.floating))
                assert_(issubclass(a.dtype.type, np.floating))

    def test_int_inputs(self):
        # Using integer frequency arguments and large N should not produce
        # np.ints that wraparound to negative numbers
        k = iirfilter(24, 100, btype='low', analog=True, ftype='bessel',
                      output='zpk')[2]
        k2 = 9.999999999999989e+47
        assert_allclose(k, k2)

    def test_invalid_wn_size(self):
        # low and high have 1 Wn, band and stop have 2 Wn
        assert_raises(ValueError, iirfilter, 1, [0.1, 0.9], btype='low')
        assert_raises(ValueError, iirfilter, 1, [0.2, 0.5], btype='high')
        assert_raises(ValueError, iirfilter, 1, 0.2, btype='bp')
        assert_raises(ValueError, iirfilter, 1, 400, btype='bs', analog=True)

    def test_invalid_wn_range(self):
        # For digital filters, 0 <= Wn <= 1
        assert_raises(ValueError, iirfilter, 1, 2, btype='low')
        assert_raises(ValueError, iirfilter, 1, -1, btype='high')
        assert_raises(ValueError, iirfilter, 1, [1, 2], btype='band')
        assert_raises(ValueError, iirfilter, 1, [10, 20], btype='stop')


class TestGroupDelay(TestCase):
    def test_identity_filter(self):
        w, gd = group_delay((1, 1))
        assert_array_almost_equal(w, pi * np.arange(512) / 512)
        assert_array_almost_equal(gd, np.zeros(512))
        w, gd = group_delay((1, 1), whole=True)
        assert_array_almost_equal(w, 2 * pi * np.arange(512) / 512)
        assert_array_almost_equal(gd, np.zeros(512))

    def test_fir(self):
        # Let's design linear phase FIR and check that the group delay
        # is constant.
        N = 100
        b = firwin(N + 1, 0.1)
        w, gd = group_delay((b, 1))
        assert_allclose(gd, 0.5 * N)

    def test_iir(self):
        # Let's design Butterworth filter and test the group delay at
        # some points against MATLAB answer.
        b, a = butter(4, 0.1)
        w = np.linspace(0, pi, num=10, endpoint=False)
        w, gd = group_delay((b, a), w=w)
        matlab_gd = np.array([8.249313898506037, 11.958947880907104,
                              2.452325615326005, 1.048918665702008,
                              0.611382575635897, 0.418293269460578,
                              0.317932917836572, 0.261371844762525,
                              0.229038045801298, 0.212185774208521])
        assert_array_almost_equal(gd, matlab_gd)

    def test_singular(self):
        # Let's create a filter with zeros and poles on the unit circle and
        # check if warning is raised and the group delay is set to zero at
        # these frequencies.
        z1 = np.exp(1j * 0.1 * pi)
        z2 = np.exp(1j * 0.25 * pi)
        p1 = np.exp(1j * 0.5 * pi)
        p2 = np.exp(1j * 0.8 * pi)
        b = np.convolve([1, -z1], [1, -z2])
        a = np.convolve([1, -p1], [1, -p2])
        w = np.array([0.1 * pi, 0.25 * pi, -0.5 * pi, -0.8 * pi])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            assert_warns(UserWarning, group_delay, (b, a), w=w)
            w, gd = group_delay((b, a), w=w)
            assert_allclose(gd, 0)


if __name__ == "__main__":
    run_module_suite()
