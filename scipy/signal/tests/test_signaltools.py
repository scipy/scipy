#this program corresponds to special.py
from decimal import Decimal

from numpy.testing import *

import scipy.signal as signal
from scipy.signal import lfilter


from numpy import array, arange
import numpy as np

# Use this to test for object arrays filtering - numpy 1.2
# assert_array_almost_equal does not handle object arrays
def assert_array_almost_equal(x, y, decimal=6, err_msg='', verbose=True):
    from numpy.core import around, number, float_
    from numpy.lib import issubdtype
    from numpy.testing.utils import assert_array_compare
    def compare(x, y):
        z = abs(x-y)
        if not issubdtype(z.dtype, number):
            z = z.astype(float_) # handle object arrays
        return around(z, decimal) <= 10.0**(-decimal)
    assert_array_compare(compare, x, y, err_msg=err_msg, verbose=verbose,
                         header='Arrays are not almost equal')

class TestConvolve(TestCase):
    def test_basic(self):
        a = [3,4,5,6,5,4]
        b = [1,2,3]
        c = signal.convolve(a,b)
        assert_array_equal(c,array([3,10,22,28,32,32,23,12]))

class TestFFTConvolve(TestCase):
    def test_real(self):
        x = array([1,2,3])
        assert_array_almost_equal(signal.fftconvolve(x,x), [1,4,10,12,9.])

    def test_complex(self):
        x = array([1+1j,2+2j,3+3j])
        assert_array_almost_equal(signal.fftconvolve(x,x),
                                  [0+2.0j, 0+8j, 0+20j, 0+24j, 0+18j])

class TestMedFilt(TestCase):
    def test_basic(self):
        f = [[50, 50, 50, 50, 50, 92, 18, 27, 65, 46],
             [50, 50, 50, 50, 50,  0, 72, 77, 68, 66],
             [50, 50, 50, 50, 50, 46, 47, 19, 64, 77],
             [50, 50, 50, 50, 50, 42, 15, 29, 95, 35],
             [50, 50, 50, 50, 50, 46, 34,  9, 21, 66],
             [70, 97, 28, 68, 78, 77, 61, 58, 71, 42],
             [64, 53, 44, 29, 68, 32, 19, 68, 24, 84],
             [ 3, 33, 53, 67,  1, 78, 74, 55, 12, 83],
             [ 7, 11, 46, 70, 60, 47, 24, 43, 61, 26],
             [32, 61, 88,  7, 39,  4, 92, 64, 45, 61]]

        d = signal.medfilt(f, [7, 3])
        e = signal.medfilt2d(np.array(f, np.float), [7, 3])
        assert_array_equal(d, [[ 0, 50, 50, 50, 42, 15, 15, 18, 27,  0],
                               [ 0, 50, 50, 50, 50, 42, 19, 21, 29,  0],
                               [50, 50, 50, 50, 50, 47, 34, 34, 46, 35],
                               [50, 50, 50, 50, 50, 50, 42, 47, 64, 42],
                               [50, 50, 50, 50, 50, 50, 46, 55, 64, 35],
                               [33, 50, 50, 50, 50, 47, 46, 43, 55, 26],
                               [32, 50, 50, 50, 50, 47, 46, 45, 55, 26],
                               [ 7, 46, 50, 50, 47, 46, 46, 43, 45, 21],
                               [ 0, 32, 33, 39, 32, 32, 43, 43, 43,  0],
                               [ 0,  7, 11,  7,  4,  4, 19, 19, 24,  0]])
        assert_array_equal(d, e)

class TestWiener(TestCase):
    def test_basic(self):
        g = array([[5,6,4,3],[3,5,6,2],[2,3,5,6],[1,6,9,7]],'d')
        correct = array([[2.16374269,3.2222222222, 2.8888888889, 1.6666666667],[2.666666667, 4.33333333333, 4.44444444444, 2.8888888888],[2.222222222, 4.4444444444, 5.4444444444, 4.801066874837],[1.33333333333, 3.92735042735, 6.0712560386, 5.0404040404]])
        h = signal.wiener(g)
        assert_array_almost_equal(h,correct,decimal=6)

class TestCSpline1DEval(TestCase):
    def test_basic(self):
        y=array([1,2,3,4,3,2,1,2,3.0])
        x=arange(len(y))
        dx=x[1]-x[0]
        cj = signal.cspline1d(y)

        x2=arange(len(y)*10.0)/10.0
        y2=signal.cspline1d_eval(cj, x2, dx=dx,x0=x[0])

        # make sure interpolated values are on knot points
        assert_array_almost_equal(y2[::10], y, decimal=5)

class TestOrderFilt(TestCase):
    def test_basic(self):
        assert_array_equal(signal.order_filter([1,2,3],[1,0,1],1),
                           [2,3,2])

class TestChebWin:
    def test_cheb_odd(self):
        cheb_odd_true = array([0.200938, 0.107729, 0.134941, 0.165348,
                               0.198891, 0.235450, 0.274846, 0.316836,
                               0.361119, 0.407338, 0.455079, 0.503883,
                               0.553248, 0.602637, 0.651489, 0.699227,
                               0.745266, 0.789028, 0.829947, 0.867485,
                               0.901138, 0.930448, 0.955010, 0.974482,
                               0.988591, 0.997138, 1.000000, 0.997138,
                               0.988591, 0.974482, 0.955010, 0.930448,
                               0.901138, 0.867485, 0.829947, 0.789028,
                               0.745266, 0.699227, 0.651489, 0.602637,
                               0.553248, 0.503883, 0.455079, 0.407338,
                               0.361119, 0.316836, 0.274846, 0.235450,
                               0.198891, 0.165348, 0.134941, 0.107729,
                               0.200938])

        cheb_odd = signal.chebwin(53, at=-40)
        assert_array_almost_equal(cheb_odd, cheb_odd_true, decimal=4)

    def test_cheb_even(self):
        cheb_even_true = array([0.203894, 0.107279, 0.133904,
                                0.163608, 0.196338, 0.231986,
                                0.270385, 0.311313, 0.354493,
                                0.399594, 0.446233, 0.493983,
                                0.542378, 0.590916, 0.639071,
                                0.686302, 0.732055, 0.775783,
                                0.816944, 0.855021, 0.889525,
                                0.920006, 0.946060, 0.967339,
                                0.983557, 0.994494, 1.000000,
                                1.000000, 0.994494, 0.983557,
                                0.967339, 0.946060, 0.920006,
                                0.889525, 0.855021, 0.816944,
                                0.775783, 0.732055, 0.686302,
                                0.639071, 0.590916, 0.542378,
                                0.493983, 0.446233, 0.399594,
                                0.354493, 0.311313, 0.270385,
                                0.231986, 0.196338, 0.163608,
                                0.133904, 0.107279, 0.203894])

        cheb_even = signal.chebwin(54, at=-40)
        assert_array_almost_equal(cheb_even, cheb_even_true, decimal=4)

class _TestLinearFilter(TestCase):
    dt = None
    def test_rank1(self):
        x = np.linspace(0, 5, 6).astype(self.dt)
        b = np.array([1, -1]).astype(self.dt)
        a = np.array([0.5, -0.5]).astype(self.dt)

        # Test simple IIR
        y_r = np.array([0, 2, 4, 6, 8, 10.]).astype(self.dt)
        assert_array_almost_equal(lfilter(b, a, x), y_r)

        # Test simple FIR
        b = np.array([1, 1]).astype(self.dt)
        a = np.array([1]).astype(self.dt)
        y_r = np.array([0, 1, 3, 5, 7, 9.]).astype(self.dt)
        assert_array_almost_equal(lfilter(b, a, x), y_r)

        # Test IIR with initial conditions
        b = np.array([1, 1]).astype(self.dt)
        a = np.array([1]).astype(self.dt)
        zi = np.array([1]).astype(self.dt)
        y_r = np.array([1, 1, 3, 5, 7, 9.]).astype(self.dt)
        zf_r = np.array([5]).astype(self.dt)
        y, zf = lfilter(b, a, x, zi=zi)
        assert_array_almost_equal(y, y_r)
        assert_array_almost_equal(zf, zf_r)

        b = np.array([1, 1, 1]).astype(self.dt)
        a = np.array([1]).astype(self.dt)
        zi = np.array([1, 1]).astype(self.dt)
        y_r = np.array([1, 2, 3, 6, 9, 12.]).astype(self.dt)
        zf_r = np.array([9, 5]).astype(self.dt)
        y, zf = lfilter(b, a, x, zi=zi)
        assert_array_almost_equal(y, y_r)
        assert_array_almost_equal(zf, zf_r)

    def test_rank2(self):
        shape = (4, 3)
        x = np.linspace(0, np.prod(shape) - 1, np.prod(shape)).reshape(shape)
        x = x.astype(self.dt)

        b = np.array([1, -1]).astype(self.dt)
        a = np.array([0.5, 0.5]).astype(self.dt)

        y_r2_a0 = np.array([[0, 2, 4], [6, 4, 2], [0, 2, 4], [6 ,4 ,2]],
                           dtype=self.dt)

        y_r2_a1 = np.array([[0, 2, 0], [6, -4, 6], [12, -10, 12],
                            [18, -16, 18]], dtype=self.dt)

        y = lfilter(b, a, x, axis = 0)
        assert_array_almost_equal(y_r2_a0, y)

        y = lfilter(b, a, x, axis = 1)
        assert_array_almost_equal(y_r2_a1, y)

    def test_rank2_init_cond_a1(self):
        # Test initial condition handling along axis 1
        shape = (4, 3)
        x = np.linspace(0, np.prod(shape) - 1, np.prod(shape)).reshape(shape)
        x = x.astype(self.dt)

        b = np.array([1, -1]).astype(self.dt)
        a = np.array([0.5, 0.5]).astype(self.dt)

        y_r2_a0_1 = np.array([[1, 1, 1], [7, -5, 7], [13, -11, 13],
                              [19, -17, 19]], dtype=self.dt)
        zf_r = np.array([-5, -17, -29, -41])[:, np.newaxis].astype(self.dt)
        y, zf = lfilter(b, a, x, axis = 1, zi = np.ones((4, 1)))
        assert_array_almost_equal(y_r2_a0_1, y)
        assert_array_almost_equal(zf, zf_r)

    #@dec.skipif(True, "Skipping lfilter test with initial condition along "\
    #                  "axis 0: it segfaults ATM")
    def test_rank2_init_cond_a0(self):
        # Test initial condition handling along axis 0
        shape = (4, 3)
        x = np.linspace(0, np.prod(shape) - 1, np.prod(shape)).reshape(shape)
        x = x.astype(self.dt)

        b = np.array([1, -1]).astype(self.dt)
        a = np.array([0.5, 0.5]).astype(self.dt)

        y_r2_a0_0 = np.array([[1, 3, 5], [5, 3, 1], [1, 3, 5], [5 ,3 ,1]], 
                             dtype=self.dt)
        zf_r = np.array([[-23, -23, -23]], dtype=self.dt)
        y, zf = lfilter(b, a, x, axis = 0, zi = np.ones((1, 3)))
        assert_array_almost_equal(y_r2_a0_0, y)
        assert_array_almost_equal(zf, zf_r)

    #@dec.skipif(True, "Skipping rank > 2 test for lfilter because its segfaults ATM")
    def test_rank3(self):
        shape = (4, 3, 2)
        x = np.linspace(0, np.prod(shape) - 1, np.prod(shape)).reshape(shape)

        b = np.array([1, -1]).astype(self.dt)
        a = np.array([0.5, 0.5]).astype(self.dt)

        # Test last axis
        y = lfilter(b, a, x)
        for i in range(x.shape[0]):
            for j in range(x.shape[1]):
                assert_array_almost_equal(y[i, j], lfilter(b, a, x[i, j]))

    def test_empty_zi(self):
        """Regression test for #880: empty array for zi crashes."""
        a = np.ones(1).astype(self.dt)
        b = np.ones(1).astype(self.dt)
        x = np.arange(5).astype(self.dt)
        zi = np.ones(0).astype(self.dt)
        lfilter(b, a, x, zi=zi)

class TestLinearFilterFloat32(_TestLinearFilter):
    dt = np.float32

class TestLinearFilterFloat64(_TestLinearFilter):
    dt = np.float64

class TestLinearFilterComplex64(_TestLinearFilter):
    dt = np.complex64

class TestLinearFilterComplex128(_TestLinearFilter):
    dt = np.complex128

class TestLinearFilterDecimal(_TestLinearFilter):
    dt = np.dtype(Decimal)

if __name__ == "__main__":
    run_module_suite()
