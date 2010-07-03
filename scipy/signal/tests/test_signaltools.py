#this program corresponds to special.py
from decimal import Decimal
import types

from numpy.testing import *

import scipy.signal as signal
from scipy.signal import lfilter, correlate, convolve, convolve2d, hilbert


from numpy import array, arange
import numpy as np

class _TestConvolve(TestCase):
    def test_basic(self):
        a = [3,4,5,6,5,4]
        b = [1,2,3]
        c = convolve(a,b, old_behavior=self.old_behavior)
        assert_array_equal(c,array([3,10,22,28,32,32,23,12]))

    def test_complex(self):
        x = array([1+1j, 2+1j, 3+1j])
        y = array([1+1j, 2+1j])
        z = convolve(x, y,old_behavior=self.old_behavior)
        assert_array_equal(z, array([2j, 2+6j, 5+8j, 5+5j]))

    def test_zero_order(self):
        a = 1289
        b = 4567
        c = convolve(a,b,old_behavior=self.old_behavior)
        assert_array_equal(c,a*b)

    def test_2d_arrays(self):
        a = [[1,2,3],[3,4,5]]
        b = [[2,3,4],[4,5,6]]
        c = convolve(a,b,old_behavior=self.old_behavior)
        d = array(  [[2 ,7 ,16,17,12],\
                     [10,30,62,58,38],\
                     [12,31,58,49,30]])
        assert_array_equal(c,d)

    def test_valid_mode(self):
        a = [1,2,3,6,5,3]
        b = [2,3,4,5,3,4,2,2,1]
        c = convolve(a,b,'valid',old_behavior=self.old_behavior)
        assert_array_equal(c,array([70,78,73,65]))

class OldTestConvolve(_TestConvolve):
    old_behavior = True
    @dec.deprecated()
    def test_basic(self):
        _TestConvolve.test_basic(self)

    @dec.deprecated()
    def test_complex(self):
        _TestConvolve.test_complex(self)

    @dec.deprecated()
    def test_2d_arrays(self):
        _TestConvolve.test_2d_arrays(self)

    @dec.deprecated()
    def test_same_mode(self):
        _TestConvolve.test_same_mode(self)

    @dec.deprecated()
    def test_valid_mode(self):
        a = [1,2,3,6,5,3]
        b = [2,3,4,5,3,4,2,2,1]
        c = convolve(a,b,'valid',old_behavior=self.old_behavior)
        assert_array_equal(c,array([70,78,73,65]))

    @dec.deprecated()
    def test_same_mode(self):
        a = [1,2,3,3,1,2]
        b = [1,4,3,4,5,6,7,4,3,2,1,1,3]
        c = convolve(a,b,'same',old_behavior=self.old_behavior)
        d = array([14,25,35,43,57,61,63,57,45,36,25,20,17])
        assert_array_equal(c,d)

class TestConvolve(_TestConvolve):
    old_behavior = False
    def test_valid_mode(self):
        # 'valid' mode if b.size > a.size does not make sense with the new
        # behavior
        a = [1,2,3,6,5,3]
        b = [2,3,4,5,3,4,2,2,1]
        def _test():
            convolve(a,b,'valid',old_behavior=self.old_behavior)
        self.failUnlessRaises(ValueError, _test)

    def test_same_mode(self):
        a = [1,2,3,3,1,2]
        b = [1,4,3,4,5,6,7,4,3,2,1,1,3]
        c = convolve(a,b,'same',old_behavior=self.old_behavior)
        d = array([57,61,63,57,45,36])
        assert_array_equal(c,d)

class _TestConvolve2d(TestCase):
    def test_2d_arrays(self):
        a = [[1,2,3],[3,4,5]]
        b = [[2,3,4],[4,5,6]]
        d = array(  [[2 ,7 ,16,17,12],\
                     [10,30,62,58,38],\
                     [12,31,58,49,30]])
        e = convolve2d(a,b,old_behavior=self.old_behavior)
        assert_array_equal(e,d)

    def test_valid_mode(self):
        e = [[2,3,4,5,6,7,8],[4,5,6,7,8,9,10]]
        f = [[1,2,3],[3,4,5]]
        g = convolve2d(e,f,'valid',old_behavior=self.old_behavior)
        h = array([[62,80,98,116,134]])
        assert_array_equal(g,h)

    def test_fillvalue(self):
        a = [[1,2,3],[3,4,5]]
        b = [[2,3,4],[4,5,6]]
        fillval = 1
        c = convolve2d(a,b,'full','fill',fillval,old_behavior=self.old_behavior)
        d = array([[24,26,31,34,32],\
                   [28,40,62,64,52],\
                   [32,46,67,62,48]])
        assert_array_equal(c,d)

    def test_wrap_boundary(self):
        a = [[1,2,3],[3,4,5]]
        b = [[2,3,4],[4,5,6]]
        c = convolve2d(a,b,'full','wrap',old_behavior=self.old_behavior)
        d = array([[80,80,74,80,80],\
                   [68,68,62,68,68],\
                   [80,80,74,80,80]])
        assert_array_equal(c,d)

    def test_sym_boundary(self):
        a = [[1,2,3],[3,4,5]]
        b = [[2,3,4],[4,5,6]]
        c = convolve2d(a,b,'full','symm',old_behavior=self.old_behavior)
        d = array([[34,30,44, 62, 66],\
                   [52,48,62, 80, 84],\
                   [82,78,92,110,114]])
        assert_array_equal(c,d)


class OldTestConvolve2d(_TestConvolve2d):
    old_behavior = True
    @dec.deprecated()
    def test_2d_arrays(self):
        _TestConvolve2d.test_2d_arrays(self)

    @dec.deprecated()
    def test_same_mode(self):
        e = [[1,2,3],[3,4,5]]
        f = [[2,3,4,5,6,7,8],[4,5,6,7,8,9,10]]
        g = convolve2d(e,f,'same',old_behavior=self.old_behavior)
        h = array([[ 7,16,22,28, 34, 40, 37],\
                   [30,62,80,98,116,134,114]])
        assert_array_equal(g,h)

    @dec.deprecated()
    def test_valid_mode(self):
        _TestConvolve2d.test_valid_mode(self)

    @dec.deprecated()
    def test_fillvalue(self):
        _TestConvolve2d.test_fillvalue(self)

    @dec.deprecated()
    def test_wrap_boundary(self):
        _TestConvolve2d.test_wrap_boundary(self)

    @dec.deprecated()
    def test_sym_boundary(self):
        _TestConvolve2d.test_sym_boundary(self)

    @dec.deprecated()
    def test_valid_mode2(self):
        # Test when in2.size > in1.size: old behavior is to do so that
        # convolve2d(in2, in1) == convolve2d(in1, in2)
        e = [[1,2,3],[3,4,5]]
        f = [[2,3,4,5,6,7,8],[4,5,6,7,8,9,10]]
        g = convolve2d(e,f,'valid',old_behavior=self.old_behavior)
        h = array([[62,80,98,116,134]])
        assert_array_equal(g,h)

#class TestConvolve2d(_TestConvolve2d):
#    old_behavior = False
#    def test_same_mode(self):
#        e = [[1,2,3],[3,4,5]]
#        f = [[2,3,4,5,6,7,8],[4,5,6,7,8,9,10]]
#        g = convolve2d(e,f,'same',old_behavior=self.old_behavior)
#        h = array([[80,98,116],\
#                   [70,82,94]])
#        assert_array_equal(g,h)
#
#    def test_valid_mode2(self):
#        # Test when in2.size > in1.size
#        e = [[1,2,3],[3,4,5]]
#        f = [[2,3,4,5,6,7,8],[4,5,6,7,8,9,10]]
#        def _test():
#            convolve2d(e,f,'valid',old_behavior=self.old_behavior)
#        self.failUnlessRaises(ValueError, _test)

class TestFFTConvolve(TestCase):
    def test_real(self):
        x = array([1,2,3])
        assert_array_almost_equal(signal.fftconvolve(x,x), [1,4,10,12,9.])

    def test_complex(self):
        x = array([1+1j,2+2j,3+3j])
        assert_array_almost_equal(signal.fftconvolve(x,x),
                                  [0+2.0j, 0+8j, 0+20j, 0+24j, 0+18j])

    def test_2d_real_same(self):
        a = array([[1,2,3],[4,5,6]])
        assert_array_almost_equal(signal.fftconvolve(a,a),\
                                                array([[1,4,10,12,9],\
                                                       [8,26,56,54,36],\
                                                       [16,40,73,60,36]]))

    def test_2d_complex_same(self):
        a = array([[1+2j,3+4j,5+6j],[2+1j,4+3j,6+5j]])
        c = signal.fftconvolve(a,a)
        d = array([[-3+4j,-10+20j,-21+56j,-18+76j,-11+60j],\
                   [10j,44j,118j,156j,122j],\
                   [3+4j,10+20j,21+56j,18+76j,11+60j]])
        assert_array_almost_equal(c,d)

    def test_real_same_mode(self):
        a = array([1,2,3])
        b = array([3,3,5,6,8,7,9,0,1])
        c = signal.fftconvolve(a,b,'same')
        d = array([9.,20.,25.,35.,41.,47.,39.,28.,2.])
        assert_array_almost_equal(c,d)

    def test_real_valid_mode(self):
        a = array([3,2,1])
        b = array([3,3,5,6,8,7,9,0,1])
        c = signal.fftconvolve(a,b,'valid')
        d = array([24.,31.,41.,43.,49.,25.,12.])
        assert_array_almost_equal(c,d)

    def test_zero_order(self):
        a = array([4967])
        b = array([3920])
        c = signal.fftconvolve(a,b)
        d = a*b
        assert_equal(c,d)

    def test_random_data(self):
        np.random.seed(1234)
        a = np.random.rand(1233) + 1j*np.random.rand(1233)
        b = np.random.rand(1321) + 1j*np.random.rand(1321)
        c = signal.fftconvolve(a, b, 'full')
        d = np.convolve(a, b, 'full')
        assert np.allclose(c, d, rtol=1e-10)

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

    def test_none(self):
        """Ticket #1124. Ensure this does not segfault."""
        try:
            signal.medfilt(None)
        except:
            pass

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
        y, zf = lfilter(b, a, x, zi=zi)
        assert_array_almost_equal(y, x)
        self.failUnless(zf.dtype == self.dt)
        self.failUnless(zf.size == 0)

class TestLinearFilterFloat32(_TestLinearFilter):
    dt = np.float32

class TestLinearFilterFloat64(_TestLinearFilter):
    dt = np.float64

class TestLinearFilterFloatExtended(_TestLinearFilter):
    dt = np.longdouble

class TestLinearFilterComplex64(_TestLinearFilter):
    dt = np.complex64

class TestLinearFilterComplex128(_TestLinearFilter):
    dt = np.complex128

class TestLinearFilterComplexxxiExtended28(_TestLinearFilter):
    dt = np.longcomplex

class TestLinearFilterDecimal(_TestLinearFilter):
    dt = np.dtype(Decimal)

class _TestCorrelateReal(TestCase):
    dt = None
    def _setup_rank1(self):
        # a.size should be greated than b.size for the tests
        a = np.linspace(0, 3, 4).astype(self.dt)
        b = np.linspace(1, 2, 2).astype(self.dt)

        y_r = np.array([0, 2, 5, 8, 3]).astype(self.dt)
        return a, b, y_r

    def test_rank1_valid(self):
        a, b, y_r = self._setup_rank1()
        y = correlate(a, b, 'valid', old_behavior=False)
        assert_array_almost_equal(y, y_r[1:4])
        self.failUnless(y.dtype == self.dt)

    def test_rank1_same(self):
        a, b, y_r = self._setup_rank1()
        y = correlate(a, b, 'same', old_behavior=False)
        assert_array_almost_equal(y, y_r[:-1])
        self.failUnless(y.dtype == self.dt)

    def test_rank1_full(self):
        a, b, y_r = self._setup_rank1()
        y = correlate(a, b, 'full', old_behavior=False)
        assert_array_almost_equal(y, y_r)
        self.failUnless(y.dtype == self.dt)

    @dec.deprecated()
    def test_rank1_valid_old(self):
        # This test assume a.size > b.size
        a, b, y_r = self._setup_rank1()
        y = correlate(b, a, 'valid')
        assert_array_almost_equal(y, y_r[1:4])
        self.failUnless(y.dtype == self.dt)

    @dec.deprecated()
    def test_rank1_same_old(self):
        # This test assume a.size > b.size
        a, b, y_r = self._setup_rank1()
        y = correlate(b, a, 'same')
        assert_array_almost_equal(y, y_r[:-1])
        self.failUnless(y.dtype == self.dt)

    @dec.deprecated()
    def test_rank1_full_old(self):
        # This test assume a.size > b.size
        a, b, y_r = self._setup_rank1()
        y = correlate(b, a, 'full')
        assert_array_almost_equal(y, y_r)
        self.failUnless(y.dtype == self.dt)

    def _setup_rank3(self):
        a = np.linspace(0, 39, 40).reshape((2, 4, 5), order='F').astype(self.dt)
        b = np.linspace(0, 23, 24).reshape((2, 3, 4), order='F').astype(self.dt)

        y_r = array([[[    0.,   184.,   504.,   912.,  1360.,   888.,   472.,   160.,],
            [   46.,   432.,  1062.,  1840.,  2672.,  1698.,   864.,   266.,],
            [  134.,   736.,  1662.,  2768.,  3920.,  2418.,  1168.,   314.,],
            [  260.,   952.,  1932.,  3056.,  4208.,  2580.,  1240.,   332.,] ,
            [  202.,   664.,  1290.,  1984.,  2688.,  1590.,   712.,   150.,] ,
            [  114.,   344.,   642.,   960.,  1280.,   726.,   296.,    38.,]],

            [[   23.,   400.,  1035.,  1832.,  2696.,  1737.,   904.,   293.,],
             [  134.,   920.,  2166.,  3680.,  5280.,  3306.,  1640.,   474.,],
             [  325.,  1544.,  3369.,  5512.,  7720.,  4683.,  2192.,   535.,],
             [  571.,  1964.,  3891.,  6064.,  8272.,  4989.,  2324.,   565.,],
             [  434.,  1360.,  2586.,  3920.,  5264.,  3054.,  1312.,   230.,],
             [  241.,   700.,  1281.,  1888.,  2496.,  1383.,   532.,    39.,]],

            [[   22.,   214.,   528.,   916.,  1332.,   846.,   430.,   132.,],
             [   86.,   484.,  1098.,  1832.,  2600.,  1602.,   772.,   206.,],
             [  188.,   802.,  1698.,  2732.,  3788.,  2256.,  1018.,   218.,],
             [  308.,  1006.,  1950.,  2996.,  4052.,  2400.,  1078.,   230.,],
             [  230.,   692.,  1290.,  1928.,  2568.,  1458.,   596.,    78.,],
             [  126.,   354.,   636.,   924.,  1212.,   654.,   234.,     0.,]]],
            dtype=self.dt)

        return a, b, y_r

    def test_rank3_valid(self):
        a, b, y_r = self._setup_rank3()
        y = correlate(a, b, "valid", old_behavior=False)
        assert_array_almost_equal(y, y_r[1:2,2:4,3:5])
        self.failUnless(y.dtype == self.dt)

    def test_rank3_same(self):
        a, b, y_r = self._setup_rank3()
        y = correlate(a, b, "same", old_behavior=False)
        assert_array_almost_equal(y, y_r[0:-1,1:-1,1:-2])
        self.failUnless(y.dtype == self.dt)

    def test_rank3_all(self):
        a, b, y_r = self._setup_rank3()
        y = correlate(a, b, old_behavior=False)
        assert_array_almost_equal(y, y_r)
        self.failUnless(y.dtype == self.dt)

    @dec.deprecated()
    def test_rank3_valid_old(self):
        a, b, y_r = self._setup_rank3()
        y = correlate(b, a, "valid")
        assert_array_almost_equal(y, y_r[1:2,2:4,3:5])
        self.failUnless(y.dtype == self.dt)

    @dec.deprecated()
    def test_rank3_same_old(self):
        a, b, y_r = self._setup_rank3()
        y = correlate(b, a, "same")
        assert_array_almost_equal(y, y_r[0:-1,1:-1,1:-2])
        self.failUnless(y.dtype == self.dt)

    @dec.deprecated()
    def test_rank3_all_old(self):
        a, b, y_r = self._setup_rank3()
        y = correlate(b, a)
        assert_array_almost_equal(y, y_r)
        self.failUnless(y.dtype == self.dt)

for i in [np.ubyte, np.byte, np.ushort, np.short, np.uint, np.int,
        np.ulonglong, np.ulonglong, np.float32, np.float64, np.longdouble,
        Decimal]:
    name = "TestCorrelate%s" % i.__name__.title()
    globals()[name] = types.ClassType(name, (_TestCorrelateReal,), {"dt": i})

class _TestCorrelateComplex(TestCase):
    dt = None
    def _setup_rank1(self, mode):
        a = np.random.randn(10).astype(self.dt)
        a += 1j * np.random.randn(10).astype(self.dt)
        b = np.random.randn(8).astype(self.dt)
        b += 1j * np.random.randn(8).astype(self.dt)

        y_r = (correlate(a.real, b.real, mode=mode, old_behavior=False) +
               correlate(a.imag, b.imag, mode=mode, old_behavior=False)).astype(self.dt)
        y_r += 1j * (-correlate(a.real, b.imag, mode=mode, old_behavior=False) +
                correlate(a.imag, b.real, mode=mode, old_behavior=False))
        return a, b, y_r

    def test_rank1_valid(self):
        a, b, y_r = self._setup_rank1('valid')
        y = correlate(a, b, 'valid', old_behavior=False)
        assert_array_almost_equal(y, y_r)
        self.failUnless(y.dtype == self.dt)

    def test_rank1_same(self):
        a, b, y_r = self._setup_rank1('same')
        y = correlate(a, b, 'same', old_behavior=False)
        assert_array_almost_equal(y, y_r)
        self.failUnless(y.dtype == self.dt)

    def test_rank1_full(self):
        a, b, y_r = self._setup_rank1('full')
        y = correlate(a, b, 'full', old_behavior=False)
        assert_array_almost_equal(y, y_r)
        self.failUnless(y.dtype == self.dt)

    def test_rank3(self):
        a = np.random.randn(10, 8, 6).astype(self.dt)
        a += 1j * np.random.randn(10, 8, 6).astype(self.dt)
        b = np.random.randn(8, 6, 4).astype(self.dt)
        b += 1j * np.random.randn(8, 6, 4).astype(self.dt)

        y_r = (correlate(a.real, b.real, old_behavior=False)
                + correlate(a.imag, b.imag, old_behavior=False)).astype(self.dt)
        y_r += 1j * (-correlate(a.real, b.imag, old_behavior=False) +
                correlate(a.imag, b.real, old_behavior=False))

        y = correlate(a, b, 'full', old_behavior=False)
        assert_array_almost_equal(y, y_r, decimal=4)
        self.failUnless(y.dtype == self.dt)

    @dec.deprecated()
    def test_rank1_valid_old(self):
        a, b, y_r = self._setup_rank1('valid')
        y = correlate(b, a.conj(), 'valid')
        assert_array_almost_equal(y, y_r)
        self.failUnless(y.dtype == self.dt)

    @dec.deprecated()
    def test_rank1_same_old(self):
        a, b, y_r = self._setup_rank1('same')
        y = correlate(b, a.conj(), 'same')
        assert_array_almost_equal(y, y_r)
        self.failUnless(y.dtype == self.dt)

    @dec.deprecated()
    def test_rank1_full_old(self):
        a, b, y_r = self._setup_rank1('full')
        y = correlate(b, a.conj(), 'full')
        assert_array_almost_equal(y, y_r)
        self.failUnless(y.dtype == self.dt)

    @dec.deprecated()
    def test_rank3_old(self):
        a = np.random.randn(10, 8, 6).astype(self.dt)
        a += 1j * np.random.randn(10, 8, 6).astype(self.dt)
        b = np.random.randn(8, 6, 4).astype(self.dt)
        b += 1j * np.random.randn(8, 6, 4).astype(self.dt)

        y_r = (correlate(a.real, b.real, old_behavior=False)
                + correlate(a.imag, b.imag, old_behavior=False)).astype(self.dt)
        y_r += 1j * (-correlate(a.real, b.imag, old_behavior=False) +
                correlate(a.imag, b.real, old_behavior=False))

        y = correlate(b, a.conj(), 'full')
        assert_array_almost_equal(y, y_r, decimal=4)
        self.failUnless(y.dtype == self.dt)

for i in [np.csingle, np.cdouble, np.clongdouble]:
    name = "TestCorrelate%s" % i.__name__.title()
    globals()[name] = types.ClassType(name, (_TestCorrelateComplex,), {"dt": i})

class TestFiltFilt:
    def test_basic(self):
        out = signal.filtfilt([1,2,3], [1,2,3], np.arange(12))
        assert_equal(out, arange(12))

class TestDecimate:
    def test_basic(self):
        x = np.arange(6)
        assert_array_equal(signal.decimate(x, 2, n=1).round(), x[::2])


class TestHilbert:
    def test_hilbert_theoretical(self):
        #test cases by Ariel Rokem
        decimal = 14

        pi = np.pi
        t = np.arange(0, 2*pi, pi/256)
        a0 = np.sin(t)
        a1 = np.cos(t)
        a2 = np.sin(2*t)
        a3 = np.cos(2*t)
        a = np.vstack([a0,a1,a2,a3])

        h = hilbert(a)
        h_abs = np.abs(h)
        h_angle = np.angle(h)
        h_real = np.real(h)

        #The real part should be equal to the original signals:
        assert_almost_equal(h_real, a, decimal)
        #The absolute value should be one everywhere, for this input:
        assert_almost_equal(h_abs, np.ones(a.shape), decimal)
        #For the 'slow' sine - the phase should go from -pi/2 to pi/2 in
        #the first 256 bins:
        assert_almost_equal(h_angle[0,:256], np.arange(-pi/2,pi/2,pi/256),
                            decimal)
        #For the 'slow' cosine - the phase should go from 0 to pi in the
        #same interval:
        assert_almost_equal(h_angle[1,:256], np.arange(0,pi,pi/256), decimal)
        #The 'fast' sine should make this phase transition in half the time:
        assert_almost_equal(h_angle[2,:128], np.arange(-pi/2,pi/2,pi/128),
                            decimal)
        #Ditto for the 'fast' cosine:
        assert_almost_equal(h_angle[3,:128], np.arange(0,pi,pi/128), decimal)

        #The imaginary part of hilbert(cos(t)) = sin(t) Wikipedia
        assert_almost_equal(h[1].imag, a0, decimal)

    def test_hilbert_axisN(self):
        # tests for axis and N arguments
        a = np.arange(18).reshape(3,6)
        # test axis
        aa = hilbert(a, axis=-1)
        yield assert_equal, hilbert(a.T, axis=0), aa.T
        # test 1d
        yield assert_equal, hilbert(a[0]), aa[0]

        # test N
        aan = hilbert(a, N=20, axis=-1)
        yield assert_equal, aan.shape, [3,20]
        yield assert_equal, hilbert(a.T, N=20, axis=0).shape, [20,3]
        #the next test is just a regression test,
        #no idea whether numbers make sense
        a0hilb = np.array(
                          [  0.000000000000000e+00-1.72015830311905j ,
                             1.000000000000000e+00-2.047794505137069j,
                             1.999999999999999e+00-2.244055555687583j,
                             3.000000000000000e+00-1.262750302935009j,
                             4.000000000000000e+00-1.066489252384493j,
                             5.000000000000000e+00+2.918022706971047j,
                             8.881784197001253e-17+3.845658908989067j,
                            -9.444121133484362e-17+0.985044202202061j,
                            -1.776356839400251e-16+1.332257797702019j,
                            -3.996802888650564e-16+0.501905089898885j,
                             1.332267629550188e-16+0.668696078880782j,
                            -1.192678053963799e-16+0.235487067862679j,
                            -1.776356839400251e-16+0.286439612812121j,
                             3.108624468950438e-16+0.031676888064907j,
                             1.332267629550188e-16-0.019275656884536j,
                            -2.360035624836702e-16-0.1652588660287j  ,
                             0.000000000000000e+00-0.332049855010597j,
                             3.552713678800501e-16-0.403810179797771j,
                             8.881784197001253e-17-0.751023775297729j,
                             9.444121133484362e-17-0.79252210110103j ])
        yield assert_almost_equal, aan[0], a0hilb, 14, 'N regression'


if __name__ == "__main__":
    run_module_suite()
