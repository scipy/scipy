#this program corresponds to special.py
from decimal import Decimal
import types

from numpy.testing import *

import scipy.signal as signal
from scipy.signal import lfilter, correlate


from numpy import array, arange
import numpy as np

class TestConvolve(TestCase):
    def test_basic(self):
        a = [3,4,5,6,5,4]
        b = [1,2,3]
        c = signal.convolve(a,b)
        assert_array_equal(c,array([3,10,22,28,32,32,23,12]))

    def test_complex(self):
        x = array([1+1j, 2+1j, 3+1j])
        y = array([1+1j, 2+1j])
        z = signal.convolve(x, y)
        assert_array_equal(z, array([2j, 2+6j, 5+8j, 5+5j]))

    def test_zero_order(self):
        a = 1289
        b = 4567
        c = signal.convolve(a,b)
        assert_array_equal(c,a*b)

    def test_2d_arrays(self):
        a = [[1,2,3],[3,4,5]]
        b = [[2,3,4],[4,5,6]]
        c = signal.convolve(a,b)
        d = array(  [[2 ,7 ,16,17,12],\
                     [10,30,62,58,38],\
                     [12,31,58,49,30]])
        e = signal.convolve2d(a,b)
        assert_array_equal(c,d)
        assert_array_equal(e,d)


    def test_same_mode(self):
        a = [1,2,3,3,1,2]
        b = [1,4,3,4,5,6,7,4,3,2,1,1,3]
        c = signal.convolve(a,b,'same')
        d = array([14,25,35,43,57,61,63,57,45,36,25,20,17])
        assert_array_equal(c,d)
        #for the 2d function
        e = [[1,2,3],[3,4,5]]
        f = [[2,3,4,5,6,7,8],[4,5,6,7,8,9,10]]
        g = signal.convolve2d(e,f,'same')
        h = array([[ 7,16,22,28, 34, 40, 37],\
                   [30,62,80,98,116,134,114]])
        assert_array_equal(g,h)

    def test_valid_mode(self):
        a = [1,2,3,6,5,3]
        b = [2,3,4,5,3,4,2,2,1]
        c = signal.convolve(a,b,'valid')
        assert_array_equal(c,array([70,78,73,65]))
        #2d function
        e = [[1,2,3],[3,4,5]]
        f = [[2,3,4,5,6,7,8],[4,5,6,7,8,9,10]]
        g = signal.convolve2d(e,f,'valid')
        h = array([[62,80,98,116,134]])
        assert_array_equal(g,h)

    def test_fillvalue(self):
        a = [[1,2,3],[3,4,5]]
        b = [[2,3,4],[4,5,6]]
        fillval = 1
        c = signal.convolve2d(a,b,'full','fill',fillval)
        d = array([[24,26,31,34,32],\
                   [28,40,62,64,52],\
                   [32,46,67,62,48]])
        assert_array_equal(c,d)

    def test_wrap_boundary(self):
        a = [[1,2,3],[3,4,5]]
        b = [[2,3,4],[4,5,6]]
        c = signal.convolve2d(a,b,'full','wrap')
        d = array([[80,80,74,80,80],\
                   [68,68,62,68,68],\
                   [80,80,74,80,80]])
        assert_array_equal(c,d)

    def test_sym_boundary(self):
        a = [[1,2,3],[3,4,5]]
        b = [[2,3,4],[4,5,6]]
        c = signal.convolve2d(a,b,'full','symm')
        d = array([[34,30,44, 62, 66],\
                   [52,48,62, 80, 84],\
                   [82,78,92,110,114]])

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

if __name__ == "__main__":
    run_module_suite()
