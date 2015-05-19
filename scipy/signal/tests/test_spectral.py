from __future__ import division, print_function, absolute_import

import warnings
import numpy as np
from numpy.testing import assert_raises, assert_approx_equal, \
                          assert_, run_module_suite, TestCase,\
                          assert_allclose, assert_array_equal,\
                          assert_array_almost_equal_nulp, dec
from scipy import signal, fftpack
from scipy._lib._version import NumpyVersion
from scipy.signal import (periodogram, welch, lombscargle, csd, coherence,
                          spectrogram)


class TestPeriodogram(TestCase):
    def test_real_onesided_even(self):
        x = np.zeros(16)
        x[0] = 1
        f, p = periodogram(x)
        assert_allclose(f, np.linspace(0, 0.5, 9))
        q = np.ones(9)
        q[0] = 0
        q[-1] /= 2.0
        q /= 8
        assert_allclose(p, q)

    def test_real_onesided_odd(self):
        x = np.zeros(15)
        x[0] = 1
        f, p = periodogram(x)
        assert_allclose(f, np.arange(8.0)/15.0)
        q = np.ones(8)
        q[0] = 0
        q *= 2.0/15.0
        assert_allclose(p, q, atol=1e-15)

    def test_real_twosided(self):
        x = np.zeros(16)
        x[0] = 1
        f, p = periodogram(x, return_onesided=False)
        assert_allclose(f, fftpack.fftfreq(16, 1.0))
        q = np.ones(16)/16.0
        q[0] = 0
        assert_allclose(p, q)

    def test_real_spectrum(self):
        x = np.zeros(16)
        x[0] = 1
        f, p = periodogram(x, scaling='spectrum')
        g, q = periodogram(x, scaling='density')
        assert_allclose(f, np.linspace(0, 0.5, 9))
        assert_allclose(p, q/16.0)

    def test_integer_even(self):
        x = np.zeros(16, dtype=np.int)
        x[0] = 1
        f, p = periodogram(x)
        assert_allclose(f, np.linspace(0, 0.5, 9))
        q = np.ones(9)
        q[0] = 0
        q[-1] /= 2.0
        q /= 8
        assert_allclose(p, q)

    def test_integer_odd(self):
        x = np.zeros(15, dtype=np.int)
        x[0] = 1
        f, p = periodogram(x)
        assert_allclose(f, np.arange(8.0)/15.0)
        q = np.ones(8)
        q[0] = 0
        q *= 2.0/15.0
        assert_allclose(p, q, atol=1e-15)

    def test_integer_twosided(self):
        x = np.zeros(16, dtype=np.int)
        x[0] = 1
        f, p = periodogram(x, return_onesided=False)
        assert_allclose(f, fftpack.fftfreq(16, 1.0))
        q = np.ones(16)/16.0
        q[0] = 0
        assert_allclose(p, q)

    def test_complex(self):
        x = np.zeros(16, np.complex128)
        x[0] = 1.0 + 2.0j
        f, p = periodogram(x)
        assert_allclose(f, fftpack.fftfreq(16, 1.0))
        q = 5.0*np.ones(16)/16.0
        q[0] = 0
        assert_allclose(p, q)

    def test_unk_scaling(self):
        assert_raises(ValueError, periodogram, np.zeros(4, np.complex128),
                scaling='foo')

    def test_nd_axis_m1(self):
        x = np.zeros(20, dtype=np.float64)
        x = x.reshape((2,1,10))
        x[:,:,0] = 1.0
        f, p = periodogram(x)
        assert_array_equal(p.shape, (2, 1, 6))
        assert_array_almost_equal_nulp(p[0,0,:], p[1,0,:], 60)
        f0, p0 = periodogram(x[0,0,:])
        assert_array_almost_equal_nulp(p0[np.newaxis,:], p[1,:], 60)

    def test_nd_axis_0(self):
        x = np.zeros(20, dtype=np.float64)
        x = x.reshape((10,2,1))
        x[0,:,:] = 1.0
        f, p = periodogram(x, axis=0)
        assert_array_equal(p.shape, (6,2,1))
        assert_array_almost_equal_nulp(p[:,0,0], p[:,1,0], 60)
        f0, p0 = periodogram(x[:,0,0])
        assert_array_almost_equal_nulp(p0, p[:,1,0])

    def test_window_external(self):
        x = np.zeros(16)
        x[0] = 1
        f, p = periodogram(x, 10, 'hanning')
        win = signal.get_window('hanning', 16)
        fe, pe = periodogram(x, 10, win)
        assert_array_almost_equal_nulp(p, pe)
        assert_array_almost_equal_nulp(f, fe)

    def test_padded_fft(self):
        x = np.zeros(16)
        x[0] = 1
        f, p = periodogram(x)
        fp, pp = periodogram(x, nfft=32)
        assert_allclose(f, fp[::2])
        assert_allclose(p, pp[::2])
        assert_array_equal(pp.shape, (17,))

    def test_empty_input(self):
        f, p = periodogram([])
        assert_array_equal(f.shape, (0,))
        assert_array_equal(p.shape, (0,))
        for shape in [(0,), (3,0), (0,5,2)]:
            f, p = periodogram(np.empty(shape))
            assert_array_equal(f.shape, shape)
            assert_array_equal(p.shape, shape)

    def test_empty_input_other_axis(self):
        for shape in [(3,0), (0,5,2)]:
            f, p = periodogram(np.empty(shape), axis=1)
            assert_array_equal(f.shape, shape)
            assert_array_equal(p.shape, shape)

    def test_short_nfft(self):
        x = np.zeros(18)
        x[0] = 1
        f, p = periodogram(x, nfft=16)
        assert_allclose(f, np.linspace(0, 0.5, 9))
        q = np.ones(9)
        q[0] = 0
        q[-1] /= 2.0
        q /= 8
        assert_allclose(p, q)

    def test_nfft_is_xshape(self):
        x = np.zeros(16)
        x[0] = 1
        f, p = periodogram(x, nfft=16)
        assert_allclose(f, np.linspace(0, 0.5, 9))
        q = np.ones(9)
        q[0] = 0
        q[-1] /= 2.0
        q /= 8
        assert_allclose(p, q)

    def test_real_onesided_even_32(self):
        x = np.zeros(16, 'f')
        x[0] = 1
        f, p = periodogram(x)
        assert_allclose(f, np.linspace(0, 0.5, 9))
        q = np.ones(9, 'f')
        q[0] = 0
        q[-1] /= 2.0
        q /= 8
        assert_allclose(p, q)
        assert_(p.dtype == q.dtype)

    def test_real_onesided_odd_32(self):
        x = np.zeros(15, 'f')
        x[0] = 1
        f, p = periodogram(x)
        assert_allclose(f, np.arange(8.0)/15.0)
        q = np.ones(8, 'f')
        q[0] = 0
        q *= 2.0/15.0
        assert_allclose(p, q, atol=1e-7)
        assert_(p.dtype == q.dtype)

    @dec.skipif(NumpyVersion(np.__version__) < '1.8.0')
    def test_real_twosided_32(self):
        x = np.zeros(16, 'f')
        x[0] = 1
        f, p = periodogram(x, return_onesided=False)
        assert_allclose(f, fftpack.fftfreq(16, 1.0))
        q = np.ones(16, 'f')/16.0
        q[0] = 0
        assert_allclose(p, q)
        assert_(p.dtype == q.dtype)

    @dec.skipif(NumpyVersion(np.__version__) < '1.8.0')
    def test_complex_32(self):
        x = np.zeros(16, 'F')
        x[0] = 1.0 + 2.0j
        f, p = periodogram(x)
        assert_allclose(f, fftpack.fftfreq(16, 1.0))
        q = 5.0*np.ones(16, 'f')/16.0
        q[0] = 0
        assert_allclose(p, q)
        assert_(p.dtype == q.dtype)


class TestWelch(TestCase):
    def test_real_onesided_even(self):
        x = np.zeros(16)
        x[0] = 1
        x[8] = 1
        f, p = welch(x, nperseg=8)
        assert_allclose(f, np.linspace(0, 0.5, 5))
        q = np.array([0.08333333, 0.15277778, 0.22222222, 0.22222222,
                      0.11111111])
        assert_allclose(p, q, atol=1e-7, rtol=1e-7)

    def test_real_onesided_odd(self):
        x = np.zeros(16)
        x[0] = 1
        x[8] = 1
        f, p = welch(x, nperseg=9)
        assert_allclose(f, np.arange(5.0)/9.0)
        q = np.array([0.15958227, 0.24193957, 0.24145224, 0.24100919,
                      0.24377353])
        assert_allclose(p, q, atol=1e-7, rtol=1e-7)

    def test_real_twosided(self):
        x = np.zeros(16)
        x[0] = 1
        x[8] = 1
        f, p = welch(x, nperseg=8, return_onesided=False)
        assert_allclose(f, fftpack.fftfreq(8, 1.0))
        q = np.array([0.08333333, 0.07638889, 0.11111111, 0.11111111,
                      0.11111111, 0.11111111, 0.11111111, 0.07638889])
        assert_allclose(p, q, atol=1e-7, rtol=1e-7)

    def test_real_spectrum(self):
        x = np.zeros(16)
        x[0] = 1
        x[8] = 1
        f, p = welch(x, nperseg=8, scaling='spectrum')
        assert_allclose(f, np.linspace(0, 0.5, 5))
        q = np.array([0.015625, 0.02864583, 0.04166667, 0.04166667,
                      0.02083333])
        assert_allclose(p, q, atol=1e-7, rtol=1e-7)

    def test_integer_onesided_even(self):
        x = np.zeros(16, dtype=np.int)
        x[0] = 1
        x[8] = 1
        f, p = welch(x, nperseg=8)
        assert_allclose(f, np.linspace(0, 0.5, 5))
        q = np.array([0.08333333, 0.15277778, 0.22222222, 0.22222222,
                      0.11111111])
        assert_allclose(p, q, atol=1e-7, rtol=1e-7)

    def test_integer_onesided_odd(self):
        x = np.zeros(16, dtype=np.int)
        x[0] = 1
        x[8] = 1
        f, p = welch(x, nperseg=9)
        assert_allclose(f, np.arange(5.0)/9.0)
        q = np.array([0.15958227, 0.24193957, 0.24145224, 0.24100919,
                      0.24377353])
        assert_allclose(p, q, atol=1e-7, rtol=1e-7)

    def test_integer_twosided(self):
        x = np.zeros(16, dtype=np.int)
        x[0] = 1
        x[8] = 1
        f, p = welch(x, nperseg=8, return_onesided=False)
        assert_allclose(f, fftpack.fftfreq(8, 1.0))
        q = np.array([0.08333333, 0.07638889, 0.11111111, 0.11111111,
                      0.11111111, 0.11111111, 0.11111111, 0.07638889])
        assert_allclose(p, q, atol=1e-7, rtol=1e-7)

    def test_complex(self):
        x = np.zeros(16, np.complex128)
        x[0] = 1.0 + 2.0j
        x[8] = 1.0 + 2.0j
        f, p = welch(x, nperseg=8)
        assert_allclose(f, fftpack.fftfreq(8, 1.0))
        q = np.array([0.41666667, 0.38194444, 0.55555556, 0.55555556,
                      0.55555556, 0.55555556, 0.55555556, 0.38194444])
        assert_allclose(p, q, atol=1e-7, rtol=1e-7)

    def test_unk_scaling(self):
        assert_raises(ValueError, welch, np.zeros(4, np.complex128),
                      scaling='foo', nperseg=4)

    def test_detrend_linear(self):
        x = np.arange(10, dtype=np.float64) + 0.04
        f, p = welch(x, nperseg=10, detrend='linear')
        assert_allclose(p, np.zeros_like(p), atol=1e-15)

    def test_no_detrending(self):
        x = np.arange(10, dtype=np.float64) + 0.04
        f1, p1 = welch(x, nperseg=10, detrend=False)
        f2, p2 = welch(x, nperseg=10, detrend=lambda x: x)
        assert_allclose(f1, f2, atol=1e-15)
        assert_allclose(p1, p2, atol=1e-15)

    def test_detrend_external(self):
        x = np.arange(10, dtype=np.float64) + 0.04
        f, p = welch(x, nperseg=10,
                     detrend=lambda seg: signal.detrend(seg, type='l'))
        assert_allclose(p, np.zeros_like(p), atol=1e-15)

    def test_detrend_external_nd_m1(self):
        x = np.arange(40, dtype=np.float64) + 0.04
        x = x.reshape((2,2,10))
        f, p = welch(x, nperseg=10,
                     detrend=lambda seg: signal.detrend(seg, type='l'))
        assert_allclose(p, np.zeros_like(p), atol=1e-15)

    def test_detrend_external_nd_0(self):
        x = np.arange(20, dtype=np.float64) + 0.04
        x = x.reshape((2,1,10))
        x = np.rollaxis(x, 2, 0)
        f, p = welch(x, nperseg=10, axis=0,
                     detrend=lambda seg: signal.detrend(seg, axis=0, type='l'))
        assert_allclose(p, np.zeros_like(p), atol=1e-15)

    def test_nd_axis_m1(self):
        x = np.arange(20, dtype=np.float64) + 0.04
        x = x.reshape((2,1,10))
        f, p = welch(x, nperseg=10)
        assert_array_equal(p.shape, (2, 1, 6))
        assert_allclose(p[0,0,:], p[1,0,:], atol=1e-13, rtol=1e-13)
        f0, p0 = welch(x[0,0,:], nperseg=10)
        assert_allclose(p0[np.newaxis,:], p[1,:], atol=1e-13, rtol=1e-13)

    def test_nd_axis_0(self):
        x = np.arange(20, dtype=np.float64) + 0.04
        x = x.reshape((10,2,1))
        f, p = welch(x, nperseg=10, axis=0)
        assert_array_equal(p.shape, (6,2,1))
        assert_allclose(p[:,0,0], p[:,1,0], atol=1e-13, rtol=1e-13)
        f0, p0 = welch(x[:,0,0], nperseg=10)
        assert_allclose(p0, p[:,1,0], atol=1e-13, rtol=1e-13)

    def test_window_external(self):
        x = np.zeros(16)
        x[0] = 1
        x[8] = 1
        f, p = welch(x, 10, 'hanning', 8)
        win = signal.get_window('hanning', 8)
        fe, pe = welch(x, 10, win, 8)
        assert_array_almost_equal_nulp(p, pe)
        assert_array_almost_equal_nulp(f, fe)

    def test_empty_input(self):
        f, p = welch([])
        assert_array_equal(f.shape, (0,))
        assert_array_equal(p.shape, (0,))
        for shape in [(0,), (3,0), (0,5,2)]:
            f, p = welch(np.empty(shape))
            assert_array_equal(f.shape, shape)
            assert_array_equal(p.shape, shape)

    def test_empty_input_other_axis(self):
        for shape in [(3,0), (0,5,2)]:
            f, p = welch(np.empty(shape), axis=1)
            assert_array_equal(f.shape, shape)
            assert_array_equal(p.shape, shape)

    def test_short_data(self):
        x = np.zeros(8)
        x[0] = 1
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            f, p = welch(x)

        f1, p1 = welch(x, nperseg=8)
        assert_allclose(f, f1)
        assert_allclose(p, p1)

    def test_window_long_or_nd(self):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            assert_raises(ValueError, welch, np.zeros(4), 1,
                          np.array([1,1,1,1,1]))
            assert_raises(ValueError, welch, np.zeros(4), 1,
                          np.arange(6).reshape((2,3)))

    def test_nondefault_noverlap(self):
        x = np.zeros(64)
        x[::8] = 1
        f, p = welch(x, nperseg=16, noverlap=4)
        q = np.array([0, 1./12., 1./3., 1./5., 1./3., 1./5., 1./3., 1./5.,
                      1./6.])
        assert_allclose(p, q, atol=1e-12)

    def test_bad_noverlap(self):
        assert_raises(ValueError, welch, np.zeros(4), 1, 'hanning', 2, 7)

    def test_nfft_too_short(self):
        assert_raises(ValueError, welch, np.ones(12), nfft=3, nperseg=4)

    def test_real_onesided_even_32(self):
        x = np.zeros(16, 'f')
        x[0] = 1
        x[8] = 1
        f, p = welch(x, nperseg=8)
        assert_allclose(f, np.linspace(0, 0.5, 5))
        q = np.array([0.08333333, 0.15277778, 0.22222222, 0.22222222,
                      0.11111111], 'f')
        assert_allclose(p, q, atol=1e-7, rtol=1e-7)
        assert_(p.dtype == q.dtype)

    def test_real_onesided_odd_32(self):
        x = np.zeros(16, 'f')
        x[0] = 1
        x[8] = 1
        f, p = welch(x, nperseg=9)
        assert_allclose(f, np.arange(5.0)/9.0)
        q = np.array([0.15958227, 0.24193957, 0.24145224, 0.24100919,
                      0.24377353], 'f')
        assert_allclose(p, q, atol=1e-7, rtol=1e-7)
        assert_(p.dtype == q.dtype)

    @dec.skipif(NumpyVersion(np.__version__) < '1.8.0')
    def test_real_twosided_32(self):
        x = np.zeros(16, 'f')
        x[0] = 1
        x[8] = 1
        f, p = welch(x, nperseg=8, return_onesided=False)
        assert_allclose(f, fftpack.fftfreq(8, 1.0))
        q = np.array([0.08333333, 0.07638889, 0.11111111,
                      0.11111111, 0.11111111, 0.11111111, 0.11111111,
                      0.07638889], 'f')
        assert_allclose(p, q, atol=1e-7, rtol=1e-7)
        assert_(p.dtype == q.dtype)

    @dec.skipif(NumpyVersion(np.__version__) < '1.8.0')
    def test_complex_32(self):
        x = np.zeros(16, 'F')
        x[0] = 1.0 + 2.0j
        x[8] = 1.0 + 2.0j
        f, p = welch(x, nperseg=8)
        assert_allclose(f, fftpack.fftfreq(8, 1.0))
        q = np.array([0.41666666, 0.38194442, 0.55555552, 0.55555552,
                      0.55555558, 0.55555552, 0.55555552, 0.38194442], 'f')
        assert_allclose(p, q, atol=1e-7, rtol=1e-7)
        assert_(p.dtype == q.dtype,
                'dtype mismatch, %s, %s' % (p.dtype, q.dtype))

class TestCSD:
    def test_pad_shorter_x(self):
        x = np.zeros(8)
        y = np.zeros(12)

        f = np.linspace(0, 0.5, 7)
        c = np.zeros(7,dtype=np.complex128)
        f1, c1 = csd(x, y, nperseg=12)

        assert_allclose(f, f1)
        assert_allclose(c, c1)

    def test_pad_shorter_y(self):
        x = np.zeros(12)
        y = np.zeros(8)

        f = np.linspace(0, 0.5, 7)
        c = np.zeros(7,dtype=np.complex128)
        f1, c1 = csd(x, y, nperseg=12)

        assert_allclose(f, f1)
        assert_allclose(c, c1)

    def test_real_onesided_even(self):
        x = np.zeros(16)
        x[0] = 1
        x[8] = 1
        f, p = csd(x, x, nperseg=8)
        assert_allclose(f, np.linspace(0, 0.5, 5))
        q = np.array([0.08333333, 0.15277778, 0.22222222, 0.22222222,
                      0.11111111])
        assert_allclose(p, q, atol=1e-7, rtol=1e-7)

    def test_real_onesided_odd(self):
        x = np.zeros(16)
        x[0] = 1
        x[8] = 1
        f, p = csd(x, x, nperseg=9)
        assert_allclose(f, np.arange(5.0)/9.0)
        q = np.array([0.15958227, 0.24193957, 0.24145224, 0.24100919,
                      0.24377353])
        assert_allclose(p, q, atol=1e-7, rtol=1e-7)

    def test_real_twosided(self):
        x = np.zeros(16)
        x[0] = 1
        x[8] = 1
        f, p = csd(x, x, nperseg=8, return_onesided=False)
        assert_allclose(f, fftpack.fftfreq(8, 1.0))
        q = np.array([0.08333333, 0.07638889, 0.11111111, 0.11111111,
                      0.11111111, 0.11111111, 0.11111111, 0.07638889])
        assert_allclose(p, q, atol=1e-7, rtol=1e-7)

    def test_real_spectrum(self):
        x = np.zeros(16)
        x[0] = 1
        x[8] = 1
        f, p = csd(x, x, nperseg=8, scaling='spectrum')
        assert_allclose(f, np.linspace(0, 0.5, 5))
        q = np.array([0.015625, 0.02864583, 0.04166667, 0.04166667,
                      0.02083333])
        assert_allclose(p, q, atol=1e-7, rtol=1e-7)

    def test_integer_onesided_even(self):
        x = np.zeros(16, dtype=np.int)
        x[0] = 1
        x[8] = 1
        f, p = csd(x, x, nperseg=8)
        assert_allclose(f, np.linspace(0, 0.5, 5))
        q = np.array([0.08333333, 0.15277778, 0.22222222, 0.22222222,
                      0.11111111])
        assert_allclose(p, q, atol=1e-7, rtol=1e-7)

    def test_integer_onesided_odd(self):
        x = np.zeros(16, dtype=np.int)
        x[0] = 1
        x[8] = 1
        f, p = csd(x, x, nperseg=9)
        assert_allclose(f, np.arange(5.0)/9.0)
        q = np.array([0.15958227, 0.24193957, 0.24145224, 0.24100919,
                      0.24377353])
        assert_allclose(p, q, atol=1e-7, rtol=1e-7)

    def test_integer_twosided(self):
        x = np.zeros(16, dtype=np.int)
        x[0] = 1
        x[8] = 1
        f, p = csd(x, x, nperseg=8, return_onesided=False)
        assert_allclose(f, fftpack.fftfreq(8, 1.0))
        q = np.array([0.08333333, 0.07638889, 0.11111111, 0.11111111,
                      0.11111111, 0.11111111, 0.11111111, 0.07638889])
        assert_allclose(p, q, atol=1e-7, rtol=1e-7)

    def test_complex(self):
        x = np.zeros(16, np.complex128)
        x[0] = 1.0 + 2.0j
        x[8] = 1.0 + 2.0j
        f, p = csd(x, x, nperseg=8)
        assert_allclose(f, fftpack.fftfreq(8, 1.0))
        q = np.array([0.41666667, 0.38194444, 0.55555556, 0.55555556,
                      0.55555556, 0.55555556, 0.55555556, 0.38194444])
        assert_allclose(p, q, atol=1e-7, rtol=1e-7)

    def test_unk_scaling(self):
        assert_raises(ValueError, csd, np.zeros(4, np.complex128),
                      np.ones(4, np.complex128), scaling='foo', nperseg=4)

    def test_detrend_linear(self):
        x = np.arange(10, dtype=np.float64) + 0.04
        f, p = csd(x, x, nperseg=10, detrend='linear')
        assert_allclose(p, np.zeros_like(p), atol=1e-15)

    def test_no_detrending(self):
        x = np.arange(10, dtype=np.float64) + 0.04
        f1, p1 = csd(x, x, nperseg=10, detrend=False)
        f2, p2 = csd(x, x, nperseg=10, detrend=lambda x: x)
        assert_allclose(f1, f2, atol=1e-15)
        assert_allclose(p1, p2, atol=1e-15)

    def test_detrend_external(self):
        x = np.arange(10, dtype=np.float64) + 0.04
        f, p = csd(x, x, nperseg=10,
                   detrend=lambda seg: signal.detrend(seg, type='l'))
        assert_allclose(p, np.zeros_like(p), atol=1e-15)

    def test_detrend_external_nd_m1(self):
        x = np.arange(40, dtype=np.float64) + 0.04
        x = x.reshape((2,2,10))
        f, p = csd(x, x, nperseg=10,
                   detrend=lambda seg: signal.detrend(seg, type='l'))
        assert_allclose(p, np.zeros_like(p), atol=1e-15)

    def test_detrend_external_nd_0(self):
        x = np.arange(20, dtype=np.float64) + 0.04
        x = x.reshape((2,1,10))
        x = np.rollaxis(x, 2, 0)
        f, p = csd(x, x, nperseg=10, axis=0,
                   detrend=lambda seg: signal.detrend(seg, axis=0, type='l'))
        assert_allclose(p, np.zeros_like(p), atol=1e-15)

    def test_nd_axis_m1(self):
        x = np.arange(20, dtype=np.float64) + 0.04
        x = x.reshape((2,1,10))
        f, p = csd(x, x, nperseg=10)
        assert_array_equal(p.shape, (2, 1, 6))
        assert_allclose(p[0,0,:], p[1,0,:], atol=1e-13, rtol=1e-13)
        f0, p0 = csd(x[0,0,:], x[0,0,:], nperseg=10)
        assert_allclose(p0[np.newaxis,:], p[1,:], atol=1e-13, rtol=1e-13)

    def test_nd_axis_0(self):
        x = np.arange(20, dtype=np.float64) + 0.04
        x = x.reshape((10,2,1))
        f, p = csd(x, x, nperseg=10, axis=0)
        assert_array_equal(p.shape, (6,2,1))
        assert_allclose(p[:,0,0], p[:,1,0], atol=1e-13, rtol=1e-13)
        f0, p0 = csd(x[:,0,0], x[:,0,0], nperseg=10)
        assert_allclose(p0, p[:,1,0], atol=1e-13, rtol=1e-13)

    def test_window_external(self):
        x = np.zeros(16)
        x[0] = 1
        x[8] = 1
        f, p = csd(x, x, 10, 'hanning', 8)
        win = signal.get_window('hanning', 8)
        fe, pe = csd(x, x, 10, win, 8)
        assert_array_almost_equal_nulp(p, pe)
        assert_array_almost_equal_nulp(f, fe)

    def test_empty_input(self):
        f, p = csd([],np.zeros(10))
        assert_array_equal(f.shape, (0,))
        assert_array_equal(p.shape, (0,))

        f, p = csd(np.zeros(10),[])
        assert_array_equal(f.shape, (0,))
        assert_array_equal(p.shape, (0,))

        for shape in [(0,), (3,0), (0,5,2)]:
            f, p = csd(np.empty(shape), np.empty(shape))
            assert_array_equal(f.shape, shape)
            assert_array_equal(p.shape, shape)

        f, p = csd(np.ones(10), np.empty((5,0)))
        assert_array_equal(f.shape, (5,0))
        assert_array_equal(p.shape, (5,0))

        f, p = csd(np.empty((5,0)), np.ones(10))
        assert_array_equal(f.shape, (5,0))
        assert_array_equal(p.shape, (5,0))

    def test_empty_input_other_axis(self):
        for shape in [(3,0), (0,5,2)]:
            f, p = csd(np.empty(shape), np.empty(shape), axis=1)
            assert_array_equal(f.shape, shape)
            assert_array_equal(p.shape, shape)

        f, p = csd(np.empty((10,10,3)), np.zeros((10,0,1)), axis=1)
        assert_array_equal(f.shape, (10,0,3))
        assert_array_equal(p.shape, (10,0,3))

        f, p = csd(np.empty((10,0,1)), np.zeros((10,10,3)), axis=1)
        assert_array_equal(f.shape, (10,0,3))
        assert_array_equal(p.shape, (10,0,3))

    def test_short_data(self):
        x = np.zeros(8)
        x[0] = 1
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            f, p = csd(x, x)

        f1, p1 = csd(x, x, nperseg=8)
        assert_allclose(f, f1)
        assert_allclose(p, p1)

    def test_window_long_or_nd(self):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            assert_raises(ValueError, csd, np.zeros(4), np.ones(4), 1,
                          np.array([1,1,1,1,1]))
            assert_raises(ValueError, csd, np.zeros(4), np.ones(4), 1,
                          np.arange(6).reshape((2,3)))

    def test_nondefault_noverlap(self):
        x = np.zeros(64)
        x[::8] = 1
        f, p = csd(x, x, nperseg=16, noverlap=4)
        q = np.array([0, 1./12., 1./3., 1./5., 1./3., 1./5., 1./3., 1./5.,
                      1./6.])
        assert_allclose(p, q, atol=1e-12)

    def test_bad_noverlap(self):
        assert_raises(ValueError, csd, np.zeros(4), np.ones(4), 1, 'hanning',
                      2, 7)

    def test_nfft_too_short(self):
        assert_raises(ValueError, csd, np.ones(12), np.zeros(12), nfft=3,
                      nperseg=4)

    def test_real_onesided_even_32(self):
        x = np.zeros(16, 'f')
        x[0] = 1
        x[8] = 1
        f, p = csd(x, x, nperseg=8)
        assert_allclose(f, np.linspace(0, 0.5, 5))
        q = np.array([0.08333333, 0.15277778, 0.22222222, 0.22222222,
                      0.11111111], 'f')
        assert_allclose(p, q, atol=1e-7, rtol=1e-7)
        assert_(p.dtype == q.dtype)

    def test_real_onesided_odd_32(self):
        x = np.zeros(16, 'f')
        x[0] = 1
        x[8] = 1
        f, p = csd(x, x, nperseg=9)
        assert_allclose(f, np.arange(5.0)/9.0)
        q = np.array([0.15958227, 0.24193957, 0.24145224, 0.24100919,
                      0.24377353], 'f')
        assert_allclose(p, q, atol=1e-7, rtol=1e-7)
        assert_(p.dtype == q.dtype)

    @dec.skipif(NumpyVersion(np.__version__) < '1.8.0')
    def test_real_twosided_32(self):
        x = np.zeros(16, 'f')
        x[0] = 1
        x[8] = 1
        f, p = csd(x, x, nperseg=8, return_onesided=False)
        assert_allclose(f, fftpack.fftfreq(8, 1.0))
        q = np.array([0.08333333, 0.07638889, 0.11111111,
                      0.11111111, 0.11111111, 0.11111111, 0.11111111,
                      0.07638889], 'f')
        assert_allclose(p, q, atol=1e-7, rtol=1e-7)
        assert_(p.dtype == q.dtype)

    @dec.skipif(NumpyVersion(np.__version__) < '1.8.0')
    def test_complex_32(self):
        x = np.zeros(16, 'F')
        x[0] = 1.0 + 2.0j
        x[8] = 1.0 + 2.0j
        f, p = csd(x, x, nperseg=8)
        assert_allclose(f, fftpack.fftfreq(8, 1.0))
        q = np.array([0.41666666, 0.38194442, 0.55555552, 0.55555552,
                      0.55555558, 0.55555552, 0.55555552, 0.38194442], 'f')
        assert_allclose(p, q, atol=1e-7, rtol=1e-7)
        assert_(p.dtype == q.dtype,
                'dtype mismatch, %s, %s' % (p.dtype, q.dtype))


class TestCoherence:
    def test_identical_input(self):
        x = np.random.randn(20)
        y = np.copy(x)  # So `y is x` -> False

        f = np.linspace(0, 0.5, 6)
        C = np.ones(6)
        f1, C1 = coherence(x, y, nperseg=10)

        assert_allclose(f, f1)
        assert_allclose(C, C1)

    def test_phase_shifted_input(self):
        x = np.random.randn(20)
        y = -x

        f = np.linspace(0, 0.5, 6)
        C = np.ones(6)
        f1, C1 = coherence(x, y, nperseg=10)

        assert_allclose(f, f1)
        assert_allclose(C, C1)


class TestSpectrogram:
    def test_average_all_segments(self):
        x = np.random.randn(1024)

        fs = 1.0
        window = ('tukey', 0.25)
        nperseg = 16
        noverlap = 2

        f, _, P = spectrogram(x, fs, window, nperseg, noverlap)
        fw, Pw = welch(x, fs, window, nperseg, noverlap)
        assert_allclose(f, fw)
        assert_allclose(np.mean(P, axis=-1), Pw)

class TestLombscargle:
    def test_frequency(self):
        """Test if frequency location of peak corresponds to frequency of
        generated input signal.
        """

        # Input parameters
        ampl = 2.
        w = 1.
        phi = 0.5 * np.pi
        nin = 100
        nout = 1000
        p = 0.7  # Fraction of points to select

        # Randomly select a fraction of an array with timesteps
        np.random.seed(2353425)
        r = np.random.rand(nin)
        t = np.linspace(0.01*np.pi, 10.*np.pi, nin)[r >= p]

        # Plot a sine wave for the selected times
        x = ampl * np.sin(w*t + phi)

        # Define the array of frequencies for which to compute the periodogram
        f = np.linspace(0.01, 10., nout)

        # Calculate Lomb-Scargle periodogram
        P = lombscargle(t, x, f)

        # Check if difference between found frequency maximum and input
        # frequency is less than accuracy
        delta = f[1] - f[0]
        assert_(w - f[np.argmax(P)] < (delta/2.))

    def test_amplitude(self):
        """Test if height of peak in normalized Lomb-Scargle periodogram
        corresponds to amplitude of the generated input signal.
        """

        # Input parameters
        ampl = 2.
        w = 1.
        phi = 0.5 * np.pi
        nin = 100
        nout = 1000
        p = 0.7  # Fraction of points to select

        # Randomly select a fraction of an array with timesteps
        np.random.seed(2353425)
        r = np.random.rand(nin)
        t = np.linspace(0.01*np.pi, 10.*np.pi, nin)[r >= p]

        # Plot a sine wave for the selected times
        x = ampl * np.sin(w*t + phi)

        # Define the array of frequencies for which to compute the periodogram
        f = np.linspace(0.01, 10., nout)

        # Calculate Lomb-Scargle periodogram
        pgram = lombscargle(t, x, f)

        # Normalize
        pgram = np.sqrt(4 * pgram / t.shape[0])

        # Check if difference between found frequency maximum and input
        # frequency is less than accuracy
        assert_approx_equal(np.max(pgram), ampl, significant=2)

    def test_wrong_shape(self):
        t = np.linspace(0, 1, 1)
        x = np.linspace(0, 1, 2)
        f = np.linspace(0, 1, 3)
        assert_raises(ValueError, lombscargle, t, x, f)

    def test_zero_division(self):
        t = np.zeros(1)
        x = np.zeros(1)
        f = np.zeros(1)
        assert_raises(ZeroDivisionError, lombscargle, t, x, f)

    def test_lombscargle_atan_vs_atan2(self):
        # https://github.com/scipy/scipy/issues/3787
        # This raised a ZeroDivisionError.
        t = np.linspace(0, 10, 1000, endpoint=False)
        x = np.sin(4*t)
        f = np.linspace(0, 50, 500, endpoint=False) + 0.1
        q = lombscargle(t, x, f*2*np.pi)


if __name__ == "__main__":
    run_module_suite()
