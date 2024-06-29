import sys

import numpy as np
from numpy.testing import (assert_,
                           assert_allclose, assert_array_equal,
                           assert_array_almost_equal_nulp, suppress_warnings)
import pytest
from pytest import raises as assert_raises

from scipy import signal
from scipy.fft import fftfreq, rfftfreq, irfft, fft
from scipy.signal import (periodogram, welch, lombscargle, coherence,
                          spectrogram, check_COLA, check_NOLA)
from scipy.signal.windows import hann
from scipy.signal._spectral_py import _spectral_helper

# Compare ShortTimeFFT.stft() / ShortTimeFFT.istft() with stft() / istft():
from scipy.signal.tests._scipy_spectral_test_shim import stft_compare as stft
from scipy.signal.tests._scipy_spectral_test_shim import istft_compare as istft
from scipy.signal.tests._scipy_spectral_test_shim import csd_compare as csd
from scipy.conftest import array_api_compatible
from scipy._lib._array_api import xp_assert_close, xp_assert_equal, copy, size

pytestmark = [array_api_compatible, pytest.mark.usefixtures("skip_xp_backends")]
skip_xp_backends = pytest.mark.skip_xp_backends

class TestPeriodogram:
    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis not available in array_api_strict"])
    def test_real_onesided_even(self, xp):
        x = xp.zeros(16)
        x[0] = 1
        f, p = periodogram(x)
        xp_assert_close(f, xp.linspace(0, 0.5, 9))
        q = xp.ones(9)
        q[0] = 0
        q[-1] /= 2.0
        q /= 8
        xp_assert_close(p, q)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis not available in array_api_strict"])
    def test_real_onesided_odd(self, xp):
        x = xp.zeros(15)
        x[0] = 1
        f, p = periodogram(x)
        xp_assert_close(f, xp.arange(8.0)/15.0)
        q = xp.ones(8)
        q[0] = 0
        q *= 2.0/15.0
        xp_assert_close(p, q, atol=1e-15)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis not available in array_api_strict"])
    def test_real_twosided(self, xp):
        x = xp.zeros(16)
        x[0] = 1
        f, p = periodogram(x, return_onesided=False)
        xp_assert_close(f, fftfreq(16, 1.0), check_namespace=False, check_dtype=False)
        q = xp.full((16,), 1/16.0)
        q[0] = 0
        xp_assert_close(p, q)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis not available in array_api_strict"])
    def test_real_spectrum(self, xp):
        x = xp.zeros(16)
        x[0] = 1
        f, p = periodogram(x, scaling='spectrum')
        g, q = periodogram(x, scaling='density')
        xp_assert_close(f, xp.linspace(0, 0.5, 9))
        xp_assert_close(p, q/16.0)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "dtype casting issue for array_api_strict"])
    def test_integer_even(self, xp):
        x = xp.zeros(16, dtype=xp.int64)
        x[0] = 1
        f, p = periodogram(x)
        xp_assert_close(f, xp.linspace(0, 0.5, 9))
        q = xp.ones(9)
        q[0] = 0
        q[-1] /= 2.0
        q /= 8
        xp_assert_close(p, q)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "dtype casting issue for array_api_strict"])
    def test_integer_odd(self, xp):
        x = xp.zeros(15, dtype=xp.int64)
        x[0] = 1
        f, p = periodogram(x)
        xp_assert_close(f, xp.arange(8.0)/15.0)
        q = xp.ones(8)
        q[0] = 0
        q *= 2.0/15.0
        xp_assert_close(p, q, atol=1e-15)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "dtype casting issue for array_api_strict"])
    def test_integer_twosided(self, xp):
        x = xp.zeros(16, dtype=xp.int64)
        x[0] = 1
        f, p = periodogram(x, return_onesided=False)
        xp_assert_close(f, fftfreq(16, 1.0), check_namespace=False, check_dtype=False)
        q = xp.full((16,), 1/16.0)
        q[0] = 0
        xp_assert_close(p, q)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "xp.mean() requires real types for array_api_strict"])
    def test_complex(self, xp):
        x = xp.zeros(16, dtype=xp.complex128)
        x[0] = 1.0 + 2.0j
        f, p = periodogram(x, return_onesided=False)
        xp_assert_close(f, fftfreq(16, 1.0), check_namespace=False, check_dtype=False)
        q = xp.full((16,), 5.0/16.0)
        q[0] = 0
        xp_assert_close(p, q, check_dtype=False)

    def test_unk_scaling(self, xp):
        assert_raises(ValueError, periodogram, xp.zeros(4, dtype=xp.complex128),
                scaling='foo')

    @pytest.mark.skipif(
        sys.maxsize <= 2**32,
        reason="On some 32-bit tolerance issue"
    )
    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "no moveaxis in array_api_strict",
                               "torch hits nulp device coercion"])
    def test_nd_axis_m1(self, xp):
        x = xp.zeros(20, dtype=xp.float64)
        x = xp.reshape(x, (2, 1, 10))
        x[:,:,0] = 1.0
        f, p = periodogram(x)
        assert p.shape == (2, 1, 6)
        assert_array_almost_equal_nulp(p[0,0,:], p[1,0,:], 60)
        f0, p0 = periodogram(x[0,0,:])
        assert_array_almost_equal_nulp(p0[xp.newaxis,:], p[1,:], 60)

    @pytest.mark.skipif(
        sys.maxsize <= 2**32,
        reason="On some 32-bit tolerance issue"
    )
    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "no moveaxis in array_api_strict",
                               "torch hits nulp device coercion"])
    def test_nd_axis_0(self, xp):
        x = xp.zeros(20, dtype=xp.float64)
        x = xp.reshape(x, (10, 2, 1))
        x[0,:,:] = 1.0
        f, p = periodogram(x, axis=0)
        assert p.shape == (6, 2, 1)
        assert_array_almost_equal_nulp(p[:,0,0], p[:,1,0], 60)
        f0, p0 = periodogram(x[:,0,0])
        assert_array_almost_equal_nulp(p0, p[:,1,0])

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "no moveaxis in array_api_strict",
                               "torch hits nulp device coercion"])
    def test_window_external(self, xp):
        x = xp.zeros(16)
        x[0] = 1
        f, p = periodogram(x, 10, 'hann')
        win = signal.get_window('hann', 16)
        fe, pe = periodogram(x, 10, win)
        assert_array_almost_equal_nulp(p, pe)
        assert_array_almost_equal_nulp(f, fe)
        win_err = signal.get_window('hann', 32)
        assert_raises(ValueError, periodogram, x,
                      10, win_err)  # win longer than signal

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "no moveaxis in array_api_strict"])
    def test_padded_fft(self, xp):
        x = xp.zeros(16)
        x[0] = 1
        f, p = periodogram(x)
        fp, pp = periodogram(x, nfft=32)
        xp_assert_close(f, fp[::2])
        xp_assert_close(p, pp[::2])
        assert pp.shape == (17,)

    def test_empty_input(self, xp):
        f, p = periodogram([])
        assert_array_equal(f.shape, (0,))
        assert_array_equal(p.shape, (0,))
        for shape in [(0,), (3,0), (0,5,2)]:
            f, p = periodogram(xp.empty(shape))
            assert_array_equal(f.shape, shape)
            assert_array_equal(p.shape, shape)

    def test_empty_input_other_axis(self, xp):
        for shape in [(3,0), (0,5,2)]:
            f, p = periodogram(xp.empty(shape), axis=1)
            assert_array_equal(f.shape, shape)
            assert_array_equal(p.shape, shape)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "no moveaxis in array_api_strict"])
    def test_short_nfft(self, xp):
        x = xp.zeros(18)
        x[0] = 1
        f, p = periodogram(x, nfft=16)
        xp_assert_close(f, xp.linspace(0, 0.5, 9))
        q = xp.ones(9)
        q[0] = 0
        q[-1] /= 2.0
        q /= 8
        xp_assert_close(p, q)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "no moveaxis in array_api_strict"])
    def test_nfft_is_xshape(self, xp):
        x = xp.zeros(16)
        x[0] = 1
        f, p = periodogram(x, nfft=16)
        xp_assert_close(f, xp.linspace(0, 0.5, 9))
        q = xp.ones(9)
        q[0] = 0
        q[-1] /= 2.0
        q /= 8
        xp_assert_close(p, q)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "no moveaxis in array_api_strict"])
    def test_real_onesided_even_32(self, xp):
        x = xp.zeros(16, dtype=xp.float32)
        x[0] = 1
        f, p = periodogram(x)
        xp_assert_close(f, xp.linspace(0, 0.5, 9))
        q = xp.ones(9, dtype=xp.float32)
        q[0] = 0
        q[-1] /= 2.0
        q /= 8
        xp_assert_close(p, q)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "no moveaxis in array_api_strict"])
    def test_real_onesided_odd_32(self, xp):
        x = xp.zeros(15, dtype=xp.float32)
        x[0] = 1
        f, p = periodogram(x)
        xp_assert_close(f, xp.arange(8.0)/15.0)
        q = xp.ones(8, dtype=xp.float32)
        q[0] = 0
        q *= 2.0/15.0
        xp_assert_close(p, q, atol=1e-7)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "no moveaxis in array_api_strict"])
    def test_real_twosided_32(self, xp):
        x = xp.zeros(16, dtype=xp.float32)
        x[0] = 1
        f, p = periodogram(x, return_onesided=False)
        xp_assert_close(f, fftfreq(16, 1.0), check_namespace=False, check_dtype=False)
        q = xp.full((16,), 1/16.0, dtype=xp.float32)
        q[0] = 0
        xp_assert_close(p, q)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "xp.mean() rqeuires real types for array_api_strict"])
    def test_complex_32(self, xp):
        x = xp.zeros(16, dtype=xp.complex64)
        x[0] = 1.0 + 2.0j
        f, p = periodogram(x, return_onesided=False)
        xp_assert_close(f, fftfreq(16, 1.0), check_dtype=False, check_namespace=False)
        q = xp.full((16,), 5.0/16.0, dtype=xp.float32)
        q[0] = 0
        xp_assert_close(p, q)

    def test_shorter_window_error(self, xp):
        x = xp.zeros(16)
        x[0] = 1
        win = signal.get_window('hann', 10)
        expected_msg = ('the size of the window must be the same size '
                        'of the input on the specified axis')
        with assert_raises(ValueError, match=expected_msg):
            periodogram(x, window=win)


class TestWelch:
    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis absent from array_api_strict"])
    def test_real_onesided_even(self, xp):
        x = xp.zeros(16)
        x[0] = 1
        x[8] = 1
        f, p = welch(x, nperseg=8)
        q = xp.asarray([0.08333333, 0.15277778, 0.22222222, 0.22222222,
                      0.11111111])
        xp_assert_close(p, q, atol=1e-7, rtol=1e-7)
        xp_assert_close(f, xp.linspace(0, 0.5, 5))

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis absent from array_api_strict"])
    def test_real_onesided_odd(self, xp):
        x = xp.zeros(16)
        x[0] = 1
        x[8] = 1
        f, p = welch(x, nperseg=9)
        xp_assert_close(f, xp.arange(5.0)/9.0)
        q = xp.asarray([0.12477455, 0.23430933, 0.17072113, 0.17072113,
                      0.17072113])
        xp_assert_close(p, q, atol=1e-7, rtol=1e-7)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis absent from array_api_strict"])
    def test_real_twosided(self, xp):
        x = xp.zeros(16, dtype=xp.float64)
        x[0] = 1
        x[8] = 1
        f, p = welch(x, nperseg=8, return_onesided=False)
        xp_assert_close(f, fftfreq(8, 1.0), check_namespace=False, check_dtype=False)
        q = xp.asarray([0.08333333, 0.07638889, 0.11111111, 0.11111111,
                      0.11111111, 0.11111111, 0.11111111, 0.07638889],
                      dtype=xp.float64)
        xp_assert_close(p, q, atol=1e-7, rtol=1e-7)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis absent from array_api_strict"])
    def test_real_spectrum(self, xp):
        x = xp.zeros(16, dtype=xp.float64)
        x[0] = 1
        x[8] = 1
        f, p = welch(x, nperseg=8, scaling='spectrum')
        xp_assert_close(f, xp.linspace(0, 0.5, 5))
        q = xp.asarray([0.015625, 0.02864583, 0.04166667, 0.04166667,
                      0.02083333], dtype=xp.float64)
        xp_assert_close(p, q, atol=1e-7, rtol=1e-7)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "dtype casting issue for array_api_strict"])
    def test_integer_onesided_even(self, xp):
        x = xp.zeros(16, dtype=xp.int64)
        x[0] = 1
        x[8] = 1
        f, p = welch(x, nperseg=8)
        xp_assert_close(f, xp.linspace(0, 0.5, 5))
        q = xp.asarray([0.08333333, 0.15277778, 0.22222222, 0.22222222,
                      0.11111111])
        xp_assert_close(p, q, atol=1e-7, rtol=1e-7)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "dtype casting issue for array_api_strict"])
    def test_integer_onesided_odd(self, xp):
        x = xp.zeros(16, dtype=xp.int64)
        x[0] = 1
        x[8] = 1
        f, p = welch(x, nperseg=9)
        xp_assert_close(f, xp.arange(5.0)/9.0)
        q = xp.asarray([0.12477455, 0.23430933, 0.17072113, 0.17072113,
                      0.17072113])
        xp_assert_close(p, q, atol=1e-7, rtol=1e-7)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "dtype casting issue for array_api_strict"])
    def test_integer_twosided(self, xp):
        x = xp.zeros(16, dtype=xp.int64)
        x[0] = 1
        x[8] = 1
        f, p = welch(x, nperseg=8, return_onesided=False)
        xp_assert_close(f, fftfreq(8, 1.0), check_namespace=False, check_dtype=False)
        q = xp.asarray([0.08333333, 0.07638889, 0.11111111, 0.11111111,
                      0.11111111, 0.11111111, 0.11111111, 0.07638889])
        xp_assert_close(p, q, atol=1e-7, rtol=1e-7)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "xp.mean() requires real types for array_api_strict"])
    def test_complex(self, xp):
        x = xp.zeros(16, dtype=xp.complex128)
        x[0] = 1.0 + 2.0j
        x[8] = 1.0 + 2.0j
        f, p = welch(x, nperseg=8, return_onesided=False)
        xp_assert_close(f, fftfreq(8, 1.0), check_namespace=False, check_dtype=False)
        q = xp.asarray([0.41666667, 0.38194444, 0.55555556, 0.55555556,
                      0.55555556, 0.55555556, 0.55555556, 0.38194444], dtype=xp.float64)
        xp_assert_close(p, q, atol=1e-7, rtol=1e-7)

    def test_unk_scaling(self, xp):
        assert_raises(ValueError, welch, xp.zeros(4, dtype=xp.complex128),
                      scaling='foo', nperseg=4)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis not available for array_api_strict"])
    def test_detrend_linear(self, xp):
        x = xp.arange(10, dtype=xp.float64) + 0.04
        f, p = welch(x, nperseg=10, detrend='linear')
        xp_assert_close(p, xp.zeros_like(p), atol=1e-15)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis not available for array_api_strict"])
    def test_no_detrending(self, xp):
        x = xp.arange(10, dtype=xp.float64) + 0.04
        f1, p1 = welch(x, nperseg=10, detrend=False)
        f2, p2 = welch(x, nperseg=10, detrend=lambda x: x)
        xp_assert_close(f1, f2, atol=1e-15)
        xp_assert_close(p1, p2, atol=1e-15)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis not available for array_api_strict"])
    def test_detrend_external(self, xp):
        x = xp.arange(10, dtype=xp.float64) + 0.04
        f, p = welch(x, nperseg=10,
                     detrend=lambda seg: signal.detrend(seg, type='l'))
        xp_assert_close(p, xp.zeros_like(p), atol=1e-15)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis not available for array_api_strict"])
    def test_detrend_external_nd_m1(self, xp):
        x = xp.arange(40, dtype=xp.float64) + 0.04
        x = xp.reshape(x, (2,2,10))
        f, p = welch(x, nperseg=10,
                     detrend=lambda seg: signal.detrend(seg, type='l'))
        xp_assert_close(p, xp.zeros_like(p), atol=1e-15)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis not available for array_api_strict"])
    def test_detrend_external_nd_0(self, xp):
        x = xp.arange(20, dtype=xp.float64) + 0.04
        x = xp.reshape(x, (2, 1, 10))
        x = xp.moveaxis(x, 2, 0)
        f, p = welch(x, nperseg=10, axis=0,
                     detrend=lambda seg: signal.detrend(seg, axis=0, type='l'))
        xp_assert_close(p, xp.zeros_like(p), atol=1e-15)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis not available for array_api_strict"])
    def test_nd_axis_m1(self, xp):
        x = xp.arange(20, dtype=xp.float64) + 0.04
        x = xp.reshape(x, (2,1,10))
        f, p = welch(x, nperseg=10)
        assert_array_equal(p.shape, (2, 1, 6))
        xp_assert_close(p[0,0,:], p[1,0,:], atol=1e-13, rtol=1e-13)
        f0, p0 = welch(x[0,0,:], nperseg=10)
        xp_assert_close(p0[None,:], p[1,:], atol=1e-13, rtol=1e-13)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis not available for array_api_strict"])
    def test_nd_axis_0(self, xp):
        x = xp.arange(20, dtype=xp.float64) + 0.04
        x = xp.reshape(x, (10,2,1))
        f, p = welch(x, nperseg=10, axis=0)
        assert_array_equal(p.shape, (6,2,1))
        xp_assert_close(p[:,0,0], p[:,1,0], atol=1e-13, rtol=1e-13)
        f0, p0 = welch(x[:,0,0], nperseg=10)
        xp_assert_close(p0, p[:,1,0], atol=1e-13, rtol=1e-13)

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis not available for array_api_strict",
                               "TODO: skip torch"])
    def test_window_external(self, xp):
        x = xp.zeros(16)
        x[0] = 1
        x[8] = 1
        f, p = welch(x, 10, 'hann', nperseg=8)
        win = signal.get_window('hann', 8)
        fe, pe = welch(x, 10, win, nperseg=None)
        assert_array_almost_equal_nulp(p, pe)
        assert_array_almost_equal_nulp(f, fe)
        assert_array_equal(fe.shape, (5,))  # because win length used as nperseg
        assert_array_equal(pe.shape, (5,))
        assert_raises(ValueError, welch, x,
                      10, win, nperseg=4)  # because nperseg != win.shape[-1]
        win_err = signal.get_window('hann', 32)
        assert_raises(ValueError, welch, x,
                      10, win_err, nperseg=None)  # win longer than signal

    def test_empty_input(self, xp):
        val = xp.asarray([])
        f, p = welch(val)
        assert_array_equal(f.shape, (0,))
        assert_array_equal(p.shape, (0,))
        for shape in [(0,), (3,0), (0,5,2)]:
            f, p = welch(np.empty(shape))
            assert_array_equal(f.shape, shape)
            assert_array_equal(p.shape, shape)

    def test_empty_input_other_axis(self, xp):
        for shape in [(3,0), (0,5,2)]:
            f, p = welch(xp.empty(shape), axis=1)
            assert_array_equal(f.shape, shape)
            assert_array_equal(p.shape, shape)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis not available for array_api_strict"])
    def test_short_data(self, xp):
        x = xp.zeros(8)
        x[0] = 1
        #for string-like window, input signal length < nperseg value gives
        #UserWarning, sets nperseg to x.shape[-1]
        with suppress_warnings() as sup:
            msg = "nperseg = 256 is greater than input length  = 8, using nperseg = 8"
            sup.filter(UserWarning, msg)
            f, p = welch(x,window='hann')  # default nperseg
            f1, p1 = welch(x,window='hann', nperseg=256)  # user-specified nperseg
        f2, p2 = welch(x, nperseg=8)  # valid nperseg, doesn't give warning
        xp_assert_close(f, f2)
        xp_assert_close(p, p2)
        xp_assert_close(f1, f2)
        xp_assert_close(p1, p2)

    def test_window_long_or_nd(self, xp):
        assert_raises(ValueError, welch, xp.zeros(4), 1, xp.asarray([1,1,1,1,1]))
        assert_raises(ValueError, welch, xp.zeros(4), 1,
                      xp.reshape(xp.arange(6), (2,3)))

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis not available for array_api_strict"])
    def test_nondefault_noverlap(self, xp):
        x = xp.zeros(64)
        x[::8] = 1
        f, p = welch(x, nperseg=16, noverlap=4)
        q = xp.asarray([0, 1./12., 1./3., 1./5., 1./3., 1./5., 1./3., 1./5.,
                      1./6.])
        xp_assert_close(p, q, atol=1e-12)

    def test_bad_noverlap(self, xp):
        assert_raises(ValueError, welch, xp.zeros(4), 1, 'hann', 2, 7)

    def test_nfft_too_short(self, xp):
        assert_raises(ValueError, welch, xp.ones(12), nfft=3, nperseg=4)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis not available for array_api_strict"])
    def test_real_onesided_even_32(self, xp):
        x = xp.zeros(16, dtype=xp.float32)
        x[0] = 1
        x[8] = 1
        f, p = welch(x, nperseg=8)
        xp_assert_close(f, xp.linspace(0, 0.5, 5))
        q = xp.asarray([0.08333333, 0.15277778, 0.22222222, 0.22222222,
                        0.11111111], dtype=xp.float32)
        xp_assert_close(p, q, atol=1e-7, rtol=1e-7)
        assert p.dtype == q.dtype

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis not available for array_api_strict"])
    def test_real_onesided_odd_32(self, xp):
        x = xp.zeros(16, dtype=xp.float32)
        x[0] = 1
        x[8] = 1
        f, p = welch(x, nperseg=9)
        xp_assert_close(f, xp.arange(5.0)/9.0)
        q = xp.asarray([0.12477458, 0.23430935, 0.17072113, 0.17072116,
                        0.17072113], dtype=xp.float32)
        xp_assert_close(p, q, atol=1e-7, rtol=1e-7)
        assert p.dtype == q.dtype

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis not available for array_api_strict"])
    def test_real_twosided_32(self, xp):
        x = xp.zeros(16, dtype=xp.float32)
        x[0] = 1
        x[8] = 1
        f, p = welch(x, nperseg=8, return_onesided=False)
        xp_assert_close(f, fftfreq(8, 1.0), check_namespace=False, check_dtype=False)
        q = xp.asarray([0.08333333, 0.07638889, 0.11111111,
                        0.11111111, 0.11111111, 0.11111111, 0.11111111,
                        0.07638889], dtype=xp.float32)
        xp_assert_close(p, q, atol=1e-7, rtol=1e-7)
        assert p.dtype == q.dtype

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "xp.mean() requires real types for array_api_strict"])
    def test_complex_32(self, xp):
        x = xp.zeros(16, dtype=xp.complex64)
        x[0] = 1.0 + 2.0j
        x[8] = 1.0 + 2.0j
        f, p = welch(x, nperseg=8, return_onesided=False)
        xp_assert_close(f, fftfreq(8, 1.0), check_namespace=False, check_dtype=False)
        q = xp.asarray([0.41666666, 0.38194442, 0.55555552, 0.55555552,
                      0.55555558, 0.55555552, 0.55555552, 0.38194442], dtype=xp.float32)
        xp_assert_close(p, q, atol=1e-7, rtol=1e-7)
        assert_(p.dtype == q.dtype,
                f'dtype mismatch, {p.dtype}, {q.dtype}')

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis absent from array_api_strict"])
    def test_padded_freqs(self, xp):
        x = xp.zeros(12)

        nfft = 24
        f = fftfreq(nfft, 1.0)[:nfft//2+1]
        f[-1] *= -1
        fodd, _ = welch(x, nperseg=5, nfft=nfft)
        feven, _ = welch(x, nperseg=6, nfft=nfft)
        xp_assert_close(fodd, f, check_namespace=False, check_dtype=False)
        xp_assert_close(feven, f, check_namespace=False, check_dtype=False)

        nfft = 25
        f = fftfreq(nfft, 1.0)[:(nfft + 1)//2]
        fodd, _ = welch(x, nperseg=5, nfft=nfft)
        feven, _ = welch(x, nperseg=6, nfft=nfft)
        xp_assert_close(fodd, f, check_namespace=False, check_dtype=False)
        xp_assert_close(feven, f, check_namespace=False, check_dtype=False)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis absent from array_api_strict"],
                      cpu_only=True)
    def test_window_correction(self, xp):
        A = 20
        fs = 1e4
        nperseg = int(fs//10)
        fsig = 300
        ii = int(fsig*nperseg//fs)  # Freq index of fsig

        tt = xp.arange(fs)/fs
        x = A*xp.sin(2*xp.pi*fsig*tt)

        for window in ['hann', 'bartlett', ('tukey', 0.1), 'flattop']:
            _, p_spec = welch(x, fs=fs, nperseg=nperseg, window=window,
                              scaling='spectrum')
            freq, p_dens = welch(x, fs=fs, nperseg=nperseg, window=window,
                                 scaling='density')

            # Check peak height at signal frequency for 'spectrum'
            xp_assert_close(p_spec[ii], A**2/2.0, rtol=5e-7, check_namespace=False)
            # Check integrated spectrum RMS for 'density'
            if np.lib.NumpyVersion(np.__version__) >= "2.0.0rc1":
                trapezoid = np.trapezoid
            else:
                trapezoid = np.trapz
            assert_allclose(np.sqrt(trapezoid(p_dens, freq)),
                            A*np.sqrt(2)/2, rtol=1e-3)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis absent from array_api_strict"])
    def test_axis_rolling(self, xp):
        np.random.seed(1234)

        x_flat = xp.asarray(np.random.randn(1024))
        _, p_flat = welch(x_flat)

        for a in range(3):
            newshape = [1,]*3
            newshape[a] = -1
            x = x_flat.reshape(newshape)

            _, p_plus = welch(x, axis=a)  # Positive axis index
            _, p_minus = welch(x, axis=a-x.ndim)  # Negative axis index

            xp_assert_equal(p_flat, p_plus.squeeze(), err_msg=a)
            xp_assert_equal(p_flat, p_minus.squeeze(), err_msg=a-x.ndim)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis absent from array_api_strict"])
    def test_average(self, xp):
        x = xp.zeros(16)
        x[0] = 1
        x[8] = 1
        f, p = welch(x, nperseg=8, average='median')
        q = xp.asarray([.1, .05, 0., 1.54074396e-33, 0.])
        xp_assert_close(f, xp.linspace(0, 0.5, 5))
        xp_assert_close(p, q, atol=1e-7, rtol=1e-7)

        assert_raises(ValueError, welch, x, nperseg=8,
                      average='unrecognised-average')


class TestCSD:
    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis absent from array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_pad_shorter_x(self, xp):
        x = xp.zeros(8)
        y = xp.zeros(12)

        f = xp.linspace(0, 0.5, 7)
        c = xp.zeros(7, dtype=xp.complex128)
        f1, c1 = csd(x, y, nperseg=12)

        xp_assert_close(f1, f)
        xp_assert_close(c1, c)

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis absent from array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_pad_shorter_y(self, xp):
        x = xp.zeros(12)
        y = xp.zeros(8)

        f = xp.linspace(0, 0.5, 7)
        c = xp.zeros(7, dtype=xp.complex128)
        f1, c1 = csd(x, y, nperseg=12)

        xp_assert_close(f1, f)
        xp_assert_close(c1, c)

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis absent from array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_real_onesided_even(self, xp):
        x = xp.zeros(16)
        x[0] = 1
        x[8] = 1
        f, p = csd(x, x, nperseg=8)
        xp_assert_close(f, xp.linspace(0, 0.5, 5))
        q = xp.asarray([0.08333333, 0.15277778, 0.22222222, 0.22222222,
                        0.11111111])
        xp_assert_close(p, q, atol=1e-7, rtol=1e-7)

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis absent from array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_real_onesided_odd(self, xp):
        x = xp.zeros(16)
        x[0] = 1
        x[8] = 1
        f, p = csd(x, x, nperseg=9)
        xp_assert_close(f, np.arange(5.0)/9.0)
        q = xp.asarray([0.12477455, 0.23430933, 0.17072113, 0.17072113,
                        0.17072113])
        xp_assert_close(p, q, atol=1e-7, rtol=1e-7)

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis absent from array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_real_twosided(self, xp):
        x = xp.zeros(16)
        x[0] = 1
        x[8] = 1
        f, p = csd(x, x, nperseg=8, return_onesided=False)
        xp_assert_close(f, fftfreq(8, 1.0))
        q = xp.asarray([0.08333333, 0.07638889, 0.11111111, 0.11111111,
                        0.11111111, 0.11111111, 0.11111111, 0.07638889])
        xp_assert_close(p, q, atol=1e-7, rtol=1e-7)

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis absent from array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_real_spectrum(self, xp):
        x = xp.zeros(16)
        x[0] = 1
        x[8] = 1
        f, p = csd(x, x, nperseg=8, scaling='spectrum')
        xp_assert_close(f, xp.linspace(0, 0.5, 5))
        q = xp.asarray([0.015625, 0.02864583, 0.04166667, 0.04166667,
                        0.02083333])
        xp_assert_close(p, q, atol=1e-7, rtol=1e-7)

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "type casting error with array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_integer_onesided_even(self, xp):
        x = xp.zeros(16, dtype=xp.int64)
        x[0] = 1
        x[8] = 1
        f, p = csd(x, x, nperseg=8)
        xp_assert_close(f, xp.linspace(0, 0.5, 5))
        q = xp.asarray([0.08333333, 0.15277778, 0.22222222, 0.22222222,
                        0.11111111])
        xp_assert_close(p, q, atol=1e-7, rtol=1e-7)

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "type casting error with array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_integer_onesided_odd(self, xp):
        x = xp.zeros(16, dtype=xp.int64)
        x[0] = 1
        x[8] = 1
        f, p = csd(x, x, nperseg=9)
        xp_assert_close(f, xp.arange(5.0)/9.0)
        q = xp.asarray([0.12477455, 0.23430933, 0.17072113, 0.17072113,
                        0.17072113])
        xp_assert_close(p, q, atol=1e-7, rtol=1e-7)

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "type casting error with array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_integer_twosided(self, xp):
        x = xp.zeros(16, dtype=xp.int64)
        x[0] = 1
        x[8] = 1
        f, p = csd(x, x, nperseg=8, return_onesided=False)
        xp_assert_close(f, fftfreq(8, 1.0))
        q = xp.asarray([0.08333333, 0.07638889, 0.11111111, 0.11111111,
                        0.11111111, 0.11111111, 0.11111111, 0.07638889])
        xp_assert_close(p, q, atol=1e-7, rtol=1e-7)

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "xp.mean() requires real input for array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_complex(self, xp):
        x = xp.zeros(16, dtype=xp.complex128)
        x[0] = 1.0 + 2.0j
        x[8] = 1.0 + 2.0j
        f, p = csd(x, x, nperseg=8, return_onesided=False)
        xp_assert_close(f, fftfreq(8, 1.0))
        q = xp.asarray([0.41666667, 0.38194444, 0.55555556, 0.55555556,
                        0.55555556, 0.55555556, 0.55555556, 0.38194444])
        xp_assert_close(p, q, atol=1e-7, rtol=1e-7)

    def test_unk_scaling(self, xp):
        assert_raises(ValueError, csd, xp.zeros(4, dtype=xp.complex128),
                      xp.ones(4, dtype=xp.complex128), scaling='foo', nperseg=4)

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis not available in array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_detrend_linear(self, xp):
        x = xp.arange(10, dtype=xp.float64) + 0.04
        f, p = csd(x, x, nperseg=10, detrend='linear')
        xp_assert_close(p, xp.zeros_like(p), atol=1e-15)

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis not available in array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_no_detrending(self, xp):
        x = xp.arange(10, dtype=xp.float64) + 0.04
        f1, p1 = csd(x, x, nperseg=10, detrend=False)
        f2, p2 = csd(x, x, nperseg=10, detrend=lambda x: x)
        xp_assert_close(f1, f2, atol=1e-15)
        xp_assert_close(p1, p2, atol=1e-15)

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis not available in array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_detrend_external(self, xp):
        x = xp.arange(10, dtype=xp.float64) + 0.04
        f, p = csd(x, x, nperseg=10,
                   detrend=lambda seg: signal.detrend(seg, type='l'))
        xp_assert_close(p, xp.zeros_like(p), atol=1e-15)

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis not available in array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_detrend_external_nd_m1(self, xp):
        x = xp.arange(40, dtype=xp.float64) + 0.04
        x = xp.reshape(x, (2, 2, 10))
        f, p = csd(x, x, nperseg=10,
                   detrend=lambda seg: signal.detrend(seg, type='l'))
        xp_assert_close(p, xp.zeros_like(p), atol=1e-15)

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis not available in array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_detrend_external_nd_0(self, xp):
        x = xp.arange(20, dtype=xp.float64) + 0.04
        x = xp.reshape(x, (2, 1, 10))
        x = xp.moveaxis(x, 2, 0)
        f, p = csd(x, x, nperseg=10, axis=0,
                   detrend=lambda seg: signal.detrend(seg, axis=0, type='l'))
        xp_assert_close(p, xp.zeros_like(p), atol=1e-15)

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis not available in array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_nd_axis_m1(self, xp):
        x = xp.arange(20, dtype=xp.float64) + 0.04
        x = xp.reshape(x, (2, 1, 10))
        f, p = csd(x, x, nperseg=10)
        assert p.shape == (2, 1, 6)
        xp_assert_close(p[0,0,:], p[1,0,:], atol=1e-13, rtol=1e-13)
        f0, p0 = csd(x[0,0,:], x[0,0,:], nperseg=10)
        xp_assert_close(p0[xp.newaxis,:],
                        p[1,:],
                        atol=1e-13,
                        rtol=1e-13,
                        check_dtype=False)

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis not available in array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_nd_axis_0(self, xp):
        x = xp.arange(20, dtype=xp.float64) + 0.04
        x = xp.reshape(x, (10, 2, 1))
        f, p = csd(x, x, nperseg=10, axis=0)
        assert p.shape == (6, 2, 1)
        xp_assert_close(p[:,0,0], p[:,1,0], atol=1e-13, rtol=1e-13)
        f0, p0 = csd(x[:,0,0], x[:,0,0], nperseg=10)
        xp_assert_close(p0, p[:,1,0], atol=1e-13, rtol=1e-13, check_dtype=False)

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["array_api_compat cupy doesn't support fft",
                               "moveaxis not available in array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_window_external(self, xp):
        x = xp.zeros(16)
        x[0] = 1
        x[8] = 1
        f, p = csd(x, x, 10, 'hann', 8)
        win = signal.get_window('hann', 8)
        fe, pe = csd(x, x, 10, win, nperseg=None)
        # TODO: no nulp array API testing funcs yet
        assert_array_almost_equal_nulp(p, pe)
        assert_array_almost_equal_nulp(f, fe)
        assert fe.shape == (5,)  # because win length used as nperseg
        assert pe.shape == (5,)
        assert_raises(ValueError, csd, x, x,
                      10, win, nperseg=256)  # because nperseg != win.shape[-1]
        win_err = signal.get_window('hann', 32)
        assert_raises(ValueError, csd, x, x,
              10, win_err, nperseg=None)  # because win longer than signal

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["array-like support issues for CuPy",
                               "moveaxis not available in array_api_strict",
                               "torch max() expects reduction dim to be specified"])
    def test_empty_input(self, xp):
        f, p = csd([], xp.zeros(10))
        assert f.shape == (0,)
        assert p.shape == (0,)

        f, p = csd(xp.zeros(10), [])
        assert f.shape == (0,)
        assert p.shape == (0,)

        for shape in [(0,), (3,0), (0,5,2)]:
            f, p = csd(xp.empty(shape), xp.empty(shape))
            assert f.shape == shape
            assert p.shape == shape

        f, p = csd(xp.ones(10), xp.empty((5,0)))
        assert f.shape == (5,0)
        assert p.shape == (5,0)

        f, p = csd(xp.empty((5,0)), xp.ones(10))
        assert f.shape == (5,0)
        assert p.shape == (5,0)

    @skip_xp_backends("array_api_strict", "torch",
                      reasons=["moveaxis not available in array_api_strict",
                               "torch max() expects reduction dim to be specified"])
    def test_empty_input_other_axis(self, xp):
        for shape in [(3,0), (0,5,2)]:
            f, p = csd(xp.empty(shape), xp.empty(shape), axis=1)
            assert f.shape == shape
            assert p.shape == shape

        f, p = csd(xp.empty((10,10,3)), xp.zeros((10,0,1)), axis=1)
        assert f.shape == (10, 0, 3)
        assert p.shape == (10, 0, 3)

        f, p = csd(xp.empty((10,0,1)), xp.zeros((10,10,3)), axis=1)
        assert f.shape == (10, 0, 3)
        assert p.shape == (10, 0, 3)

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["lack of fft support in array_api_compat cupy",
                               "moveaxis not available in array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_short_data(self, xp):
        x = xp.zeros(8)
        x[0] = 1

        #for string-like window, input signal length < nperseg value gives
        #UserWarning, sets nperseg to x.shape[-1]
        with suppress_warnings() as sup:
            msg = "nperseg = 256 is greater than input length  = 8, using nperseg = 8"
            sup.filter(UserWarning, msg)
            f, p = csd(x, x, window='hann')  # default nperseg
            f1, p1 = csd(x, x, window='hann', nperseg=256)  # user-specified nperseg
        f2, p2 = csd(x, x, nperseg=8)  # valid nperseg, doesn't give warning
        xp_assert_close(f, f2)
        xp_assert_close(p, p2)
        xp_assert_close(f1, f2)
        xp_assert_close(p1, p2)

    def test_window_long_or_nd(self, xp):
        assert_raises(ValueError, csd, xp.zeros(4), xp.ones(4), 1,
                      xp.asarray([1,1,1,1,1]))
        assert_raises(ValueError, csd, xp.zeros(4), xp.ones(4), 1,
                      xp.reshape(xp.arange(6), (2, 3)))

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["lack of fft support in array_api_compat cupy",
                               "moveaxis not available in array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_nondefault_noverlap(self, xp):
        x = xp.zeros(64)
        x[::8] = 1
        f, p = csd(x, x, nperseg=16, noverlap=4)
        q = xp.asarray([0, 1./12., 1./3., 1./5., 1./3., 1./5., 1./3., 1./5.,
                        1./6.])
        xp_assert_close(p, q, atol=1e-12)

    def test_bad_noverlap(self, xp):
        assert_raises(ValueError, csd, xp.zeros(4), xp.ones(4), 1, 'hann',
                      2, 7)

    def test_nfft_too_short(self, xp):
        assert_raises(ValueError, csd, xp.ones(12), xp.zeros(12), nfft=3,
                      nperseg=4)

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["lack of fft support in array_api_compat cupy",
                               "moveaxis not available in array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_real_onesided_even_32(self, xp):
        x = xp.zeros(16, dtype=xp.float32)
        x[0] = 1
        x[8] = 1
        f, p = csd(x, x, nperseg=8)
        xp_assert_close(f, xp.linspace(0, 0.5, 5))
        q = xp.asarray([0.08333333, 0.15277778, 0.22222222, 0.22222222,
                        0.11111111], xp.float32)
        xp_assert_close(p, q, atol=1e-7, rtol=1e-7)

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["lack of fft support in array_api_compat cupy",
                               "moveaxis not available in array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_real_onesided_odd_32(self, xp):
        x = xp.zeros(16, dtype=xp.float32)
        x[0] = 1
        x[8] = 1
        f, p = csd(x, x, nperseg=9)
        xp_assert_close(f, xp.arange(5.0)/9.0)
        q = xp.asarray([0.12477458, 0.23430935, 0.17072113, 0.17072116,
                        0.17072113], xp.float32)
        xp_assert_close(p, q, atol=1e-7, rtol=1e-7)

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["lack of fft support in array_api_compat cupy",
                               "moveaxis not available in array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_real_twosided_32(self, xp):
        x = xp.zeros(16, dtype=xp.float32)
        x[0] = 1
        x[8] = 1
        f, p = csd(x, x, nperseg=8, return_onesided=False)
        xp_assert_close(f, fftfreq(8, 1.0))
        q = xp.asarray([0.08333333, 0.07638889, 0.11111111,
                        0.11111111, 0.11111111, 0.11111111, 0.11111111,
                        0.07638889], xp.float32)
        xp_assert_close(p, q, atol=1e-7, rtol=1e-7)

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["lack of fft support in array_api_compat cupy",
                               "xp.mean requires real input for array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_complex_32(self, xp):
        x = xp.zeros(16, dtype=xp.complex64)
        x[0] = 1.0 + 2.0j
        x[8] = 1.0 + 2.0j
        f, p = csd(x, x, nperseg=8, return_onesided=False)
        xp_assert_close(f, fftfreq(8, 1.0))
        q = xp.asarray([0.41666666, 0.38194442, 0.55555552, 0.55555552,
                        0.55555558, 0.55555552, 0.55555552, 0.38194442], xp.float32)
        xp_assert_close(p, q, atol=1e-7, rtol=1e-7)

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["lack of fft support in array_api_compat cupy",
                               "moveaxis not availabe in array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_padded_freqs(self, xp):
        x = xp.zeros(12)
        y = xp.ones(12)

        nfft = 24
        f = fftfreq(nfft, 1.0)[:nfft//2+1]
        f[-1] *= -1
        fodd, _ = csd(x, y, nperseg=5, nfft=nfft)
        feven, _ = csd(x, y, nperseg=6, nfft=nfft)
        xp_assert_close(fodd, f)
        xp_assert_close(feven, f)

        nfft = 25
        f = fftfreq(nfft, 1.0)[:(nfft + 1)//2]
        fodd, _ = csd(x, y, nperseg=5, nfft=nfft)
        feven, _ = csd(x, y, nperseg=6, nfft=nfft)
        xp_assert_close(fodd, f)
        xp_assert_close(feven, f)

    @skip_xp_backends("cupy", "array_api_strict", "torch",
                      reasons=["lack of fft support in array_api_compat cupy",
                               "moveaxis not availabe in array_api_strict",
                               "torch hits messy np.pad codepath"])
    def test_copied_data(self, xp):
        x = np.random.randn(64)
        x = xp.asarray(x)
        y = copy(x)

        _, p_same = csd(x, x, nperseg=8, average='mean',
                        return_onesided=False)
        _, p_copied = csd(x, y, nperseg=8, average='mean',
                          return_onesided=False)
        xp_assert_close(p_same, p_copied, check_dtype=False)

        _, p_same = csd(x, x, nperseg=8, average='median',
                        return_onesided=False)
        _, p_copied = csd(x, y, nperseg=8, average='median',
                          return_onesided=False)
        xp_assert_close(p_same, p_copied, check_dtype=False)


class TestCoherence:
    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["lack of fft support in array_api_compat cupy",
                               "moveaxis not availabe in array_api_strict"])
    def test_identical_input(self, xp):
        x = np.random.randn(20)
        x = xp.asarray(x)
        y = copy(x)  # So `y is x` -> False

        f = xp.linspace(0, 0.5, 6)
        C = xp.ones(6)
        f1, C1 = coherence(x, y, nperseg=10)

        xp_assert_close(f, f1)
        xp_assert_close(C, C1, check_dtype=False)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["lack of fft support in array_api_compat cupy",
                               "moveaxis not availabe in array_api_strict"])
    def test_phase_shifted_input(self, xp):
        x = np.random.randn(20)
        x = xp.asarray(x)
        y = -x

        f = xp.linspace(0, 0.5, 6)
        C = xp.ones(6)
        f1, C1 = coherence(x, y, nperseg=10)

        xp_assert_close(f, f1)
        xp_assert_close(C, C1, check_dtype=False)


class TestSpectrogram:
    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["lack of fft support in array_api_compat cupy",
                               "moveaxis not availabe in array_api_strict"])
    def test_average_all_segments(self, xp):
        x = np.random.randn(1024)
        x = xp.asarray(x)

        fs = 1.0
        window = ('tukey', 0.25)
        nperseg = 16
        noverlap = 2

        f, _, P = spectrogram(x, fs, window, nperseg, noverlap)
        fw, Pw = welch(x, fs, window, nperseg, noverlap)
        xp_assert_close(f, fw)
        xp_assert_close(xp.mean(P, axis=-1), Pw)

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["lack of fft support in array_api_compat cupy",
                               "moveaxis not availabe in array_api_strict"])
    def test_window_external(self, xp):
        x = np.random.randn(1024)
        x = xp.asarray(x)

        fs = 1.0
        window = ('tukey', 0.25)
        nperseg = 16
        noverlap = 2
        f, _, P = spectrogram(x, fs, window, nperseg, noverlap)

        win = signal.get_window(('tukey', 0.25), 16)
        fe, _, Pe = spectrogram(x, fs, win, nperseg=None, noverlap=2)
        assert_array_equal(fe.shape, (9,))  # because win length used as nperseg
        assert_array_equal(Pe.shape, (9,73))
        assert_raises(ValueError, spectrogram, x,
                      fs, win, nperseg=8)  # because nperseg != win.shape[-1]
        win_err = signal.get_window(('tukey', 0.25), 2048)
        assert_raises(ValueError, spectrogram, x,
                      fs, win_err, nperseg=None)  # win longer than signal

    @skip_xp_backends("cupy", "array_api_strict",
                      reasons=["lack of fft support in array_api_compat cupy",
                               "moveaxis not availabe in array_api_strict"])
    def test_short_data(self, xp):
        x = np.random.randn(1024)
        x = xp.asarray(x)
        fs = 1.0

        #for string-like window, input signal length < nperseg value gives
        #UserWarning, sets nperseg to x.shape[-1]
        f, _, p = spectrogram(x, fs, window=('tukey',0.25))  # default nperseg
        with suppress_warnings() as sup:
            sup.filter(UserWarning,
                       "nperseg = 1025 is greater than input length  = 1024, "
                       "using nperseg = 1024",)
            f1, _, p1 = spectrogram(x, fs, window=('tukey',0.25),
                                    nperseg=1025)  # user-specified nperseg
        f2, _, p2 = spectrogram(x, fs, nperseg=256)  # to compare w/default
        f3, _, p3 = spectrogram(x, fs, nperseg=1024)  # compare w/user-spec'd
        xp_assert_close(f, f2)
        xp_assert_close(p, p2)
        xp_assert_close(f1, f3)
        xp_assert_close(p1, p3)

class TestLombscargle:
    @skip_xp_backends(np_only=True,
                      reasons=["_lombscargle is a Cython function"])
    def test_frequency(self, xp):
        """Test if frequency location of peak corresponds to frequency of
        generated input signal.
        """

        # Input parameters
        ampl = 2.
        w = 1.
        phi = 0.5 * xp.pi
        nin = 100
        nout = 1000
        p = 0.7  # Fraction of points to select

        # Randomly select a fraction of an array with timesteps
        np.random.seed(2353425)
        r = np.random.rand(nin)
        r = xp.asarray(r)
        t = xp.linspace(0.01*xp.pi, 10.*xp.pi, nin)[r >= p]

        # Plot a sine wave for the selected times
        x = ampl * xp.sin(w*t + phi)

        # Define the array of frequencies for which to compute the periodogram
        f = xp.linspace(0.01, 10., nout)

        # Calculate Lomb-Scargle periodogram
        P = lombscargle(t, x, f)

        # Check if difference between found frequency maximum and input
        # frequency is less than accuracy
        delta = f[1] - f[0]
        assert_(w - f[xp.argmax(P)] < (delta/2.))

    @skip_xp_backends(np_only=True,
                      reasons=["_lombscargle is a Cython function"])
    def test_amplitude(self, xp):
        # Test if height of peak in normalized Lomb-Scargle periodogram
        # corresponds to amplitude of the generated input signal.

        # Input parameters
        ampl = 2.
        w = 1.
        phi = 0.5 * xp.pi
        nin = 100
        nout = 1000
        p = 0.7  # Fraction of points to select

        # Randomly select a fraction of an array with timesteps
        np.random.seed(2353425)
        r = np.random.rand(nin)
        r = xp.asarray(r)
        t = xp.linspace(0.01*xp.pi, 10.*xp.pi, nin)[r >= p]

        # Plot a sine wave for the selected times
        x = ampl * xp.sin(w*t + phi)

        # Define the array of frequencies for which to compute the periodogram
        f = xp.linspace(0.01, 10., nout)

        # Calculate Lomb-Scargle periodogram
        pgram = lombscargle(t, x, f)

        # Normalize
        pgram = xp.sqrt(4 * pgram / t.shape[0])

        # Check if difference between found frequency maximum and input
        # frequency is less than accuracy
        xp_assert_close(xp.max(pgram), ampl, rtol=0.035)

    @skip_xp_backends("cupy", "torch",
                      reasons=["_lombscargle is a Cython function",
                               "_lombscargle is a Cython function"])
    def test_precenter(self, xp):
        # Test if precenter gives the same result as manually precentering.

        # Input parameters
        ampl = 2.
        w = 1.
        phi = 0.5 * xp.pi
        nin = 100
        nout = 1000
        p = 0.7  # Fraction of points to select
        offset = 0.15  # Offset to be subtracted in pre-centering

        # Randomly select a fraction of an array with timesteps
        np.random.seed(2353425)
        r = np.random.rand(nin)
        r = xp.asarray(r)
        t = xp.linspace(0.01*xp.pi, 10.*xp.pi, nin)[r >= p]

        # Plot a sine wave for the selected times
        x = ampl * xp.sin(w*t + phi) + offset

        # Define the array of frequencies for which to compute the periodogram
        f = xp.linspace(0.01, 10., nout)

        # Calculate Lomb-Scargle periodogram
        pgram = lombscargle(t, x, f, precenter=True)
        pgram2 = lombscargle(t, x - xp.mean(x), f, precenter=False)

        # check if centering worked
        xp_assert_close(pgram, pgram2)

    @skip_xp_backends(np_only=True,
                      reasons=["_lombscargle is a Cython function"])
    def test_normalize(self, xp):
        # Test normalize option of Lomb-Scarge.

        # Input parameters
        ampl = 2.
        w = 1.
        phi = 0.5 * xp.pi
        nin = 100
        nout = 1000
        p = 0.7  # Fraction of points to select

        # Randomly select a fraction of an array with timesteps
        np.random.seed(2353425)
        r = np.random.rand(nin)
        r = xp.asarray(r)
        t = xp.linspace(0.01*xp.pi, 10.*xp.pi, nin)[r >= p]

        # Plot a sine wave for the selected times
        x = ampl * xp.sin(w*t + phi)

        # Define the array of frequencies for which to compute the periodogram
        f = xp.linspace(0.01, 10., nout)

        # Calculate Lomb-Scargle periodogram
        pgram = lombscargle(t, x, f)
        pgram2 = lombscargle(t, x, f, normalize=True)

        # check if normalization works as expected
        xp_assert_close(pgram * 2 / (x @ x), pgram2)
        xp_assert_close(xp.max(pgram2), 1.0, rtol=0.35)

    @skip_xp_backends("torch", "cupy",
                      reasons=["_lombscargle is a Cython function",
                               "_lombscargle is a Cython function"])
    def test_wrong_shape(self, xp):
        t = xp.linspace(0, 1, 1)
        x = xp.linspace(0, 1, 2)
        f = xp.linspace(0, 1, 3)
        assert_raises(ValueError, lombscargle, t, x, f)

    @skip_xp_backends("torch", "cupy",
                      reasons=["_lombscargle is a Cython function",
                               "_lombscargle is a Cython function"])
    def test_zero_division(self, xp):
        t = xp.zeros(1)
        x = xp.zeros(1)
        f = xp.zeros(1)
        assert_raises(ZeroDivisionError, lombscargle, t, x, f)

    @skip_xp_backends("torch", "cupy",
                      reasons=["endpoint usage in linspace",
                               "_lombscargle is a Cython function"])
    def test_lombscargle_atan_vs_atan2(self, xp):
        # https://github.com/scipy/scipy/issues/3787
        # This raised a ZeroDivisionError.
        t = xp.linspace(0, 10, 1000, endpoint=False)
        x = xp.sin(4*t)
        f = xp.linspace(0, 50, 500, endpoint=False) + 0.1
        lombscargle(t, x, f*2*xp.pi)


class TestSTFT:
    @skip_xp_backends("torch", "cupy", "array_api_strict",
                      reasons=["torch hits messy np.pad codepath",
                               "lack of fft support in array_api_compat",
                               "moveaxis not available in array_api_strict"])
    def test_input_validation(self, xp):

        def chk_VE(match):
            """Assert for a ValueError matching regexp `match`.

            This little wrapper allows a more concise code layout.
            """
            return pytest.raises(ValueError, match=match)

        # Checks for check_COLA():
        with chk_VE('nperseg must be a positive integer'):
            check_COLA('hann', -10, 0)
        with chk_VE('noverlap must be less than nperseg.'):
            check_COLA('hann', 10, 20)
        with chk_VE('window must be 1-D'):
            check_COLA(xp.ones((2, 2)), 10, 0)
        with chk_VE('window must have length of nperseg'):
            check_COLA(xp.ones(20), 10, 0)

        # Checks for check_NOLA():
        with chk_VE('nperseg must be a positive integer'):
            check_NOLA('hann', -10, 0)
        with chk_VE('noverlap must be less than nperseg'):
            check_NOLA('hann', 10, 20)
        with chk_VE('window must be 1-D'):
            check_NOLA(xp.ones((2, 2)), 10, 0)
        with chk_VE('window must have length of nperseg'):
            check_NOLA(xp.ones(20), 10, 0)
        with chk_VE('noverlap must be a nonnegative integer'):
            check_NOLA('hann', 64, -32)

        x = xp.zeros(1024)
        z = stft(x)[2]

        # Checks for stft():
        with chk_VE('window must be 1-D'):
            stft(x, window=xp.ones((2, 2)))
        with chk_VE('value specified for nperseg is different ' +
                    'from length of window'):
            stft(x, window=xp.ones(10), nperseg=256)
        with chk_VE('nperseg must be a positive integer'):
            stft(x, nperseg=-256)
        with chk_VE('noverlap must be less than nperseg.'):
            stft(x, nperseg=256, noverlap=1024)
        with chk_VE('nfft must be greater than or equal to nperseg.'):
            stft(x, nperseg=256, nfft=8)

        # Checks for istft():
        with chk_VE('Input stft must be at least 2d!'):
            istft(x)
        with chk_VE('window must be 1-D'):
            istft(z, window=xp.ones((2, 2)))
        with chk_VE('window must have length of 256'):
            istft(z, window=xp.ones(10), nperseg=256)
        with chk_VE('nperseg must be a positive integer'):
            istft(z, nperseg=-256)
        with chk_VE('noverlap must be less than nperseg.'):
            istft(z, nperseg=256, noverlap=1024)
        with chk_VE('nfft must be greater than or equal to nperseg.'):
            istft(z, nperseg=256, nfft=8)
        with pytest.warns(UserWarning, match="NOLA condition failed, " +
                          "STFT may not be invertible"):
            istft(z, nperseg=256, noverlap=0, window='hann')
        with chk_VE('Must specify differing time and frequency axes!'):
            istft(z, time_axis=0, freq_axis=0)

        # Checks for _spectral_helper():
        with chk_VE("Unknown value for mode foo, must be one of: " +
                    r"\{'psd', 'stft'\}"):
            _spectral_helper(x, x, mode='foo')
        with chk_VE("x and y must be equal if mode is 'stft'"):
            _spectral_helper(x[:512], x[512:], mode='stft')
        with chk_VE("Unknown boundary option 'foo', must be one of: " +
                    r"\['even', 'odd', 'constant', 'zeros', None\]"):
            _spectral_helper(x, x, boundary='foo')

        scaling = "not_valid"
        with chk_VE(fr"Parameter {scaling=} not in \['spectrum', 'psd'\]!"):
            stft(x, scaling=scaling)
        with chk_VE(fr"Parameter {scaling=} not in \['spectrum', 'psd'\]!"):
            istft(z, scaling=scaling)

    def test_check_COLA(self, xp):
        settings = [
                    ('boxcar', 10, 0),
                    ('boxcar', 10, 9),
                    ('bartlett', 51, 26),
                    ('hann', 256, 128),
                    ('hann', 256, 192),
                    ('blackman', 300, 200),
                    (('tukey', 0.5), 256, 64),
                    ('hann', 256, 255),
                    ]

        for setting in settings:
            msg = '{}, {}, {}'.format(*setting)
            # NOTE: there is no array input here with string
            # window type--I think we have to assume NumPy return
            # type with no input array type to key off of?
            xp_assert_equal(xp.asarray(True),
                            check_COLA(*setting),
                            err_msg=msg,
                            check_dtype=False,
                            check_namespace=False,
                            check_shape=False)

    def test_check_NOLA(self, xp):
        settings_pass = [
                    ('boxcar', 10, 0),
                    ('boxcar', 10, 9),
                    ('boxcar', 10, 7),
                    ('bartlett', 51, 26),
                    ('bartlett', 51, 10),
                    ('hann', 256, 128),
                    ('hann', 256, 192),
                    ('hann', 256, 37),
                    ('blackman', 300, 200),
                    ('blackman', 300, 123),
                    (('tukey', 0.5), 256, 64),
                    (('tukey', 0.5), 256, 38),
                    ('hann', 256, 255),
                    ('hann', 256, 39),
                    ]
        for setting in settings_pass:
            msg = '{}, {}, {}'.format(*setting)
            # NOTE: there is no array input here with string
            # window type--I think we have to assume NumPy return
            # type with no input array type to key off of?
            xp_assert_equal(xp.asarray(True), check_NOLA(*setting),
                            err_msg=msg,
                            check_dtype=False,
                            check_namespace=False,
                            check_shape=False)

        w_fail = xp.ones(16)
        w_fail[::2] = 0
        settings_fail = [
                    (w_fail, size(w_fail), size(w_fail) // 2),
                    ('hann', 64, 0),
        ]
        for setting in settings_fail:
            msg = '{}, {}, {}'.format(*setting)
            if isinstance(setting[0], str):
                check_namespace = False
                check_dtype = False
            else:
                # we can do a stricter check when there is an
                # actual array type pass-through to key off of
                check_namespace = True
                check_dtype = True
            xp_assert_equal(xp.asarray(False), check_NOLA(*setting),
                            err_msg=msg,
                            check_dtype=check_dtype,
                            check_namespace=check_namespace,
                            check_shape=False)

    @skip_xp_backends("torch", "cupy", "array_api_strict",
                      reasons=["torch hits messy np.pad codepath",
                               "lack of fft support in array_api_compat",
                               "moveaxis not available in array_api_strict"])
    def test_average_all_segments(self, xp):
        np.random.seed(1234)
        x = np.random.randn(1024)
        x = xp.asarray(x)

        fs = 1.0
        window = 'hann'
        nperseg = 16
        noverlap = 8

        # Compare twosided, because onesided welch doubles non-DC terms to
        # account for power at negative frequencies. stft doesn't do this,
        # because it breaks invertibility.
        f, _, Z = stft(x, fs, window, nperseg, noverlap, padded=False,
                       return_onesided=False, boundary=None)
        fw, Pw = welch(x, fs, window, nperseg, noverlap, return_onesided=False,
                       scaling='spectrum', detrend=False)

        xp_assert_close(f, fw)
        xp_assert_close(xp.mean(xp.abs(Z)**2, axis=-1), Pw)

    @skip_xp_backends("torch", "cupy", "array_api_strict",
                      reasons=["torch hits messy np.pad codepath",
                               "lack of fft support in array_api_compat",
                               "moveaxis not available in array_api_strict"])
    def test_permute_axes(self, xp):
        np.random.seed(1234)
        x = np.random.randn(1024)
        x = xp.asarray(x)

        fs = 1.0
        window = 'hann'
        nperseg = 16
        noverlap = 8

        f1, t1, Z1 = stft(x, fs, window, nperseg, noverlap)
        f2, t2, Z2 = stft(x.reshape((-1, 1, 1)), fs, window, nperseg, noverlap,
                          axis=0)

        t3, x1 = istft(Z1, fs, window, nperseg, noverlap)
        t4, x2 = istft(Z2.T, fs, window, nperseg, noverlap, time_axis=0,
                       freq_axis=-1)

        xp_assert_close(f1, f2)
        xp_assert_close(t1, t2)
        xp_assert_close(t3, t4)
        xp_assert_close(Z1, Z2[:, 0, 0, :])
        xp_assert_close(x1, x2[:, 0, 0])

    @skip_xp_backends("torch", "cupy", "array_api_strict",
                      reasons=["torch hits messy np.pad codepath",
                               "lack of fft support in array_api_compat",
                               "moveaxis not available in array_api_strict"])
    @pytest.mark.parametrize('scaling', ['spectrum', 'psd'])
    def test_roundtrip_real(self, scaling, xp):
        np.random.seed(1234)

        settings = [
                    ('boxcar', 100, 10, 0),           # Test no overlap
                    ('boxcar', 100, 10, 9),           # Test high overlap
                    ('bartlett', 101, 51, 26),        # Test odd nperseg
                    ('hann', 1024, 256, 128),         # Test defaults
                    (('tukey', 0.5), 1152, 256, 64),  # Test Tukey
                    ('hann', 1024, 256, 255),         # Test overlapped hann
                    ]

        for window, N, nperseg, noverlap in settings:
            t = xp.arange(N)
            x = 10*np.random.randn(size(t))
            x = xp.asarray(x)

            _, _, zz = stft(x, nperseg=nperseg, noverlap=noverlap,
                            window=window, detrend=None, padded=False,
                            scaling=scaling)

            tr, xr = istft(zz, nperseg=nperseg, noverlap=noverlap,
                           window=window, scaling=scaling)

            msg = f'{window}, {noverlap}'
            # NOTE: when the skips above can be removed, it seems
            # unlikely we'll be able to enforce namespace matches here
            # because there are no array type inputs (string window type)
            xp_assert_close(t, tr, err_msg=msg, check_dtype=False)
            xp_assert_close(x, xr, err_msg=msg)

    @skip_xp_backends("torch", "cupy", "array_api_strict",
            reasons=["torch error with x.imag called but x is real",
                     "lack of fft support in array_api_compat",
                     "moveaxis not available in array_api_strict"])
    def test_roundtrip_not_nola(self, xp):
        np.random.seed(1234)

        w_fail = xp.ones(16)
        w_fail[::2] = 0
        settings = [
                    (w_fail, 256, size(w_fail), size(w_fail) // 2),
                    ('hann', 256, 64, 0),
        ]

        for window, N, nperseg, noverlap in settings:
            msg = f'{window}, {N}, {nperseg}, {noverlap}'
            assert not check_NOLA(window, nperseg, noverlap), msg

            t = xp.arange(N)
            x = 10 * np.random.randn(size(t))
            x = xp.asarray(x)

            _, _, zz = stft(x, nperseg=nperseg, noverlap=noverlap,
                            window=window, detrend=None, padded=True,
                            boundary='zeros')
            with pytest.warns(UserWarning, match='NOLA'):
                tr, xr = istft(zz, nperseg=nperseg, noverlap=noverlap,
                               window=window, boundary=True)

            xp_assert_close(t, tr[:len(t)], err_msg=msg, check_dtype=False)
            with pytest.raises(AssertionError):
                xp_assert_close(x, xr[:len(x)], err_msg=msg)

    @skip_xp_backends("torch", "cupy", "array_api_strict",
            reasons=["torch hits messy np.pad codepath",
                     "lack of fft support in array_api_compat",
                     "moveaxis not available in array_api_strict"])
    def test_roundtrip_nola_not_cola(self, xp):
        np.random.seed(1234)

        settings = [
                    ('boxcar', 100, 10, 3),           # NOLA True, COLA False
                    ('bartlett', 101, 51, 37),        # NOLA True, COLA False
                    ('hann', 1024, 256, 127),         # NOLA True, COLA False
                    (('tukey', 0.5), 1152, 256, 14),  # NOLA True, COLA False
                    ('hann', 1024, 256, 5),           # NOLA True, COLA False
                    ]

        for window, N, nperseg, noverlap in settings:
            msg = f'{window}, {nperseg}, {noverlap}'
            assert check_NOLA(window, nperseg, noverlap), msg
            assert not check_COLA(window, nperseg, noverlap), msg

            t = xp.arange(N)
            x = 10 * np.random.randn(size(t))
            x = xp.asarray(x)

            _, _, zz = stft(x, nperseg=nperseg, noverlap=noverlap,
                            window=window, detrend=None, padded=True,
                            boundary='zeros')

            tr, xr = istft(zz, nperseg=nperseg, noverlap=noverlap,
                           window=window, boundary=True)

            msg = f'{window}, {noverlap}'
            xp_assert_close(t, tr[:len(t)], err_msg=msg, check_dtype=False)
            xp_assert_close(x, xr[:len(x)], err_msg=msg)

    @skip_xp_backends("torch", "cupy", "array_api_strict",
            reasons=["torch hits messy np.pad codepath",
                     "lack of fft support in array_api_compat",
                     "moveaxis not available in array_api_strict"])
    def test_roundtrip_float32(self, xp):
        np.random.seed(1234)

        settings = [('hann', 1024, 256, 128)]

        for window, N, nperseg, noverlap in settings:
            t = xp.arange(N)
            x = 10*np.random.randn(size(t))
            x = xp.asarray(x, dtype=xp.float32)

            _, _, zz = stft(x, nperseg=nperseg, noverlap=noverlap,
                            window=window, detrend=None, padded=False)

            tr, xr = istft(zz, nperseg=nperseg, noverlap=noverlap,
                           window=window)

            msg = f'{window}, {noverlap}'
            xp_assert_close(t, t, err_msg=msg)
            xp_assert_close(x, xr, err_msg=msg, rtol=1e-4, atol=1e-5)

    @pytest.mark.parametrize('scaling', ['spectrum', 'psd'])
    def test_roundtrip_complex(self, scaling, xp):
        np.random.seed(1234)

        settings = [
                    ('boxcar', 100, 10, 0),           # Test no overlap
                    ('boxcar', 100, 10, 9),           # Test high overlap
                    ('bartlett', 101, 51, 26),        # Test odd nperseg
                    ('hann', 1024, 256, 128),         # Test defaults
                    (('tukey', 0.5), 1152, 256, 64),  # Test Tukey
                    ('hann', 1024, 256, 255),         # Test overlapped hann
                    ]

        for window, N, nperseg, noverlap in settings:
            t = xp.arange(N)
            x = 10*np.random.randn(size(t)) + 10j*np.random.randn(size(t))

            _, _, zz = stft(x, nperseg=nperseg, noverlap=noverlap,
                            window=window, detrend=None, padded=False,
                            return_onesided=False, scaling=scaling)

            tr, xr = istft(zz, nperseg=nperseg, noverlap=noverlap,
                           window=window, input_onesided=False,
                           scaling=scaling)

            msg = f'{window}, {nperseg}, {noverlap}'
            # NOTE: here and below it may not be surprising that namespace checks
            # fail because there are no input array types to key off of
            xp_assert_close(t, tr, err_msg=msg,
                            check_namespace=False, check_dtype=False)
            xp_assert_close(x, xr, err_msg=msg)

        # Check that asking for onesided switches to twosided
        with suppress_warnings() as sup:
            sup.filter(UserWarning,
                       "Input data is complex, switching to return_onesided=False")
            _, _, zz = stft(x, nperseg=nperseg, noverlap=noverlap,
                            window=window, detrend=None, padded=False,
                            return_onesided=True, scaling=scaling)

        tr, xr = istft(zz, nperseg=nperseg, noverlap=noverlap,
                       window=window, input_onesided=False, scaling=scaling)

        msg = f'{window}, {nperseg}, {noverlap}'
        xp_assert_close(t, tr, err_msg=msg, check_namespace=False, check_dtype=False)
        xp_assert_close(x, xr, err_msg=msg)

    @skip_xp_backends("torch", "cupy", "array_api_strict",
            reasons=["torch hits messy np.pad codepath",
                     "lack of fft support in array_api_compat",
                     "moveaxis not available in array_api_strict"])
    def test_roundtrip_boundary_extension(self, xp):
        np.random.seed(1234)

        # Test against boxcar, since window is all ones, and thus can be fully
        # recovered with no boundary extension

        settings = [
                    ('boxcar', 100, 10, 0),           # Test no overlap
                    ('boxcar', 100, 10, 9),           # Test high overlap
                    ]

        for window, N, nperseg, noverlap in settings:
            t = xp.arange(N)
            x = 10*np.random.randn(size(t))
            x = xp.asarray(x)

            _, _, zz = stft(x, nperseg=nperseg, noverlap=noverlap,
                           window=window, detrend=None, padded=True,
                           boundary=None)

            _, xr = istft(zz, noverlap=noverlap, window=window, boundary=False)

            for boundary in ['even', 'odd', 'constant', 'zeros']:
                _, _, zz_ext = stft(x, nperseg=nperseg, noverlap=noverlap,
                                window=window, detrend=None, padded=True,
                                boundary=boundary)

                _, xr_ext = istft(zz_ext, noverlap=noverlap, window=window,
                                boundary=True)

                msg = f'{window}, {noverlap}, {boundary}'
                xp_assert_close(x, xr, err_msg=msg)
                xp_assert_close(x, xr_ext, err_msg=msg)

    @skip_xp_backends("torch", "cupy", "array_api_strict",
            reasons=["torch hits messy np.pad codepath",
                     "lack of fft support in array_api_compat",
                     "moveaxis not available in array_api_strict"])
    def test_roundtrip_padded_signal(self, xp):
        np.random.seed(1234)

        settings = [
                    ('boxcar', 101, 10, 0),
                    ('hann', 1000, 256, 128),
                    ]

        for window, N, nperseg, noverlap in settings:
            t = xp.arange(N)
            x = 10*np.random.randn(size(t))
            x = xp.asarray(x)

            _, _, zz = stft(x, nperseg=nperseg, noverlap=noverlap,
                            window=window, detrend=None, padded=True)

            tr, xr = istft(zz, noverlap=noverlap, window=window)

            msg = f'{window}, {noverlap}'
            # Account for possible zero-padding at the end
            xp_assert_close(t, tr[:t.size], err_msg=msg, check_dtype=False)
            xp_assert_close(x, xr[:x.size], err_msg=msg)

    @skip_xp_backends("torch", "cupy", "array_api_strict",
            reasons=["torch hits messy np.pad codepath",
                     "lack of fft support in array_api_compat",
                     "moveaxis not available in array_api_strict"])
    def test_roundtrip_padded_FFT(self, xp):
        np.random.seed(1234)

        settings = [
                    ('hann', 1024, 256, 128, 512),
                    ('hann', 1024, 256, 128, 501),
                    ('boxcar', 100, 10, 0, 33),
                    (('tukey', 0.5), 1152, 256, 64, 1024),
                    ]

        for window, N, nperseg, noverlap, nfft in settings:
            t = xp.arange(N)
            x = 10*np.random.randn(size(t))
            x = xp.asarray(x)
            xc = x*xp.exp(xp.asarray(1j*xp.pi/4))

            # real signal
            _, _, z = stft(x, nperseg=nperseg, noverlap=noverlap, nfft=nfft,
                            window=window, detrend=None, padded=True)

            # complex signal
            _, _, zc = stft(xc, nperseg=nperseg, noverlap=noverlap, nfft=nfft,
                            window=window, detrend=None, padded=True,
                            return_onesided=False)

            tr, xr = istft(z, nperseg=nperseg, noverlap=noverlap, nfft=nfft,
                           window=window)

            tr, xcr = istft(zc, nperseg=nperseg, noverlap=noverlap, nfft=nfft,
                            window=window, input_onesided=False)

            msg = f'{window}, {noverlap}'
            xp_assert_close(t, tr, err_msg=msg, check_dtype=False)
            xp_assert_close(x, xr, err_msg=msg)
            xp_assert_close(xc, xcr, err_msg=msg)

    @skip_xp_backends("torch", "cupy", "array_api_strict",
            reasons=["torch hits messy np.pad codepath",
                     "lack of fft support in array_api_compat",
                     "moveaxis not available in array_api_strict"])
    def test_axis_rolling(self, xp):
        np.random.seed(1234)

        x_flat = np.random.randn(1024)
        x_flat = xp.asarray(x_flat)
        _, _, z_flat = stft(x_flat)

        for a in range(3):
            newshape = [1,]*3
            newshape[a] = -1
            x = x_flat.reshape(newshape)

            _, _, z_plus = stft(x, axis=a)  # Positive axis index
            _, _, z_minus = stft(x, axis=a-x.ndim)  # Negative axis index

            xp_assert_equal(z_flat, z_plus.squeeze(), err_msg=a)
            xp_assert_equal(z_flat, z_minus.squeeze(), err_msg=a-x.ndim)

        # z_flat has shape [n_freq, n_time]

        # Test vs. transpose
        _, x_transpose_m = istft(z_flat.T, time_axis=-2, freq_axis=-1)
        _, x_transpose_p = istft(z_flat.T, time_axis=0, freq_axis=1)

        xp_assert_close(x_flat, x_transpose_m, err_msg='istft transpose minus')
        xp_assert_close(x_flat, x_transpose_p, err_msg='istft transpose plus')

    @skip_xp_backends("torch", "cupy", "array_api_strict",
            reasons=["torch has issue with slice with negative step",
                     "lack of fft support in array_api_compat",
                     "moveaxis not available in array_api_strict"])
    def test_roundtrip_scaling(self, xp):
        """Verify behavior of scaling parameter. """
        # Create 1024 sample cosine signal with amplitude 2:
        # NOTE: don't use xp here, coerce expected value
        # after
        X = np.zeros(513, dtype=np.complex128)
        X[256] = 1024
        x = np.fft.irfft(X)
        x = xp.asarray(x)
        power_x = sum(x**2) / size(x)  # power of signal x is 2

        # Calculate magnitude-scaled STFT:
        Zs = stft(x, boundary='even', scaling='spectrum')[2]

        # Test round trip:
        x1 = istft(Zs, boundary=True, scaling='spectrum')[1]
        xp_assert_close(x1, x)

        # For a Hann-windowed 256 sample length FFT, we expect a peak at
        # frequency 64 (since it is 1/4 the length of X) with a height of 1
        # (half the amplitude). A Hann window of a perfectly centered sine has
        # the magnitude [..., 0, 0, 0.5, 1, 0.5, 0, 0, ...].
        # Note that in this case the 'even' padding works for the beginning
        # but not for the end of the STFT.
        xp_assert_close(abs(Zs[63, :-1]), 0.5, check_shape=False)
        xp_assert_close(abs(Zs[64, :-1]), 1, check_dtype=False, check_shape=False)
        xp_assert_close(abs(Zs[65, :-1]), 0.5, check_shape=False)
        # All other values should be zero:
        Zs[63:66, :-1] = 0
        # Note since 'rtol' does not have influence here, atol needs to be set:
        xp_assert_close(Zs[:, :-1], 0, atol=xp.finfo(Zs.dtype).resolution,
                        check_dtype=False,
                        check_shape=False)

        # Calculate two-sided psd-scaled STFT:
        #  - using 'even' padding since signal is axis symmetric - this ensures
        #    stationary behavior on the boundaries
        #  - using the two-sided transform allows determining the spectral
        #    power by `sum(abs(Zp[:, k])**2) / len(f)` for the k-th time slot.
        Zp = stft(x, return_onesided=False, boundary='even', scaling='psd')[2]

        # Calculate spectral power of Zd by summing over the frequency axis:
        psd_Zp = xp.sum(Zp.real**2 + Zp.imag**2, axis=0) / Zp.shape[0]
        # Spectral power of Zp should be equal to the signal's power:
        xp_assert_close(psd_Zp, power_x, check_shape=False)

        # Test round trip:
        x1 = istft(Zp, input_onesided=False, boundary=True, scaling='psd')[1]
        xp_assert_close(x1, x, check_dtype=False)

        # The power of the one-sided psd-scaled STFT can be determined
        # analogously (note that the two sides are not of equal shape):
        Zp0 = stft(x, return_onesided=True, boundary='even', scaling='psd')[2]

        # Since x is real, its Fourier transform is conjugate symmetric, i.e.,
        # the missing 'second side' can be expressed through the 'first side':
        Zp1 = xp.conj(Zp0[-2:0:-1, :])  # 'second side' is conjugate reversed
        xp_assert_close(Zp[:129, :], Zp0, atol=9e-16)
        xp_assert_close(Zp[129:, :], Zp1, atol=9e-16)

        # Calculate the spectral power:
        s2 = (xp.sum(Zp0.real ** 2 + Zp0.imag ** 2, axis=0) +
              xp.sum(Zp1.real ** 2 + Zp1.imag ** 2, axis=0))
        psd_Zp01 = s2 / (Zp0.shape[0] + Zp1.shape[0])
        xp_assert_close(psd_Zp01, power_x, check_shape=False)

        # Test round trip:
        x1 = istft(Zp0, input_onesided=True, boundary=True, scaling='psd')[1]
        assert_allclose(x1, x)


class TestSampledSpectralRepresentations:
    """Check energy/power relations from `Spectral Analysis` section in the user guide.

    A 32 sample cosine signal is used to compare the numerical to the expected results
    stated in :ref:`tutorial_SpectralAnalysis` in
    file ``doc/source/tutorial/signal.rst``
    """
    n: int = 32  #: number of samples
    T: float = 1/16  #: sampling interval
    a_ref: float = 3  #: amplitude of reference
    l_a: int = 3  #: index in fft for defining frequency of test signal

    x_ref: np.ndarray  #: reference signal
    X_ref: np.ndarray  #: two-sided FFT of x_ref
    E_ref: float  #: energy of signal
    P_ref: float  #: power of signal

    def setup_method(self):
        """Create Cosine signal with amplitude a from spectrum. """
        f = rfftfreq(self.n, self.T)
        X_ref = np.zeros_like(f)
        self.l_a = 3
        X_ref[self.l_a] = self.a_ref/2 * self.n  # set amplitude
        self.x_ref = irfft(X_ref)
        self.X_ref = fft(self.x_ref)

        # Closed form expression for continuous-time signal:
        self.E_ref = self.tau * self.a_ref**2 / 2  # energy of signal
        self.P_ref = self.a_ref**2 / 2  # power of signal

    @property
    def tau(self) -> float:
        """Duration of signal. """
        return self.n * self.T

    @property
    def delta_f(self) -> float:
        """Bin width """
        return 1 / (self.n * self.T)

    def test_reference_signal(self):
        """Test energy and power formulas. """
        # Verify that amplitude is a:
        assert_allclose(2*self.a_ref, np.ptp(self.x_ref), rtol=0.1)
        # Verify that energy expression for sampled signal:
        assert_allclose(self.T * sum(self.x_ref ** 2), self.E_ref)

        # Verify that spectral energy and power formulas are correct:
        sum_X_ref_squared = sum(self.X_ref.real**2 + self.X_ref.imag**2)
        assert_allclose(self.T/self.n * sum_X_ref_squared, self.E_ref)
        assert_allclose(1/self.n**2 * sum_X_ref_squared, self.P_ref)

    def test_windowed_DFT(self):
        """Verify spectral representations of windowed DFT.

        Furthermore, the scalings of `periodogram` and `welch` are verified.
        """
        w = hann(self.n, sym=False)
        c_amp, c_rms = abs(sum(w)), np.sqrt(sum(w.real**2 + w.imag**2))
        Xw = fft(self.x_ref*w)  # unnormalized windowed DFT

        # Verify that the *spectrum* peak is consistent:
        assert_allclose(self.tau * Xw[self.l_a] / c_amp, self.a_ref * self.tau / 2)
        # Verify that the *amplitude spectrum* peak is consistent:
        assert_allclose(Xw[self.l_a] / c_amp, self.a_ref/2)

        # Verify spectral power/energy equals signal's power/energy:
        X_ESD = self.tau * self.T * abs(Xw / c_rms)**2  # Energy Spectral Density
        X_PSD = self.T * abs(Xw / c_rms)**2  # Power Spectral Density
        assert_allclose(self.delta_f * sum(X_ESD), self.E_ref)
        assert_allclose(self.delta_f * sum(X_PSD), self.P_ref)

        # Verify scalings of periodogram:
        kw = dict(fs=1/self.T, window=w, detrend=False, return_onesided=False)
        _, P_mag = periodogram(self.x_ref, scaling='spectrum', **kw)
        _, P_psd = periodogram(self.x_ref, scaling='density', **kw)

        # Verify that periodogram calculates a squared magnitude spectrum:
        float_res = np.finfo(P_mag.dtype).resolution
        assert_allclose(P_mag, abs(Xw/c_amp)**2, atol=float_res*max(P_mag))
        # Verify that periodogram calculates a PSD:
        assert_allclose(P_psd, X_PSD, atol=float_res*max(P_psd))

        # Ensure that scaling of welch is the same as of periodogram:
        kw = dict(nperseg=len(self.x_ref), noverlap=0, **kw)
        assert_allclose(welch(self.x_ref, scaling='spectrum', **kw)[1], P_mag,
                        atol=float_res*max(P_mag))
        assert_allclose(welch(self.x_ref, scaling='density', **kw)[1], P_psd,
                        atol=float_res*max(P_psd))
