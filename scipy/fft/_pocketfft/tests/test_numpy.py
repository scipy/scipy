from __future__ import division, absolute_import, print_function

import numpy as np
import pytest
from numpy.random import random
from numpy.testing import (
        assert_array_almost_equal, assert_array_equal, assert_raises,
        )
from scipy.fft import _pocketfft as pfft


def fft1(x):
    L = len(x)
    phase = -2j*np.pi*(np.arange(L)/float(L))
    phase = np.arange(L).reshape(-1, 1) * phase
    return np.sum(x*np.exp(phase), axis=1)


class TestFFTShift(object):

    def test_fft_n(self):
        assert_raises(ValueError, pfft.fft, [1, 2, 3], 0)


class TestFFT1D(object):

    def test_identity(self):
        maxlen = 512
        x = random(maxlen) + 1j*random(maxlen)
        xr = random(maxlen)
        for i in range(1,maxlen):
            assert_array_almost_equal(pfft.ifft(pfft.fft(x[0:i])), x[0:i],
                                      decimal=12)
            assert_array_almost_equal(pfft.irfft(pfft.rfft(xr[0:i]),i),
                                      xr[0:i], decimal=12)

    def test_fft(self):
        x = random(30) + 1j*random(30)
        assert_array_almost_equal(fft1(x), pfft.fft(x))
        assert_array_almost_equal(fft1(x) / np.sqrt(30),
                                  pfft.fft(x, norm="ortho"))

    def test_ifft(self):
        x = random(30) + 1j*random(30)
        assert_array_almost_equal(x, pfft.ifft(pfft.fft(x)))
        assert_array_almost_equal(
            x, pfft.ifft(pfft.fft(x, norm="ortho"), norm="ortho"))

    def test_fft2(self):
        x = random((30, 20)) + 1j*random((30, 20))
        assert_array_almost_equal(pfft.fft(pfft.fft(x, axis=1), axis=0),
                                  pfft.fft2(x))
        assert_array_almost_equal(pfft.fft2(x) / np.sqrt(30 * 20),
                                  pfft.fft2(x, norm="ortho"))

    def test_ifft2(self):
        x = random((30, 20)) + 1j*random((30, 20))
        assert_array_almost_equal(pfft.ifft(pfft.ifft(x, axis=1), axis=0),
                                  pfft.ifft2(x))
        assert_array_almost_equal(pfft.ifft2(x) * np.sqrt(30 * 20),
                                  pfft.ifft2(x, norm="ortho"))

    def test_fftn(self):
        x = random((30, 20, 10)) + 1j*random((30, 20, 10))
        assert_array_almost_equal(
            pfft.fft(pfft.fft(pfft.fft(x, axis=2), axis=1), axis=0),
            pfft.fftn(x))
        assert_array_almost_equal(pfft.fftn(x) / np.sqrt(30 * 20 * 10),
                                  pfft.fftn(x, norm="ortho"))

    def test_ifftn(self):
        x = random((30, 20, 10)) + 1j*random((30, 20, 10))
        assert_array_almost_equal(
            pfft.ifft(pfft.ifft(pfft.ifft(x, axis=2), axis=1), axis=0),
            pfft.ifftn(x))
        assert_array_almost_equal(pfft.ifftn(x) * np.sqrt(30 * 20 * 10),
                                  pfft.ifftn(x, norm="ortho"))

    def test_rfft(self):
        x = random(30)
        for n in [x.size, 2*x.size]:
            for norm in [None, 'ortho']:
                assert_array_almost_equal(
                    pfft.fft(x, n=n, norm=norm)[:(n//2 + 1)],
                    pfft.rfft(x, n=n, norm=norm))
            assert_array_almost_equal(pfft.rfft(x, n=n) / np.sqrt(n),
                                      pfft.rfft(x, n=n, norm="ortho"))

    def test_irfft(self):
        x = random(30)
        assert_array_almost_equal(x, pfft.irfft(pfft.rfft(x)))
        assert_array_almost_equal(
            x, pfft.irfft(pfft.rfft(x, norm="ortho"), norm="ortho"))

    def test_rfft2(self):
        x = random((30, 20))
        assert_array_almost_equal(pfft.fft2(x)[:, :11], pfft.rfft2(x))
        assert_array_almost_equal(pfft.rfft2(x) / np.sqrt(30 * 20),
                                  pfft.rfft2(x, norm="ortho"))

    def test_irfft2(self):
        x = random((30, 20))
        assert_array_almost_equal(x, pfft.irfft2(pfft.rfft2(x)))
        assert_array_almost_equal(
            x, pfft.irfft2(pfft.rfft2(x, norm="ortho"), norm="ortho"))

    def test_rfftn(self):
        x = random((30, 20, 10))
        assert_array_almost_equal(pfft.fftn(x)[:, :, :6], pfft.rfftn(x))
        assert_array_almost_equal(pfft.rfftn(x) / np.sqrt(30 * 20 * 10),
                                  pfft.rfftn(x, norm="ortho"))

    def test_irfftn(self):
        x = random((30, 20, 10))
        assert_array_almost_equal(x, pfft.irfftn(pfft.rfftn(x)))
        assert_array_almost_equal(
            x, pfft.irfftn(pfft.rfftn(x, norm="ortho"), norm="ortho"))

    @pytest.mark.skip(reason="hfft not currently implemented")
    def test_hfft(self):
        x = random(14) + 1j*random(14)
        x_herm = np.concatenate((random(1), x, random(1)))
        x = np.concatenate((x_herm, x[::-1].conj()))
        assert_array_almost_equal(pfft.fft(x), pfft.hfft(x_herm))
        assert_array_almost_equal(pfft.hfft(x_herm) / np.sqrt(30),
                                  pfft.hfft(x_herm, norm="ortho"))

    @pytest.mark.skip(reason="ihfft not currently implemented")
    def test_ihfft(self):
        x = random(14) + 1j*random(14)
        x_herm = np.concatenate((random(1), x, random(1)))
        x = np.concatenate((x_herm, x[::-1].conj()))
        assert_array_almost_equal(x_herm, pfft.ihfft(pfft.hfft(x_herm)))
        assert_array_almost_equal(
            x_herm, pfft.ihfft(pfft.hfft(x_herm, norm="ortho"),
                                 norm="ortho"))

    @pytest.mark.parametrize("op", [pfft.fftn, pfft.ifftn,
                                    pfft.rfftn, pfft.irfftn])
    def test_axes(self, op):
        x = random((30, 20, 10))
        axes = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
        for a in axes:
            op_tr = op(np.transpose(x, a))
            tr_op = np.transpose(op(x, axes=a), a)
            assert_array_almost_equal(op_tr, tr_op)

    def test_all_1d_norm_preserving(self):
        # verify that round-trip transforms are norm-preserving
        x = random(30)
        x_norm = np.linalg.norm(x)
        n = x.size * 2
        func_pairs = [(pfft.fft, pfft.ifft),
                      (pfft.rfft, pfft.irfft),
                      # hfft: order so the first function takes x.size samples
                      #       (necessary for comparison to x_norm above)
                      #(pfft.ihfft, pfft.hfft),
                      ]
        for forw, back in func_pairs:
            for n in [x.size, 2*x.size]:
                for norm in [None, 'ortho']:
                    tmp = forw(x, n=n, norm=norm)
                    tmp = back(tmp, n=n, norm=norm)
                    assert_array_almost_equal(x_norm,
                                              np.linalg.norm(tmp))

    @pytest.mark.parametrize("dtype", [np.half, np.single, np.double,
                                       np.longdouble])
    def test_dtypes(self, dtype):
        # make sure that all input precisions are accepted
        x = random(30).astype(dtype)
        assert_array_almost_equal(pfft.ifft(pfft.fft(x)), x)
        assert_array_almost_equal(pfft.irfft(pfft.rfft(x)), x)


@pytest.mark.parametrize(
        "dtype",
        [np.float32, np.float64, np.longfloat,
         np.complex64, np.complex128, np.longcomplex])
@pytest.mark.parametrize("order", ["F", 'non-contiguous'])
@pytest.mark.parametrize(
        "fft",
        [pfft.fft, pfft.fft2, pfft.fftn,
         pfft.ifft, pfft.ifft2, pfft.ifftn])
def test_fft_with_order(dtype, order, fft):
    # Check that FFT/IFFT produces identical results for C, Fortran and
    # non contiguous arrays
    rng = np.random.RandomState(42)
    X = rng.rand(8, 7, 13).astype(dtype, copy=False)
    if order == 'F':
        Y = np.asfortranarray(X)
    else:
        # Make a non contiguous array
        Y = X[::-1]
        X = np.ascontiguousarray(X[::-1])

    if fft.__name__.endswith('fft'):
        for axis in range(3):
            X_res = fft(X, axis=axis)
            Y_res = fft(Y, axis=axis)
            assert_array_almost_equal(X_res, Y_res)
    elif fft.__name__.endswith(('fft2', 'fftn')):
        axes = [(0, 1), (1, 2), (0, 2)]
        if fft.__name__.endswith('fftn'):
            axes.extend([(0,), (1,), (2,), None])
        for ax in axes:
            X_res = fft(X, axes=ax)
            Y_res = fft(Y, axes=ax)
            assert_array_almost_equal(X_res, Y_res)
    else:
        raise ValueError
