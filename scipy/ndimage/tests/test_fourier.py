import numpy as np
#from numpy import fft
from scipy._lib._array_api import (
    xp_assert_equal,
    assert_array_almost_equal,
    assert_almost_equal,
)

import pytest

from scipy import ndimage

from scipy.conftest import array_api_compatible
pytestmark = [array_api_compatible, pytest.mark.usefixtures("skip_xp_backends")]
skip_xp_backends = pytest.mark.skip_xp_backends


class TestNdimageFourier:

    @pytest.mark.parametrize('shape', [(32, 16), (31, 15), (1, 10)])
    @pytest.mark.parametrize('dtype, dec', [("float32", 6), ("float64", 14)])
    def test_fourier_gaussian_real01(self, shape, dtype, dec, xp):
        fft = getattr(xp, 'fft')
        dtype = getattr(xp, dtype)

        a = xp.zeros(shape, dtype=dtype)
        a[0, 0] = 1.0
        a = fft.rfft(a, shape[0], 0)
        a = fft.fft(a, shape[1], 1)
        a = ndimage.fourier_gaussian(a, [5.0, 2.5], shape[0], 0)
        a = fft.ifft(a, shape[1], 1)
        a = fft.irfft(a, shape[0], 0)
        assert_almost_equal(ndimage.sum(a), xp.asarray(1), decimal=dec)

    @pytest.mark.parametrize('shape', [(32, 16), (31, 15)])
    @pytest.mark.parametrize('dtype, dec', [("complex64", 6), ("complex128", 14)])
    def test_fourier_gaussian_complex01(self, shape, dtype, dec, xp):
        fft = getattr(xp, 'fft')
        dtype = getattr(xp, dtype)

        a = xp.zeros(shape, dtype=dtype)
        a[0, 0] = 1.0
        a = fft.fft(a, shape[0], 0)
        a = fft.fft(a, shape[1], 1)
        a = ndimage.fourier_gaussian(a, [5.0, 2.5], -1, 0)
        a = fft.ifft(a, shape[1], 1)
        a = fft.ifft(a, shape[0], 0)
        assert_almost_equal(ndimage.sum(a.real), xp.asarray(1.0), decimal=dec)

    @pytest.mark.parametrize('shape', [(32, 16), (31, 15), (1, 10)])
    @pytest.mark.parametrize('dtype, dec', [("float32", 6), ("float64", 14)])
    def test_fourier_uniform_real01(self, shape, dtype, dec, xp):
        fft = getattr(xp, 'fft')
        dtype = getattr(xp, dtype)

        a = xp.zeros(shape, dtype=dtype)
        a[0, 0] = 1.0
        a = fft.rfft(a, shape[0], 0)
        a = fft.fft(a, shape[1], 1)
        a = ndimage.fourier_uniform(a, [5.0, 2.5], shape[0], 0)
        a = fft.ifft(a, shape[1], 1)
        a = fft.irfft(a, shape[0], 0)
        assert_almost_equal(ndimage.sum(a), xp.asarray(1.0), decimal=dec)

    @pytest.mark.parametrize('shape', [(32, 16), (31, 15)])
    @pytest.mark.parametrize('dtype, dec', [("complex64", 6), ("complex128", 14)])
    def test_fourier_uniform_complex01(self, shape, dtype, dec, xp):
        fft = getattr(xp, 'fft')
        dtype = getattr(xp, dtype)

        a = xp.zeros(shape, dtype=dtype)
        a[0, 0] = 1.0
        a = fft.fft(a, shape[0], 0)
        a = fft.fft(a, shape[1], 1)
        a = ndimage.fourier_uniform(a, [5.0, 2.5], -1, 0)
        a = fft.ifft(a, shape[1], 1)
        a = fft.ifft(a, shape[0], 0)
        assert_almost_equal(ndimage.sum(a.real), xp.asarray(1.0), decimal=dec)

    @pytest.mark.parametrize('shape', [(32, 16), (31, 15)])
    @pytest.mark.parametrize('dtype, dec', [("float32", 4), ("float64", 11)])
    def test_fourier_shift_real01(self, shape, dtype, dec, xp):
        fft = getattr(xp, 'fft')
        dtype = getattr(xp, dtype)

        expected = xp.arange(shape[0] * shape[1], dtype=dtype)
        expected.shape = shape
        a = fft.rfft(expected, shape[0], 0)
        a = fft.fft(a, shape[1], 1)
        a = ndimage.fourier_shift(a, [1, 1], shape[0], 0)
        a = fft.ifft(a, shape[1], 1)
        a = fft.irfft(a, shape[0], 0)
        assert_array_almost_equal(a[1:, 1:], expected[:-1, :-1], decimal=dec)
        assert_array_almost_equal(a.imag, xp.zeros(shape), decimal=dec)

    @pytest.mark.parametrize('shape', [(32, 16), (31, 15)])
    @pytest.mark.parametrize('dtype, dec', [("complex64", 4), ("complex128", 11)])
    def test_fourier_shift_complex01(self, shape, dtype, dec, xp):
        fft = getattr(xp, 'fft')
        dtype = getattr(xp, dtype)

        expected = xp.arange(shape[0] * shape[1], dtype=dtype)
        expected.shape = shape
        a = fft.fft(expected, shape[0], 0)
        a = fft.fft(a, shape[1], 1)
        a = ndimage.fourier_shift(a, [1, 1], -1, 0)
        a = fft.ifft(a, shape[1], 1)
        a = fft.ifft(a, shape[0], 0)
        assert_array_almost_equal(a.real[1:, 1:], expected[:-1, :-1], decimal=dec)
        assert_array_almost_equal(a.imag, xp.zeros(shape), decimal=dec)

    @pytest.mark.parametrize('shape', [(32, 16), (31, 15), (1, 10)])
    @pytest.mark.parametrize('dtype, dec', [("float32", 5), ("float64", 14)])
    def test_fourier_ellipsoid_real01(self, shape, dtype, dec, xp):
        fft = getattr(xp, 'fft')
        dtype = getattr(xp, dtype)

        a = xp.zeros(shape, dtype=dtype)
        a[0, 0] = 1.0
        a = fft.rfft(a, shape[0], 0)
        a = fft.fft(a, shape[1], 1)
        a = ndimage.fourier_ellipsoid(a, [5.0, 2.5],
                                      shape[0], 0)
        a = fft.ifft(a, shape[1], 1)
        a = fft.irfft(a, shape[0], 0)
        assert_almost_equal(ndimage.sum(a), xp.asarray(1.0), decimal=dec)

    @pytest.mark.parametrize('shape', [(32, 16), (31, 15)])
    @pytest.mark.parametrize('dtype, dec', [("complex64", 5), ("complex128", 14)])
    def test_fourier_ellipsoid_complex01(self, shape, dtype, dec, xp):
        fft = getattr(xp, 'fft')
        dtype = getattr(xp, dtype)

        a = xp.zeros(shape, dtype=dtype)
        a[0, 0] = 1.0
        a = fft.fft(a, shape[0], 0)
        a = fft.fft(a, shape[1], 1)
        a = ndimage.fourier_ellipsoid(a, [5.0, 2.5], -1, 0)
        a = fft.ifft(a, shape[1], 1)
        a = fft.ifft(a, shape[0], 0)
        assert_almost_equal(ndimage.sum(a.real), xp.asarray(1.0), decimal=dec)

    def test_fourier_ellipsoid_unimplemented_ndim(self):
        # arrays with ndim > 3 raise NotImplementedError
        x = np.ones((4, 6, 8, 10), dtype="complex128")
        with pytest.raises(NotImplementedError):
            ndimage.fourier_ellipsoid(x, 3)

    def test_fourier_ellipsoid_1d_complex(self, xp):
        # expected result of 1d ellipsoid is the same as for fourier_uniform
        for shape in [(32, ), (31, )]:
            for type_, dec in zip(["complex64", "complex128"], [5, 14]):
                x = np.ones(shape, dtype=type_)
                a = ndimage.fourier_ellipsoid(x, 5, -1, 0)
                b = ndimage.fourier_uniform(x, 5, -1, 0)
                assert_array_almost_equal(a, b, decimal=dec)

    @pytest.mark.parametrize('shape', [(0, ), (0, 10), (10, 0)])
    @pytest.mark.parametrize('dtype', ["float32", "float64",
                                       "complex64", "complex128"])
    @pytest.mark.parametrize('test_func',
                             [ndimage.fourier_ellipsoid,
                              ndimage.fourier_gaussian,
                              ndimage.fourier_uniform])
    def test_fourier_zero_length_dims(self, shape, dtype, test_func, xp):
        dtype = getattr(xp, dtype)
        a = np.ones(shape, dtype=dtype)
        b = test_func(a, 3)
        xp_assert_equal(a, b)
