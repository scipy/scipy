import numpy as np
import scipy.fft
from scipy.fft import set_backend
import scipy.fftpack as fftpack
from  scipy.fft.tests import mock_backend

from numpy.testing import assert_allclose, assert_equal, assert_
import pytest

fnames = ('fft', 'fft2', 'fftn',
          'ifft', 'ifft2', 'ifftn',
          'rfft', 'rfft2', 'rfftn',
          'irfft', 'irfft2', 'irfftn',
          'dct', 'idct', 'dctn', 'idctn',
          'dst', 'idst', 'dstn', 'idstn')

np_funcs = (np.fft.fft, np.fft.fft2, np.fft.fftn,
            np.fft.ifft, np.fft.ifft2, np.fft.ifftn,
            np.fft.rfft, np.fft.rfft2, np.fft.rfftn,
            np.fft.irfft, np.fft.irfft2, np.fft.irfftn,
            fftpack.dct, fftpack.idct, fftpack.dctn, fftpack.idctn,
            fftpack.dst, fftpack.idst, fftpack.dstn, fftpack.idstn)

funcs = (scipy.fft.fft, scipy.fft.fft2, scipy.fft.fftn,
         scipy.fft.ifft, scipy.fft.ifft2, scipy.fft.ifftn,
         scipy.fft.rfft, scipy.fft.rfft2, scipy.fft.rfftn,
         scipy.fft.irfft, scipy.fft.irfft2, scipy.fft.irfftn,
         scipy.fft.dct, scipy.fft.idct, scipy.fft.dctn, scipy.fft.idctn,
         scipy.fft.dst, scipy.fft.idst, scipy.fft.dstn, scipy.fft.idstn)

mocks = (mock_backend.fft, mock_backend.fft2, mock_backend.fftn,
         mock_backend.ifft, mock_backend.ifft2, mock_backend.ifftn,
         mock_backend.rfft, mock_backend.rfft2, mock_backend.rfftn,
         mock_backend.irfft, mock_backend.irfft2, mock_backend.irfftn,
         mock_backend.dct, mock_backend.idct, mock_backend.dctn, mock_backend.idctn,
         mock_backend.dst, mock_backend.idst, mock_backend.dstn, mock_backend.idstn)


@pytest.mark.parametrize("func, np_func, mock", zip(funcs, np_funcs, mocks))
def test_backend_call(func, np_func, mock):
    x = np.arange(20).reshape((10,2))
    answer = np_func(x)
    assert_allclose(func(x), answer, atol=1e-10)

    with set_backend(mock_backend.MockBackend(), only=True):
        mock.number_calls = 0
        y = func(x)
        assert_equal(y, mock.return_value)
        assert_equal(mock.number_calls, 1)

    assert_allclose(func(x), answer, atol=1e-10)
