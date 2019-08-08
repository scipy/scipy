from scipy import fft
import numpy as np
import pytest
from numpy.testing import assert_allclose
import multiprocessing

@pytest.fixture(scope='module')
def x():
    return np.random.randn(512, 128)  # Must be large enough to qualify for mt


@pytest.mark.slow
@pytest.mark.parametrize("func", [
    fft.fft, fft.ifft, fft.fft2, fft.ifft2, fft.fftn, fft.ifftn,
    fft.rfft, fft.irfft, fft.rfft2, fft.irfft2, fft.rfftn, fft.irfftn,
    fft.hfft, fft.ihfft, fft.hfft2, fft.ihfft2, fft.hfftn, fft.ihfftn,
    fft.dct, fft.idct, fft.dctn, fft.idctn,
    fft.dst, fft.idst, fft.dstn, fft.idstn,
])
@pytest.mark.parametrize("workers", [2, -1])
def test_threaded_same(x, func, workers):
    expected = func(x, workers=1)
    actual = func(x, workers=workers)
    assert_allclose(actual, expected)


def _mt_fft(x):
    return fft.fft(x, workers=2)


def test_mixed_threads_processes(x):
    # Test that the fft threadpool is safe to use before & after fork
    expect = fft.fft(x, workers=2)

    with multiprocessing.Pool(2) as p:
        res = p.map(_mt_fft, [x for _ in range(4)])

    for r in res:
        assert_allclose(r, expect)

    fft.fft(x, workers=2)

def test_invalid_workers(x):

    with pytest.raises(ValueError, match='workers must be'):
        fft.fft(x, workers=0)

    with pytest.raises(ValueError, match='workers must be'):
        fft.ifft(x, workers=-2)
