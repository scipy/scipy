from numpy.testing import assert_allclose

def test_fft_function():
    # Many NumPy symbols are imported into the scipy namespace, including
    # numpy.fft.fft as scipy.fft, conflicting with this module (gh-10253)
    import scipy
    scipy.random.seed(1234)

    # Callable before scipy.fft is imported
    x = scipy.randn(10) + 1j * scipy.randn(10)
    y = scipy.ifft(scipy.fft(x))
    assert_allclose(y, x)

    import scipy.fft
    # Callable after scipy.fft is imported
    z = scipy.ifft(scipy.fft(x))
    assert_allclose(z, x)
    assert_allclose(scipy.fft(x), scipy.fft.fft(x))
