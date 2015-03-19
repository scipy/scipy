# This program is public domain
# Authors: Paul Kienzle, Nadav Horesh
'''
A unit test module for czt.py
'''
from __future__ import division, absolute_import, print_function

from numpy.testing import run_module_suite, assert_, assert_allclose, assert_raises
from czt import czt, zoomfft, scaledfft, czt_points
import numpy as np

fft = np.fft.fft
fftshift = np.fft.fftshift
norm = np.linalg.norm


def check_czt(x):
    # Check that czt is the equivalent of normal fft
    y = fft(x)
    y1 = czt(x)
    err = norm(y-y1)/norm(y)
    assert_(err < 1e-10, "error for direct transform is %g" % (err,))

    # Check that interpolated czt is the equivalent of normal fft
    y = fft(x, 100*len(x))
    y1 = czt(x, 100*len(x))
    err = norm(y-y1)/norm(y)
    assert_(err < 3e-10, "error for direct transform is %g" % (err,))


def check_zoomfft(x):
    # Check that zoomfft is the equivalent of normal fft
    y = fft(x)
    y1 = zoomfft(x, [0, 2-2./len(y)])
    err = norm(y-y1)/norm(y)
    assert_(err < 1e-10, "error for direct transform is %g" % (err,))

    # Check that zoomfft with oversampling is equivalent to zero padding
    over = 10
    yover = fft(x, over*len(x))
    y2 = zoomfft(x, [0, 2-2./len(yover)], m=len(yover))
    err = norm(yover-y2)/norm(yover)
    assert_(err < 1e-10, "error for oversampling is %g" % (err,))

    # Check that zoomfft works on a subrange
    w = np.linspace(0, 2-2./len(x), len(x))
    f1, f2 = w[3], w[6]
    y3 = zoomfft(x, [f1, f2], m=3*over+1)
    idx3 = slice(3*over, 6*over+1)
    err = norm(yover[idx3]-y3)/norm(yover[idx3])
    assert_(err < 1e-10, "error for subrange is %g" % (err,))


def check_scaledfft(x):
    # Test for unscaled equivalence with FFT
    n = len(x)
    assert_(norm(fft(x)-scaledfft(x)) < 2E-10)

    # Check that scaled FFT matches
    n = n//4 * 4
    x = x[:n]
    st = n//4
    m = n//2
    sliced_fft = fftshift(fft(x))[st:st+m]
    scl_fft = fftshift(scaledfft(x, scale=0.5, m=m))
    assert_(norm(sliced_fft-scl_fft) < 2E-10)


def test_1D():
    # Test of 1D version of the transforms

    # Deterministic randomness
    np.random.seed(0)

    # Random signals
    lengths = np.random.randint(8, 200, 20)
    np.append(lengths, 1)
    for length in lengths:
        x = np.random.random(length)
        yield check_zoomfft, x
        yield check_scaledfft, x
        yield check_czt, x

    # Gauss
    t = np.linspace(-2, 2, 128)
    x = np.exp(-t**2/0.01)
    check_zoomfft(x)

    # Linear
    x = [1, 2, 3, 4, 5, 6, 7]
    check_zoomfft(x)

    # Check near powers of two
    check_zoomfft(range(126-31))
    check_zoomfft(range(127-31))
    check_zoomfft(range(128-31))
    check_zoomfft(range(129-31))
    check_zoomfft(range(130-31))

    # Check transform on n-D array input
    x = np.reshape(np.arange(3*2*28), (3, 2, 28))
    y1 = zoomfft(x, [0, 2-2./28])
    y2 = zoomfft(x[2, 0, :], [0, 2-2./28])
    err = np.linalg.norm(y2-y1[2, 0])
    assert err < 1e-15, "error for n-D array is %g" % (err,)

    # Random (not a test condition)
    x = np.random.rand(101)
    check_zoomfft(x)

    # Spikes
    t = np.linspace(0, 1, 128)
    x = np.sin(2*np.pi*t*5)+np.sin(2*np.pi*t*13)
    check_zoomfft(x)

    # Sines
    x = np.zeros(100, dtype=complex)
    x[[1, 5, 21]] = 1
    check_zoomfft(x)

    # Sines plus complex component
    x += 1j*np.linspace(0, 0.5, x.shape[0])
    check_zoomfft(x)

    # Scaled FFT on complex sines
    x += 1j*np.linspace(0, 0.5, x.shape[0])
    check_scaledfft(x)


def test_large_prime():
    np.random.seed(0)
    x = np.random.rand(10007)
    y = fft(x)
    y1 = czt(x)
    assert_allclose(y, y1, rtol=1e-6)


def test_empty_input():
    assert_raises(ValueError, czt, [])
    assert_raises(ValueError, scaledfft, [])
    assert_raises(ValueError, zoomfft, [], 0.5)


def test_0_rank_input():
    assert_raises(IndexError, czt, 5)
    assert_raises(IndexError, scaledfft, 5)
    assert_raises(IndexError, zoomfft, 5, 0.5)


def test_int_args():
    # Integer argument `a` was producing all 0s
    assert_allclose(abs(czt([0, 1], m=10, a=2)), 0.5*np.ones(10))
    assert_allclose(czt_points(11, w=2), 1/(2**np.arange(11)))

def test_czt_points():
    for N in (1, 2, 3, 8, 11, 100, 101, 10007):
        assert_allclose(czt_points(N), np.exp(2j*np.pi*np.arange(N)/N))

    assert_allclose(czt_points(7, w=1), np.ones(7))
    assert_allclose(czt_points(11, w=2.), 1/(2**np.arange(11)))

    # Test every branch with every combination of args

if __name__ == '__main__':
    np.random.seed()
    run_module_suite()
