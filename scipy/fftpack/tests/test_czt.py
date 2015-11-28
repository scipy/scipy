# This program is public domain
# Authors: Paul Kienzle, Nadav Horesh
'''
A unit test module for czt.py
'''
from __future__ import division, absolute_import, print_function

from numpy.testing import (run_module_suite, assert_, assert_allclose,
                           assert_raises, dec)
from scipy.fftpack.czt import (czt, zoomfft, czt_points, CZT, ZoomFFT)
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
    assert_(err < 1e-10, "error for direct transform is %g" % (err,))


def check_zoomfft(x):
    # Check that zoomfft is the equivalent of normal fft
    y = fft(x)
    y1 = zoomfft(x, [0, 2-2./len(y)])
    err = norm(y-y1)/norm(y)
    assert_(err < 1e-10, "error for direct transform is %g" % (err,))

    # Test fn scalar
    y1 = zoomfft(x, 2-2./len(y))
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


def test_1D():
    # Test of 1D version of the transforms

    np.random.seed(0)  # Deterministic randomness

    # Random signals
    lengths = np.random.randint(8, 200, 20)
    np.append(lengths, 1)
    for length in lengths:
        x = np.random.random(length)
        yield check_zoomfft, x
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


def test_large_prime_lengths():
    np.random.seed(0)  # Deterministic randomness
    for N in (101, 1009, 10007):
        x = np.random.rand(N)
        y = fft(x)
        y1 = czt(x)
        assert_allclose(y, y1, rtol=1e-9)


@dec.slow
def test_czt_vs_fft():
    np.random.seed(123)
    random_lengths = np.random.exponential(100000, size=10).astype('int')
    for n in random_lengths:
        a = np.random.randn(n)
        assert_allclose(czt(a), fft(a), rtol=1e-7, atol=1e-8)


def test_empty_input():
    assert_raises(ValueError, czt, [])
    assert_raises(ValueError, zoomfft, [], 0.5)


def test_0_rank_input():
    assert_raises(IndexError, czt, 5)
    assert_raises(IndexError, zoomfft, 5, 0.5)


def test_czt_math():
    for impulse in ([0, 0, 1],
                    [0, 0, 1, 0, 0],
                    np.concatenate(([0, 0, 1], np.zeros(100)))):
        for m in (1, 3, 5, 8, 101, 1021):
            for a in (1, 2, 0.5, 1.1):
                for w in (None, 0.7+0.7j):
                    # z-transform of an impulse is 1 everywhere
                    assert_allclose(czt(impulse[2:], m=m, a=a),
                                    np.ones(m))

                    # z-transform of a delayed impulse is z**-1
                    assert_allclose(czt(impulse[1:], m=m, a=a),
                                    czt_points(m=m, a=a)**-1)

                    # z-transform of a 2-delayed impulse is z**-2
                    assert_allclose(czt(impulse, m=m, a=a),
                                    czt_points(m=m, a=a)**-2)


def test_int_args():
    # Integer argument `a` was producing all 0s
    assert_allclose(abs(czt([0, 1], m=10, a=2)), 0.5*np.ones(10))
    assert_allclose(czt_points(11, w=2), 1/(2**np.arange(11)))


def test_conflicting_args():
    # Cannot specify scale and w at the same time
    assert_raises(ValueError, czt, x=np.ones(8),
                  w=0.70710678118654746+1j*0.70710678118654746,
                  scale=1)


def test_czt_points():
    for N in (1, 2, 3, 8, 11, 100, 101, 10007):
        assert_allclose(czt_points(N), np.exp(2j*np.pi*np.arange(N)/N))

    assert_allclose(czt_points(7, w=1), np.ones(7))
    assert_allclose(czt_points(11, w=2.), 1/(2**np.arange(11)))

    func = CZT(12, m=11, w=2., a=1)
    assert_allclose(func.points(), 1/(2**np.arange(11)))


def test_czt_points_errors():
    # Invalid number of points
    assert_raises(ValueError, czt_points, 0)
    assert_raises(ValueError, czt_points, -11)
    assert_raises(ValueError, czt_points, 5.5)

    # Cannot specify scale and w at the same time
    assert_raises(ValueError, czt_points, 5,
                  0.70710678118654746+1j*0.70710678118654746, 1, 1)


def test_invalid_size():
    # Data size doesn't match function's expected size
    for myfunc in (CZT(100), ZoomFFT(100, 0.2)):
        assert_raises(ValueError, myfunc, np.arange(5))

    # Nonsense input and output sizes
    # Numpy and Scipy fft() give ValueError for 0 output size, so we do, too
    for size in (0, -5, 3.5):
        assert_raises(ValueError, CZT, size, 3)
        assert_raises(ValueError, ZoomFFT, size, 0.2, 3)
        assert_raises(ValueError, CZT, 3, size)
        assert_raises(ValueError, ZoomFFT, 3, 0.2, size)
        assert_raises(ValueError, czt, [1, 2, 3], size)
        assert_raises(ValueError, zoomfft, [1, 2, 3], 0.2, size)


def test_invalid_range():
    assert_raises(ValueError, ZoomFFT, 100, [1, 2, 3])


if __name__ == '__main__':
    np.random.seed()

    old_settings = np.seterr(all='ignore')  # seterr to known value
    np.seterr(over='raise')  # Test should fail if internal overflow occurs

    run_module_suite()

    np.seterr(**old_settings)  # reset to default
