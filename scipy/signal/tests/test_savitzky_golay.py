
import numpy as np
from numpy.testing import run_module_suite, assert_allclose, assert_equal

from scipy.ndimage import convolve1d

from scipy.signal import savgol_coeffs, savgol_filter


def alt_sg_coeffs(window_length, polyorder):
    """This is an alternative implementation of the SG coefficients.

    It uses numpy.polyfit and numpy.polyval.  The results should be
    equivalent to those of savgol_coeffs(), but this implementation
    is slower.

    window_length should be odd.
    """

    ind = np.arange(window_length) - window_length // 2
    a = np.polyfit(ind, ind == 0, polyorder)
    h = np.polyval(a, ind)
    return h


def test_sg_coeffs_trivial():
    """Test some trivial edge cases of savgolay_coeffs()"""
    h = savgol_coeffs(1, 0)
    assert_equal(h, [1])

    h = savgol_coeffs(3, 2)
    assert_equal(h, [0, 1, 0])

    h = savgol_coeffs(5, 4)
    assert_equal(h, [0, 0, 1, 0, 0])


def test_sg_coeffs_compare():
    """ Compare savitzky_golay_fir() to alt_sg_coeffs()."""
    for window_length in [11, 21, 31]:
        for order in [2, 4, 6]:
            h1 = savgol_coeffs(window_length, order)
            h2 = alt_sg_coeffs(window_length, order)
            assert_allclose(h1, h2, atol=1e-14,
                            err_msg=("window_length = %d, order = %d" %
                                     (window_length, order)))


def test_sg_coeffs_exact():
    # For this test, window_length must be odd.
    window_length = 9
    half_wind = window_length // 2
    x = np.arange(21.0)
    # The data is a cubic polynomial.  We'll use an order 4
    # SG filter, so the filtered values should equal the input data
    # (except within half window_length of the egdes).
    y = 0.5 * x ** 3 - x
    h = savgol_coeffs(7, 4)
    y0 = convolve1d(y, h)
    assert_allclose(y0[half_wind:-half_wind], y[half_wind:-half_wind])

    # Check the same input, but use deriv=1.  dy is the exact result.
    dy = 1.5 * x ** 2 - 1
    h = savgol_coeffs(7, 4, deriv=1)
    y1 = convolve1d(y, h)
    assert_allclose(y1[half_wind:-half_wind], dy[half_wind:-half_wind])

    # Check the same input, but use deriv=2. d2y is the exact result.
    d2y = 3.0 * x
    h = savgol_coeffs(7, 4, deriv=2)
    y2 = convolve1d(y, h)
    assert_allclose(y2[half_wind:-half_wind], d2y[half_wind:-half_wind])


def test_sg_filter_trivial():
    """ Test some trivial edge cases for savgol_filter()."""
    x = np.array([1.0])
    y = savgol_filter(x, 1, 0)
    assert_equal(y, [1.0])

    # Input is a single value.  With a window length of 3 and polyorder 1,
    # the value in y is from the straight-line fit of (-1,0), (0,3) and
    # (1, 0) at 0. This is just the average of the three values, hence 1.0.
    x = np.array([3.0])
    y = savgol_filter(x, 3, 1)
    assert_equal(y, [1.0])

    x = np.array([3.0])
    y = savgol_filter(x, 3, 1, mode='nearest')
    assert_equal(y, [3.0])


def test_sg_filter_basic():
    """ Some basic test cases for savgol_filter()."""
    x = np.array([1.0, 2.0, 1.0])
    y = savgol_filter(x, 3, 1)
    assert_allclose(y, [1.0, 4.0 / 3, 1.0])

    y = savgol_filter(x, 3, 1, mode='mirror')
    assert_allclose(y, [5.0 / 3, 4.0 / 3, 5.0 / 3])


def test_sg_filter_2d():
    x = np.array([[1.0, 2.0, 1.0],
                  [2.0, 4.0, 2.0]])
    expected = np.array([[1.0, 4.0 / 3, 1.0],
                         [2.0, 8.0 / 3, 2.0]])
    y = savgol_filter(x, 3, 1)
    assert_allclose(y, expected)

    y = savgol_filter(x.T, 3, 1, axis=0)
    assert_allclose(y, expected.T)


if __name__ == "__main__":
    run_module_suite()
