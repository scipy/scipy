from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import (run_module_suite, assert_allclose,
                           assert_array_equal)
from scipy.optimize import fmin
from scipy.signal.signaltools import _filtfilt_gust
from scipy.signal import filtfilt, lfilter, lfilter_zi, ellip, tf2zpk


def func(ics, b, a, x):
    """Objective function used in filtfilt_opt."""
    m = max(len(a), len(b)) - 1
    z0f = ics[:m]
    z0b = ics[m:]
    y_f = lfilter(b, a, x, zi=z0f)[0]
    y_fb = lfilter(b, a, y_f[::-1], zi=z0b)[0][::-1]

    y_b = lfilter(b, a, x[::-1], zi=z0b)[0][::-1]
    y_bf = lfilter(b, a, y_b, zi=z0f)[0]
    value = np.sum((y_fb - y_bf)**2)
    return value


def filtfilt_zi(x, b, a, z0f, z0b):
    """Forward-backward filter using the given initial conditions."""
    y_b = lfilter(b, a, x[::-1], zi=z0b)[0][::-1]
    y_bf = lfilter(b, a, y_b, zi=z0f)[0]
    return y_bf


def filtfilt_opt(b, a, x):
    """
    An alternative implementation of filtfilt with Gustafsson edges.

    This function computes the same result as `_filtfilt_gust`, but only
    1-d arrays are accepted.  The problem is solved using `fmin` from
    `scipy.optimize`.  `_filtfilt_gust` is significanly faster than this
    implementation.
    """
    m = max(len(a), len(b)) - 1
    zi = lfilter_zi(b, a)
    ics = np.concatenate((x[:m].mean()*zi, x[-m:].mean()*zi))
    result = fmin(func, ics, args=(b, a, x),
                  xtol=1e-10, ftol=1e-12,
                  maxfun=10000, maxiter=10000,
                  full_output=True, disp=False)
    opt, fopt, niter, funcalls, warnflag = result
    if warnflag > 0:
        raise RuntimeError("minimization failed in filtfilt_opt: "
                           "warnflag=%d" % warnflag)
    z0f = opt[:m]
    z0b = opt[m:]

    y = filtfilt_zi(x, b, a, z0f, z0b)
    return y, z0f, z0b


def check_filtfilt_gust(b, a, x, irlen=None, result=None):
    """
    A generator that yields tests for comparing the results of
    _filtfilt_gust and filtfilt_opt.  `x` must be a 1-d array.
    """
    if result is None:
        y = filtfilt(b, a, x, method="gust", irlen=irlen)
        # Also call the private function so we can test the ICs.
        yg, zg1, zg2 = _filtfilt_gust(b, a, x, irlen=irlen)
    else:
        # The caller provided the values to compare to the "expected" result.
        yg, zg1, zg2 = result
        y = None

    # filtfilt_opt gives the "expected" result.
    yo, zo1, zo2 = filtfilt_opt(b, a, x)

    if y is not None:
        assert_allclose(y, yo)
    assert_allclose(yg, yo)
    assert_allclose(zg1, zo1)
    assert_allclose(zg2, zo2)


def test_simple():
    # The input array has length 2.  The exact solution for this case
    # was computed "by hand".
    x = np.array([1.0, 2.0])
    b = np.array([0.5])
    a = np.array([1.0, -0.5])
    y, z1, z2 = _filtfilt_gust(b, a, x)
    assert_allclose([z1[0], z2[0]],
                    [0.3*x[0] + 0.2*x[1], 0.2*x[0] + 0.3*x[1]])
    assert_allclose(y, [z1[0] + 0.25*z2[0] + 0.25*x[0] + 0.125*x[1],
                        0.25*z1[0] + z2[0] + 0.125*x[0] + 0.25*x[1]])


def test_scalars():
    x = np.arange(12)
    b = 3.0
    a = 2.0
    y = filtfilt(b, a, x, method="gust")
    expected = (b/a)**2 * x
    assert_allclose(y, expected)


def test_basic():
    b, a = ellip(3, 0.01, 120, 0.0875)
    z, p, k = tf2zpk(b, a)
    eps = 1e-10
    r = np.max(np.abs(p))
    # Approximate impulse response length.
    approx_impulse_len = int(np.ceil(np.log(eps) / np.log(r)))

    np.random.seed(123)

    for irlen in [None, approx_impulse_len]:
        signal_len = 5 * approx_impulse_len

        # 1-d test case
        x = np.random.randn(signal_len)
        yield check_filtfilt_gust, b, a, x, irlen

        # 3-d test case using irlen; test each axis.
        x = np.random.randn(4 * signal_len)
        for axis in range(3):
            shape = [2, 2, 2]
            shape[axis] = signal_len
            x = x.reshape(*shape)
            yg, zg1, zg2 = _filtfilt_gust(b, a, x, axis=axis, irlen=irlen)
            for i in range(2):
                for j in range(2):
                    slc = [i, j]
                    slc.insert(axis, slice(None, None))
                    result = (yg[slc], zg1[slc], zg2[slc])
                    yield (check_filtfilt_gust, b, a, x[slc],
                           irlen, result)
            # The above code tested the private function _filtfilt_gust.
            # Because the `result` argument was passed to check_filtfilt_gust,
            # the public function was not tested, so we do that here.
            y = filtfilt(b, a, x, method="gust", axis=axis, irlen=irlen)
            # Since the public function calls the private function,
            # they should return exactly the same array (although there is
            # some risk of non-determinacy with floating point calculations).
            yield assert_array_equal, y, yg

    # Test case with length less than 2*approx_impulse_len.
    # In this case, `filtfilt_gust` should behave the same as if
    # `irlen=None` was given.
    length = 2*approx_impulse_len - 50
    x = np.random.randn(length)
    yield check_filtfilt_gust, b, a, x, approx_impulse_len


if __name__ == "__main__":
    run_module_suite()
