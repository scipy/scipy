''' Some tests for filters '''
from __future__ import division, print_function, absolute_import

import sys
import numpy as np

from numpy.testing import (assert_equal, assert_raises, assert_allclose,
                           assert_array_equal, TestCase, run_module_suite)

import scipy.ndimage as sndi


def test_ticket_701():
    # Test generic filter sizes
    arr = np.arange(4).reshape((2,2))
    func = lambda x: np.min(x)
    res = sndi.generic_filter(arr, func, size=(1,1))
    # The following raises an error unless ticket 701 is fixed
    res2 = sndi.generic_filter(arr, func, size=1)
    assert_equal(res, res2)


def test_gh_5430():
    # At least one of these raises an error unless gh-5430 is
    # fixed. In py2k an int is implemented using a C long, so
    # which one fails depends on your system. In py3k there is only
    # one arbitrary precision integer type, so both should fail.
    x = np.zeros(1)
    sigma = np.int32(1)
    y = sndi.gaussian_filter(x, sigma)
    assert_allclose(x, y)
    sigma = np.int64(1)
    y = sndi.gaussian_filter(x, sigma)
    assert_allclose(x, y)
    sigma = 1
    # This worked before; make sure it still works
    y = sndi.gaussian_filter(x, sigma)
    assert_allclose(x, y)
    # This worked before; make sure it still works
    x = np.zeros((2, 2))
    sigma = [1, 1]
    y = sndi.gaussian_filter(x, sigma)
    assert_allclose(x, y)


def test_orders_gauss():
    # Check order inputs to Gaussians
    arr = np.zeros((1,))
    yield assert_equal, 0, sndi.gaussian_filter(arr, 1, order=0)
    yield assert_equal, 0, sndi.gaussian_filter(arr, 1, order=3)
    yield assert_raises, ValueError, sndi.gaussian_filter, arr, 1, -1
    yield assert_raises, ValueError, sndi.gaussian_filter, arr, 1, 4
    yield assert_equal, 0, sndi.gaussian_filter1d(arr, 1, axis=-1, order=0)
    yield assert_equal, 0, sndi.gaussian_filter1d(arr, 1, axis=-1, order=3)
    yield assert_raises, ValueError, sndi.gaussian_filter1d, arr, 1, -1, -1
    yield assert_raises, ValueError, sndi.gaussian_filter1d, arr, 1, -1, 4


def test_valid_origins():
    """Regression test for #1311."""
    func = lambda x: np.mean(x)
    data = np.array([1,2,3,4,5], dtype=np.float64)
    assert_raises(ValueError, sndi.generic_filter, data, func, size=3,
                  origin=2)
    func2 = lambda x, y: np.mean(x + y)
    assert_raises(ValueError, sndi.generic_filter1d, data, func,
                  filter_size=3, origin=2)
    assert_raises(ValueError, sndi.percentile_filter, data, 0.2, size=3,
                  origin=2)

    for filter in [sndi.uniform_filter, sndi.minimum_filter,
                   sndi.maximum_filter, sndi.maximum_filter1d,
                   sndi.median_filter, sndi.minimum_filter1d]:
        # This should work, since for size == 3, the valid range for origin is
        # -1 to 1.
        list(filter(data, 3, origin=-1))
        list(filter(data, 3, origin=1))
        # Just check this raises an error instead of silently accepting or
        # segfaulting.
        assert_raises(ValueError, filter, data, 3, origin=2)


def test_gaussian_truncate():
    # Test that Gaussian filters can be truncated at different widths.
    # These tests only check that the result has the expected number
    # of nonzero elements.
    arr = np.zeros((100, 100), float)
    arr[50, 50] = 1
    num_nonzeros_2 = (sndi.gaussian_filter(arr, 5, truncate=2) > 0).sum()
    assert_equal(num_nonzeros_2, 21**2)
    num_nonzeros_5 = (sndi.gaussian_filter(arr, 5, truncate=5) > 0).sum()
    assert_equal(num_nonzeros_5, 51**2)

    # Test truncate when sigma is a sequence.
    f = sndi.gaussian_filter(arr, [0.5, 2.5], truncate=3.5)
    fpos = f > 0
    n0 = fpos.any(axis=0).sum()
    # n0 should be 2*int(2.5*3.5 + 0.5) + 1
    assert_equal(n0, 19)
    n1 = fpos.any(axis=1).sum()
    # n1 should be 2*int(0.5*3.5 + 0.5) + 1
    assert_equal(n1, 5)

    # Test gaussian_filter1d.
    x = np.zeros(51)
    x[25] = 1
    f = sndi.gaussian_filter1d(x, sigma=2, truncate=3.5)
    n = (f > 0).sum()
    assert_equal(n, 15)

    # Test gaussian_laplace
    y = sndi.gaussian_laplace(x, sigma=2, truncate=3.5)
    nonzero_indices = np.where(y != 0)[0]
    n = nonzero_indices.ptp() + 1
    assert_equal(n, 15)

    # Test gaussian_gradient_magnitude
    y = sndi.gaussian_gradient_magnitude(x, sigma=2, truncate=3.5)
    nonzero_indices = np.where(y != 0)[0]
    n = nonzero_indices.ptp() + 1
    assert_equal(n, 15)


class TestThreading(TestCase):
    def check_func_thread(self, n, fun, args, out):
        from threading import Thread
        thrds = [Thread(target=fun, args=args, kwargs={'output': out[x]}) for x in range(n)]
        [t.start() for t in thrds]
        [t.join() for t in thrds]

    def check_func_serial(self, n, fun, args, out):
        for i in range(n):
            fun(*args, output=out[i])

    def test_correlate1d(self):
        d = np.random.randn(5000)
        os = np.empty((4, d.size))
        ot = np.empty_like(os)
        self.check_func_serial(4, sndi.correlate1d, (d, np.arange(5)), os)
        self.check_func_thread(4, sndi.correlate1d, (d, np.arange(5)), ot)
        assert_array_equal(os, ot)

    def test_correlate(self):
        d = np.random.randn(500, 500)
        k = np.random.randn(10, 10)
        os = np.empty([4] + list(d.shape))
        ot = np.empty_like(os)
        self.check_func_serial(4, sndi.correlate, (d, k), os)
        self.check_func_thread(4, sndi.correlate, (d, k), ot)
        assert_array_equal(os, ot)

    def test_median_filter(self):
        d = np.random.randn(500, 500)
        os = np.empty([4] + list(d.shape))
        ot = np.empty_like(os)
        self.check_func_serial(4, sndi.median_filter, (d, 3), os)
        self.check_func_thread(4, sndi.median_filter, (d, 3), ot)
        assert_array_equal(os, ot)

    def test_uniform_filter1d(self):
        d = np.random.randn(5000)
        os = np.empty((4, d.size))
        ot = np.empty_like(os)
        self.check_func_serial(4, sndi.uniform_filter1d, (d, 5), os)
        self.check_func_thread(4, sndi.uniform_filter1d, (d, 5), ot)
        assert_array_equal(os, ot)

    def test_minmax_filter(self):
        d = np.random.randn(500, 500)
        os = np.empty([4] + list(d.shape))
        ot = np.empty_like(os)
        self.check_func_serial(4, sndi.maximum_filter, (d, 3), os)
        self.check_func_thread(4, sndi.maximum_filter, (d, 3), ot)
        assert_array_equal(os, ot)
        self.check_func_serial(4, sndi.minimum_filter, (d, 3), os)
        self.check_func_thread(4, sndi.minimum_filter, (d, 3), ot)
        assert_array_equal(os, ot)


def test_minmaximum_filter1d():
    # Regression gh-3898
    in_ = np.arange(10)
    out = sndi.minimum_filter1d(in_, 1)
    assert_equal(in_, out)
    out = sndi.maximum_filter1d(in_, 1)
    assert_equal(in_, out)
    # Test reflect
    out = sndi.minimum_filter1d(in_, 5, mode='reflect')
    assert_equal([0, 0, 0, 1, 2, 3, 4, 5, 6, 7], out)
    out = sndi.maximum_filter1d(in_, 5, mode='reflect')
    assert_equal([2, 3, 4, 5, 6, 7, 8, 9, 9, 9], out)
    #Test constant
    out = sndi.minimum_filter1d(in_, 5, mode='constant', cval=-1)
    assert_equal([-1, -1, 0, 1, 2, 3, 4, 5, -1, -1], out)
    out = sndi.maximum_filter1d(in_, 5, mode='constant', cval=10)
    assert_equal([10, 10, 4, 5, 6, 7, 8, 9, 10, 10], out)
    # Test nearest
    out = sndi.minimum_filter1d(in_, 5, mode='nearest')
    assert_equal([0, 0, 0, 1, 2, 3, 4, 5, 6, 7], out)
    out = sndi.maximum_filter1d(in_, 5, mode='nearest')
    assert_equal([2, 3, 4, 5, 6, 7, 8, 9, 9, 9], out)
    # Test wrap
    out = sndi.minimum_filter1d(in_, 5, mode='wrap')
    assert_equal([0, 0, 0, 1, 2, 3, 4, 5, 0, 0], out)
    out = sndi.maximum_filter1d(in_, 5, mode='wrap')
    assert_equal([9, 9, 4, 5, 6, 7, 8, 9, 9, 9], out)


if __name__ == "__main__":
    run_module_suite(argv=sys.argv)
