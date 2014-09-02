''' Some tests for filters '''
from __future__ import division, print_function, absolute_import

import numpy as np

from numpy.testing import assert_equal, assert_raises

import scipy.ndimage as sndi


def test_ticket_701():
    # Test generic filter sizes
    arr = np.arange(4).reshape((2,2))
    func = lambda x: np.min(x)
    res = sndi.generic_filter(arr, func, size=(1,1))
    # The following raises an error unless ticket 701 is fixed
    res2 = sndi.generic_filter(arr, func, size=1)
    assert_equal(res, res2)


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
    arr = np.zeros((100, 100), np.float)
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

def test_minimum_filter1d():
    # Regression gh-3898
    out = sndi.minimum_filter1d(np.arange(10), 1)
    assert_equal(np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]), out)
