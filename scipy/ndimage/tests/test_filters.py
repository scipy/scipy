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
