''' Some tests for filters '''

import numpy as np

from numpy.testing import assert_equal, assert_raises

from nose.tools import assert_true

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
