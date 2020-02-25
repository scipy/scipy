from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import (assert_almost_equal, assert_equal, assert_allclose,
                           assert_array_almost_equal, assert_)

from scipy.special import log_softmax

def test_log_softmax():
    assert_allclose(log_softmax([1000, 1]), np.array([0, -999]), rtol=1e-13)
    
    # Expected value computed using mpmath (with mpmath.mp.dps = 200) and then
    # converted to float.
    x = np.arange(4)
    expected = np.array([-3.4401896985611953,
                         -2.4401896985611953,
                         -1.4401896985611953,
                         -0.44018969856119533])

    assert_allclose(log_softmax(x), expected, rtol=1e-13)

    # Translation property.  If all the values are changed by the same amount,
    # the softmax result does not change.
    assert_allclose(log_softmax(x + 100), expected, rtol=1e-13)

    # When axis=None, softmax operates on the entire array, and preserves
    # the shape.
    assert_allclose(log_softmax(x.reshape(2, 2)), expected.reshape(2, 2),
                    rtol=1e-13)
    
    
def test_log_softmax_multi_axes():
    assert_allclose(log_softmax([[1000, 1], [1000, 1]], axis=0), 
                    np.log(0.5)*np.ones((2,2)), rtol=1e-13)
    assert_allclose(log_softmax([[1000, 1], [1000, 1]], axis=1), 
                    np.array([[0, -999],[0, -999]]), rtol=1e-13)
    
    # Expected value computed using mpmath (with mpmath.mp.dps = 200) and then
    # converted to float.
    x = np.arange(8).reshape(2, 4)
    expected = np.array([[-3.4401896985611953,
                         -2.4401896985611953,
                         -1.4401896985611953,
                         -0.44018969856119533],
                        [-3.4401896985611953,
                         -2.4401896985611953,
                         -1.4401896985611953,
                         -0.44018969856119533]])

    assert_allclose(log_softmax(x, axis=1), expected, rtol=1e-13)
    assert_allclose(log_softmax(x.T, axis=0), expected.T, rtol=1e-13)

    # 3-d input, with a tuple for the axis.
    x3d = x.reshape(2, 2, 2)
    assert_allclose(log_softmax(x3d, axis=(1, 2)), expected.reshape(2, 2, 2),
                    rtol=1e-13)
    
    
def test_log_softmax_scalar():
    assert_allclose(log_softmax(1.0), 0.0, rtol=1e-13)
