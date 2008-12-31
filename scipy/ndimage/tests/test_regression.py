import numpy as np
from numpy.testing import *

import scipy.ndimage as ndimage

def test_byte_order_median():
    """Regression test for #413: median_filter does not handle bytes orders."""
    a = np.arange(9, dtype='<f4').reshape(3, 3)
    ref = ndimage.filters.median_filter(a,(3, 3))
    b = np.arange(9, dtype='>f4').reshape(3, 3)
    t = ndimage.filters.median_filter(b, (3, 3))
    assert_array_almost_equal(ref, t)

def test_zoom_output_shape():
    """Ticket #643"""
    x = np.arange(12).reshape((3,4))
    ndimage.zoom(x, 2, output=np.zeros((6,8)))

if __name__ == "__main__":
    run_module_suite()
