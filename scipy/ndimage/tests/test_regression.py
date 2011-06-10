import numpy as np
from numpy.testing import assert_array_almost_equal, run_module_suite

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

def test_ticket_742():
    def SE(img, thresh=.7, size=4):
        mask = img > thresh
        rank = len(mask.shape)
        la, co = ndimage.label(mask,
                               ndimage.generate_binary_structure(rank, rank))
        slices = ndimage.find_objects(la)

    if np.dtype(np.intp) != np.dtype('i'):
        shape=(3,1240,1240)
        a = np.random.rand(np.product(shape)).reshape(shape)
        # shouldn't crash
        SE(a)

if __name__ == "__main__":
    run_module_suite()
