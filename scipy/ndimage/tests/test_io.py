from numpy.testing import *
import scipy.ndimage as ndi

import os

def test_imread():
    lp = os.path.join(os.path.dirname(__file__), 'dots.png')
    img = ndi.imread(lp)
    assert_array_equal(img.shape, (300, 420, 3))

    img = ndi.imread(lp, flatten=True)
    assert_array_equal(img.shape, (300, 420))

if __name__ == "__main__":
    run_module_suite()
