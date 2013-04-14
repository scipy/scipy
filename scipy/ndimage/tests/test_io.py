from __future__ import division, print_function, absolute_import

from numpy.testing import assert_array_equal, dec, run_module_suite
import scipy.ndimage as ndi

import os

try:
    from PIL import Image
    pil_missing = False
except ImportError:
    pil_missing = True


@dec.skipif(pil_missing, msg="The Python Image Library could not be found.")
def test_imread():
    lp = os.path.join(os.path.dirname(__file__), 'dots.png')
    img = ndi.imread(lp, mode="RGB")
    assert_array_equal(img.shape, (300, 420, 3))

    img = ndi.imread(lp, flatten=True)
    assert_array_equal(img.shape, (300, 420))

if __name__ == "__main__":
    run_module_suite()
