import numpy as np
from numpy.testing import assert_allclose
from scipy.misc import derest
import os


class TestDerest(object):

    file_path = os.path.join(os.path.dirname(__file__), "data/derestData.npz")
    data = np.load(file_path)

    x = data['x']
    outtest = data['out']
    tBasic = data['tBasic']
    tDenseSmall = data['tDenseSmall']
    tDenseLarge = data['tDenseLarge']
    denseSmalltest = data['denseSmall']
    denseLargetest = data['denseLarge']

    def test_basic(self):
        out = derest.derest(tBasic, x)
        assert_allclose(out, outtest)

    def test_dense_small(self):
        __, out = derest(tBasic, x, tDenseSmall)
        assert_allclose(out, denseSmalltest)

    def test_dense_large(self):
        __, out = derest(tBasic, x, tDenseLarge)
        assert_allclose(out, denseLargetest)
