import numpy as np
from numpy.testing import assert_allclose
from scipy.misc import derest
import os


class TestDerest(object):

    def test_basic(self):
        file_path = os.path.join(os.path.dirname(__file__), "data/derestData.npz")
        data = np.load(file_path)

        x = data['x']
        outtest = data['out']
        tbasic = data['tBasic']

        out = derest.derest(tbasic, x)
        assert_allclose(out, outtest)

    def test_dense_small(self):
        file_path = os.path.join(os.path.dirname(__file__), "data/derestData.npz")
        data = np.load(file_path)

        x = data['x']
        tbasic = data['tBasic']
        tdensesmall = data['tDenseSmall']
        densesmalltest = data['denseSmall']

        __, out = derest(tbasic, x, tdensesmall)
        assert_allclose(out, densesmalltest)

    def test_dense_large(self):
        file_path = os.path.join(os.path.dirname(__file__), "data/derestData.npz")
        data = np.load(file_path)

        x = data['x']
        tbasic = data['tBasic']
        tdenselarge = data['tDenseLarge']
        denselargetest = data['denseLarge']

        __, out = derest(tbasic, x, tdenselarge)
        assert_allclose(out, denselargetest)
