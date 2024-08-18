import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_equal
from scipy.special import pow1p


class TestPow1p:

    @pytest.mark.parametrize('x, y', [
        (1.2, 3.4),        # "ordinary" value
        (-2.0, 3.0),       # (1+x) < 0
        (1e-37, 1e38),     # normal float32, but < 2e-16
        (-1e-45, 3.4e38),  # subnormal float32, and < 2e-16
    ])
    def test_float32(self, x, y):
        # float32 inputs should give the same answer as float64 inputs
        x32, y32 = np.float32(x), np.float32(y)
        x64, y64 = np.float64(x32), np.float64(y32)
        res32 = pow1p(x32, y32)
        assert isinstance(res32, np.float32)
        res64 = pow1p(x64, y64)
        assert isinstance(res64, np.float64)
        assert_equal(res32, np.float32(res64))
