import numpy as np
import scipy.special as sc
import pytest
from numpy.testing import assert_allclose


class Testbetainc:
    def test(self):
        val = sc.betainc(0.5, 1, 1)
        assert_allclose(val, 1)

    @pytest.mark.parametrize('a, b, x', [
        (np.inf, 2, 0.5),
        (1.0, np.inf, 0.5),
        (1.0, 2, np.inf),
        (-1.1, 1, 0.5)
    ])
    def test_inf(self, a, b, x):
        val = sc.betainc(a, b, x)
        assert np.isnan(val)

    @pytest.mark.parametrize('a, b, x, result', [
        (0.0, 3, 0.5, 1.0),
        (3.0, 0, 0.5, 1.0),
        (0.0, 0, 0.5, np.nan),
    ])
    def test_domain_edge(self, a, b, x, result):
        val = sc.betainc(a, b, x)
        assert_allclose(val, result)
