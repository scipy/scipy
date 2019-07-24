import numpy as np

import scipy.special as sc


class TestHyperu:

    def test_negative_x(self):
        a, b, x = np.meshgrid(
            [-1, -0.5, 0, 0.5, 1],
            [-1, -0.5, 0, 0.5, 1],
            np.linspace(-100, -1, 10),
        )
        assert np.all(np.isnan(sc.hyperu(a, b, x)))

    def test_special_cases(self):
        assert sc.hyperu(0, 1, 1) == 1.0
