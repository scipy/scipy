import numpy as np
from numpy.testing import assert_allclose
import pytest

from scipy.stats import f_ishigami, sobol_indices


@pytest.fixture(scope='session')
def ishigami_ref_indices():
    a = 7.
    b = 0.1

    var = 0.5 + a ** 2 / 8 + b * np.pi ** 4 / 5\
        + b ** 2 * np.pi ** 8 / 18
    v1 = 0.5 + b * np.pi ** 4 / 5 + b ** 2 * np.pi ** 8 / 50
    v2 = a ** 2 / 8
    v3 = 0
    v12 = 0
    v13 = b ** 2 * np.pi ** 8 * 8 / 225
    v23 = 0

    s_first = np.array([v1 / var, v2 / var, v3 / var])
    s_second = np.array([
        [0., 0., v13 / var],
        [v12 / var, 0., v23 / var],
        [v13 / var, v23 / var, 0.]
    ])
    s_total = s_first + s_second.sum(axis=1)

    return s_first.reshape(-1, 1), s_total.reshape(-1, 1)


class TestSobolIndices:

    def test_ishigami(self, ishigami_ref_indices):
        indices = sobol_indices(
            func=f_ishigami, n=4096, d=3,
            l_bounds=[-np.pi, -np.pi, -np.pi], u_bounds=[np.pi, np.pi, np.pi]
        )

        assert_allclose(indices[0], ishigami_ref_indices[0], atol=1e-2)
        assert_allclose(indices[1], ishigami_ref_indices[1], atol=1e-2)
