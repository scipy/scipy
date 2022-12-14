import numpy as np
from numpy.testing import assert_allclose
import pytest

from scipy.stats import f_ishigami, sobol_indices, uniform


@pytest.fixture(scope='session')
def ishigami_ref_indices():
    """Reference values for Ishigami from Saltelli2007.

    Chapter 4, exercise 5 pages 179-182.
    """
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

    return [s_first.reshape(1, -1), s_total.reshape(1, -1)]


def f_ishigami_vec(x):
    res = f_ishigami(x)
    return np.column_stack([res, res])


class TestSobolIndices:

    @pytest.mark.parametrize(
        'func',
        [f_ishigami, f_ishigami_vec],
        ids=['scalar', 'vector']
    )
    def test_ishigami(self, ishigami_ref_indices, func):
        rng = np.random.default_rng(28631265345463262246170309650372465332)
        res = sobol_indices(
            func=func, n=4096,
            dists=[
                uniform(loc=-np.pi, scale=2*np.pi),
                uniform(loc=-np.pi, scale=2*np.pi),
                uniform(loc=-np.pi, scale=2*np.pi)
            ],
            random_state=rng
        )

        if func.__name__ == 'f_ishigami_vec':
            ishigami_ref_indices[0] = np.concatenate(
                [ishigami_ref_indices[0], ishigami_ref_indices[0]]
            )
            ishigami_ref_indices[1] = np.concatenate(
                [ishigami_ref_indices[1], ishigami_ref_indices[1]]
            )

        assert_allclose(res.first_order, ishigami_ref_indices[0], atol=1e-2)
        assert_allclose(res.total_order, ishigami_ref_indices[1], atol=1e-2)

    def test_raises(self):

        message = r"Either of the methods `pdf` or `logpdf` must be specified"
        with pytest.raises(ValueError, match=message):
            sobol_indices(n=0, func=f_ishigami, dists="uniform")

        with pytest.raises(ValueError, match=message):
            sobol_indices(n=0, func=f_ishigami, dists=[lambda x: x])

        message = r"The balance properties of Sobol'"
        with pytest.raises(ValueError, match=message):
            sobol_indices(n=7, func=f_ishigami, dists=[uniform()])

        def func_wrong_shape_output(x):
            return x.reshape(1, -1)

        message = r"'func' output should have a shape"
        with pytest.raises(ValueError, match=message):
            sobol_indices(n=2, func=func_wrong_shape_output, dists=[uniform()])
