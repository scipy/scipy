
import numpy as np
from numpy.testing import (TestCase, run_module_suite, assert_allclose,
        assert_, assert_equal, decorators)

import scipy.linalg
from scipy.sparse.linalg import _expm_action


def less_than_or_close(a, b):
    return np.allclose(a, b) or (a < b)


class TestExpmAction(TestCase):

    def test_p_max_default(self):
        m_max = 55
        expected_p_max = 8
        observed_p_max = _expm_action._compute_p_max(m_max)
        assert_equal(observed_p_max, expected_p_max)

    def test_p_max_range(self):
        for m_max in range(1, 55+1):
            p_max = _expm_action._compute_p_max(m_max)
            assert_(p_max*(p_max - 1) <= m_max + 1)
            p_too_big = p_max + 1
            assert_(p_too_big*(p_too_big - 1) > m_max + 1)

    @decorators.slow
    def test_onenormest_matrix_power(self):
        np.random.seed(1234)
        n = 40
        nsamples = 10
        for i in range(nsamples):
            A = scipy.linalg.inv(np.random.randn(n, n))
            for p in range(4):
                if not p:
                    M = np.identity(n)
                else:
                    M = np.dot(M, A)
                estimated = _expm_action._onenormest_matrix_power(A, p)
                exact = np.linalg.norm(M, 1)
                assert_(less_than_or_close(estimated, exact))
                assert_(less_than_or_close(exact, 3*estimated))


if __name__ == '__main__':
    run_module_suite()

