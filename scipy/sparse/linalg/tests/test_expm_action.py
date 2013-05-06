"""Test functions for the sparse.linalg._expm_action module
"""

from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import (TestCase, run_module_suite, assert_allclose,
        assert_, assert_equal, decorators)

import scipy.linalg
from scipy.sparse.linalg import _expm_action


def less_than_or_close(a, b):
    return np.allclose(a, b) or (a < b)


class TestExpmAction(TestCase):

    def test_theta_monotonicity(self):
        pairs = sorted(_expm_action._theta.items())
        for (m_a, theta_a), (m_b, theta_b) in zip(pairs[:-1], pairs[1:]):
            assert_(theta_a < theta_b)

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

    def test_expm_action(self):
        np.random.seed(1234)
        n = 40
        k = 3
        nsamples = 10
        for i in range(nsamples):
            A = scipy.linalg.inv(np.random.randn(n, n))
            B = np.random.randn(n, k)
            observed = _expm_action.expm_action(A, B)
            expected = np.dot(scipy.linalg.expm(A), B)
            assert_allclose(observed, expected)

    def test_matrix_vector_action(self):
        np.random.seed(1234)
        n = 40
        nsamples = 10
        for i in range(nsamples):
            A = scipy.linalg.inv(np.random.randn(n, n))
            v = np.random.randn(n)
            observed = _expm_action.expm_action(A, v)
            expected = np.dot(scipy.linalg.expm(A), v)
            assert_allclose(observed, expected)

    def test_scaled_expm_action(self):
        np.random.seed(1234)
        n = 40
        k = 3
        nsamples = 10
        for i in range(nsamples):
            for t in (0.2, 1.0, 1.5):
                A = scipy.linalg.inv(np.random.randn(n, n))
                B = np.random.randn(n, k)
                observed = _expm_action._expm_action_simple(A, B, t=t)
                expected = np.dot(scipy.linalg.expm(t*A), B)
                assert_allclose(observed, expected)

    def test_sparse_expm_action(self):
        n = 40
        k = 3
        nsamples = 10
        for i in range(nsamples):
            A = scipy.sparse.rand(n, n, density=0.05)
            B_sparse = scipy.sparse.rand(n, k, density=0.4)
            B_dense = np.array(B_sparse.todense())
            for B in (B_sparse, B_dense):
                observed = _expm_action.expm_action(A, B)
                expected = np.dot(scipy.linalg.expm(A), B)
                assert_allclose(observed, expected)


if __name__ == '__main__':
    run_module_suite()

