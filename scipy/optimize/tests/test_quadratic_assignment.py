import pytest
import numpy as np
from numpy.testing import assert_array_equal
from scipy.optimize import quadratic_assignment


@pytest.mark.parametrize("n", [10, 100])
def test_weighted_graph_isomorphism(n):
    rng = np.random.RandomState(seed=0)
    A = np.abs(rng.randn(n, n))
    A = A + A.T

    P = rng.permutation(n)
    B = A[P][:, P]

    result = quadratic_assignment(A, B, maximize=True)
    assert_array_equal(result.col_ind, P)


class TestQAP():

    def setup_method(self):
        self.A = np.array([[0, 5, 8, 6],
                           [5, 0, 5, 1],
                           [8, 5, 0, 2],
                           [6, 1, 2, 0]])

        self.B = np.array([[0, 1, 8, 4],
                           [1, 0, 5, 2],
                           [8, 5, 0, 5],
                           [4, 2, 5, 0]])

    # Example in section IIIB of Umeyama: "Weighted graph matching problems".
    @pytest.mark.parametrize("maximize,permutation,score",
                             [(False, [3, 2, 1, 0], 200),
                              (True, [2, 3, 0, 1], 286)])
    def test_weighted_graph_matching(self, maximize, permutation, score):
        result = quadratic_assignment(self.A, self.B, maximize=maximize)
        assert_array_equal(result.col_ind, permutation)
        assert result.score == score

        n = len(self.A)
        P = np.zeros((n, n), dtype=int)
        P[np.arange(n), result.col_ind] = 1
        assert score == np.trace(self.A @ P @ self.B.T @ P.T)

    def test_square_A(self):
        with pytest.raises(ValueError, match="must be a square matrix"):
            quadratic_assignment(self.A[0], self.B)

    def test_square_B(self):
        with pytest.raises(ValueError, match="must be a square matrix"):
            quadratic_assignment(self.A, self.B[0])

    def test_equal_size(self):
        with pytest.raises(ValueError, match="must be of equal size"):
            quadratic_assignment(self.A, self.B[:-1, :-1])

    def test_negative_A(self):
        with pytest.raises(ValueError, match="contains negative entries"):
            quadratic_assignment(-self.A, self.B)

    def test_negative_B(self):
        with pytest.raises(ValueError, match="contains negative entries"):
            quadratic_assignment(self.A, -self.B)
