import pytest
import numpy as np
from numpy.testing import assert_array_equal
from scipy.optimize import quadratic_assignment


@pytest.mark.parametrize("n", [10, 100])
@pytest.mark.parametrize("directed", [False, True])
def test_weighted_graph_isomorphism(n, directed):
    rng = np.random.RandomState(seed=0)
    A = np.abs(rng.randn(n, n))
    if not directed:
        A = (A + A.T) / 2

    P = rng.permutation(n)
    B = A[P][:, P]

    result = quadratic_assignment(A, B, maximize=True)
    assert_array_equal(result.col_ind, P)


def test_directed_weighted_graph_matching():
    # Example in section IVB of Umeyama:
    # "An eigendecomposition approach to weighted graph matching"
    AG = np.array([[0, 3, 4, 2],
                  [0, 0, 1, 2],
                  [1, 0, 0, 1],
                  [0, 0, 1, 0]])

    AH = np.array([[0, 4, 2, 4],
                  [0, 0, 1, 0],
                  [0, 2, 0, 2],
                  [0, 1, 2, 0]])

    result = quadratic_assignment(AG, AH, maximize=True)
    assert_array_equal(result.col_ind, [0, 2, 3, 1])


class TestQAP():

    def setup_method(self):
        self.AG = np.array([[0, 5, 8, 6],
                           [5, 0, 5, 1],
                           [8, 5, 0, 2],
                           [6, 1, 2, 0]])

        self.AH = np.array([[0, 1, 8, 4],
                           [1, 0, 5, 2],
                           [8, 5, 0, 5],
                           [4, 2, 5, 0]])

    # Example in section IIIB of Umeyama:
    # "An eigendecomposition approach to weighted graph matching"
    @pytest.mark.parametrize("maximize,permutation,score",
                             [(False, [3, 2, 1, 0], 200),
                              (True, [2, 3, 0, 1], 286)])
    def test_weighted_graph_matching(self, maximize, permutation, score):
        result = quadratic_assignment(self.AG, self.AH, maximize=maximize)
        assert_array_equal(result.col_ind, permutation)
        assert result.score == score

        n = len(self.AG)
        P = np.zeros((n, n), dtype=int)
        P[np.arange(n), result.col_ind] = 1
        assert score == np.trace(self.AG @ P @ self.AH.T @ P.T)

    def test_square_A(self):
        with pytest.raises(ValueError, match="must be a square matrix"):
            quadratic_assignment(self.AG[0], self.AH)

    def test_square_B(self):
        with pytest.raises(ValueError, match="must be a square matrix"):
            quadratic_assignment(self.AG, self.AH[0])

    def test_equal_size(self):
        with pytest.raises(ValueError, match="must be of equal size"):
            quadratic_assignment(self.AG, self.AH[:-1, :-1])

    def test_negative_A(self):
        with pytest.raises(ValueError, match="contains negative entries"):
            quadratic_assignment(-self.AG, self.AH)

    def test_negative_B(self):
        with pytest.raises(ValueError, match="contains negative entries"):
            quadratic_assignment(self.AG, -self.AH)
