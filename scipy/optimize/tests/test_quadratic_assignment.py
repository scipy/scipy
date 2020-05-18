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

    res = quadratic_assignment(A, B, maximize=True)
    assert_array_equal(res[0], np.arange(n))
    assert_array_equal(res[1], P)


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
    @pytest.mark.parametrize("maximize,expected",
                             [(False, [3, 2, 1, 0]), (True, [2, 3, 0, 1])])
    def test_weighted_graph_matching(self, maximize, expected):
        res = quadratic_assignment(self.A, self.B, maximize=maximize)
        assert_array_equal(res[0], np.arange(len(self.A)))
        assert_array_equal(res[1], expected)

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
