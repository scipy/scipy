import pytest
import numpy as np
from scipy.optimize import quadratic_assignment


def test_quadratic_assignment():
    # cost and distance matrices of QAPLIB instance chr12c

    cost_matrix = [
        [0, 90, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [90, 0, 0, 23, 0, 0, 0, 0, 0, 0, 0, 0],
        [10, 0, 0, 0, 43, 0, 0, 0, 0, 0, 0, 0],
        [0, 23, 0, 0, 0, 88, 0, 0, 0, 0, 0, 0],
        [0, 0, 43, 0, 0, 0, 26, 0, 0, 0, 0, 0],
        [0, 0, 0, 88, 0, 0, 0, 16, 0, 0, 0, 0],
        [0, 0, 0, 0, 26, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 16, 0, 0, 0, 96, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 29, 0],
        [0, 0, 0, 0, 0, 0, 0, 96, 0, 0, 0, 37],
        [0, 0, 0, 0, 0, 0, 0, 0, 29, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 37, 0, 0],
    ]
    dist_matrix = [
        [0, 36, 54, 26, 59, 72, 9, 34, 79, 17, 46, 95],
        [36, 0, 73, 35, 90, 58, 30, 78, 35, 44, 79, 36],
        [54, 73, 0, 21, 10, 97, 58, 66, 69, 61, 54, 63],
        [26, 35, 21, 0, 93, 12, 46, 40, 37, 48, 68, 85],
        [59, 90, 10, 93, 0, 64, 5, 29, 76, 16, 5, 76],
        [72, 58, 97, 12, 64, 0, 96, 55, 38, 54, 0, 34],
        [9, 30, 58, 46, 5, 96, 0, 83, 35, 11, 56, 37],
        [34, 78, 66, 40, 29, 55, 83, 0, 44, 12, 15, 80],
        [79, 35, 69, 37, 76, 38, 35, 44, 0, 64, 39, 33],
        [17, 44, 61, 48, 16, 54, 11, 12, 64, 0, 70, 86],
        [46, 79, 54, 68, 5, 0, 56, 15, 39, 70, 0, 18],
        [95, 36, 63, 85, 76, 34, 37, 80, 33, 86, 18, 0],
    ]
    n = 12
    opt_perm = np.array([7, 5, 1, 3, 10, 4, 8, 6, 9, 11, 2, 12]) - [1] * n
    seed_cost = [4, 8, 10]
    seed_dist = [opt_perm[z] for z in seed_cost]

    row, col = quadratic_assignment(cost_matrix, dist_matrix)

    assert 11156 <= _score(cost_matrix, dist_matrix, col) < 21000

    row, col = quadratic_assignment(
        cost_matrix, dist_matrix, init_method="rand", n_init=100
    )

    assert 11156 <= _score(cost_matrix, dist_matrix, col) < 13500

    row, col = quadratic_assignment(cost_matrix, dist_matrix, seed_cost, seed_dist)

    assert 11156 <= _score(cost_matrix, dist_matrix, col) < 21000


def test_linear_sum_assignment_input_validation():
    A = np.identity(2)
    B = A
    with pytest.raises(TypeError):
        quadratic_assignment(A, B, n_init=-1.5)
    with pytest.raises(ValueError):
        quadratic_assignment(A, B, init_method="random")
    with pytest.raises(TypeError):
        quadratic_assignment(A, B, max_iter=-1.5)
    with pytest.raises(TypeError):
        quadratic_assignment(A, B, shuffle_input="hey")
    with pytest.raises(TypeError):
        quadratic_assignment(A, B, eps=-1)
    with pytest.raises(TypeError):
        quadratic_assignment(A, B, gmp="hey")
    with pytest.raises(ValueError):
        quadratic_assignment(
            np.random.random((3, 3)),
            np.random.random((4, 4)),
            np.arange(2),
            np.arange(2),
        )
    with pytest.raises(ValueError):
        quadratic_assignment(
            np.random.random((3, 4)),
            np.random.random((3, 4)),
            np.arange(2),
            np.arange(2),
        )
    with pytest.raises(ValueError):
        quadratic_assignment(
            np.identity(3), np.identity(3), np.identity(3), np.arange(2)
        )
    with pytest.raises(ValueError):
        quadratic_assignment(np.identity(3), np.identity(3), np.arange(1), np.arange(2))
    with pytest.raises(ValueError):
        quadratic_assignment(np.identity(3), np.identity(3), np.arange(5), np.arange(5))
    with pytest.raises(ValueError):
        quadratic_assignment(
            np.identity(3), np.identity(3), -1 * np.arange(2), -1 * np.arange(2)
        )


def _score(A, B, col):
    A = np.asarray(A)
    B = np.asarray(B)
    col = np.asarray(col)
    return np.trace(np.transpose(A) @ B[np.ix_(col, col)])
