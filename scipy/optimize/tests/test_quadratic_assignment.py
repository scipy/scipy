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
    seed = np.asarray([seed_cost, [opt_perm[z] for z in seed_cost]]).T
    np.random.seed(0)
    res = quadratic_assignment(cost_matrix, dist_matrix)

    # check ofv with barycenter initialization
    assert 11156 <= res['score'] < 21000
    assert res['score'] == _score(cost_matrix, dist_matrix, res['col_ind'])

    # check ofv with seeds
    res = quadratic_assignment(cost_matrix, dist_matrix,
                               options={'seed': seed})

    assert 11156 <= res['score'] < 21000

    # check performance when seeds are the global optimum
    seed = np.asarray([np.arange(n), [opt_perm[z] for z in np.arange(n)]]).T
    res = quadratic_assignment(cost_matrix, dist_matrix,
                               options={'seed': seed})

    assert 11156 == res['score']

    # check that max_iter is obeying with low input value
    iter = 5
    res = quadratic_assignment(cost_matrix, dist_matrix,
                               options={'maxiter': iter})

    assert iter >= res['nit']


@pytest.mark.slow
def test_rand_qap():
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

    # check ofv with 100 random initializations
    res = quadratic_assignment(
        cost_matrix, dist_matrix, options={'init': "rand", 'n_init': 100}
    )

    assert 11156 <= res['score'] < 14500


def test_quadratic_assignment_input_validation():
    A = np.identity(2)
    B = A

    # Type Error Checks: making sure single value parameters are of
    # correct type

    with pytest.raises(TypeError):
        quadratic_assignment(A, B, options={'n_init': 1.5})
    with pytest.raises(ValueError):
        quadratic_assignment(A, B, options={'n_init': -1})
    with pytest.raises(ValueError):
        quadratic_assignment(A, B, options={'init': "random"})
    with pytest.raises(TypeError):
        quadratic_assignment(A, B, options={'maxiter': 1.5})
    with pytest.raises(ValueError):
        quadratic_assignment(A, B, options={'maxiter': -1})
    with pytest.raises(TypeError):
        quadratic_assignment(A, B, options={'shuffle_input': "hey"})
    with pytest.raises(TypeError):
        quadratic_assignment(A, B, options={'eps': -1})
    with pytest.raises(TypeError):
        quadratic_assignment(A, B, maximize="hey")
    # test that cost and dist matrices of different sizes return error
    with pytest.raises(ValueError):
        quadratic_assignment(
            np.random.random((3, 3)),
            np.random.random((4, 4)),
        )
    # test that cost/dist matrices must be nonnegative
    with pytest.raises(ValueError):
        quadratic_assignment(
            -1 * np.random.random((3, 3)),
            np.random.random((3, 3))
        )
    # test that non square matrices return error
    with pytest.raises(ValueError):
        quadratic_assignment(
            np.random.random((3, 4)),
            np.random.random((3, 3)),
        )
    with pytest.raises(ValueError):
        quadratic_assignment(
            np.random.random((3, 3)),
            np.random.random((3, 4)),
        )
    # test for only two seed columns
    with pytest.raises(ValueError):
        quadratic_assignment(
            np.identity(3), np.identity(3),
            options={'seed': _range_matrix(2, 3)}
        )
    # can't have more seed nodes than cost/dist nodes
    with pytest.raises(ValueError):
        quadratic_assignment(np.identity(3), np.identity(3),
                             options={'seed': _range_matrix(5, 2)})
    # seeds cannot be negative valued
    with pytest.raises(ValueError):
        quadratic_assignment(
            np.identity(3), np.identity(3),
            options={'seed': -1 * _range_matrix(2, 2)}
        )
    # seeds can't have values greater than number of nodes
    with pytest.raises(ValueError):
        quadratic_assignment(
            np.identity(5), np.identity(5),
            options={'seed': 2 * _range_matrix(4, 2)}
        )
    # columns of seed matrix must be unique
    with pytest.raises(ValueError):
        quadratic_assignment(
            np.identity(3), np.identity(3),
            options={'seed': np.ones((2, 2))}
        )
    # test init matrix not doubly stochastic
    with pytest.raises(ValueError):
        quadratic_assignment(
            np.identity(3), np.identity(3), options={'init': np.ones((3, 3))}
        )


def _range_matrix(a, b):
    mat = np.zeros((a, b))
    for i in range(b):
        mat[:, i] = np.arange(a)
    return mat


def _score(A, B, col):
    A = np.asarray(A)
    B = np.asarray(B)
    col = np.asarray(col)
    return np.trace(np.transpose(A) @ B[np.ix_(col, col)])
