import numpy as np
from numpy.testing import run_module_suite, assert_equal, assert_raises
from scipy.sparse import bsr_matrix, coo_matrix, csr_matrix, csc_matrix


def test_argmax():
    D1 = np.array([
        [-1, 5, 2, 3],
        [0, 0, -1, -2],
        [-1, -2, -3, -4],
        [1, 2, 3, 4],
        [1, 2, 0, 0],
    ])

    D2 = D1.transpose()

    classes = [bsr_matrix, coo_matrix, csr_matrix, csc_matrix]

    for D in [D1, D2]:
        for spmatrix in classes:
            argmin = np.argmin(D)
            argmax = np.argmax(D)
            argmin_0 = np.argmin(D, axis=0)
            argmax_0 = np.argmax(D, axis=0)
            argmax_1 = np.argmax(D, axis=1)
            argmin_1 = np.argmin(D, axis=1)

            mat = spmatrix(D)

            assert_equal(mat.argmax(), argmax)
            assert_equal(mat.argmin(), argmin)

            assert_equal(mat.argmax(axis=0), argmax_0)
            assert_equal(mat.argmin(axis=0), argmin_0)

            assert_equal(mat.argmax(axis=1), argmax_1)
            assert_equal(mat.argmin(axis=1), argmin_1)

    D1 = np.empty((0, 5))
    D2 = np.empty((5, 0))
    for spmatrix in classes:
        for axis in [None, 0]:
            mat = spmatrix(D1)
            assert_raises(ValueError, mat.argmax, axis=axis)
            assert_raises(ValueError, mat.argmin, axis=axis)

        for axis in [None, 1]:
            mat = spmatrix(D2)
            assert_raises(ValueError, mat.argmax, axis=axis)
            assert_raises(ValueError, mat.argmin, axis=axis)


if __name__ == '__main__':
    run_module_suite()
