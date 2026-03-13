from pyprima.cobyla.cobyla import get_lincon
from pyprima.common.consts import BOUNDMAX
import numpy as np

def test_get_lincon():
    Aeq = np.array([[1, 2], [3, 4]])
    Aineq = np.array([[5, 6], [7, 8]])
    beq = np.array([9, 10])
    bineq = np.array([11, 12])
    xl = np.array([0, -1])
    xu = np.array([13, 14])
    amat, bvec = get_lincon(Aeq, Aineq, beq, bineq, xl, xu)
    assert np.allclose(amat, np.array([
        [-1, 0],
        [0, -1],
        [1, 0],
        [0, 1],
        [-1, -2],
        [-3, -4],
        [1, 2],
        [3, 4],
        [5, 6],
        [7, 8],
    ]))
    assert np.allclose(bvec, np.array([0, 1, 13, 14, -9, -10, 9, 10, 11, 12]))


def test_get_lincon_boundmax():
    Aeq = np.array([[1, 2], [3, 4]])
    Aineq = np.array([[5, 6], [7, 8]])
    beq = np.array([9, 10])
    bineq = np.array([11, 12])
    # Since the first element is below BOUNDMAX, we should ultimately
    # see only 3 bounds in the resultant matrix/vector.
    xl = np.array([-BOUNDMAX - 1, -1])
    xu = np.array([13, 14])
    amat, bvec = get_lincon(Aeq, Aineq, beq, bineq, xl, xu)
    assert np.allclose(amat, np.array([
        # Note that the first row is missing because the first element of xl is below BOUNDMAX.
        [0, -1],
        [1, 0],
        [0, 1],
        [-1, -2],
        [-3, -4],
        [1, 2],
        [3, 4],
        [5, 6],
        [7, 8],
    ]))
    assert np.allclose(bvec, np.array([1, 13, 14, -9, -10, 9, 10, 11, 12]))


def test_none():
    Aeq = None
    Aineq = None
    beq = None
    bineq = None
    xl = None
    xu = None
    amat, bvec = get_lincon(Aeq, Aineq, beq, bineq, xl, xu)
    assert amat is None
    assert bvec is None
