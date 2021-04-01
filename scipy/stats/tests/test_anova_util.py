
import numpy as np
from numpy.testing import assert_array_equal, assert_allclose
from scipy.stats._anova_util import (_unpack_2d_data_grid, _encode,
                                     _encode_nested, _interaction, _proj)


def test_unpack_2d_data_grid1():
    # data has shape (2, 3, 2)
    data = np.array([[[10, 11], [12, 13], [14, 15]],
                     [[16, 17], [18, 19], [20, 21]]])
    a, b, values = _unpack_2d_data_grid(data)
    assert_array_equal(a, [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1])
    assert_array_equal(b, [0, 0, 1, 1, 2, 2, 0, 0, 1, 1, 2, 2])
    assert_array_equal(values, data.ravel())


def test_unpack_2d_data_grid2():
    # data is a nested list.  The lengths of the innermost lists are
    # not all the same.
    data = [[[10], [11, 12, 13], [14, 15]],
            [[16, 17, 18], [19, 20], [20, 21, 22]]]
    a, b, values = _unpack_2d_data_grid(data)
    assert_array_equal(a, [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1])
    assert_array_equal(b, [0, 1, 1, 1, 2, 2, 0, 0, 0, 1, 1, 2, 2, 2])
    assert_array_equal(values, sum(sum(data, []), []))


def test_encode1():
    a = [10, 12, 10, 11, 12]
    v, u, inv = _encode(a)
    assert_array_equal(u, [10, 11, 12])
    assert_array_equal(inv, [0, 2, 0, 1, 2])
    # There are three levels; with the "sum to zero" contrasts, the
    # encoding for each level should be
    #   10:  [-1, -1]
    #   11:  [ 1,  0]
    #   12:  [ 0,  1]
    assert_array_equal(v, [[-1, -1], [0, 1], [-1, -1], [1, 0], [0, 1]])


def test_encode2():
    a = ['a', 'b', 'b', 'a', 'a']
    v, u, inv = _encode(a)
    assert_array_equal(u, ['a', 'b'])
    assert_array_equal(inv, [0, 1, 1, 0, 0])
    # There are two levels; with the "sum to zero" contrasts, the
    # encoding for each level should be
    #   'a':  [-1]
    #   'b':  [ 1]
    assert_array_equal(v, [[-1], [1], [1], [-1], [-1]])


def test_encode_nested1():
    a = np.array([100, 100, 100, 101, 101, 102, 102, 102, 102, 100, 101])
    b = np.array(['x', 'y', 'x', 'x', 'y', 'x', 'y', 'x', 'y', 'y', 'y'])
    # Note that we're nesting `b` in `a`, so the levels 'x' and 'y' for
    # each distinct value in `a` are from different categories.  That is,
    # the 'x' in the group corresponding to `a` = 100 is not the same
    # level as the 'x' in the group corresponding to `a` = 101.
    v, sizes = _encode_nested(a, b)
    assert len(sizes) == 3
    assert_array_equal(sizes[0], [2, 2])
    assert_array_equal(sizes[1], [1, 2])
    assert_array_equal(sizes[2], [2, 2])
    # The number of levels within each nested group is just 2, so the
    # encoding for the values in each nested group is -1 or 1.
    assert_array_equal(v, [[-1,  0,  0],
                           [ 1,  0,  0],
                           [-1,  0,  0],
                           [ 0, -1,  0],
                           [ 0,  1,  0],
                           [ 0,  0, -1],
                           [ 0,  0,  1],
                           [ 0,  0, -1],
                           [ 0,  0,  1],
                           [ 1,  0,  0],
                           [ 0,  1,  0]])


def test_encode_nested2():
    a = np.array([2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2])
    b = np.array([0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 0, 0, 0, 1, 1, 2])
    #
    #        b  sum-to-zero encoding
    # -----  -  --------------------
    # a=1    0    [-1, -1]
    #        1    [ 1,  0]
    #        2    [ 0,  1]
    # a=2    0    [-1, -1]
    #        1    [ 1,  0]
    #        2    [ 0,  1]
    v, sizes = _encode_nested(a, b)
    assert len(sizes) == 2
    assert_array_equal(sizes[0], [3, 3, 2])
    assert_array_equal(sizes[1], [4, 3, 3])
    assert_array_equal(v, [[ 0,  0, -1, -1],
                           [ 0,  0, -1, -1],
                           [-1, -1,  0,  0],
                           [-1, -1,  0,  0],
                           [ 0,  0,  1,  0],
                           [ 0,  0,  1,  0],
                           [ 1,  0,  0,  0],
                           [ 1,  0,  0,  0],
                           [ 0,  0,  0,  1],
                           [ 0,  0,  0,  1],
                           [ 0,  1,  0,  0],
                           [ 0,  1,  0,  0],
                           [ 0,  0, -1, -1],
                           [ 0,  0, -1, -1],
                           [-1, -1,  0,  0],
                           [ 1,  0,  0,  0],
                           [ 0,  0,  1,  0],
                           [ 0,  0,  0,  1]])


def test_encode_nested3():
    # In this test, the number of levels in the nested subgroups is not
    # the same for each subgroup.
    a = np.array([1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3])
    b = np.array([0, 1, 2, 0, 1, 2, 0, 0, 0, 1, 1, 1, 0, 1, 2, 3, 0, 1, 2, 3])
    #  a level   distinct b levels
    #  -------   -----------------
    #      1         0 1 2
    #      2         0 1
    #      3         0 1 2 3
    #
    # Encoding the sublevels for each level of a:
    #
    #        b  sum-to-zero encoding
    # -----  -  --------------------
    # a=1    0     [-1, -1]
    #        1     [ 1,  0]
    #        2     [ 0,  1]
    # a=2    0     [-1]
    #        1     [ 1]
    # a=3    0     [-1, -1, -1]
    #        1     [ 1,  0,  0]
    #        2     [ 0,  1,  0]
    #        3     [ 0,  0,  1]
    v, sizes = _encode_nested(a, b)
    assert len(sizes) == 3
    assert_array_equal(sizes[0], [2, 2, 2])
    assert_array_equal(sizes[1], [3, 3])
    assert_array_equal(sizes[2], [2, 2, 2, 2])
    # The first two columns of v hold the encoding for a == 1.
    # The third column is the encoding for a == 2.
    # The last three columns hold the encoding for a == 3.
    assert_array_equal(v, [[-1, -1,  0,  0,  0,  0],
                           [ 1,  0,  0,  0,  0,  0],
                           [ 0,  1,  0,  0,  0,  0],
                           [-1, -1,  0,  0,  0,  0],
                           [ 1,  0,  0,  0,  0,  0],
                           [ 0,  1,  0,  0,  0,  0],
                           [ 0,  0, -1,  0,  0,  0],
                           [ 0,  0, -1,  0,  0,  0],
                           [ 0,  0, -1,  0,  0,  0],
                           [ 0,  0,  1,  0,  0,  0],
                           [ 0,  0,  1,  0,  0,  0],
                           [ 0,  0,  1,  0,  0,  0],
                           [ 0,  0,  0, -1, -1, -1],
                           [ 0,  0,  0,  1,  0,  0],
                           [ 0,  0,  0,  0,  1,  0],
                           [ 0,  0,  0,  0,  0,  1],
                           [ 0,  0,  0, -1, -1, -1],
                           [ 0,  0,  0,  1,  0,  0],
                           [ 0,  0,  0,  0,  1,  0],
                           [ 0,  0,  0,  0,  0,  1]])


def test_interaction1():
    X1 = np.array([[-1, -1],
                   [ 0,  1],
                   [ 0,  1],
                   [ 1,  0]])
    X2 = np.array([[-1],
                   [-1],
                   [ 1],
                   [ 1]])
    X = _interaction(X1, X2)
    assert_array_equal(X, np.column_stack((X1[:, 0] * X2[:, 0],
                                           X1[:, 1] * X2[:, 0])))


def test_interaction2():
    X1 = np.array([[-1, -1],
                   [ 0,  1],
                   [ 0,  1],
                   [ 1,  0]])
    X2 = np.array([[-1, -1, -1],
                   [-1,  0,  1],
                   [ 1, -1, -1],
                   [ 1,  0,  1]])
    X = _interaction(X1, X2)
    assert_array_equal(X, np.column_stack((X1[:, 0] * X2[:, 0],
                                           X1[:, 0] * X2[:, 1],
                                           X1[:, 0] * X2[:, 2],
                                           X1[:, 1] * X2[:, 0],
                                           X1[:, 1] * X2[:, 1],
                                           X1[:, 1] * X2[:, 2])))


def test_proj():
    X = np.array([[1, 0],
                  [1, 1],
                  [2, 1]])
    P = _proj(X)
    assert_allclose(P @ X, X, rtol=1e-12, atol=1e-14)
    # Xperp is perpendicular to the columns of X.
    Xperp = np.array([1, 1, -1])
    assert_allclose(P @ Xperp, np.zeros(3), atol=1e-14)
