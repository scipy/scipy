import numpy as np
from numpy.testing import assert_array_equal
from scipy.stats.contingency import crosstab


def test_crosstab_basic():
    a = [0, 0, 9, 9, 0, 0, 9]
    b = [2, 1, 3, 1, 2, 3, 3]
    expected_avals = [0, 9]
    expected_bvals = [1, 2, 3]
    expected_count = np.array([[1, 2, 1],
                               [1, 0, 2]])
    (avals, bvals), count = crosstab(a, b)
    assert_array_equal(avals, expected_avals)
    assert_array_equal(bvals, expected_bvals)
    assert_array_equal(count, expected_count)


def test_crosstab_basic_1d():
    # Verify that a single input sequence works as expected.
    x = [1, 2, 3, 1, 2, 3, 3]
    expected_xvals = [1, 2, 3]
    expected_count = np.array([2, 2, 3])
    (xvals,), count = crosstab(x)
    assert_array_equal(xvals, expected_xvals)
    assert_array_equal(count, expected_count)


def test_crosstab_basic_3d():
    # Verify the function for three input sequences.
    a = 'a'
    b = 'b'
    x = [0, 0, 9, 9, 0, 0, 9, 9]
    y = [a, a, a, a, b, b, b, a]
    z = [1, 2, 3, 1, 2, 3, 3, 1]
    expected_xvals = [0, 9]
    expected_yvals = [a, b]
    expected_zvals = [1, 2, 3]
    expected_count = np.array([[[1, 1, 0],
                                [0, 1, 1]],
                               [[2, 0, 1],
                                [0, 0, 1]]])
    (xvals, yvals, zvals), count = crosstab(x, y, z)
    assert_array_equal(xvals, expected_xvals)
    assert_array_equal(yvals, expected_yvals)
    assert_array_equal(zvals, expected_zvals)
    assert_array_equal(count, expected_count)


def test_crosstab_levels():
    a = [0, 0, 9, 9, 0, 0, 9]
    b = [1, 2, 3, 1, 2, 3, 3]
    expected_avals = [0, 9]
    expected_bvals = [0, 1, 2, 3]
    expected_count = np.array([[0, 1, 2, 1],
                               [0, 1, 0, 2]])
    (avals, bvals), count = crosstab(a, b, levels=[None, [0, 1, 2, 3]])
    assert_array_equal(avals, expected_avals)
    assert_array_equal(bvals, expected_bvals)
    assert_array_equal(count, expected_count)


def test_crosstab_extra_levels():
    # The pair of values (-1, 3) will be ignored, because we explicitly
    # request the counted `a` values to be [0, 9].
    a = [0, 0, 9, 9, 0, 0, 9, -1]
    b = [1, 2, 3, 1, 2, 3, 3, 3]
    expected_avals = [0, 9]
    expected_bvals = [0, 1, 2, 3]
    expected_count = np.array([[0, 1, 2, 1],
                               [0, 1, 0, 2]])
    (avals, bvals), count = crosstab(a, b, levels=[[0, 9], [0, 1, 2, 3]])
    assert_array_equal(avals, expected_avals)
    assert_array_equal(bvals, expected_bvals)
    assert_array_equal(count, expected_count)
