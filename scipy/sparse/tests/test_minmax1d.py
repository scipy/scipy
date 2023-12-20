"""Test of min-max 1D features of sparse array classes"""

import pytest

import numpy as np

from numpy.testing import assert_equal, assert_array_equal

from scipy.sparse import coo_array  # , dok_array
from scipy.sparse._sputils import isscalarlike


def toarray(a):
    if isinstance(a, np.ndarray) or isscalarlike(a):
        return a
    return a.toarray()


formats_for_minmax = [coo_array]  # , dok_array]


@pytest.mark.parametrize("spcreator", formats_for_minmax)
class Test_MinMaxMixin1D:
    def test_minmax(self, spcreator):
        D = np.arange(5)
        X = spcreator(D)

        assert_equal(X.min(), 0)
        assert_equal(X.max(), 4)
        assert_equal((-X).min(), -4)
        assert_equal((-X).max(), 0)


    def test_minmax_axis(self, spcreator):
        D = np.arange(50)
        X = spcreator(D)

        for axis in [0, -1]:
            assert_array_equal(
                toarray(X.max(axis=axis)), D.max(axis=axis, keepdims=True)
            )
            assert_array_equal(
                toarray(X.min(axis=axis)), D.min(axis=axis, keepdims=True)
            )
        for axis in [-2, 1]:
            pytest.raises(ValueError, X.min, axis=axis)
            pytest.raises(ValueError, X.max, axis=axis)


    def test_numpy_minmax(self, spcreator):
        dat = np.array([0, 1, 2])
        datsp = spcreator(dat)
        assert_array_equal(np.min(datsp), np.min(dat))
        assert_array_equal(np.max(datsp), np.max(dat))


    def test_argmax(self, spcreator):
        D1 = np.array([-1, 5, 2, 3])
        D2 = np.array([0, 0, -1, -2])
        D3 = np.array([-1, -2, -3, -4])
        D4 = np.array([1, 2, 3, 4])
        D5 = np.array([1, 2, 0, 0])

        for D in [D1, D2, D3, D4, D5]:
            mat = spcreator(D)

            assert_equal(mat.argmax(), np.argmax(D))
            assert_equal(mat.argmin(), np.argmin(D))

            assert_equal(mat.argmax(axis=0), np.argmax(D, axis=0))
            assert_equal(mat.argmin(axis=0), np.argmin(D, axis=0))

        D6 = np.empty((0,))

        for axis in [None, 0]:
            mat = spcreator(D6)
            pytest.raises(ValueError, mat.argmax, axis=axis)
            pytest.raises(ValueError, mat.argmin, axis=axis)
