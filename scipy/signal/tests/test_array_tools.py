from numpy.testing import assert_array_equal
from pytest import raises as assert_raises
import pytest

from scipy.signal._arraytools import (axis_slice, axis_reverse,
     odd_ext, even_ext, const_ext, zero_ext)

from scipy.signal.tests._utils import pytest_enable_array_api


@pytest_enable_array_api
class TestArrayTools:
    def test_axis_slice(self, xp):
        a = xp.reshape(xp.arange(12), (3, 4))

        s = axis_slice(a, start=0, stop=1, axis=0)
        assert_array_equal(s, a[0:1, :])

        s = axis_slice(a, start=-1, axis=0)
        assert_array_equal(s, a[-1:, :])

        s = axis_slice(a, start=0, stop=1, axis=1)
        assert_array_equal(s, a[:, 0:1])

        s = axis_slice(a, start=-1, axis=1)
        assert_array_equal(s, a[:, -1:])

        s = axis_slice(a, start=0, step=2, axis=0)
        assert_array_equal(s, a[::2, :])

        s = axis_slice(a, start=0, step=2, axis=1)
        assert_array_equal(s, a[:, ::2])

    def test_axis_reverse(self, xp):
        a = xp.reshape(xp.arange(12), (3, 4))

        r = axis_reverse(a, axis=0)
        assert_array_equal(r, a[::-1, :])

        r = axis_reverse(a, axis=1)
        assert_array_equal(r, a[:, ::-1])

    def test_odd_ext(self, xp):
        a = xp.asarray([[1, 2, 3, 4, 5], [9, 8, 7, 6, 5]])

        odd = odd_ext(a, 2, axis=1)
        expected = xp.asarray(
            [[-1, 0, 1, 2, 3, 4, 5, 6, 7], [11, 10, 9, 8, 7, 6, 5, 4, 3]]
        )
        assert_array_equal(odd, expected)

        odd = odd_ext(a, 1, axis=0)
        expected = xp.asarray(
            [
                [-7, -4, -1, 2, 5],
                [1, 2, 3, 4, 5],
                [9, 8, 7, 6, 5],
                [17, 14, 11, 8, 5],
            ]
        )
        assert_array_equal(odd, expected)

        assert_raises(ValueError, odd_ext, a, 2, axis=0)
        assert_raises(ValueError, odd_ext, a, 5, axis=1)

    def test_even_ext(self, xp):
        a = xp.asarray([[1, 2, 3, 4, 5], [9, 8, 7, 6, 5]])

        even = even_ext(a, 2, axis=1)
        expected = xp.asarray(
            [[3, 2, 1, 2, 3, 4, 5, 4, 3], [7, 8, 9, 8, 7, 6, 5, 6, 7]]
        )
        assert_array_equal(even, expected)

        even = even_ext(a, 1, axis=0)
        expected = xp.asarray(
            [
                [9, 8, 7, 6, 5],
                [1, 2, 3, 4, 5],
                [9, 8, 7, 6, 5],
                [1, 2, 3, 4, 5],
            ]
        )
        assert_array_equal(even, expected)

        assert_raises(ValueError, even_ext, a, 2, axis=0)
        assert_raises(ValueError, even_ext, a, 5, axis=1)

    def test_const_ext(self, xp):
        a = xp.asarray([[1, 2, 3, 4, 5], [9, 8, 7, 6, 5]])

        const = const_ext(a, 2, axis=1)
        expected = xp.asarray(
            [[1, 1, 1, 2, 3, 4, 5, 5, 5], [9, 9, 9, 8, 7, 6, 5, 5, 5]]
        )
        assert_array_equal(const, expected)

        const = const_ext(a, 1, axis=0)
        expected = xp.asarray(
            [
                [1, 2, 3, 4, 5],
                [1, 2, 3, 4, 5],
                [9, 8, 7, 6, 5],
                [9, 8, 7, 6, 5],
            ]
        )
        assert_array_equal(const, expected)

    def test_zero_ext(self, xp):
        a = xp.asarray([[1, 2, 3, 4, 5], [9, 8, 7, 6, 5]])

        zero = zero_ext(a, 2, axis=1)
        expected = xp.asarray(
            [[0, 0, 1, 2, 3, 4, 5, 0, 0], [0, 0, 9, 8, 7, 6, 5, 0, 0]]
        )
        assert_array_equal(zero, expected)

        zero = zero_ext(a, 1, axis=0)
        expected = xp.asarray(
            [
                [0, 0, 0, 0, 0],
                [1, 2, 3, 4, 5],
                [9, 8, 7, 6, 5],
                [0, 0, 0, 0, 0],
            ]
        )
        assert_array_equal(zero, expected)


@pytest_enable_array_api
class TestArrayToolsReturnTypes:
    @pytest.mark.parametrize("func", [axis_slice, axis_reverse])
    def test_type_axis_funcs(self, xp, func):
        x = xp.asarray([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        result = func(x)
        assert type(x) == type(result)

    @pytest.mark.parametrize("func", [odd_ext, even_ext, const_ext, zero_ext])
    def test_type_ext_funcs(self, xp, func):
        x = xp.asarray([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        result = func(x, 2)
        assert type(x) == type(result)
