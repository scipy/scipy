from __future__ import division, print_function, absolute_import

import copy

import numpy as np
from numpy.testing import (
    assert_,
    assert_equal,
    assert_allclose,
    assert_array_equal
)
import pytest
from pytest import raises, warns

from scipy.signal._crossing_finding import (
    _select_cross_comparator,
    _boolcross
)


class TestSelectCrossComparator(object):

    def test_select_up(self):
        """Select up-crossing comparators"""
        comp_a, comp_b = _select_cross_comparator('up')
        comp_a_expected, comp_b_expected = np.less_equal, np.greater
        assert_((comp_a, comp_b) == (comp_a_expected, comp_b_expected))

    def test_select_down(self):
        """Select down-crossing comparators"""
        comp_a, comp_b = _select_cross_comparator('down')
        comp_a_expected, comp_b_expected = np.greater_equal, np.less
        assert_((comp_a, comp_b) == (comp_a_expected, comp_b_expected))

    def test_select_invalid(self):
        """Pass invalid arg"""
        with raises(ValueError):
            comp_a, comp_b = _select_cross_comparator('sideways')


class TestBoolCross(object):
    def test_no_upcross_1d(self):
        "Test up-crossing on 1d-array - no crossings."
        false_array = np.zeros(5, dtype=bool)
        z1 = np.zeros(5)

        i = _boolcross(z1, -1, 'up')
        assert_equal(len(i), 5)
        assert_array_equal(i, false_array)

    def test_no_upcross_nd(self):
        "Test up-crossing on ndarray - no crossings."
        false_array = np.zeros((3, 5, 8), dtype=bool)
        z1 = np.zeros((3, 5, 8))

        i = _boolcross(z1, -1, 'up', axis=0)
        assert_equal(i.shape, (3, 5, 8))
        assert_array_equal(i, false_array)

        j = _boolcross(z1, -1, 'up', axis=1)
        assert_equal(j.shape, (3, 5, 8))
        assert_array_equal(j, false_array)

    def test_no_downcross_1d(self):
        "Test down-crossing on 1d-array - no crossings."
        false_array = np.zeros(5, dtype=bool)
        z1 = np.zeros(5)

        i = _boolcross(z1, 1, 'down')
        assert_equal(len(i), 5)
        assert_array_equal(i, false_array)

    def test_no_downcross_nd(self):
        "Test down-crossing on ndarray - no crossings."
        false_array = np.zeros((3, 5, 8), dtype=bool)
        z1 = np.zeros((3, 5, 8))

        i = _boolcross(z1, 1, 'down', axis=0)
        assert_equal(i.shape, (3, 5, 8))
        assert_array_equal(i, false_array)

        j = _boolcross(z1, 1, 'down', axis=1)
        assert_equal(j.shape, (3, 5, 8))
        assert_array_equal(j, false_array)

    def test_upcross_1d(self):
        "Test up-crossing on 1d-array."
        result_array = np.array([False, True, False, False, False, True, False])

        z1 = np.array([-1., 0., 1., 0., -1., 0., 1.])
        i = _boolcross(z1, 0., 'up')
        assert_equal(len(i), 7)
        assert_array_equal(i, result_array)

        z2 = np.array([-1., 0., 1., 0., -1., 0., 1.]) + 4.
        j = _boolcross(z2, None, 'up')
        assert_equal(len(j), 7)
        assert_array_equal(j, result_array)

    def test_upcross_nd_axis_0(self):
        "Test up-crossing on ndarray along axis 0 - threshold set to 0."
        result_array = np.array([
            [False, True, False, False, False, True, False],
            [True, False, False, False, True, False, False],
            [False, False, False, False, False, False, False]
        ])

        z1 = np.array([
            [-1., 0., 1., 0., -1., 0., 1.],
            [0., 1., 0., -1., 0., 1., 0.],
            [1., 0., -1., 0., 1., 0., -1.]
        ])

        i = _boolcross(z1, 0., 'up', axis=0)
        assert_equal(i.shape, (3, 7))
        assert_array_equal(i, result_array)

    def test_upcross_nd_axis_1(self):
        "Test up-crossing on ndarray along axis 1 - threshold set to None."
        result_array = np.array([
            [False, True, False, False, False, True, False],
            [True, False, False, False, True, False, False],
            [False, False, False, True, False, False, False]
        ])

        z1 = np.array([
            [-1., 0., 1., 0., -1., 0., 1.],
            [0., 1., 0., -1., 0., 1., 0.],
            [1., 0., -1., 0., 1., 0., -1.]
        ])

        i = _boolcross(z1, None, 'up', axis=1)
        assert_equal(i.shape, (3, 7))
        assert_array_equal(i, result_array)

    def test_downcross_nd_axis_0(self):
        "Test down-crossing on ndarray along axis 0 - threshold set to 0."
        result_array = np.array([
            [False, False, False, True, False, False, False],
            [False, False, True, False, False, False, True],
            [False, False, False, False, False, False, False]
        ])

        z1 = np.array([
            [-1., 0., 1., 0., -1., 0., 1.],
            [0., 1., 0., -1., 0., 1., 0.],
            [1., 0., -1., 0., 1., 0., -1.]
        ])

        i = _boolcross(z1, 0., 'down', axis=0)
        assert_equal(i.shape, (3, 7))
        assert_array_equal(i, result_array)

    def test_downcross_nd_axis_1(self):
        "Test down-crossing on ndarray along axis 1 - threshold set to None."
        result_array = np.array([
            [False, False, False, True, False, False, False],
            [False, True, False, False, False, True, False],
            [False, True, False, False, False, True, False]
        ])

        z1 = np.array([
            [-1., 0., 1., 0., -1., 0., 1.],
            [0., 1., 0., -1., 0., 1., 0.],
            [1., 0., -1., 0., 1., 0., -1.]
        ])

        i = _boolcross(z1, None, 'down', axis=1)
        assert_equal(i.shape, (3, 7))
        assert_array_equal(i, result_array)
