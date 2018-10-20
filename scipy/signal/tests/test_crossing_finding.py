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
    _select_cross_comparator
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
        comp_a_expected, comp_b_expected = np.greater, np.less_equal
        assert_((comp_a, comp_b) == (comp_a_expected, comp_b_expected))

    def test_select_invalid(self):
        """Pass invalid arg"""
        with raises(ValueError):
            comp_a, comp_b = _select_cross_comparator('sideways')
