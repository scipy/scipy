#!/usr/bin/env python
#
"""
Test functions for linalg._levinson_durbin module

"""
from __future__ import division, print_function, absolute_import

import numpy as np
import scipy.linalg

from numpy.testing import (
        run_module_suite, assert_equal, assert_allclose, assert_raises)


def test_solve_equivalence():
    # For toeplitz matrices, solve_toeplitz() should be equivalent to solve().
    np.random.seed(1234)
    for n in (1, 2, 3, 10):
        c = np.random.randn(n)
        if np.random.rand() < 0.5:
            c = c + 1j * np.random.randn(n)
        r = np.random.randn(n)
        if np.random.rand() < 0.5:
            r = r + 1j * np.random.randn(n)
        y = np.random.randn(n)
        if np.random.rand() < 0.5:
            y = y + 1j * np.random.randn(n)

        # Check equivalence when both the column and row are provided.
        actual = scipy.linalg.solve_toeplitz(c, r=r, y=y)
        desired = scipy.linalg.solve(scipy.linalg.toeplitz(c, r=r), y)
        assert_allclose(actual, desired)

        # Check equivalence when the column is provided but not the row.
        actual = scipy.linalg.solve_toeplitz(c, y=y)
        desired = scipy.linalg.solve(scipy.linalg.toeplitz(c), y)
        assert_allclose(actual, desired)


def test_multiple_rhs():
    np.random.seed(1234)
    c = np.random.randn(4)
    r = np.random.randn(4)
    for yshape in ((4,), (4, 3), (4, 3, 2)):
        y = np.random.randn(*yshape)
        actual = scipy.linalg.solve_toeplitz(c, r=r, y=y)
        desired = scipy.linalg.solve(scipy.linalg.toeplitz(c, r=r), y)
        assert_equal(actual.shape, yshape)
        assert_equal(desired.shape, yshape)
        assert_allclose(actual, desired)


def test_zero_diag_error():
    # The Levinson-Durbin implementation fails when the diagonal is zero.
    np.random.seed(1234)
    n = 4
    c = np.random.randn(n)
    r = np.random.randn(n)
    y = np.random.randn(n)
    c[0] = 0
    assert_raises(np.linalg.LinAlgError,
            scipy.linalg.solve_toeplitz, c, r=r, y=y)


def test_wikipedia_counterexample():
    # The Levinson-Durbin implementation also fails in other cases.
    # This example is from the talk page of the wikipedia article.
    np.random.seed(1234)
    c = [2, 2, 1]
    y = np.random.randn(3)
    assert_raises(np.linalg.LinAlgError, scipy.linalg.solve_toeplitz, c, y=y)


if __name__ == '__main__':
    run_module_suite()
