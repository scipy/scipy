"""Tests for _sketches.py."""

from __future__ import division, print_function, absolute_import
import numpy as np
from scipy.linalg import clarkson_woodruff_transform

from numpy.testing import (assert_equal, assert_almost_equal,
                           assert_array_almost_equal, assert_array_equal,
                           assert_raises, assert_, assert_allclose)

import pytest


def make_random_dense_gaussian_matrix(n_rows, n_columns, mu=0, sigma=0.01):
    """
    Make some random data with Gaussian distributed values
    """
    res = np.random.normal(mu, sigma, n_rows*n_columns)
    return np.reshape(res, (n_rows, n_columns))

class TestClarksonWoodruffTransform(object):
    """
    Testing the Clarkson Woodruff Transform
    """
    # Big dense matrix dimensions
    n_matrix_rows = 2000
    n_matrix_columns = 100

    # Sketch matrix dimensions
    n_sketch_rows = 100

    # Repetitions/max_errors per test
    repetitions_per_test = 10
    n_max_errors = 3

    # Error threshold
    threshold = 0.1

    dense_big_matrix = make_random_dense_gaussian_matrix(n_matrix_rows,
        n_matrix_columns)
    
    def test_sketch_dimensions(self):
        sketch = clarkson_woodruff_transform(
            self.dense_big_matrix,
            self.n_sketch_rows)

        assert_(sketch.shape == (
            self.n_sketch_rows,
            self.dense_big_matrix.shape[1]))

    def test_sketch_rows_norm(self):
        # Given the probabilistic nature of the sketches
        # we run the 'test' multiple times and check that
        # we pass all/almost all the tries
        n_errors = 0
        for _ in range(self.repetitions_per_test):
            sketch = clarkson_woodruff_transform(
                self.dense_big_matrix,
                self.n_sketch_rows)

            #we could use other norms (like L2)
            err = np.linalg.norm(self.dense_big_matrix) - np.linalg.norm(sketch)
            if err > self.threshold:
                n_errors += 1
        assert_(n_errors <= self.n_max_errors)
