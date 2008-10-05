#!/usr/bin/env python
"""Tests for parsing full arff files."""
import os

import numpy as np

from numpy.testing import *

from scipy.io.arff.arffread import loadarff

data_path = os.path.join(os.path.dirname(__file__), 'data')

test4 = os.path.join(data_path, 'test4.arff')
test5 = os.path.join(data_path, 'test5.arff')
expect4_data = [(0.1, 0.2, 0.3, 0.4, 'class1'),
        (-0.1, -0.2, -0.3, -0.4, 'class2'),
        (1, 2, 3, 4, 'class3')]

missing = os.path.join(data_path, 'missing.arff')
expect_missing_raw = np.array([[1, 5], [2, 4], [np.nan, np.nan]])
expect_missing = np.empty(3, [('yop', np.float), ('yap', np.float)])
expect_missing['yop'] = expect_missing_raw[:, 0]
expect_missing['yap'] = expect_missing_raw[:, 1]

class DataTest(TestCase):
    def test1(self):
        """Parsing trivial file with nothing."""
        self._test(test4)

    def test2(self):
        """Parsing trivial file with some comments in the data section."""
        self._test(test5)

    def _test(self, test_file):
        data, meta = loadarff(test_file)
        for i in range(len(data)):
            for j in range(4):
                assert_array_almost_equal(expect4_data[i][j], data[i][j])

class MissingDataTest(TestCase):
    def test_missing(self):
        data, meta = loadarff(missing)
        for i in ['yop', 'yap']:
            assert_array_almost_equal(data[i], expect_missing[i])

if __name__ == "__main__":
    nose.run(argv=['', __file__])
