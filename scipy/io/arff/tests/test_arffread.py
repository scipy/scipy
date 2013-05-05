#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

import os
import sys
from os.path import join as pjoin

if sys.version_info[0] >= 3:
    from io import StringIO
else:
    from cStringIO import StringIO

import numpy as np

from numpy.testing import TestCase, assert_array_almost_equal, assert_equal, \
        assert_, assert_raises

from scipy.io.arff.arffread import loadarff
from scipy.io.arff.arffread import read_header, parse_type, ParseArffError

data_path = pjoin(os.path.dirname(__file__), 'data')

test1 = os.path.join(data_path, 'test1.arff')
test2 = os.path.join(data_path, 'test2.arff')
test3 = os.path.join(data_path, 'test3.arff')

test4 = pjoin(data_path, 'test4.arff')
test5 = pjoin(data_path, 'test5.arff')
expect4_data = [(0.1, 0.2, 0.3, 0.4, 'class1'),
        (-0.1, -0.2, -0.3, -0.4, 'class2'),
        (1, 2, 3, 4, 'class3')]
expected_types = ['numeric', 'numeric', 'numeric', 'numeric', 'nominal']

missing = pjoin(data_path, 'missing.arff')
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
        assert_equal(meta.types(), expected_types)

    def test_filelike(self):
        """Test reading from file-like object (StringIO)"""
        f1 = open(test1)
        data1, meta1 = loadarff(f1)
        f1.close()
        f2 = open(test1)
        data2, meta2 = loadarff(StringIO(f2.read()))
        f2.close()
        assert_(data1 == data2)
        assert_(repr(meta1) == repr(meta2))


class MissingDataTest(TestCase):
    def test_missing(self):
        data, meta = loadarff(missing)
        for i in ['yop', 'yap']:
            assert_array_almost_equal(data[i], expect_missing[i])


class HeaderTest(TestCase):
    def test_type_parsing(self):
        """Test parsing type of attribute from their value."""
        ofile = open(test2)
        rel, attrs = read_header(ofile)
        ofile.close()

        expected = ['numeric', 'numeric', 'numeric', 'numeric', 'numeric',
                    'numeric', 'string', 'string', 'nominal', 'nominal']

        for i in range(len(attrs)):
            assert_(parse_type(attrs[i][1]) == expected[i])

    def test_badtype_parsing(self):
        """Test parsing wrong type of attribute from their value."""
        ofile = open(test3)
        rel, attrs = read_header(ofile)
        ofile.close()

        for name, value in attrs:
            assert_raises(ParseArffError, parse_type, value)

    def test_fullheader1(self):
        """Parsing trivial header with nothing."""
        ofile = open(test1)
        rel, attrs = read_header(ofile)
        ofile.close()

        # Test relation
        assert_(rel == 'test1')

        # Test numerical attributes
        assert_(len(attrs) == 5)
        for i in range(4):
            assert_(attrs[i][0] == 'attr%d' % i)
            assert_(attrs[i][1] == 'REAL')
        classes = attrs[4][1]

        # Test nominal attribute
        assert_(attrs[4][0] == 'class')
        assert_(attrs[4][1] == '{class0, class1, class2, class3}')

if __name__ == "__main__":
    import nose
    nose.run(argv=['', __file__])
