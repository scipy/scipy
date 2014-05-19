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

from numpy.testing import (TestCase, assert_array_almost_equal, assert_equal,
        assert_, assert_raises, run_module_suite)

from scipy.io.arff.arffread import loadarff
from scipy.io.arff.arffread import read_header, parse_type, ParseArffError

import scipy.sparse as sp
import scipy

data_path = pjoin(os.path.dirname(__file__), 'data')

test1 = os.path.join(data_path, 'test1.arff')
test2 = os.path.join(data_path, 'test2.arff')
test3 = os.path.join(data_path, 'test3.arff')

test4 = pjoin(data_path, 'test4.arff')
test5 = pjoin(data_path, 'test5.arff')
test6 = pjoin(data_path, 'test6.arff')
test7 = pjoin(data_path, 'test7.arff')
expect4_data = [(0.1, 0.2, 0.3, 0.4, 'class1'),
        (-0.1, -0.2, -0.3, -0.4, 'class2'),
        (1, 2, 3, 4, 'class3')]
expected_types = ['numeric', 'numeric', 'numeric', 'numeric', 'nominal']

missing = pjoin(data_path, 'missing.arff')
expect_missing_raw = np.array([[1, 5], [2, 4], [np.nan, np.nan]])
expect_missing = np.empty(3, [('yop', np.float), ('yap', np.float)])
expect_missing['yop'] = expect_missing_raw[:, 0]
expect_missing['yap'] = expect_missing_raw[:, 1]
expect7_data = np.matrix([[0.488144,0.198068,0.,0.,0.],[0.,0., 0.247897,0.153404,0.],[-0.0621701,0.184011,-0.19966,0.0513624,0.133107]])
expected_types7 = ['nominal', 'nominal', 'nominal', 'nominal', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric']


class DataTest(TestCase):
    def test1(self):
        # Parsing trivial file with nothing.
        self._test(test4)

    def test2(self):
        # Parsing trivial file with some comments in the data section.
        self._test(test5)

    def test3(self):
        # Parsing trivial file with nominal attribute of 1 character.
        self._test(test6)
        
    def test4(self):
        # Parsing sparse
        self._test(test7)

    def _test(self, test_file):
        data, meta, classes = loadarff(test_file,returnclasses=True)

        if sp.issparse(data):
            ld = data.shape[0]
        else:
            ld = len(data)
        if sp.issparse(data):
            assert_(scipy.absolute(data.todense()-expect7_data).sum() == 0)
        else:
            for i in range(ld):
                for j in range(4):
                    assert_array_almost_equal(expect4_data[i][j], data[i][j])
        if sp.issparse(data):
            assert_equal(meta.types(), expected_types7)
        else:
            assert_equal(meta.types(), expected_types)

    def test_filelike(self):
        # Test reading from file-like object (StringIO)
        f1 = open(test1)
        data1, meta1, classes = loadarff(f1,returnclasses=True)
        f1.close()
        f2 = open(test1)
        data2, meta2, classes = loadarff(StringIO(f2.read()),returnclasses=True)
        f2.close()
        assert_(data1 == data2)
        assert_(repr(meta1) == repr(meta2))


class MissingDataTest(TestCase):
    def test_missing(self):
        data, meta, classes = loadarff(missing,returnclasses=True)
        for i in ['yop', 'yap']:
            assert_array_almost_equal(data[i], expect_missing[i])


class HeaderTest(TestCase):
    def test_type_parsing(self):
        # Test parsing type of attribute from their value.
        ofile = open(test2)
        rel, attrs, nrclasses = read_header(ofile,returnclasses=True)
        ofile.close()

        expected = ['numeric', 'numeric', 'numeric', 'numeric', 'numeric',
                    'numeric', 'string', 'string', 'nominal', 'nominal']

        for i in range(len(attrs)):
            assert_(parse_type(attrs[i][1]) == expected[i])

    def test_badtype_parsing(self):
        # Test parsing wrong type of attribute from their value.
        ofile = open(test3)
        rel, attrs, nrclasses = read_header(ofile,returnclasses=True)
        ofile.close()

        for name, value in attrs:
            assert_raises(ParseArffError, parse_type, value)

    def test_fullheader1(self):
        # Parsing trivial header with nothing.
        ofile = open(test1)
        rel, attrs, nrclasses = read_header(ofile,returnclasses=True)
        ofile.close()

        # Test relation
        assert_(rel == 'test1')

        # Test numerical attributes
        assert_(len(attrs) == 5)
        for i in range(4):
            assert_(attrs[i][0] == 'attr%d' % i)
            assert_(attrs[i][1] == 'REAL')

        # Test nominal attribute
        assert_(attrs[4][0] == 'class')
        assert_(attrs[4][1] == '{class0, class1, class2, class3}')


if __name__ == "__main__":
    run_module_suite()
