#!/usr/bin/env python
"""Test for parsing arff headers only."""
import os

from numpy.testing import *

from scipy.io.arff.arffread import read_header, parse_type, ParseArffError

data_path = os.path.join(os.path.dirname(__file__), 'data')

test1 = os.path.join(data_path, 'test1.arff')
test2 = os.path.join(data_path, 'test2.arff')
test3 = os.path.join(data_path, 'test3.arff')

class HeaderTest(TestCase):
    def test_type_parsing(self):
        """Test parsing type of attribute from their value."""
        ofile = open(test2)
        rel, attrs = read_header(ofile)

        expected = ['numeric', 'numeric', 'numeric', 'numeric', 'numeric',
                    'numeric', 'string', 'string', 'nominal', 'nominal']

        for i in range(len(attrs)):
            assert parse_type(attrs[i][1]) == expected[i]

    def test_badtype_parsing(self):
        """Test parsing wrong type of attribute from their value."""
        ofile = open(test3)
        rel, attrs = read_header(ofile)

        for name, value in attrs:
            try:
                parse_type(value)
                raise Error("Could parse type of crap, should not happen.")
            except ParseArffError:
                pass

    def test_fullheader1(self):
        """Parsing trivial header with nothing."""
        ofile = open(test1)
        rel, attrs = read_header(ofile)

        # Test relation
        assert rel == 'test1'

        # Test numerical attributes
        assert len(attrs) == 5
        for i in range(4):
            assert attrs[i][0] == 'attr%d' % i
            assert attrs[i][1] == 'REAL'
        classes = attrs[4][1]

        # Test nominal attribute
        assert attrs[4][0] == 'class'
        assert attrs[4][1] == '{class0, class1, class2, class3}'

if __name__ == "__main__":
    nose.run(argv=['', __file__])
