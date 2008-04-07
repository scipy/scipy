#!/usr/bin/env python
"""Test for parsing arff headers only."""
import os

from scipy.testing import *

from scipy.io.arff.arffread import read_header, MetaData

data_path = os.path.join(os.path.dirname(__file__), 'data')

test1 = os.path.join(data_path, 'test1.arff')

class HeaderTest(TestCase):
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
