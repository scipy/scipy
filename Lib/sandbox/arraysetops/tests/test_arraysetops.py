#!/usr/bin/env python
#
# Created by: Travis Oliphant, April 2004
#
""" Test functions for sparse matrices

"""
__usage__ = """
Build arraysetops:
  python setup.py build
Run tests if scipy is installed:
  python -c 'import scipy;scipy.arraysetops.test(<level>)'
Run tests if arraysetops is not installed:
  python tests/test_sparse.py [<level>]
"""

import numpy
from numpy import array, rand

import sys
import random
from numpy.test.testing import *
set_package_path()
from scipy.arraysetops import *
restore_path()

class test_aso( ScipyTestCase ):

    def chech_all():
        test_unique1d()
        test_intersect1d()
        test_intersect1d_nu()
        test_setxor1d()
        test_setmember1d()
        test_union1d()
        test_setdiff1d()
        test_manyways()

if __name__ == "__main__":
    ScipyTest().run()
