#!/usr/bin/env python
#
# Created by: Pearu Peterson, September 2002
#

__usage__ = """
Build linalg:
  python setup_linalg.py build
Run tests if scipy is installed:
  python -c 'import scipy;scipy.linalg.test(<level>)'
Run tests if linalg is not installed:
  python tests/test_lapack.py [<level>]
"""

from scipy_test.testing import ScipyTestCase
import unittest
from scipy_distutils.misc_util import PostponedException

import sys
from scipy_test.testing import set_package_path
set_package_path()
try: from linalg import flapack
except: flapack = PostponedException()
try: from linalg import clapack
except: clapack = PostponedException()
del sys.path[0]

def test_suite(level=1):
    suites = []
    if level > 0:
        if isinstance(flapack,PostponedException):
            print """
****************************************************************
WARNING: Importing flapack failed with the following exception:
-----------
%s
-----------
See scipy/INSTALL.txt for troubleshooting.
****************************************************************
""" %(flapack.__doc__)
        if isinstance(clapack,PostponedException):
            print """
****************************************************************
WARNING: Importing clapack failed with the following exception:
-----------
%s
-----------
See scipy/INSTALL.txt for troubleshooting.
Notes:
* If atlas library is not found by scipy/system_info.py,
  then scipy skips building clapack and uses flapack instead.
****************************************************************
""" %(clapack.__doc__)

    total_suite = unittest.TestSuite(suites)
    return total_suite

def test(level=10):
    all_tests = test_suite(level)
    runner = unittest.TextTestRunner()
    runner.run(all_tests)
    return runner

if __name__ == "__main__":
    if len(sys.argv)>1:
        level = eval(sys.argv[1])
    else:
        level = 1
    test(level)
