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

import sys
from scipy_test.testing import *
set_package_path()
from linalg import flapack
from linalg import clapack
del sys.path[0]

class test_lapack(ScipyTestCase):

    def check_flapack(self):
        if hasattr(flapack,'empty_module'):
            print """
****************************************************************
WARNING: flapack module is empty
-----------
See scipy/INSTALL.txt for troubleshooting.
****************************************************************
"""
    def check_clapack(self):
        if hasattr(clapack,'empty_module'):
            print """
****************************************************************
WARNING: clapack module is empty
-----------
See scipy/INSTALL.txt for troubleshooting.
Notes:
* If atlas library is not found by scipy/system_info.py,
  then scipy uses flapack instead of clapack.
****************************************************************
"""

if __name__ == "__main__":
    ScipyTest('linalg.lapack').run()
