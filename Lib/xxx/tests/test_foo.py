#!/usr/bin/env python
#
# Created by: Pearu Peterson, October 2003
#
__usage__ = """
First ensure that scipy_core modules are installed.

Build xxx:
  python setup_xxx.py build
Run tests locally:
  python tests/test_foo.py [-l<int>] [-v<int>]

Run tests if scipy is installed:
  python -c 'import scipy;scipy.xxx.test(level=<int>,verbosity=<int>)'
"""

import sys
from scipy_test.testing import *

set_package_path()
from xxx.foo import *
del sys.path[0]

class test_foo_bar(ScipyTestCase):

    def check_simple(self, level=1):
        assert xxx_foo_bar()=='Hello from xxx_foo_bar'

class test_foo_gun(ScipyTestCase):

    def check_simple(self, level=1):
        assert foo_gun()=='Hello from foo_gun'

if __name__ == "__main__":
    ScipyTest('xxx.foo').run()
