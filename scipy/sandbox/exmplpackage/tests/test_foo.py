#!/usr/bin/env python
#
# Created by: Pearu Peterson, October 2003
# Adjusted to numpy.distutils: Pearu Peterson, October 2005
#
__usage__ = """
First ensure that scipy core modules are installed.

Build exmplpackage:
  python setup.py build
Run tests locally:
  python tests/test_foo.py [-l<int>] [-v<int>]

Run tests if scipy is installed:
  python -c 'import scipy;scipy.exmplpackage.test(level=<int>,verbosity=<int>)'
"""

import sys
from scipy.testing import *

from scipy.sandbox.exmplpackage.foo import *

class TestFooBar(TestCase):

    def check_simple(self, level=1):
        assert exmplpackage_foo_bar()=='Hello from exmplpackage_foo_bar'

class TestFooGun(TestCase):

    def check_simple(self, level=1):
        assert foo_gun()=='Hello from foo_gun'

if __name__ == "__main__":
    nose.run(argv=['', __file__])
