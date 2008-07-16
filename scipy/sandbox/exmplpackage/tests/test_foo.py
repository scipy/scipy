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
  python tests/test_foo.py

Run tests if scipy is installed:
  python -c 'import scipy;scipy.exmplpackage.test(label=<str>,verbose=<int>)'
"""

import sys
from numpy.testing import *

from scipy.sandbox.exmplpackage.foo import *

class TestFooBar(TestCase):

    def test_simple(self):
        assert exmplpackage_foo_bar()=='Hello from exmplpackage_foo_bar'

class TestFooGun(TestCase):

    def test_simple(self):
        assert foo_gun()=='Hello from foo_gun'

if __name__ == "__main__":
    nose.run(argv=['', __file__])
