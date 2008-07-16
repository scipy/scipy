#!/usr/bin/env python
#
__usage__ = """
First ensure that scipy core modules are installed.

Run benchmarks if scipy is installed:
  python -c 'import scipy;scipy.exmplpackage.bench(label=<str>,verbose=<int>)'
"""

import sys
from numpy.testing import *

from scipy.sandbox.exmplpackage.foo import *

class BenchFooBar(TestCase):

    def bench_simple(self):
        print 'A long timed benchmark here'

