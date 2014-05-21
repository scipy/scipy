"""
Locally Optimal Block Preconditioned Conjugate Gradient Method (LOBPCG)

LOBPCG is a preconditioned eigensolver for large symmetric positive definite
(SPD) generalized eigenproblems.

Call the function lobpcg - see help for lobpcg.lobpcg.

"""
from __future__ import division, print_function, absolute_import

from .lobpcg import *

from numpy.testing import Tester
test = Tester().test
bench = Tester().bench
