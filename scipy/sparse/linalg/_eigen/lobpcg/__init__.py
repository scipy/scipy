"""
Locally Optimal Block Preconditioned Conjugate Gradient Method (LOBPCG)

LOBPCG is a preconditioned eigensolver for large symmetric positive definite
(SPD) generalized eigenproblems.

Call the function lobpcg - see help for lobpcg.lobpcg.

"""
from . import lobpcg as _lobpcg
from .lobpcg import *

__all__ = []
__all__ += _lobpcg.__all__

from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
