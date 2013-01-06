"""
Sparse Eigenvalue Solvers
-------------------------

The submodules of sparse.linalg.eigen:
    1. lobpcg: Locally Optimal Block Preconditioned Conjugate Gradient Method

"""
from __future__ import division, print_function, absolute_import

from .arpack import *
from .lobpcg import *

__all__ = [s for s in dir() if not s.startswith('_')]
from numpy.testing import Tester
test = Tester().test
bench = Tester().bench
