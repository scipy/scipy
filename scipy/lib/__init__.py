"""
Python wrappers to external libraries
=====================================

  lapack -- wrappers for LAPACK/ATLAS libraries
  blas -- wrappers for BLAS/ATLAS libraries

"""

__all__ = ['lapack','blas']

from numpy.testing import Tester
test = Tester().test
