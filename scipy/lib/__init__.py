"""
Python wrappers to external libraries
=====================================

  lapack -- wrappers for LAPACK/ATLAS libraries
  blas -- wrappers for BLAS/ATLAS libraries

"""
from __future__ import division, print_function, absolute_import

__all__ = ['lapack','blas']

from numpy.testing import Tester
test = Tester().test
