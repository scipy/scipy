"""
Python wrappers to external libraries
=====================================

- lapack -- wrappers for `LAPACK/ATLAS <http://netlib.org/lapack/>`_
            libraries
- blas -- wrappers for `BLAS/ATLAS <http://www.netlib.org/blas/>`_
          libraries

"""
from __future__ import division, print_function, absolute_import

__all__ = ['lapack','blas']

from numpy.testing import Tester
test = Tester().test
