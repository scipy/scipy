"""
===========================================================
Finite-differences derivatives routines (:mod:`scipy.diff`)
===========================================================

.. currentmodule:: scipy.diff

Functions for finite differences using numdifftools.

.. autosummary::
   :toctree: generated/

   dea3 - Extrapolate a slowly convergent sequence
   Derivative - Estimate the nth derivative of a function at a point
   Jacobian - Estimate the Jacobian matrix of a function at a point
   Gradient - Estimate the gradient of a function at a point
   Hessian - Estimate the Hessian matrix of a function at a point
   Hessdiag - Estimate the diagonal elements of the Hessian at a point

"""
from __future__ import division, print_function, absolute_import

__all__ = []

from ._numdifftools_core import *

from . import _numdifftools_core
__all__ += _numdifftools_core.__all__
del _numdifftools_core

from numpy.testing import Tester
test = Tester().test
