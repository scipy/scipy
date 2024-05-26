"""
===================================================================
Elementwise Scalar Optimization (:mod:`scipy.optimize.elementwise`)
===================================================================

.. currentmodule:: scipy.optimize.elementwise

This module provides a collection of functions for rootfinding and
minimization of scalar, real-valued functions of one variable. Unlike their
counterparts in the base :mod:`scipy.optimize` namespace, these functions work
elementwise, enabling the solution of many related problems in an efficient,
vectorized call. Furthermore, when environment variable ``SCIPY_ARRAY_API=1``,
these functions can accept non-NumPy, array API standard compatible arrays and
perform all calculations using the corresponding array library (e.g. PyTorch,
JAX, CuPy).

Rootfinding
===========

.. autosummary::
   :toctree: generated/

   rootfind
   bracket_root

Rootfinding
===========

.. autosummary::
   :toctree: generated/

   minimize
   bracket_minimum

"""
from ._elementwise import rootfind, bracket_root, minimize, bracket_minimum  # noqa: F401, E501

__all__ = ["rootfind", "bracket_root", "minimize", "bracket_minimum"]
