"""
========================================================
Spatial Transformations (:mod:`scipy.spatial.transform`)
========================================================

.. currentmodule:: scipy.spatial.transform

This package implements various spatial transformations. For now,
only rotations are supported.

Rotations in 3 dimensions
-------------------------
.. autosummary::
    :toctree: generated/

    Rotation
    Slerp
    Spline

"""
from __future__ import division, print_function, absolute_import

from .rotation import Rotation, Slerp, Spline

__all__ = ['Rotation', 'Slerp', 'Spline']

from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
