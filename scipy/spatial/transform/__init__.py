"""
Spatial Transformations (:mod:`scipy.spatial.transform`)
========================================================

.. currentmodule:: scipy.spatial.transform

This package implements various spatial transformations. For now,
rotations and translation are individually supported.

Rotations and translation in 3 dimensions
-------------------------
.. autosummary::
   :toctree: generated/

   Rotation
   Slerp
   RotationSpline
   Translation
"""
from ._rotation import Rotation, Slerp
from ._rotation_spline import RotationSpline
from ._translation import Translation

# Deprecated namespaces, to be removed in v2.0.0
from . import rotation
from . import translation

__all__ = ['Rotation', 'Slerp', 'RotationSpline', 'Translation']

from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
