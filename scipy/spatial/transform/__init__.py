"""
========================================================
Spatial Transformations (:mod:`scipy.spatial.transform`)
========================================================

.. currentmodule:: scipy.spatial.transform

This package implements various spatial transformations. For now,
only rotations are supported.

Rotations in 3 dimensions
-------------------------
This package supports `Rotation` transforms with the following representations

* Quaternions
* Discrete Cosine Matrices
* Euler angles
* Axis angle representation

It also implements some useful algorithms such as

* Quaternion slerp
* Quaternion spline
* Rotation estimation (Wahba's problem)

Class Reference
---------------

.. autosummary::
    :toctree: generated/

    Rotation

"""
from __future__ import division, print_function, absolute_import

from .rotation import *

__all__ = [s for s in dir() if not s.startswith('_')]

from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
