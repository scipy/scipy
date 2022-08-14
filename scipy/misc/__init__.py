"""
==========================================
Miscellaneous routines (:mod:`scipy.misc`)
==========================================

.. currentmodule:: scipy.misc

Various utilities that don't have another home.

.. autosummary::
   :toctree: generated/

   ascent - Get example image for processing
   central_diff_weights - Weights for an n-point central mth derivative
   derivative - Find the nth derivative of a function at a point
   face - Get example image for processing
   electrocardiogram - Load an example of a 1-D signal

"""


from ._common import *
from . import _common

# Deprecated namespaces, to be removed in v2.0.0
from . import common, doccer

__all__ = _common.__all__

del _common

from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
