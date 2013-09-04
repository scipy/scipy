"""
===========================================================
Finite-differences derivatives routines (:mod:`scipy.diff`)
===========================================================

.. currentmodule:: scipy.diff

Functions for finite differences.

TODO replace the function names below

.. autosummary::
   :toctree: generated/

   bytescale - Byte scales an array (image)
   central_diff_weights - Weights for an n-point central m-th derivative
   comb - Combinations of N things taken k at a time, "N choose k"
   derivative - Find the n-th derivative of a function at a point
   factorial  - The factorial function, n! = special.gamma(n+1)
   factorial2 - Double factorial, (n!)!
   factorialk - (...((n!)!)!...)! where there are k '!'
   fromimage - Return a copy of a PIL image as a numpy array
   imfilter - Simple filtering of an image
   imread - Read an image file from a filename
   imresize - Resize an image
   imrotate - Rotate an image counter-clockwise
   imsave - Save an array to an image file
   imshow - Simple showing of an image through an external viewer
   info - Get help information for a function, class, or module
   lena - Get classic image processing example image Lena
   logsumexp - Compute the log of the sum of exponentials of input elements
   pade - Pade approximation to function as the ratio of two polynomials
   toimage - Takes a numpy array and returns a PIL image
   who - Print the Numpy arrays in the given dictionary

"""
from __future__ import division, print_function, absolute_import

__all__ = []

from ._numdifftools_core import *

from . import _numdifftools_core
__all__ += _numdifftools_core.__all__
del _numdifftools_core

from numpy.testing import Tester
test = Tester().test
