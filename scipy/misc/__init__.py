"""
==========================================
Miscellaneous routines (:mod:`scipy.misc`)
==========================================

.. currentmodule:: scipy.misc

Various utilities that don't have another home.

Note that the Python Imaging Library (PIL) is not a dependency
of SciPy and therefore the `pilutil` module is not available on
systems that don't have PIL installed.

.. autosummary::
   :toctree: generated/

   bytescale - Byte scales an array (image)
   central_diff_weights - Weights for an n-point central m-th derivative
   comb - Combinations of N things taken k at a time, "N choose k"
   derivative -\tFind the n-th derivative of a function at a point
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
   pade - Pade approximation to function as the ratio of two polynomials
   radon -
   toimage - Takes a numpy array and returns a PIL image

"""

__all__ = ['who', 'source', 'info', 'doccer']

import doccer
from common import *
from numpy import who, source, info as _info

import sys
def info(object=None,maxwidth=76,output=sys.stdout,toplevel='scipy'):
    return _info(object, maxwidth, output, toplevel)
info.__doc__ = _info.__doc__
del sys

try:
    from pilutil import *
    import pilutil
    __all__ += pilutil.__all__
except ImportError:
    pass

import common
__all__ += common.__all__

from numpy.testing import Tester
test = Tester().test
