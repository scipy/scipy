"""
==========================================
Miscellaneous routines (:mod:`scipy.misc`)
==========================================

.. currentmodule:: scipy.misc

Various utilities that don't have another home.

Note that Pillow (https://python-pillow.org/) is not a dependency
of SciPy, but the image manipulation functions indicated in the list
below are not available without it.

.. autosummary::
   :toctree: generated/

   ascent - Get example image for processing
   bytescale - Byte scales an array (image) [requires Pillow]
   central_diff_weights - Weights for an n-point central m-th derivative
   comb - Combinations of N things taken k at a time, "N choose k" (imported from scipy.special)
   derivative - Find the n-th derivative of a function at a point
   face - Get example image for processing
   factorial  - The factorial function, n! = special.gamma(n+1) (imported from scipy.special)
   factorial2 - Double factorial, (n!)! (imported from scipy.special)
   factorialk - (...((n!)!)!...)! where there are k '!' (imported from scipy.special)
   fromimage - Return a copy of a PIL image as a numpy array [requires Pillow]
   imfilter - Simple filtering of an image [requires Pillow]
   imread - Read an image file from a filename [requires Pillow]
   imresize - Resize an image [requires Pillow]
   imrotate - Rotate an image counter-clockwise [requires Pillow]
   imsave - Save an array to an image file [requires Pillow] 
   imshow - Simple showing of an image through an external viewer [requires Pillow]
   info - Get help information for a function, class, or module
   lena - Get classic image processing example image Lena
   logsumexp - Compute the log of the sum of exponentials of input elements
               (imported from scipy.special)
   pade - Pade approximation to function as the ratio of two polynomials.
          (imported from scipy.interpolate)
   toimage - Takes a numpy array and returns a PIL image [requires Pillow]
   source - Print function source code
   who - Print the Numpy arrays in the given dictionary

"""

from __future__ import division, print_function, absolute_import

__all__ = ['who', 'source', 'info', 'doccer', 'pade',
           'comb', 'factorial', 'factorial2', 'factorialk', 'logsumexp']

from . import doccer
from .common import *
from numpy import who, source, info as _info
from scipy.interpolate._pade import pade
from scipy.special import comb, factorial, factorial2, factorialk, logsumexp

import sys


def info(object=None,maxwidth=76,output=sys.stdout,toplevel='scipy'):
    return _info(object, maxwidth, output, toplevel)
info.__doc__ = _info.__doc__
del sys

try:
    from .pilutil import *
    from . import pilutil
    __all__ += pilutil.__all__
    del pilutil
except ImportError:
    pass

from . import common
__all__ += common.__all__
del common

from numpy.testing import Tester
test = Tester().test
