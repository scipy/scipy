"""
N-dimensional image package
===========================

This package contains various functions for multi-dimensional image
processing.

Modules
-------

.. autosummary::
   :toctree: generated/

   filters -
   fourier -
   interpolation -
   io -
   measurements -
   morphology -

Functions (partial list)
------------------------

.. autosummary::
   :toctree: generated/

   affine_transform - Apply an affine transformation
   center_of_mass - The center of mass of the values of an array at labels
   convolve - Multi-dimensional convolution
   convolve1d - 1-D convolution along the given axis
   correlate - Multi-dimensional correlation
   correlate1d - 1-D correlation along the given axis
   extrema - Min's and max's of an array at labels, with their positions
   find_objects - Find objects in a labeled array
   generic_filter - Multi-dimensional filter using a given function
   generic_filter1d - 1-D generic filter along the given axis
   geometric_transform - Apply an arbritrary geometric transform
   histogram - Histogram of the values of an array, optionally at labels
   imread - Load an image from a file
   label - Label features in an array
   laplace - n-D Laplace filter based on approximate second derivatives
   map_coordinates - Map input array to new coordinates by interpolation
   mean - Mean of the values of an array at labels
   median_filter - Calculates a multi-dimensional median filter
   percentile_filter - Calculates a multi-dimensional percentile filter
   rank_filter - Calculates a multi-dimensional rank filter
   rotate - Rotate an array
   shift - Shift an array
   standard_deviation - Standard deviation of an n-D image array
   sum - Sum of the values of the array
   uniform_filter - Multi-dimensional uniform filter
   uniform_filter1d - 1-D uniform filter along the given axis
   variance - Variance of the values of an n-D image array
   zoom - Zoom an array

Note: the above is only roughly half the functions available in this
package

Objects
-------

.. autosummary::
   :toctree: generated/

   docdict -

"""
# Copyright (C) 2003-2005 Peter J. Verveer
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above
#    copyright notice, this list of conditions and the following
#    disclaimer in the documentation and/or other materials provided
#    with the distribution.
#
# 3. The name of the author may not be used to endorse or promote
#    products derived from this software without specific prior
#    written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
# OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import numpy
from filters import *
from fourier import *
from interpolation import *
from measurements import *
from morphology import *
from io import *

# doccer is moved to scipy.misc in scipy 0.8
from scipy.misc import doccer
doccer = numpy.deprecate(doccer, old_name='doccer',
                         new_name='scipy.misc.doccer')

from info import __doc__
__version__ = '2.0'

from numpy.testing import Tester
test = Tester().test
