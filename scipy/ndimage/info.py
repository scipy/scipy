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

postpone_import = 1
depends = []
