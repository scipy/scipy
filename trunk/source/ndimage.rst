=========================================================
Multi-dimensional image processing (:mod:`scipy.ndimage`)
=========================================================

.. module:: scipy.ndimage

Functions for multi-dimensional image processing.

Filters :mod:`scipy.ndimage.filters`
====================================

.. module:: scipy.ndimage.filters

.. autosummary::
   :toctree: generated/

   convolve
   convolve1d
   correlate
   correlate1d
   gaussian_filter
   gaussian_filter1d
   gaussian_gradient_magnitude
   gaussian_laplace
   generic_filter
   generic_filter1d
   generic_gradient_magnitude
   generic_laplace
   laplace
   maximum_filter
   maximum_filter1d
   median_filter
   minimum_filter
   minimum_filter1d
   percentile_filter
   prewitt
   rank_filter
   sobel
   uniform_filter
   uniform_filter1d

Fourier filters :mod:`scipy.ndimage.fourier`
============================================

.. module:: scipy.ndimage.fourier

.. autosummary::
   :toctree: generated/

   fourier_ellipsoid
   fourier_gaussian
   fourier_shift
   fourier_uniform

Interpolation :mod:`scipy.ndimage.interpolation`
================================================

.. module:: scipy.ndimage.interpolation

.. autosummary::
   :toctree: generated/

   affine_transform
   geometric_transform
   map_coordinates
   rotate
   shift
   spline_filter
   spline_filter1d
   zoom

Measurements :mod:`scipy.ndimage.measurements`
==============================================

.. module:: scipy.ndimage.measurements

.. autosummary::
   :toctree: generated/

   center_of_mass
   extrema
   find_objects
   histogram
   label
   maximum
   maximum_position
   mean
   minimum
   minimum_position
   standard_deviation
   sum
   variance
   watershed_ift

Morphology :mod:`scipy.ndimage.morphology`
==========================================

.. module:: scipy.ndimage.morphology

.. autosummary::
   :toctree: generated/

   binary_closing
   binary_dilation
   binary_erosion
   binary_fill_holes
   binary_hit_or_miss
   binary_opening
   binary_propagation
   black_tophat
   distance_transform_bf
   distance_transform_cdt
   distance_transform_edt
   generate_binary_structure
   grey_closing
   grey_dilation
   grey_erosion
   grey_opening
   iterate_structure
   morphological_gradient
   morphological_laplace
   white_tophat
