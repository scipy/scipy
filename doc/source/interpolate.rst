========================================
Interpolation (:mod:`scipy.interpolate`)
========================================

.. module:: scipy.interpolate

Univariate interpolation
========================

.. autosummary::
   :toctree: generated/

   interp1d
   BarycentricInterpolator
   KroghInterpolator
   PiecewisePolynomial
   barycentric_interpolate
   krogh_interpolate
   piecewise_polynomial_interpolate


Multivariate interpolation
==========================

Unstructured data:

.. autosummary::
   :toctree: generated/

   griddata
   LinearNDInterpolator
   NearestNDInterpolator
   CloughTocher2DInterpolator
   Rbf
   interp2d

For data on a grid:

.. autosummary::

   RectBivariateSpline

.. seealso:: `scipy.ndimage.map_coordinates`


1-D Splines
===========

.. autosummary::
   :toctree: generated/

   UnivariateSpline
   InterpolatedUnivariateSpline
   LSQUnivariateSpline

The above univariate spline classes have the following methods:


.. autosummary::
   :toctree: generated/

   UnivariateSpline.__call__
   UnivariateSpline.derivatives
   UnivariateSpline.integral
   UnivariateSpline.roots
   UnivariateSpline.get_coeffs
   UnivariateSpline.get_knots
   UnivariateSpline.get_residual
   UnivariateSpline.set_smoothing_factor


Low-level interface to FITPACK functions:

.. autosummary::
   :toctree: generated/

   splrep
   splprep
   splev
   splint
   sproot
   spalde
   bisplrep
   bisplev


2-D Splines
===========

For data on a grid:

.. autosummary::
   :toctree: generated/

   RectBivariateSpline

For unstructured data:

.. autosummary::
   :toctree: generated/

   BivariateSpline
   SmoothBivariateSpline
   LSQBivariateSpline

Low-level interface to FITPACK functions:

.. autosummary::
   :toctree: generated/

   bisplrep
   bisplev

Additional tools
================

.. autosummary::
   :toctree: generated/

   lagrange
   approximate_taylor_polynomial
