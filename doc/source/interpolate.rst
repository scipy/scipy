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

.. autosummary::
   :toctree: generated/

   interp2d
   Rbf


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

.. seealso:: scipy.ndimage.map_coordinates

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
