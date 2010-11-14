"""
Sub-package for objects used in interpolation.

As listed below, this sub-package contains spline functions and classes,
one-dimensional and multi-dimensional (univariate and multivariate)
interpolation classes, Lagrange and Taylor polynomial interpolators, and
wrappers for `FITPACK <http://www.cisl.ucar.edu/softlib/FITPACK.html>`_
and DFITPACK functions.

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

Data given on a regular grid:

.. autosummary::
   :toctree: generated/

   interp2d


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

.. seealso::

   ndimage.map_coordinates
   ndimage.spline_filter
   signal.resample
   signal.bspline
   signal.gauss_spline
   signal.qspline1d
   signal.cspline1d
   signal.qspline1d_eval
   signal.cspline1d_eval
   signal.qspline2d
   signal.cspline2d

"""

postpone_import = 1
