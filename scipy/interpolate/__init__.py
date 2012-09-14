"""
========================================
Interpolation (:mod:`scipy.interpolate`)
========================================

.. currentmodule:: scipy.interpolate

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
   pchip
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
   RectSphereBivariateSpline

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

.. seealso::

   `scipy.ndimage.map_coordinates`,
   `scipy.ndimage.spline_filter`,
   `scipy.signal.resample`,
   `scipy.signal.bspline`,
   `scipy.signal.gauss_spline`,
   `scipy.signal.qspline1d`,
   `scipy.signal.cspline1d`,
   `scipy.signal.qspline1d_eval`,
   `scipy.signal.cspline1d_eval`,
   `scipy.signal.qspline2d`,
   `scipy.signal.cspline2d`.

"""

from interpolate import *
from fitpack import *

# New interface to fitpack library:
from fitpack2 import *

from rbf import Rbf

from polyint import *

from ndgriddata import *

__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import Tester
test = Tester().test
