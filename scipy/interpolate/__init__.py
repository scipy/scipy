"""========================================
Interpolation (:mod:`scipy.interpolate`)
========================================

.. currentmodule:: scipy.interpolate

Sub-package for objects used in interpolation.

As listed below, this sub-package contains spline functions and classes,
one-dimensional and multi-dimensional (univariate and multivariate)
interpolation classes, Lagrange and Taylor polynomial interpolators, and
wrappers for `FITPACK <http://www.netlib.org/dierckx/>`__
and DFITPACK functions.

Univariate interpolation
========================

.. autosummary::
   :toctree: generated/

   interp1d
   BarycentricInterpolator
   KroghInterpolator
   PchipInterpolator
   barycentric_interpolate
   krogh_interpolate
   pchip_interpolate
   Akima1DInterpolator
   CubicSpline
   PPoly
   BPoly


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
   :toctree: generated/

   interpn
   RegularGridInterpolator
   RectBivariateSpline

.. seealso:: `scipy.ndimage.interpolation.map_coordinates`


1-D Splines
===========

.. autosummary::
   :toctree: generated/

   UnivariateSpline
   InterpolatedUnivariateSpline
   LSQUnivariateSpline


Functional interface to FITPACK functions:

.. autosummary::
   :toctree: generated/

   splrep
   splprep
   splev
   splint
   sproot
   spalde
   splder
   splantider
   insert


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
   SmoothSphereBivariateSpline
   LSQBivariateSpline
   LSQSphereBivariateSpline

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

   `scipy.ndimage.interpolation.map_coordinates`,
   `scipy.ndimage.interpolation.spline_filter`,
   `scipy.signal.resample`,
   `scipy.signal.bspline`,
   `scipy.signal.gauss_spline`,
   `scipy.signal.qspline1d`,
   `scipy.signal.cspline1d`,
   `scipy.signal.qspline1d_eval`,
   `scipy.signal.cspline1d_eval`,
   `scipy.signal.qspline2d`,
   `scipy.signal.cspline2d`.

Functions existing for backward compatibility (should not be used in
new code):

.. autosummary::
   :toctree: generated/

   ppform
   spleval
   spline
   splmake
   spltopp
   pchip

"""
from __future__ import division, print_function, absolute_import

from .interpolate import *
from .fitpack import *

# New interface to fitpack library:
from .fitpack2 import *

from .rbf import Rbf

from .polyint import *

from ._monotone import *

from .ndgriddata import *

__all__ = [s for s in dir() if not s.startswith('_')]
from numpy.testing import Tester
test = Tester().test
bench = Tester().bench
