"""
========================================
Interpolation (:mod:`scipy.interpolate`)
========================================

.. currentmodule:: scipy.interpolate

Sub-package for functions and objects used in interpolation.

See the :ref:`user guide <tutorial-interpolate>` for recommendations on choosing a
routine, and other usage details.


Univariate interpolation
========================

.. autosummary::
   :toctree: generated/

   make_interp_spline
   CubicSpline
   PchipInterpolator
   Akima1DInterpolator
   FloaterHormannInterpolator
   BarycentricInterpolator
   KroghInterpolator
   CubicHermiteSpline

**Low-level data structures for univariate interpolation:**

.. autosummary::
   :toctree: generated/

   PPoly
   BPoly
   BSpline


Multivariate interpolation
==========================

**Unstructured data**

.. autosummary::
   :toctree: generated/

   LinearNDInterpolator
   NearestNDInterpolator
   CloughTocher2DInterpolator
   RBFInterpolator

**For data on a grid:**

.. autosummary::
   :toctree: generated/

   RegularGridInterpolator

.. seealso::

    `scipy.ndimage.map_coordinates`,
    :ref:`An example wrapper for map_coordinates <tutorial-interpolate_cartesian-grids>`


**Low-level data structures for tensor product polynomials and splines:**


.. autosummary::
   :toctree: generated/

   NdPPoly
   NdBSpline


1-D spline smoothing and approximation
======================================

.. autosummary::
   :toctree: generated/

   make_lsq_spline
   make_smoothing_spline
   make_splrep
   make_splprep
   generate_knots

Rational Approximation
======================

.. autosummary::
   :toctree: generated/

   AAA


Interfaces to FITPACK routines for 1D and 2D spline fitting
===========================================================

This section lists wrappers for `FITPACK <http://www.netlib.org/dierckx/>`__
functionality for 1D and 2D smoothing splines. In most cases, users are better off
using higher-level routines listed in previous sections.


1D FITPACK splines
------------------

This package provides two sets of functionally equivalent wrappers: object-oriented and
functional.

**Functional FITPACK interface:**


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

**Object-oriented FITPACK interface:**

.. autosummary::
   :toctree: generated/

   UnivariateSpline
   InterpolatedUnivariateSpline
   LSQUnivariateSpline


2D FITPACK splines
------------------

**For data on a grid:**

.. autosummary::
   :toctree: generated/

   RectBivariateSpline
   RectSphereBivariateSpline

**For unstructured data (OOP interface):**

.. autosummary::
   :toctree: generated/

   BivariateSpline
   SmoothBivariateSpline
   SmoothSphereBivariateSpline
   LSQBivariateSpline
   LSQSphereBivariateSpline

**For unstructured data (functional interface):**

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
   pade

   interpn
   griddata
   barycentric_interpolate
   krogh_interpolate
   pchip_interpolate
   Rbf
   interp1d
   interp2d

.. seealso::

   `scipy.ndimage.map_coordinates`,
   `scipy.ndimage.spline_filter`,

"""  # noqa: E501
from . import _interpolate
from ._interpolate import *
from . import _fitpack_py
from ._fitpack_py import *

from . import _fitpack2
from ._fitpack2 import *

from ._rbf import Rbf

from . import _rbfinterp
from ._rbfinterp import *

from . import _polyint
from ._polyint import *

from . import _cubic
from ._cubic import *

from . import _ndgriddata
from ._ndgriddata import *

from . import _bsplines
from ._bsplines import *
from ._fitpack_repro import generate_knots, make_splrep, make_splprep

from . import _pade
from ._pade import *

from . import _rgi
from ._rgi import *

from ._ndbspline import NdBSpline

from . import _bary_rational
from ._bary_rational import *

# Deprecated namespaces, to be removed in v2.0.0
from . import fitpack, fitpack2, interpolate, ndgriddata, polyint, rbf, interpnd

__all__ = []
__all__ += _interpolate.__all__
__all__ += _fitpack_py.__all__
__all__ += _fitpack2.__all__
__all__ += ['Rbf']
__all__ += _rbfinterp.__all__
__all__ += _polyint.__all__
__all__ += _cubic.__all__
__all__ += _ndgriddata.__all__
__all__ += _bsplines.__all__
__all__ += ['generate_knots', 'make_splrep', 'make_splprep']
__all__ += _pade.__all__
__all__ += _rgi.__all__
__all__ += ['NdBSpline']
__all__ += _bary_rational.__all__
__all__ += ['fitpack', 'fitpack2', 'interpolate', 'ndgriddata', 'polyint',
            'rbf', 'interpnd']

from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester

# Backward compatibility
pchip = PchipInterpolator
