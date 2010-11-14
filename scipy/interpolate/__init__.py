"""
Sub-package for objects used in interpolation.

As listed below, this sub-package contains spline functions and classes,
one-dimensional and multi-dimensional (univariate and multivariate)
interpolation classes, Lagrange and Taylor polynomial interpolators, and
wrappers for `FITPACK <http://www.cisl.ucar.edu/softlib/FITPACK.html>`_
and DFITPACK functions.

Spline Functions
----------------

.. autosummary::
   :toctree: generated/

   bisplev
   bisplrep
   insert
   spalde
   splev
   spleval
   splint
   spline
   splmake
   splprep
   splrep
   spltopp
   sproot

Spline Classes
--------------

.. autosummary::
   :toctree: generated/

   UnivariateSpline
   InterpolatedUnivariateSpline
   LSQUnivariateSpline
   BivariateSpline
   SmoothBivariateSpline

Interpolation Classes (univariate)
----------------------------------

.. autosummary::
   :toctree: generated/

   interp1d
   BarycentricInterpolator
   barycentric_interpolate
   KroghInterpolator
   krogh_interpolate
   PiecewisePolynomial
   piecewise_polynomial_interpolate
   ppform

Interpolation Classes (multivariate)
------------------------------------

.. autosummary::
   :toctree: generated/

   interp2d
   Rbf

Additional tools
----------------

.. autosummary::
   :toctree: generated/

   lagrange
   approximate_taylor_polynomial

Wrappers around FITPACK functions
---------------------------------

.. autosummary::
   :toctree: generated/

   fitpack.bisplev
   fitpack.bisplrep
   fitpack.insert
   fitpack.spalde
   fitpack.splev
   fitpack.splint
   fitpack.splprep
   fitpack.splrep
   fitpack.sproot

Wrappers around DFITPACK functions
----------------------------------

   `dfitpack.bispeu`
   `dfitpack.bispev`
   `dfitpack.curfit`
   `dfitpack.dblint`
   `dfitpack.fpcurf0`
   `dfitpack.fpcurf1`
   `dfitpack.fpcurfm1`
   `dfitpack.parcur`
   `dfitpack.percur`
   `dfitpack.regrid_smth`
   `dfitpack.spalde`
   `dfitpack.splder`
   `dfitpack.splev`
   `dfitpack.splint`
   `dfitpack.sproot`
   `dfitpack.surfit_lsq`
   `dfitpack.surfit_smth`

See Also
--------

.. autosummary::
   :toctree: generated/

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
#
# interpolate - Interpolation Tools
#

from info import __doc__

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
