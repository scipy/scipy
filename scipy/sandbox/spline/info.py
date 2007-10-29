"""
Spline Tools
===================

This package contains wrappers around the netlib Dierckx library for curve
and surface fitting with splines.

Two interfaces are provided. The prefered interface is an object orientated
interface defined in spline.spline:

  UnivariateSpline     -- 1 dimensional, smoothing or interpolating spline
  InterpolatedUnivariateSpline -- i dimensional interpolating spline
  LSQUnivariateSpline  -- 1 dimensional least squares fitted spline
  SmoothBivariateSpline  -- 2 dimensional smoothing spline on scattered data
  LSQBivariateSpline  -- a least scquares fitted spline for scattered data
  RectBivariateSpline -- a 2 dimensional smoothing or interpolating spline
                         on regular gridded rectangular data

An alternative interface is a more direct wrapping of the fortran rountines:

  splrep    -- find smoothing spline given (x,y) points on curve.
  splprep   -- find smoothing spline given parametrically defined curve.
  splev     -- evaluate the spline or its derivatives.
  splint    -- compute definite integral of a spline.
  sproot    -- find the roots of a cubic spline.
  spalde    -- compute all derivatives of a spline at given points.
  bisplrep  -- find bivariate smoothing spline representation.
  bisplev   -- evaluate bivariate smoothing spline.

"""

postpone_import = 1
