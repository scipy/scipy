"""
Interpolation Tools
===================

Wrappers around FITPACK functions:

  splrep    -- find smoothing spline given (x,y) points on curve.
  splprep   -- find smoothing spline given parametrically defined curve.
  splev     -- evaluate the spline or its derivatives.
  splint    -- compute definite integral of a spline.
  sproot    -- find the roots of a cubic spline.
  spalde    -- compute all derivatives of a spline at given points.
  bisplrep   -- find bivariate smoothing spline representation.
  bisplev   -- evaluate bivariate smoothing spline. 

Interpolation Class

  interp1d -- Create a class whose instances can linearly interpolate
               to compute unknown values.
"""

postpone_import = 1
