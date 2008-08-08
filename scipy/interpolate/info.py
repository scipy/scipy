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

  UnivariateSpline             -- A more recent, object-oriented wrapper;
                                  finds a (possibly smoothed) interpolating
                                  spline.
  InterpolatedUnivariateSpline
  LSQUnivariateSpline
  BivariateSpline              -- A more recent, object-oriented wrapper;
                                  finds a interpolating spline for a
                                  bivariate function.

  SmoothBivariateSpline

Interpolation Classes (univariate)

  interp1d -- Create a class whose instances can linearly interpolate
               to compute unknown values of a univariate function.
  BarycentricInterpolator -- Compute with a numerically-stable version
               of the Lagrange interpolating polynomial.
  barycentric_interpolate -- procedural interface to the above
  KroghInterpolator -- Compute with the Hermite interpolating polynomial
               (allows the specification of derivatives at some points).
  krogh_interpolate -- procedural interface to the above
  PiecewisePolynomial -- Spline that is specified by giving positions and
               derivatives at every knot; allows high orders and
               efficient appending.
  piecewise_polynomial_interpolate -- procedural interface to the above

Interpolation Classes (multivariate)

  interp2d -- Create a class whose instances can interpolate
               to compute unknown values of a bivariate function.
  Rbf -- Apply Radial Basis Functions to interpolate scattered N-D data.

Additional tools

  lagrange -- Compute the Lagrange interpolating polynomial.
  approximate_taylor_polynomial -- compute an approximate Taylor polynomial for
               a function using polynomial interpolation
"""

postpone_import = 1
