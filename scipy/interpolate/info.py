"""
Interpolation Tools
===================

Wrappers around FITPACK functions
----------------------------------

  splrep
          find smoothing spline given (x,y) points on curve.
  splprep
          find smoothing spline given parametrically defined curve.
  splev
          evaluate the spline or its derivatives.
  splint
          compute definite integral of a spline.
  sproot
          find the roots of a cubic spline.
  spalde
          compute all derivatives of a spline at given points.
  bisplrep
          find bivariate smoothing spline representation.
  bisplev
          evaluate bivariate smoothing spline.

  UnivariateSpline
          A more recent, object-oriented wrapper; finds a (possibly
          smoothed) interpolating spline.

  InterpolatedUnivariateSpline

  LSQUnivariateSpline

  BivariateSpline
           A more recent, object-oriented wrapper; finds a
           interpolating spline for a bivariate function.

  SmoothBivariateSpline

Low-level Piece-wise Spline Tools
-----------------------------------
  splmake
              Create a spline representation from data-points
              where the internal knots are the data-points.
              
  spleval
              Evaluate a spline representation on a new set of
              input data values.

  spline
              Single-call interface to splmake and spleval
              
  spltopp
              Return piecewise polynomial representation from a
              spline representation. 
                
  

Interpolation Classes (univariate)
-----------------------------------

  interp1d
               Create a class whose instances can linearly interpolate
               to compute unknown values of a univariate function.

  BarycentricInterpolator
               Compute with a numerically-stable version
               of the Lagrange interpolating polynomial.

  barycentric_interpolate
               procedural interface to the above

  KroghInterpolator
               Compute with the Hermite interpolating polynomial
               (allows the specification of derivatives at some points).

  krogh_interpolate
               procedural interface to the above

  PiecewisePolynomial
               Spline that is specified by giving positions and
               derivatives at every knot; allows high orders and
               efficient appending.

  piecewise_polynomial_interpolate
               procedural interface to the above

  ppform
               Class to create a piecewise polynomial representation of
               a spline from the coefficients of the polynomial in each
               section and the break-points

Interpolation Classes (multivariate)
-------------------------------------

  interp2d
               Create a class whose instances can interpolate
               to compute unknown values of a bivariate function.

  Rbf
               Apply Radial Basis Functions to interpolate scattered N-D data.

Additional tools
-----------------

  lagrange
               Compute the Lagrange interpolating polynomial.

  approximate_taylor_polynomial
               compute an approximate Taylor polynomial for
               a function using polynomial interpolation


See Also
------------

  ndimage
     map_coordinates
               N-d interpolation from evenly-spaced data using
               fast B-splines. 

     spline_filter
               Method to pre-compute spline coefficients to make
               map_coordinates efficient for muliple calls on the same
               set of interpolated points.

  signal
     resample

              Perform sinc-interpolation using a Fourier filter.
              This function can decimate or interpolate to an evenly
              sampled grid.
          
     bspline
     gauss_spline
     qspline1d
     cspline1d
     qspline1d_eval
     cspline1d_eval
     qspline2d
     cspline2d
               Low-level spline tools for regularly spaced data using
               fast B-spline algorithms.

"""

postpone_import = 1
