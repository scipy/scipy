""" Signal Processing Tools: A collection of signal processing tools:

 Convolution:
 
    convolve      --  N-dimensional convolution.
    correlate     --  N-dimensional correlation.
    convolve2     --  2-dimensional convolution (faster).
    correlate2    --  2-dimensional correlation (faster).

 B-splines:
 
    cspline1d     --  Coefficients for 1-D cubic B-spline.
    bspline       --  B-spline basis function of order n.
    gauss_approx  --  Gaussian approximation to the B-spline basis function.
    sepfir2d      --  Convolve with a 2-D separable FIR filter.
    cspline2d     --  Coefficients for 2-D cubic B-spline.
    qspline2d     --  Coefficients for 2-D quadratic B-spline.

 Filtering:

    symiirorder1  --  2nd-order IIR filter (cascade of first-order systems).
    symiirorder2  --  4th-order IIR filter (cascade of second-order systems).
    spline_filter --  Smoothing spline (cubic) filtering of a rank-2 array.
    order_filter  --  N-dimensional order filter.
    medfilt       --  N-dimensional median filter.
    medfilt2      --  2-dimensional median filter (faster).
    wiener        --  N-dimensional wiener filter.
    remez         --  Optimal FIR filter design.
    lfilter       --  1-dimensional FIR and IIR digital linear filtering.

"""
from signaltools import *
from bsplines import *
