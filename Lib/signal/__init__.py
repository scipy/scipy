""" Signal Processing Tools: A collection of signal processing tools:

 Convolution:
 
    convolve      --  N-dimensional convolution.
    correlate     --  N-dimensional correlation.
    convolve2d    --  2-dimensional convolution (more options).
    correlate2d   --  2-dimensional correlation (more options).
    sepfir2d      --  Convolve with a 2-D separable FIR filter.

 B-splines:
    
    bspline       --  B-spline basis function of order n.
    gauss_spline  --  Gaussian approximation to the B-spline basis function.
    cspline1d     --  Coefficients for 1-D cubic (3rd order) B-spline.
    qspline1d     --  Coefficients for 1-D quadratic (2nd order) B-spline.
    cspline2d     --  Coefficients for 2-D cubic (3rd order) B-spline.
    qspline2d     --  Coefficients for 2-D quadratic (2nd order) B-spline.
    spline_filter --  Smoothing spline (cubic) filtering of a rank-2 array. 

 Filtering:

    order_filter  --  N-dimensional order filter.
    medfilt       --  N-dimensional median filter.
    medfilt2      --  2-dimensional median filter (faster).
    wiener        --  N-dimensional wiener filter.

    symiirorder1  --  2nd-order IIR filter (cascade of first-order systems).
    symiirorder2  --  4th-order IIR filter (cascade of second-order systems).
    lfilter       --  1-dimensional FIR and IIR digital linear filtering.

 Filter design:
 
    remez         --  Optimal FIR filter design.
    iirdesign     --- IIR filter design given bands and gains
    iirfilter     --- IIR filter design given order and critical frequencies
    
 
"""
_modules = ['sigtools']
_namespaces = ['signaltools', 'bsplines', 'filter_design']

__all__ = []

import scipy
scipy.modules2all(__all__, _modules, globals())
scipy.names2all(__all__, _namespaces, globals())
del scipy
