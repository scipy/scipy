# -*-cython-*-
#
# Implementation of the Lambert W function [1]. Based on the MPMath 
# implementation [2], and documentaion [3].
#
# Copyright: Yosef Meller, 2009
# Author email: mellerf@netvision.net.il
# 
# Distributed under the same license as SciPy
#
# References:
# [1] On the Lambert W function, Adv. Comp. Math. 5 (1996) 329-359,
#     available online: http://www.apmaths.uwo.ca/~djeffrey/Offprints/W-adv-cm.pdf
# [2] mpmath source code, Subversion revision 990
#     http://code.google.com/p/mpmath/source/browse/trunk/mpmath/functions.py?spec=svn994&r=992
# [3] mpmath source code, Subversion revision 994
#     http://code.google.com/p/mpmath/source/browse/trunk/mpmath/function_docs.py?spec=svn994&r=994

# NaN checking as per suggestions of the cython-users list,
# http://groups.google.com/group/cython-users/browse_thread/thread/ff03eed8221bc36d

# TODO: use a series expansion when extremely close to the branch point
# at `-1/e` and make sure that the proper branch is chosen there

import cython

cimport sf_error

cdef extern from "math.h":
    double exp(double x) nogil
    double log(double x) nogil

from _complexstuff cimport *

# Heavy lifting is here:

@cython.cdivision(True)
cdef inline double complex lambertw_scalar(double complex z, long k, double tol) nogil:
    """
    This is just the implementation of W for a single input z.
    See the docstring for lambertw() below for the full description.
    """
    # Comments copied verbatim from [2] are marked with '>'
    if zisnan(z):
        return z

    # Return value:
    cdef double complex w
    
    #> We must be extremely careful near the singularities at -1/e and 0
    cdef double u
    u = exp(-1)
    
    cdef double absz
    absz = zabs(z)
    if absz <= u:
        if z == 0:
            #> w(0,0) = 0; for all other branches we hit the pole
            if k == 0:
                return z
            sf_error.error("lambertw", sf_error.SINGULAR, NULL)
            return -inf
        
        if k == 0:
            w = z # Initial guess for iteration
        #> For small real z < 0, the -1 branch beaves roughly like log(-z)
        elif k == -1 and z.imag ==0 and z.real < 0:
            w = log(-z.real)
        #> Use a simple asymptotic approximation.
        else:
            w = zlog(z)
            #> The branches are roughly logarithmic. This approximation
            #> gets better for large |k|; need to check that this always
            #> works for k ~= -1, 0, 1.
            if k: w = w + k*2*pi*1j
    
    elif k == 0 and z.imag and zabs(z) <= 0.7:
        #> Both the W(z) ~= z and W(z) ~= ln(z) approximations break
        #> down around z ~= -0.5 (converging to the wrong branch), so patch
        #> with a constant approximation (adjusted for sign)
        if zabs(z+0.5) < 0.1:
            if z.imag > 0:
                w = 0.7 + 0.7j
            else:
                w = 0.7 - 0.7j
        else:
            w = z
    
    else:
        if z.real == inf:
            if k == 0:
                return z
            else:
                return z + 2*k*pi*1j
        
        if z.real == -inf:
            return (-z) + (2*k+1)*pi*1j
                
        #> Simple asymptotic approximation as above
        w = zlog(z)
        if k: w = w + k*2*pi*1j

    #> Use Halley iteration to solve w*exp(w) = z
    cdef double complex ew, wew, wewz, wn
    cdef int i
    for i in range(100):
        ew = zexp(w)
        wew = w*ew
        wewz = wew-z
        wn = w - wewz / (wew + ew - (w + 2)*wewz/(2*w + 2))
        if zabs(wn-w) < tol*zabs(wn):
            return wn
        else:
            w = wn

    sf_error.error("lambertw", sf_error.SLOW,
                   "iteration failed to converge: %g + %gj",
                   <double>z.real, <double>z.imag)
    return nan
