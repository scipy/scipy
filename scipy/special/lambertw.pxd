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
import warnings

cdef extern from "math.h":
    double exp(double x) nogil
    double log(double x) nogil

# Use Numpy's portable C99-compatible complex functios

cdef extern from "numpy/npy_math.h":
    ctypedef struct npy_cdouble:
        double real
        double imag

    double npy_cabs(npy_cdouble z) nogil
    npy_cdouble npy_clog(npy_cdouble z) nogil
    npy_cdouble npy_cexp(npy_cdouble z) nogil
    int npy_isnan(double x) nogil
    double NPY_INFINITY
    double NPY_PI

cdef inline bint zisnan(double complex x) nogil:
    return npy_isnan(x.real) or npy_isnan(x.imag)

cdef inline double zabs(double complex x) nogil:
    cdef double r
    r = npy_cabs((<npy_cdouble*>&x)[0])
    return r

cdef inline double complex zlog(double complex x) nogil:
    cdef npy_cdouble r
    r = npy_clog((<npy_cdouble*>&x)[0])
    return (<double complex*>&r)[0]

cdef inline double complex zexp(double complex x) nogil:
    cdef npy_cdouble r
    r = npy_cexp((<npy_cdouble*>&x)[0])
    return (<double complex*>&r)[0]

cdef inline void lambertw_raise_warning(double complex z) with gil:
    warnings.warn("Lambert W iteration failed to converge: %r" % z)

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
            return -NPY_INFINITY
        
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
            if k: w = w + k*2*NPY_PI*1j
    
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
        if z.real == NPY_INFINITY:
            if k == 0:
                return z
            else:
                return z + 2*k*NPY_PI*1j
        
        if z.real == -NPY_INFINITY:
            return (-z) + (2*k+1)*NPY_PI*1j
                
        #> Simple asymptotic approximation as above
        w = zlog(z)
        if k: w = w + k*2*NPY_PI*1j

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

    lambertw_raise_warning(z)
    return wn
