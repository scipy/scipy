"""
Evaluate orthogonal polynomial values using recurrence relations
or by calling special functions.

References
----------

.. [AMS55] Abramowitz & Stegun, Section 22.5.

.. [MH] Mason & Handscombe, Chebyshev Polynomials, CRC Press (2003).

"""
#
# Copyright (C) 2009 Pauli Virtanen
# Distributed under the same license as Scipy.
#

#------------------------------------------------------------------------------
# Direct evaluation of polynomials
#------------------------------------------------------------------------------

cdef extern from "math.h":
    double sqrt(double x) nogil

cdef double eval_poly_chebyt(long k, double x) nogil:
    # Use Chebyshev T recurrence directly, see [MH]
    cdef long m
    cdef double b2, b1, b0

    b2 = 0
    b1 = -1
    b0 = 0
    x = 2*x
    for m in range(k+1):
        b2 = b1
        b1 = b0
        b0 = x*b1 - b2
    return (b0 - b2)/2.0

#------------------------------------------------------------------------------
# Ufunc boilerplate
#------------------------------------------------------------------------------

cdef extern from "numpy/arrayobject.h":
    void import_array()
    ctypedef int npy_intp
    cdef enum NPY_TYPES:
        NPY_LONG
        NPY_DOUBLE

cdef extern from "numpy/ufuncobject.h":
    void import_ufunc()
    ctypedef void (*PyUFuncGenericFunction)(char**, npy_intp*, npy_intp*, void*)
    object PyUFunc_FromFuncAndData(PyUFuncGenericFunction* func, void** data,
                                   char* types, int ntypes, int nin, int nout,
                                   int identity, char* name, char* doc, int c)

cdef void _loop_id_d(char **args, npy_intp *dimensions, npy_intp *steps,
                     void *func) nogil:
    cdef int i
    cdef double x
    cdef char *ip1=args[0], *ip2=args[1], *op=args[2]
    for i in range(0, dimensions[0]):
        (<double*>op)[0] = (<double(*)(long,double) nogil>func)(
            (<long*>ip1)[0], (<double*>ip2)[0])
        ip1 += steps[0]; ip2 += steps[1]; op += steps[2]

cdef char _id_d_types[3]

cdef PyUFuncGenericFunction _id_d_funcs[1]

_id_d_types[0] = NPY_LONG
_id_d_types[1] = NPY_DOUBLE
_id_d_types[2] = NPY_DOUBLE

_id_d_funcs[0] = _loop_id_d

import_array()
import_ufunc()

#--

cdef void *chebyt_data[1]
chebyt_data[0] = <void*>eval_poly_chebyt
_eval_chebyt = PyUFunc_FromFuncAndData(_id_d_funcs, chebyt_data,
                                       _id_d_types, 1, 2, 1, 0, "", "", 0)


#------------------------------------------------------------------------------
# Actual evaluation functions
#------------------------------------------------------------------------------

import numpy as np
from scipy.special._cephes import gamma, hyp2f1, hyp1f1, gammaln
from numpy import exp

def binom(n, k):
    """
    binom(n, k)

    Binomial coefficient
    """
    return np.exp(gammaln(1+n) - gammaln(1+k) - gammaln(1+n-k))

def eval_jacobi(n, alpha, beta, x, out=None):
    """
    eval_jacobi(n, alpha, beta, x, out=None)

    Evaluate Jacobi polynomial at a point.
    """
    d = binom(n+alpha, n)
    a = -n
    b = n + alpha + beta + 1
    c = alpha + 1
    g = (1-x)/2.0
    return hyp2f1(a, b, c, g) * d

def eval_sh_jacobi(n, p, q, x, out=None):
    """
    eval_sh_jacobi(n, p, q, x, out=None)

    Evaluate shifted Jacobi polynomial at a point.
    """
    factor = np.exp(gammaln(1+n) + gammaln(n+p) - gammaln(2*n+p))
    return factor * eval_jacobi(n, p-q, q-1, 2*x-1)

def eval_gegenbauer(n, alpha, x, out=None):
    """
    eval_gegenbauer(n, alpha, x, out=None)

    Evaluate Gegenbauer polynomial at a point.
    """
    d = gamma(n+2*alpha)/gamma(1+n)/gamma(2*alpha)
    a = -n
    b = n + 2*alpha
    c = alpha + 0.5
    g = (1-x)/2.0
    return hyp2f1(a, b, c, g) * d

def eval_chebyt(n, x, out=None):
    """
    eval_chebyt(n, x, out=None)

    Evaluate Chebyshev T polynomial at a point.

    This routine is numerically stable for `x` in ``[-1, 1]`` at least
    up to order ``10000``.
    """
    return _eval_chebyt(n, x, out)

def eval_chebyu(n, x, out=None):
    """
    eval_chebyu(n, x, out=None)

    Evaluate Chebyshev U polynomial at a point.
    """
    d = n+1
    a = -n
    b = n+2
    c = 1.5
    g = (1-x)/2.0
    return hyp2f1(a, b, c, g) * d

def eval_chebys(n, x, out=None):
    """
    eval_chebys(n, x, out=None)

    Evaluate Chebyshev S polynomial at a point.
    """
    return eval_chebyu(n, x/2, out=out)

def eval_chebyc(n, x, out=None):
    """
    eval_chebyc(n, x, out=None)

    Evaluate Chebyshev C polynomial at a point.
    """
    return 2*eval_chebyt(n, x/2.0, out)

def eval_sh_chebyt(n, x, out=None):
    """
    eval_sh_chebyt(n, x, out=None)

    Evaluate shifted Chebyshev T polynomial at a point.
    """
    return eval_chebyt(n, 2*x-1, out=out)

def eval_sh_chebyu(n, x, out=None):
    """
    eval_sh_chebyu(n, x, out=None)

    Evaluate shifted Chebyshev U polynomial at a point.
    """
    return eval_chebyu(n, 2*x-1, out=out)

def eval_legendre(n, x, out=None):
    """
    eval_legendre(n, x, out=None)

    Evaluate Legendre polynomial at a point.
    """
    d = 1
    a = -n
    b = n+1
    c = 1
    g = (1-x)/2.0
    return hyp2f1(a, b, c, g) * d

def eval_sh_legendre(n, x, out=None):
    """
    eval_sh_legendre(n, x, out=None)

    Evaluate shifted Legendre polynomial at a point.
    """
    return eval_legendre(n, 2*x-1, out=out)

def eval_genlaguerre(n, alpha, x, out=None):
    """
    eval_genlaguerre(n, alpha, x, out=None)

    Evaluate generalized Laguerre polynomial at a point.
    """
    d = binom(n+alpha, n)
    a = -n
    b = alpha + 1
    g = x
    return hyp1f1(a, b, g) * d

def eval_laguerre(n, x, out=None):
    """
    eval_laguerre(n, x, out=None)

    Evaluate Laguerre polynomial at a point.
    """
    return eval_genlaguerre(n, 0., x, out=out)

def eval_hermite(n, x, out=None):
    """
    eval_hermite(n, x, out=None)

    Evaluate Hermite polynomial at a point.
    """
    n, x = np.broadcast_arrays(n, x)
    n, x = np.atleast_1d(n, x)

    if out is None:
        out = np.zeros_like(0*n + 0*x)
    if (n % 1 != 0).any():
        raise ValueError("Order must be integer")

    even = (n % 2 == 0)

    m = n[even]/2
    out[even] = ((-1)**m * 2**(2*m) * gamma(1+m)
                 * eval_genlaguerre(m, -0.5, x[even]**2))

    m = (n[~even]-1)/2
    out[~even] = ((-1)**m * 2**(2*m+1) * gamma(1+m)
                  * x[~even] * eval_genlaguerre(m, 0.5, x[~even]**2))

    return out

def eval_hermitenorm(n, x, out=None):
    """
    eval_hermitenorm(n, x, out=None)

    Evaluate normalized Hermite polynomial at a point.
    """
    return eval_hermite(n, x/sqrt(2)) * 2**(-n/2.0)
