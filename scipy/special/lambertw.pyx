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

cdef void lambertw_raise_warning(double complex z) with gil:
    warnings.warn("Lambert W iteration failed to converge: %r" % z)

# Heavy lifting is here:

@cython.cdivision(True)
cdef double complex lambertw_scalar(double complex z, long k, double tol) nogil:
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


# Turn the above function into a Ufunc:
#--------------------------------------
cdef extern from "numpy/arrayobject.h":
    void import_array()
    ctypedef int npy_intp
    cdef enum NPY_TYPES:
        NPY_LONG
        NPY_CDOUBLE
        NPY_DOUBLE

cdef extern from "numpy/ufuncobject.h":
    void import_ufunc()
    ctypedef void (*PyUFuncGenericFunction)(char**, npy_intp*, npy_intp*, void*)
    object PyUFunc_FromFuncAndData(PyUFuncGenericFunction* func, void** data,
        char* types, int ntypes, int nin, int nout,
        int identity, char* name, char* doc, int c)

cdef void _apply_func_to_1d_vec(char **args, npy_intp *dimensions, npy_intp *steps,
                     void *func) nogil:
    cdef npy_intp i
    cdef char *ip1=args[0], *ip2=args[1], *ip3=args[2], *op=args[3]
    for i in range(0, dimensions[0]):
        (<double complex*>op)[0] = (<double complex(*)(double complex, long, double) nogil>func)(
            (<double complex*>ip1)[0], (<long*>ip2)[0], (<double*>ip3)[0])
        ip1 += steps[0]; ip2 += steps[1]; ip3 += steps[2]; op += steps[3]

cdef PyUFuncGenericFunction _loop_funcs[1]
_loop_funcs[0] = _apply_func_to_1d_vec

cdef char _inp_outp_types[4]
_inp_outp_types[0] = NPY_CDOUBLE
_inp_outp_types[1] = NPY_LONG
_inp_outp_types[2] = NPY_DOUBLE
_inp_outp_types[3] = NPY_CDOUBLE

import_array()
import_ufunc()

# The actual ufunc declaration:
cdef void *the_func_to_apply[1]
the_func_to_apply[0] = <void*>lambertw_scalar
_lambertw = PyUFunc_FromFuncAndData(_loop_funcs, the_func_to_apply,
    _inp_outp_types, 1, 3, 1, 0, "", "", 0)

def lambertw(z, k=0, tol=1e-8):
    r"""
    lambertw(z, k=0, tol=1e-8)

    Lambert W function.

    The Lambert W function `W(z)` is defined as the inverse function
    of ``w * exp(w)``. In other words, the value of ``W(z)`` is
    such that ``z = W(z) * exp(W(z))`` for any complex number
    ``z``.

    The Lambert W function is a multivalued function with infinitely
    many branches. Each branch gives a separate solution of the
    equation ``w exp(w)``. Here, the branches are indexed by the
    integer `k`.
    
    Parameters
    ----------
    z : array_like
        Input argument.
    k : int, optional
        Branch index.
    tol : float, optional
        Evaluation tolerance.

    Notes
    -----
    All branches are supported by `lambertw`:

    * ``lambertw(z)`` gives the principal solution (branch 0)
    * ``lambertw(z, k)`` gives the solution on branch `k`

    The Lambert W function has two partially real branches: the
    principal branch (`k = 0`) is real for real ``z > -1/e``, and the
    ``k = -1`` branch is real for ``-1/e < z < 0``. All branches except
    ``k = 0`` have a logarithmic singularity at ``z = 0``.

    **Possible issues**
    
    The evaluation can become inaccurate very close to the branch point
    at ``-1/e``. In some corner cases, `lambertw` might currently
    fail to converge, or can end up on the wrong branch.

    **Algorithm**

    Halley's iteration is used to invert ``w * exp(w)``, using a first-order
    asymptotic approximation (O(log(w)) or `O(w)`) as the initial estimate.

    The definition, implementation and choice of branches is based on [1]_.
    
    TODO: use a series expansion when extremely close to the branch point
    at ``-1/e`` and make sure that the proper branch is chosen there

    References
    ----------
    .. [1] Corless et al, "On the Lambert W function", Adv. Comp. Math. 5
       (1996) 329-359.
       http://www.apmaths.uwo.ca/~djeffrey/Offprints/W-adv-cm.pdf

    Examples
    --------
    The Lambert W function is the inverse of ``w exp(w)``:

    >>> from scipy.special import lambertw
    >>> w = lambertw(1)
    >>> w
    0.56714329040978387299996866221035555
    >>> w*exp(w)
    1.0

    Any branch gives a valid inverse:

    >>> w = lambertw(1, k=3)
    >>> w
    (-2.8535817554090378072068187234910812 +
    17.113535539412145912607826671159289j)
    >>> w*exp(w)
    (1.0 + 3.5075477124212226194278700785075126e-36j)

    **Applications to equation-solving**

    The Lambert W function may be used to solve various kinds of
    equations, such as finding the value of the infinite power
    tower :math:`z^{z^{z^{\ldots}}}`:

    >>> def tower(z, n):
    ... if n == 0:
    ... return z
    ... return z ** tower(z, n-1)
    ...
    >>> tower(0.5, 100)
    0.641185744504986
    >>> -lambertw(-log(0.5))/log(0.5)
    0.6411857445049859844862004821148236665628209571911

    **Properties**

    The Lambert W function grows roughly like the natural logarithm
    for large arguments:

    >>> lambertw(1000)
    5.2496028524016
    >>> log(1000)
    6.90775527898214
    >>> lambertw(10**100)
    224.843106445119
    >>> log(10**100)
    230.258509299405
    
    The principal branch of the Lambert W function has a rational
    Taylor series expansion around `z = 0`::
    
    >>> nprint(taylor(lambertw, 0, 6), 10)
    [0.0, 1.0, -1.0, 1.5, -2.666666667, 5.208333333, -10.8]
    
    Some special values and limits are:
    
    >>> lambertw(0)
    0.0
    >>> lambertw(1)
    0.567143290409784
    >>> lambertw(e)
    1.0
    >>> lambertw(inf)
    +inf
    >>> lambertw(0, k=-1)
    -inf
    >>> lambertw(0, k=3)
    -inf
    >>> lambertw(inf, k=3)
    (+inf + 18.8495559215388j)

    The `k = 0` and `k = -1` branches join at `z = -1/e` where
    `W(z) = -1` for both branches. Since `-1/e` can only be represented
    approximately with mpmath numbers, evaluating the Lambert W function
    at this point only gives `-1` approximately:

    >>> lambertw(-1/e, 0)
    -0.999999999999837133022867
    >>> lambertw(-1/e, -1)
    -1.00000000000016286697718
    
    If `-1/e` happens to round in the negative direction, there might be
    a small imaginary part:
    
    >>> lambertw(-1/e)
    (-1.0 + 8.22007971511612e-9j)

    """
    return _lambertw(z, k, tol)

