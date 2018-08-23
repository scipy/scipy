from __future__ import division, print_function, absolute_import
import warnings
import cython
from . cimport c_zeros_array

DEF MAXARGS = 10


@cython.cdivision(True)
cdef double newton(callback_type_array func, double p0, callback_type_array fprime, int n, double* args):
    # Newton-Rapheson method
    cdef int i
    cdef double fder, fval, p
    cdef double[MAXARGS] p0args
    for iter in range(MAXITER):
        p0args[0] = p0
        for i in range(n):
            p0args[i+1] = args[i]
        fder = fprime(n, p0args)
        if fder == 0:
            msg = "derivative was zero."
            warnings.warn(msg, RuntimeWarning)
            return p0
        fval = func(n, p0args)
        # Newton step
        p = p0 - fval / fder
        if abs(p - p0) < TOL:  # np_abs(p - p0).max() < tol:
            return p
        p0 = p
    msg = "Failed to converge after %d iterations, value is %s" % (MAXITER, p)
    raise RuntimeError(msg)


# callback function wrapper that extracts function, args from params struct
cdef double scipy_zeros_functions_func(double x, void *params):
    cdef c_zeros_array.scipy_zeros_parameters *myparams = <c_zeros_array.scipy_zeros_parameters *> params
    cdef int n = myparams.n, i
    cdef double[MAXARGS] myargs
    cdef callback_type_array f = myparams.function

    myargs[0] = x
    for i in range(n):
        myargs[i+1] = myparams.args[i]

    return f(myargs)


# cythonized way to call scalar bisect
cdef double bisect(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter):
    cdef c_zeros_array.scipy_zeros_parameters myparams
    # create params struct
    myparams.n = n
    myparams.args = args
    myparams.function = f
    return c_zeros_array.bisect(scipy_zeros_functions_func, xa, xb, xtol, rtol, iter, <c_zeros_array.default_parameters *> &myparams)
