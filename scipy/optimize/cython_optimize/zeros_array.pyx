from __future__ import division, print_function, absolute_import
import warnings
import cython
from . cimport c_zeros_array

cdef double TOL = 1.48e-8
cdef int MAXITER = 50


# the new standard callback function that uses the params struct instead of tuple
@staticmethod
cdef double scipy_zeros_functions_func(double x, c_zeros_array.scipy_zeros_parameters *params):
    cdef c_zeros_array.scipy_zeros_parameters *myparams
    cdef int n
    cdef double[] args
    cdef callback_type f

    myparams = params
    n = myparams.n
    args = myparams.args
    f = myparams.function

    cdef double[n] xargs
    xargs = [x] + args

    return f(n, xargs)


@cython.cdivision(True)
cdef double newton(callback_type_array func, double p0, callback_type_array fprime, int n, double[] args):
    # Newton-Rapheson method
    cdef double fder, fval, p
    cdef double[n] p0args
    for iter in range(MAXITER):
        p0args = [p0] + args
        fder = fprime(n, p0args)
        if fder == 0:
            msg = "derivative was zero."
            warnings.warn(msg, RuntimeWarning)
            return p0
        p0args = [p0] + args
        fval = func(n, p0args)
        # Newton step
        p = p0 - fval / fder
        if abs(p - p0) < TOL:  # np_abs(p - p0).max() < tol:
            return p
        p0 = p
    msg = "Failed to converge after %d iterations, value is %s" % (MAXITER, p)
    raise RuntimeError(msg)


# cythonized way to call scalar bisect
cdef double bisect(callback_type_array f, double xa, double xb, int n, double[] args, double xtol, double rtol, int iter):
    cdef c_zeros_struct.scipy_zeros_parameters myparams
    # create params struct
    myparmas.n = n
    myparams.args = args
    myparams.function = f
    return c_zeros_array.bisect(scipy_zeros_functions_func, xa, xb, xtol, rtol, iter, <c_zeros_struct.default_parameters *> &myparams)
