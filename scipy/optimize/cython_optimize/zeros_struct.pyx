from __future__ import division, print_function, absolute_import
import warnings
import cython
from scipy.optimize.cython_optimize cimport c_zeros_struct

cdef double TOL = 1.48e-8
cdef int MAXITER = 50


# the new standard callback function that uses the params struct instead of tuple
@staticmethod
cdef double scipy_zeros_functions_func(double x, c_zeros.scipy_zeros_parameters *params):
    cdef c_zeros_struct.scipy_zeros_parameters *myparams
    cdef void* args
    cdef callback_type f

    myparams = params
    args = myparams.args
    f = myparams.function

    return f(x, args)  # recall callback_type takes a double and a tuple


@cython.cdivision(True)
cdef double newton(callback_type func, double p0, callback_type fprime, void *args):
    # Newton-Rapheson method
    cdef double fder, fval, p
    for iter in range(MAXITER):
        fder = fprime(p0, args)
        if fder == 0:
            msg = "derivative was zero."
            warnings.warn(msg, RuntimeWarning)
            return p0
        fval = func(p0, args)
        # Newton step
        p = p0 - fval / fder
        if abs(p - p0) < TOL:  # np_abs(p - p0).max() < tol:
            return p
        p0 = p
    msg = "Failed to converge after %d iterations, value is %s" % (MAXITER, p)
    raise RuntimeError(msg)


# cythonized way to call scalar bisect
cdef double bisect(callback_type f, double xa, double xb, void *args, double xtol, double rtol, int iter):
    cdef c_zeros.scipy_zeros_parameters myparams
    # create params struct
    myparams.args = args
    myparams.function = f
    return c_zeros.bisect(scipy_zeros_functions_func, xa, xb, xtol, rtol, iter, <c_zeros.default_parameters *> &myparams)
