from __future__ import division, print_function, absolute_import

import warnings

TOL = 1.48e-8
MAXITER = 50


cdef double newton(callback_type func, double p0, callback_type fprime, tuple args):
    # Newton-Rapheson method
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
