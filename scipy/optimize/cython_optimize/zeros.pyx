from __future__ import division, print_function, absolute_import

import warnings


cdef float newton(callback_type func, float x0, callback_type fprime, tuple args):
    tol=1.48e-8
    maxiter=50,
    fprime2=None
    if tol <= 0:
        raise ValueError("tol too small (%g <= 0)" % tol)
    if maxiter < 1:
        raise ValueError("maxiter must be greater than 0")
    if fprime:  # is not None:
        # Newton-Rapheson method
        # Multiply by 1.0 to convert to floating point.  We don't use float(x0)
        # so it still works if x0 is complex.
        p0 = 1.0 * x0  # asarray(x0)  # convert to ndarray
        for iter in range(maxiter):
            myargs = (p0,) + args
            fder = fprime(p0, args)  # *myargs)  # asarray(fprime(*myargs))  # convert to ndarray
            if fder == 0:  # (fder == 0).any():
                msg = "derivative was zero."
                warnings.warn(msg, RuntimeWarning)
                return p0
            fval = func(p0, args)  # *myargs)  # asarray(func(*myargs))  # convert to ndarray
            if fprime2 is None:
                # Newton step
                p = p0 - fval / fder
            else:
                fder2 = fprime2(*myargs)  # asarray(fprime2(*myargs))  # convert to ndarray
                # Halley's method
                # https://en.wikipedia.org/wiki/Halley%27s_method
                p = p0 - 2 * fval * fder / (2 * fder ** 2 - fval * fder2)
            if abs(p - p0) < tol:  # np_abs(p - p0).max() < tol:
                return p
            p0 = p
    else:
        pass
        # # Secant method
        # p0 = x0  # asarray(x0)
        # dx = 1e-4  # finfo(float).eps ** 0.33
        # dp = dx if x0 >= 0 else -dx # where(p0 >= 0, dx, -dx)
        # p1 = p0 * (1 + dx) + dp
        # q0 = func(*((p0,) + args))  # asarray(func(*((p0,) + args)))
        # q1 = func(*((p1,) + args))  # asarray(func(*((p1,) + args)))
        # for iter in range(maxiter):
        #     # divide_by_zero = (q1 == q0)
        #     if q1 == q0:  # divide_by_zero.any():
        #         # tolerance_reached = (p1 != p0)
        #         if p1 != p0:  # (divide_by_zero & tolerance_reached).any():
        #             msg = "Tolerance of %s reached" % p1 - p0  # sqrt(sum((p1 - p0) ** 2))
        #             warnings.warn(msg, RuntimeWarning)
        #         return (p1 + p0) / 2.0
        #     else:
        #         p = p1 - q1 * (p1 - p0) / (q1 - q0)
        #     if abs(p - p1) < tol:  # np_abs(p - p1).max() < tol:
        #         return p
        #     p0 = p1
        #     q0 = q1
        #     p1 = p
        #     q1 = func(*((p1,) + args))  # asarray(func(*((p1,) + args)))
    msg = "Failed to converge after %d iterations, value is %s" % (maxiter, p)
    raise RuntimeError(msg)
