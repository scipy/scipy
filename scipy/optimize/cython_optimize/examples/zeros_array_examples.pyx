"""
Use an array to hold callback args so there's no Python and functions can free
the global interpreter lock.
"""

from __future__ import division, print_function, absolute_import
from math import exp, sin

from .. cimport zeros_array

DEF N = 6  # number of arguments

NUM_OF_IRRAD = 10
IL = [sin(il) + 6.0 for il in range(NUM_OF_IRRAD)]
ARGS = (1e-09, 0.004, 10.0, 0.27456)
TOL, MAXITER = 1.48e-8, 50


# governing equations
cdef double f_solarcell(int n, double* args):
    cdef double i = args[0]
    cdef double v = args[1]
    cdef double il = args[2]
    cdef double io = args[3]
    cdef double rs = args[4]
    cdef double rsh = args[5]
    cdef double vt = args[6]
    vd = v + i * rs
    return il - io * (exp(vd / vt) - 1.0) - vd / rsh - i


cdef double fprime(int n, double* args):
    cdef double i = args[0]
    cdef double v = args[1]
    cdef double il = args[2]
    cdef double io = args[3]
    cdef double rs = args[4]
    cdef double rsh = args[5]
    cdef double vt = args[6]
    return -io * exp((v + i * rs) / vt) * rs / vt - rs / rsh - 1


cdef double fprime2(int n, double* args):
    cdef double i = args[0]
    cdef double v = args[1]
    cdef double il = args[2]
    cdef double io = args[3]
    cdef double rs = args[4]
    cdef double rsh = args[5]
    cdef double vt = args[6]
    return -io * exp((v + i * rs) / vt) * (rs / vt)**2


# cython newton solver
cdef double solarcell_newton(tuple args):
    """test newton with array"""
    cdef int n = N
    cdef double[N] myargs = args
    return zeros_array.newton(f_solarcell, 6.0, fprime, n, myargs, TOL, MAXITER)


# test cython newton solver in a loop
def test_cython_newton(v=5.25, il=IL, args=ARGS):
    return map(solarcell_newton, ((v, il_,) + args for il_ in il))


# cython secant solver
cdef double solarcell_secant(tuple args):
    """test secant with array"""
    cdef int n = N
    cdef double[N] myargs = args
    return zeros_array.secant(f_solarcell, 6.0, n, myargs, TOL, MAXITER)


# test cython secant solver in a loop
def test_cython_secant(v=5.25, il=IL, args=ARGS):
    return map(solarcell_secant, ((v, il_,) + args for il_ in il))


# cython halley solver
cdef double solarcell_halley(tuple args):
    """test halley with array"""
    cdef int n = N
    cdef double[N] myargs = args
    return zeros_array.halley(f_solarcell, 6.0, fprime, n, myargs, TOL, MAXITER, fprime2)


# test cython halley solver in a loop
def test_cython_halley(v=5.25, il=IL, args=ARGS):
    return map(solarcell_halley, ((v, il_,) + args for il_ in il))


# cython bisect solver
cdef double solarcell_bisect(tuple args):
    """test newton with array"""
    cdef int n = N
    cdef double[N] myargs = args
    return zeros_array.bisect(f_solarcell, 7, 0, n, myargs, 0.001, 0.001, 10)


# test cython bisect in a loop
def test_cython_bisect(v=5.25, il=IL, args=ARGS):
    return map(solarcell_bisect, ((v, il_,) + args for il_ in il))
