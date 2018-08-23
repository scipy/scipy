"""
Use a tuple to hold callback args so Python calls functions but can't free the
global interpreter lock.
"""

from __future__ import division, print_function, absolute_import
from math import exp, sin

from .. cimport zeros

NUM_OF_IRRAD = 10
IL = [sin(il) + 6.0 for il in range(NUM_OF_IRRAD)]
ARGS = (1e-09, 0.004, 10.0, 0.27456)
TOL, MAXITER = 1.48e-8, 50


# callback function
cdef double f_solarcell(double i, tuple args):
    v, il, io, rs, rsh, vt = args
    vd = v + i * rs
    return il - io * (exp(vd / vt) - 1.0) - vd / rsh - i


cdef double fprime(double i, tuple args):
    v, il, io, rs, rsh, vt = args
    return -io * exp((v + i * rs) / vt) * rs / vt - rs / rsh - 1


cdef double fprime2(double i, tuple args):
    v, il, io, rs, rsh, vt = args
    return -io * exp((v + i * rs) / vt) * (rs / vt)**2


# cython newton solver
cdef double solarcell_newton(tuple args):
    """test newton with array"""
    return zeros.newton(f_solarcell, 6.0, fprime, args, TOL, MAXITER)


# test cython newton solver in a loop
def test_cython_newton(v=5.25, il=IL, args=ARGS):
    return map(solarcell_newton, ((v, il_,) + args for il_ in il))


# cython secant solver
cdef double solarcell_secant(tuple args):
    """test secant with array"""
    return zeros.secant(f_solarcell, 6.0, args, TOL, MAXITER)


# test cython secant solver in a loop
def test_cython_secant(v=5.25, il=IL, args=ARGS):
    return map(solarcell_secant, ((v, il_,) + args for il_ in il))


# cython halley solver
cdef double solarcell_halley(tuple args):
    """test halley with array"""
    return zeros.halley(f_solarcell, 6.0, fprime, args, TOL, MAXITER, fprime2)


# test cython halley solver in a loop
def test_cython_halley(v=5.25, il=IL, args=ARGS):
    return map(solarcell_halley, ((v, il_,) + args for il_ in il))


# cython bisect solver
cdef double solarcell_bisect(tuple args):
    """test newton with array"""
    return zeros.bisect(f_solarcell, 7.0, 0.0, args, 0.001, 0.001, 10)


# test cython bisect in a loop
def test_cython_bisect(v=5.25, il=IL, args=ARGS):
    return map(solarcell_bisect, ((v, il_,) + args for il_ in il))
