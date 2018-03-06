from __future__ import division, print_function, absolute_import
from math import exp, sin

from scipy.optimize.cython_optimize cimport zeros

NUM_OF_IRRAD = 10
IL = [sin(il) + 6.0 for il in range(NUM_OF_IRRAD)]


# governing equations
cdef double f_solarcell(double i, tuple args):
    v, il, io, rs, rsh, vt = args
    vd = v + i * rs
    return il - io * (exp(vd / vt) - 1.0) - vd / rsh - i


cdef double fprime(double i, tuple args):
    v, il, io, rs, rsh, vt = args
    return -io * exp((v + i * rs) / vt) * rs / vt - rs / rsh - 1


#solver
cdef double solarcell_newton(tuple args):
    """test newton with array"""
    return zeros.newton(f_solarcell, 6.0, fprime, args)


# cython
def test_cython_newton(v=5.25, il=IL, args=(1e-09, 0.004, 10, 0.27456)):
    return map(solarcell_newton,
               ((v, il_,) + args for il_ in il))


#solver
cdef double solarcell_bisect(tuple args):
    """test newton with array"""
    return zeros.bisect(f_solarcell, 7, 0, args, 0.001, 0.001, 10)


# cython
def test_cython_bisect(v=5.25, il=IL, args=(1e-09, 0.004, 10, 0.27456)):
    return map(solarcell_bisect,
               ((v, il_,) + args for il_ in il))
