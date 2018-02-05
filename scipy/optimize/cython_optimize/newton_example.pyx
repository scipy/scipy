from __future__ import division, print_function, absolute_import
from math import exp, sin

from scipy.optimize.cython_optimize cimport zeros

NUM_OF_IRRAD = 100000
IL = [sin(il) + 6.0 for il in range(NUM_OF_IRRAD)]


# governing equations
cdef float f_solarcell(float i, tuple args):
    v, il, io, rs, rsh, vt = args
    vd = v + i * rs
    return il - io * (exp(vd / vt) - 1.0) - vd / rsh - i


cdef float fprime(float i, tuple args):
    v, il, io, rs, rsh, vt = args
    return -io * exp((v + i * rs) / vt) * rs / vt - rs / rsh - 1


# def fprime2(i, v, il, io, rs, rsh, vt):
#     return -io * exp((v + i * rs) / vt) * (rs / vt)**2


#solvers
cdef float solarcell_newton(tuple args):
    """test newton with array"""
    return zeros.newton(f_solarcell, 7.0, fprime, args)


# def solarcell_halley(args):
#     """test newton with array"""
#     return zeros.newton(f_solarcell, 7.0, fprime, args, fprime2=fprime2)
#
#
# def solarcell_secant(args):
#     """test newton with array"""
#     return zeros.newton(f_solarcell, 7.0, args=args)


# cython
def bench_cython_newton():
    return map(solarcell_newton,
                 ((5.25, il, 1e-09, 0.004, 10, 0.27456) for il in IL))


# def bench_cython_halley():
#     return map(solarcell_halley,
#                  ((5.25, il, 1e-09, 0.004, 10, 0.27456) for il in IL))
#
#
# def bench_cython_secant():
#     return map(solarcell_secant,
#                  ((5.25, il, 1e-09, 0.004, 10, 0.27456) for il in IL))
