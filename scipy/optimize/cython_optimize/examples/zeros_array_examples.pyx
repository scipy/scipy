from __future__ import division, print_function, absolute_import
from math import exp, sin

from .. cimport zeros_array

DEF N = 6  # number of arguments

NUM_OF_IRRAD = 10
IL = [sin(il) + 6.0 for il in range(NUM_OF_IRRAD)]


# governing equations
cdef double f_solarcell(double* args):
    cdef double i = args[0]
    cdef double v = args[1]
    cdef double il = args[2]
    cdef double io = args[3]
    cdef double rs = args[4]
    cdef double rsh = args[5]
    cdef double vt = args[6]
    vd = v + i * rs
    return il - io * (exp(vd / vt) - 1.0) - vd / rsh - i


cdef double fprime(double* args):
    cdef double i = args[0]
    cdef double v = args[1]
    cdef double il = args[2]
    cdef double io = args[3]
    cdef double rs = args[4]
    cdef double rsh = args[5]
    cdef double vt = args[6]
    return -io * exp((v + i * rs) / vt) * rs / vt - rs / rsh - 1


cdef double fprime2(double* args):
    cdef double i = args[0]
    cdef double v = args[1]
    cdef double il = args[2]
    cdef double io = args[3]
    cdef double rs = args[4]
    cdef double rsh = args[5]
    cdef double vt = args[6]
    return -io * exp((v + i * rs) / vt) * (rs / vt)**2


#solver
cdef double solarcell_newton(tuple args):
    """test newton with array"""
    cdef int n = N
    cdef double[N] myargs = args
    return zeros_array.newton(f_solarcell, 6.0, fprime, n, myargs)


# cython
def test_cython_newton(v=5.25, il=IL, args=(1e-09, 0.004, 10, 0.27456)):
    return map(solarcell_newton,
               ((v, il_,) + args for il_ in il))


#solver
cdef double solarcell_bisect(tuple args):
    """test newton with array"""
    cdef int n = N
    cdef double[N] myargs = args
    return zeros_array.bisect(f_solarcell, 7, 0, n, myargs, 0.001, 0.001, 10)


# cython
def test_cython_bisect(v=5.25, il=IL, args=(1e-09, 0.004, 10, 0.27456)):
    return map(solarcell_bisect,
               ((v, il_,) + args for il_ in il))
