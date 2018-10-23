"""
Use an array to hold callback args so there's no Python and functions can free
the global interpreter lock.
"""

from __future__ import division, print_function, absolute_import
from math import sin

import cython
from libc.math cimport exp

from ... cimport zeros_array

DEF N = 6  # number of arguments

NUM_OF_IRRAD = 10
IL = [sin(il) + 6.0 for il in range(NUM_OF_IRRAD)]
ARGS = (1e-09, 0.004, 10.0, 0.27456)
XTOL, RTOL, MITR = 0.001, 0.001, 10


# callback functions

@cython.cdivision(True)
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


# cython bisect solver
cdef double solarcell_bisect(tuple args):
    cdef int n = N
    cdef double[N] myargs = args
    return zeros_array.bisect(f_solarcell, 7, 0, n, myargs, XTOL, RTOL, MITR, NULL)


# test cython bisect in a loop
def test_cython_bisect(v=5.25, il=IL, args=ARGS):
    """test bisect with array"""
    return map(solarcell_bisect, ((v, il_,) + args for il_ in il))


# cython ridder solver
cdef double solarcell_ridder(tuple args):
    cdef int n = N
    cdef double[N] myargs = args
    return zeros_array.ridder(f_solarcell, 7, 0, n, myargs, XTOL, RTOL, MITR, NULL)


# test cython ridder in a loop
def test_cython_ridder(v=5.25, il=IL, args=ARGS):
    """test ridder with array"""
    return map(solarcell_ridder, ((v, il_,) + args for il_ in il))


# cython brenth solver
cdef double solarcell_brenth(tuple args):
    cdef int n = N
    cdef double[N] myargs = args
    return zeros_array.brenth(f_solarcell, 7, 0, n, myargs, XTOL, RTOL, MITR, NULL)


# test cython brenth in a loop
def test_cython_brenth(v=5.25, il=IL, args=ARGS):
    """test brenth with array"""
    return map(solarcell_brenth, ((v, il_,) + args for il_ in il))


# cython brentq solver
cdef double solarcell_brentq(tuple args):
    cdef int n = N
    cdef double[N] myargs = args
    return zeros_array.brentq(f_solarcell, 7, 0, n, myargs, XTOL, RTOL, MITR, NULL)


# test cython brentq in a loop
def test_cython_brentq(v=5.25, il=IL, args=ARGS):
    """test brentq with array"""
    return map(solarcell_brentq, ((v, il_,) + args for il_ in il))
