"""
Use a tuple to hold callback args so Python calls functions but can't free the
global interpreter lock.
"""

from __future__ import division, print_function, absolute_import
from math import sin

import cython
from libc.math cimport exp

from .. cimport zeros_tuple

NUM_OF_IRRAD = 10
IL = [sin(il) + 6.0 for il in range(NUM_OF_IRRAD)]
ARGS = (1e-09, 0.004, 10.0, 0.27456)
XTOL, RTOL, MITR = 0.001, 0.001, 10

DEF SIGNERR = -1
DEF CONVERR = -2

ctypedef struct scipy_zeros_full_output:
    int funcalls
    int iterations
    int error_num
    double root
    char* flag


# callback functions

@cython.cdivision(True)
cdef double f_solarcell(double i, tuple args):
    v, il, io, rs, rsh, vt = args
    vd = v + i * rs
    return il - io * (exp(vd / vt) - 1.0) - vd / rsh - i


# cython bisect solver
cdef double solarcell_bisect(tuple args):
    return zeros_tuple.bisect(f_solarcell, 7.0, 0.0, args, XTOL, RTOL, MITR, NULL)


# test cython bisect in a loop
def test_cython_bisect(v=5.25, il=IL, args=ARGS):
    """test bisect with array"""
    return map(solarcell_bisect, ((v, il_,) + args for il_ in il))


# cython bisect solver with full output
cdef scipy_zeros_full_output solarcell_bisect_full_output(tuple args, double xa, double xb, double xtol, double rtol, int mitr):
    cdef scipy_zeros_full_output full_output
    full_output.root = zeros_tuple.bisect(f_solarcell, xa, xb, args, xtol, rtol, mitr, <zeros_tuple.scipy_zeros_parameters *> &full_output)
    if full_output.error_num == SIGNERR:
        full_output.flag = "F(XA) and F(XB) must have opposite signs"
        full_output.funcalls = 0
        full_output.iterations = 0
    elif full_output.error_num == CONVERR:
        full_output.flag = "Failed to converge"
    else:
        full_output.error_num = 0
        full_output.flag = "Converged successfully"
    return full_output


def test_bisect_full_output(v=5.25, il=6.0, args=ARGS, xa=7.0, xb=0.0, xtol=XTOL, rtol=RTOL, mitr=MITR):
    """test bisect with full output"""
    return solarcell_bisect_full_output((v, il,) + args, xa, xb, xtol, rtol, mitr)


# cython ridder solver
cdef double solarcell_ridder(tuple args):
    return zeros_tuple.ridder(f_solarcell, 7.0, 0.0, args, XTOL, RTOL, MITR, NULL)


# test cython ridder in a loop
def test_cython_ridder(v=5.25, il=IL, args=ARGS):
    """test ridder with array"""
    return map(solarcell_ridder, ((v, il_,) + args for il_ in il))


# cython brenth solver
cdef double solarcell_brenth(tuple args):
    return zeros_tuple.brenth(f_solarcell, 7.0, 0.0, args, XTOL, RTOL, MITR, NULL)


# test cython brenth in a loop
def test_cython_brenth(v=5.25, il=IL, args=ARGS):
    """test brenth with array"""
    return map(solarcell_brenth, ((v, il_,) + args for il_ in il))


# cython brentq solver
cdef double solarcell_brentq(tuple args):
    return zeros_tuple.brentq(f_solarcell, 7.0, 0.0, args, XTOL, RTOL, MITR, NULL)


# test cython brentq in a loop
def test_cython_brentq(v=5.25, il=IL, args=ARGS):
    """test brentq with array"""
    return map(solarcell_brentq, ((v, il_,) + args for il_ in il))


# cython brentq solver with full output
cdef scipy_zeros_full_output solarcell_brentq_full_output(tuple args, double xa, double xb, double xtol, double rtol, int mitr):
    cdef scipy_zeros_full_output full_output
    full_output.root = zeros_tuple.brentq(f_solarcell, xa, xb, args, xtol, rtol, mitr, <zeros_tuple.scipy_zeros_parameters *> &full_output)
    if full_output.error_num == SIGNERR:
        full_output.flag = "F(XA) and F(XB) must have opposite signs"
        full_output.funcalls = 0
        full_output.iterations = 0
    elif full_output.error_num == CONVERR:
        full_output.flag = "Failed to converge"
    else:
        full_output.error_num = 0
        full_output.flag = "Converged successfully"
    return full_output


def test_brentq_full_output(v=5.25, il=6.0, args=ARGS, xa=7.0, xb=0.0, xtol=XTOL, rtol=RTOL, mitr=MITR):
    """test brentq with full output"""
    return solarcell_brentq_full_output((v, il,) + args, xa, xb, xtol, rtol, mitr)
