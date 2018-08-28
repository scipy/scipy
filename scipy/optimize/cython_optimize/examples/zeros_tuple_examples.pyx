"""
Use a tuple to hold callback args so Python calls functions but can't free the
global interpreter lock.
"""

from __future__ import division, print_function, absolute_import
from math import exp, sin

from .. cimport zeros_tuple

NUM_OF_IRRAD = 10
IL = [sin(il) + 6.0 for il in range(NUM_OF_IRRAD)]
ARGS = (1e-09, 0.004, 10.0, 0.27456)
TOL, MAXITER = 1.48e-8, 50
XTOL, RTOL, MITR = 0.001, 0.001, 10

DEF SIGNERR = -1
DEF CONVERR = -2

ctypedef struct scipy_newton_full_output:
    int funcalls
    int iterations
    int error_num
    double root
    char* flag


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
    return zeros_tuple.newton(f_solarcell, 6.0, fprime, args, TOL, MAXITER, NULL)


# test cython newton solver in a loop
def test_cython_newton(v=5.25, il=IL, args=ARGS):
    """test newton with array"""
    return map(solarcell_newton, ((v, il_,) + args for il_ in il))


# cython newton solver with full output
cdef scipy_newton_full_output solarcell_newton_full_output(tuple args, int tol, int maxiter):
    cdef scipy_newton_full_output full_output
    full_output.root = zeros_tuple.newton(f_solarcell, 6.0, fprime, args, tol, maxiter, <zeros_tuple.scipy_newton_parameters *> &full_output)
    if full_output.error_num == SIGNERR:
        full_output.flag = "TOL and MAXITER must be positive integers"
        full_output.funcalls = 0
        full_output.iterations = 0
    elif full_output.error_num == CONVERR:
        full_output.flag = "Failed to converge"
    else:
        full_output.error_num = 0
        full_output.flag = "Converged successfully"
    return full_output


def test_newton_full_output(v=5.25, il=6.0, args=ARGS, tol=TOL, maxiter=MAXITER):
    """test newton with full output"""
    return solarcell_newton_full_output((v, il,) + args, tol, maxiter)


# cython secant solver
cdef double solarcell_secant(tuple args):
    return zeros_tuple.secant(f_solarcell, 6.0, args, TOL, MAXITER, NULL)


# test cython secant solver in a loop
def test_cython_secant(v=5.25, il=IL, args=ARGS):
    """test secant with array"""
    return map(solarcell_secant, ((v, il_,) + args for il_ in il))


# cython halley solver
cdef double solarcell_halley(tuple args):
    return zeros_tuple.halley(f_solarcell, 6.0, fprime, args, TOL, MAXITER, fprime2, NULL)


# test cython halley solver in a loop
def test_cython_halley(v=5.25, il=IL, args=ARGS):
    """test halley with array"""
    return map(solarcell_halley, ((v, il_,) + args for il_ in il))


# cython bisect solver
cdef double solarcell_bisect(tuple args):
    return zeros_tuple.bisect(f_solarcell, 7.0, 0.0, args, XTOL, RTOL, MITR, NULL)


# test cython bisect in a loop
def test_cython_bisect(v=5.25, il=IL, args=ARGS):
    """test bisect with array"""
    return map(solarcell_bisect, ((v, il_,) + args for il_ in il))


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
