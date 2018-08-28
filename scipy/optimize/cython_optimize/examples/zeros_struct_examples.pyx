"""
Use a structure to hold callback args so there's no Python and functions can
free the global interpreter lock.
"""

from __future__ import division, print_function, absolute_import
from math import exp, sin

from .. cimport zeros_struct

# test parameter structure
ctypedef struct test_params:
    double voltage
    double light_current
    double dark_current
    double series_resistance
    double shunt_resistance
    double thermal_voltage

NUM_OF_IRRAD = 10
IL = [sin(il) + 6.0 for il in range(NUM_OF_IRRAD)]
ARGS = (1e-09, 0.004, 10.0, 0.27456)
TOL, MAXITER = 1.48e-8, 50
XTOL, RTOL, MITR = 0.001, 0.001, 10


# callabck function
cdef double f_solarcell(double i, void *args):
    cdef test_params *myargs = <test_params *> args
    cdef double v, il, io, rs, rsh, vt
    # unpack structure
    v = myargs.voltage
    il = myargs.light_current
    io = myargs.dark_current
    rs = myargs.series_resistance
    rsh = myargs.shunt_resistance
    vt = myargs.thermal_voltage
    vd = v + i * rs
    return il - io * (exp(vd / vt) - 1.0) - vd / rsh - i


cdef double fprime(double i, void *args):
    cdef test_params *myargs = <test_params *> args
    cdef double v, il, io, rs, rsh, vt
    # unpack structure
    v = myargs.voltage
    il = myargs.light_current
    io = myargs.dark_current
    rs = myargs.series_resistance
    rsh = myargs.shunt_resistance
    vt = myargs.thermal_voltage
    return -io * exp((v + i * rs) / vt) * rs / vt - rs / rsh - 1


cdef double fprime2(double i, void *args):
    cdef test_params *myargs = <test_params *> args
    cdef double v, il, io, rs, rsh, vt
    # unpack structure
    v = myargs.voltage
    il = myargs.light_current
    io = myargs.dark_current
    rs = myargs.series_resistance
    rsh = myargs.shunt_resistance
    vt = myargs.thermal_voltage
    return -io * exp((v + i * rs) / vt) * (rs / vt)**2


# cython newton solver
cdef double solarcell_newton(tuple args):
    """test newton with array"""
    cdef test_params myargs
    myargs.voltage = args[0]
    myargs.light_current = args[1]
    myargs.dark_current = args[2]
    myargs.series_resistance = args[3]
    myargs.shunt_resistance = args[4]
    myargs.thermal_voltage = args[5]
    return zeros_struct.newton(f_solarcell, 6.0, fprime, <test_params *> &myargs, TOL, MAXITER, NULL)


# test cython newton solver in a loop
def test_cython_newton(v=5.25, il=IL, args=ARGS):
    return map(solarcell_newton, ((v, il_,) + args for il_ in il))


# cython secant solver
cdef double solarcell_secant(tuple args):
    """test secant with array"""
    cdef test_params myargs
    myargs.voltage = args[0]
    myargs.light_current = args[1]
    myargs.dark_current = args[2]
    myargs.series_resistance = args[3]
    myargs.shunt_resistance = args[4]
    myargs.thermal_voltage = args[5]
    return zeros_struct.secant(f_solarcell, 6.0, <test_params *> &myargs, TOL, MAXITER, NULL)


# test cython secant solver in a loop
def test_cython_secant(v=5.25, il=IL, args=ARGS):
    return map(solarcell_secant, ((v, il_,) + args for il_ in il))


# cython halley solver
cdef double solarcell_halley(tuple args):
    """test halley with array"""
    cdef test_params myargs
    myargs.voltage = args[0]
    myargs.light_current = args[1]
    myargs.dark_current = args[2]
    myargs.series_resistance = args[3]
    myargs.shunt_resistance = args[4]
    myargs.thermal_voltage = args[5]
    return zeros_struct.halley(f_solarcell, 6.0, fprime, <test_params *> &myargs, TOL, MAXITER, fprime2, NULL)


# test cython halley solver in a loop
def test_cython_halley(v=5.25, il=IL, args=ARGS):
    return map(solarcell_halley, ((v, il_,) + args for il_ in il))


# cython bisect solver
cdef double solarcell_bisect(tuple args):
    """test bisect with array"""
    cdef test_params myargs
    myargs.voltage = args[0]
    myargs.light_current = args[1]
    myargs.dark_current = args[2]
    myargs.series_resistance = args[3]
    myargs.shunt_resistance = args[4]
    myargs.thermal_voltage = args[5]
    return zeros_struct.bisect(f_solarcell, 7.0, 0.0, <test_params *> &myargs, XTOL, RTOL, MITR, NULL)


# test cython bisect in a loop
def test_cython_bisect(v=5.25, il=IL, args=ARGS):
    return map(solarcell_bisect, ((v, il_,) + args for il_ in il))


# cython ridder solver
cdef double solarcell_ridder(tuple args):
    """test ridder with array"""
    cdef test_params myargs
    myargs.voltage = args[0]
    myargs.light_current = args[1]
    myargs.dark_current = args[2]
    myargs.series_resistance = args[3]
    myargs.shunt_resistance = args[4]
    myargs.thermal_voltage = args[5]
    return zeros_struct.ridder(f_solarcell, 7.0, 0.0, <test_params *> &myargs, XTOL, RTOL, MITR, NULL)


# test cython ridder in a loop
def test_cython_ridder(v=5.25, il=IL, args=ARGS):
    return map(solarcell_ridder, ((v, il_,) + args for il_ in il))


# cython brenth solver
cdef double solarcell_brenth(tuple args):
    """test brenth with array"""
    cdef test_params myargs
    myargs.voltage = args[0]
    myargs.light_current = args[1]
    myargs.dark_current = args[2]
    myargs.series_resistance = args[3]
    myargs.shunt_resistance = args[4]
    myargs.thermal_voltage = args[5]
    return zeros_struct.brenth(f_solarcell, 7.0, 0.0, <test_params *> &myargs, XTOL, RTOL, MITR, NULL)


# test cython brenth in a loop
def test_cython_brenth(v=5.25, il=IL, args=ARGS):
    return map(solarcell_brenth, ((v, il_,) + args for il_ in il))


# cython brentq solver
cdef double solarcell_brentq(tuple args):
    """test brentq with array"""
    cdef test_params myargs
    myargs.voltage = args[0]
    myargs.light_current = args[1]
    myargs.dark_current = args[2]
    myargs.series_resistance = args[3]
    myargs.shunt_resistance = args[4]
    myargs.thermal_voltage = args[5]
    return zeros_struct.brentq(f_solarcell, 7.0, 0.0, <test_params *> &myargs, XTOL, RTOL, MITR, NULL)


# test cython brentq in a loop
def test_cython_brentq(v=5.25, il=IL, args=ARGS):
    return map(solarcell_brentq, ((v, il_,) + args for il_ in il))
