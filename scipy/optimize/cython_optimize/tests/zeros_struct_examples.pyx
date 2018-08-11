"""
Use a structure to hold callback args so there's no Python and functions can
free the global interpreter lock.
"""

from __future__ import division, print_function, absolute_import
from math import exp, sin

from scipy.optimize.cython_optimize cimport zeros_struct

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


# governing equations
cdef double f_solarcell(double i, test_params *args):
    cdef test_params *myargs
    cdef double v, il, io, rs, rsh, vt
    # unpack structure
    myargs = args
    v = myargs.voltage
    il = myargs.light_current
    io = myargs.dark_current
    rs = myargs.series_resistance
    rsh = myargs.shunt_resistsance
    vt = myargs.thermal_voltage
    vd = v + i * rs
    return il - io * (exp(vd / vt) - 1.0) - vd / rsh - i


cdef double fprime(double i, test_params *args):
    cdef test_params *myargs
    cdef double v, il, io, rs, rsh, vt
    # unpack structure
    myargs = args
    v = myargs.voltage
    il = myargs.light_current
    io = myargs.dark_current
    rs = myargs.series_resistance
    rsh = myargs.shunt_resistsance
    vt = myargs.thermal_voltage
    return -io * exp((v + i * rs) / vt) * rs / vt - rs / rsh - 1


#solver
cdef double solarcell_newton(tuple args):
    """test newton with array"""
    cdef test_params *myargs
    myargs.voltage = args[0]
    myargs.light_current = args[1]
    myargs.dark_current = args[2]
    myargs.series_resistance = args[3]
    myargs.shunt_resistsance = args[4]
    myargs.thermal_voltage = args[5]
    return zeros_struct.newton(f_solarcell, 6.0, fprime, myargs)


# cython
def test_cython_newton(v=5.25, il=IL, args=(1e-09, 0.004, 10, 0.27456)):
    return map(solarcell_newton,
               ((v, il_,) + args for il_ in il))


#solver
cdef double solarcell_bisect(tuple args):
    """test newton with array"""
    cdef test_params *myargs
    myargs.voltage = args[0]
    myargs.light_current = args[1]
    myargs.dark_current = args[2]
    myargs.series_resistance = args[3]
    myargs.shunt_resistsance = args[4]
    myargs.thermal_voltage = args[5]
    return zeros_struct.bisect(f_solarcell, 7, 0, myargs, 0.001, 0.001, 10)


# cython
def test_cython_bisect(v=5.25, il=IL, args=(1e-09, 0.004, 10, 0.27456)):
    return map(solarcell_bisect,
               ((v, il_,) + args for il_ in il))
