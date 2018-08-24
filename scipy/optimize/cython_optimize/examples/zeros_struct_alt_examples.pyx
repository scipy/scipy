"""
Use a structure to hold callback args so there's no Python and functions can
free the global interpreter lock.
"""

from __future__ import division, print_function, absolute_import
from math import exp, sin

from .. cimport zeros_struct
from . cimport zeros_struct_examples
from .zeros_struct_examples import IL, TOL, MAXITER, XTOL, RTOL, MITR

ARGS = {'dark_current': 1e-09, 'series_resistance': 0.004,
    'shunt_resistance': 10.0, 'thermal_voltage': 0.27456}

ctypedef struct test_params:
    double voltage
    double light_current
    double dark_current
    double series_resistance
    double shunt_resistance
    double thermal_voltage


# cython newton solver
cdef double solarcell_newton(dict args):
    """test newton with dictionary"""
    cdef test_params myargs
    myargs = args
    return zeros_struct.newton(zeros_struct_examples.f_solarcell, 6.0, zeros_struct_examples.fprime, <test_params *> &myargs, TOL, MAXITER)


# test cython newton solver in a loop
def test_cython_newton(v=5.25, il=IL, args=ARGS):
    return map(solarcell_newton,
               (dict(voltage=v, light_current=il_, **args) for il_ in il))


# cython bisect solver
cdef double solarcell_bisect(dict args):
    """test bisect with dictionary"""
    cdef test_params myargs
    myargs = args
    return zeros_struct.bisect(zeros_struct_examples.f_solarcell, 7.0, 0.0, <test_params *> &myargs, XTOL, RTOL, MITR)


# test cython bisect solver in a loop
def test_cython_bisect(v=5.25, il=IL, args=ARGS):
    return map(solarcell_bisect,
               (dict(voltage=v, light_current=il_, **args) for il_ in il))
