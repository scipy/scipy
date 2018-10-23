"""
Use a structure to hold callback args so there's no Python and functions can
free the global interpreter lock.
"""

from __future__ import division, print_function, absolute_import

from ... cimport zeros_struct
from . cimport zeros_struct_examples
from .zeros_struct_examples import IL, XTOL, RTOL, MITR

ARGS = {'dark_current': 1e-09, 'series_resistance': 0.004,
    'shunt_resistance': 10.0, 'thermal_voltage': 0.27456}

ctypedef struct test_params:
    double voltage
    double light_current
    double dark_current
    double series_resistance
    double shunt_resistance
    double thermal_voltage


# cython bisect solver
cdef double solarcell_bisect(dict args):
    cdef test_params myargs = args
    return zeros_struct.bisect(zeros_struct_examples.f_solarcell, 7.0, 0.0,
        <test_params *> &myargs, XTOL, RTOL, MITR, NULL)


# test cython bisect solver in a loop
def test_cython_bisect(v=5.25, il=IL, args=ARGS):
    """test bisect with dictionary"""
    return map(solarcell_bisect,
               (dict(voltage=v, light_current=il_, **args) for il_ in il))
