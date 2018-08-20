"""
Use a structure to hold callback args so there's no Python and functions can
free the global interpreter lock.
"""

from __future__ import division, print_function, absolute_import
from math import exp, sin

from .. cimport zeros_struct
from . cimport zeros_struct_examples
from .zeros_struct_examples import IL


#solver
cdef double solarcell_newton(dict args):
    """test newton with dictionary"""
    cdef zeros_struct_examples.test_params myargs
    myargs = args
    return zeros_struct.newton(zeros_struct_examples.f_solarcell, 6.0, zeros_struct_examples.fprime, <zeros_struct_examples.test_params *> &myargs)


# cython
def test_cython_newton(v=5.25, il=IL, args=(1e-09, 0.004, 10.0, 0.27456)):
    return map(solarcell_newton,
               ((v, il_,) + args for il_ in il))


#solver
cdef double solarcell_bisect(dict args):
    """test bisect with dictionary"""
    cdef zeros_struct_examples.test_params myargs
    myargs = args
    return zeros_struct.bisect(zeros_struct_examples.f_solarcell, 7.0, 0.0, <zeros_struct_examples.test_params *> &myargs, 0.001, 0.001, 10)


# cython
def test_cython_bisect(v=5.25, il=IL, args=(1e-09, 0.004, 10.0, 0.27456)):
    return map(solarcell_bisect,
               ((v, il_,) + args for il_ in il))
