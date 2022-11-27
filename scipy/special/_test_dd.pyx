# This module provides Python wrappers for a few of the "double-double"
# C functions defined in cephes/dd_*.  The wrappers are not part of the
# public API; they are for use in scipy.special unit tests only.

import numpy as np
from numpy.testing import assert_allclose


cdef extern from "cephes/dd_real.h":
    cdef struct double2:
        pass
    double2 dd_create(double, double)
    double dd_hi(double2)
    double dd_lo(double2)
    double2 dd_exp(const double2 a)
    double2 dd_log(const double2 a)
    double2 dd_expm1(const double2 a)


def _dd_exp(double xhi, double xlo):
    cdef double2 x = dd_create(xhi, xlo)
    cdef double2 y = dd_exp(x)
    return dd_hi(y), dd_lo(y)


def _dd_log(double xhi, double xlo):
    cdef double2 x = dd_create(xhi, xlo)
    cdef double2 y = dd_log(x)
    return dd_hi(y), dd_lo(y)


def _dd_expm1(double xhi, double xlo):
    cdef double2 x = dd_create(xhi, xlo)
    cdef double2 y = dd_expm1(x)
    return dd_hi(y), dd_lo(y)
