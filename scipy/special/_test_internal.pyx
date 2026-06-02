# cython: cpow=True

"""
Wrappers to allow unit tests of internal C code.

This module includes wrappers for:
* Tests of the functions add_round_up() and add_round_down() from _round.h
* Several double-double functions from cephes/dd_real.c

"""

import numpy as np
from numpy.testing import assert_

from libc.math cimport isnan

cdef extern from "_round.h":
    double add_round_up(double, double) nogil
    double add_round_down(double, double) nogil
    int fesetround(int) nogil
    int fegetround() nogil
    int FE_UPWARD
    int FE_DOWNWARD


cdef extern from "dd_real_wrappers.h" nogil:
    cdef struct double2:
        double hi
        double lo
    double2 dd_create(double a, double b)
    double2 dd_exp(const double2* x)
    double2 dd_log(const double2* x)


def have_fenv():
    old_round = fegetround()
    have_getround = True if old_round >= 0 else False
    if not have_getround:
        return False

    have_setround = True
    try:
        if fesetround(FE_UPWARD) != 0:
            have_setround = False
        if fesetround(FE_DOWNWARD) != 0:
            have_setround = False
    finally:
        fesetround(old_round)
    return have_setround


def random_double(size, rng):
    # This code is a little hacky to work around some issues:
    # - randint doesn't have a dtype keyword until 1.11
    #   and the default dtype of randint is np.dtype(int)
    # - Something like
    #   >>> low = np.iinfo(np.dtype(int)).min
    #   >>> high = np.iinfo(np.dtype(int)).max + 1
    #   >>> np.random.randint(low=low, high=high)
    #   fails in NumPy 1.10.4 (note that the 'high' value in randint
    #   is exclusive); this is fixed in 1.11.
    x = rng.randint(low=0, high=2**16, size=4*size)
    return x.astype(np.uint16).view(np.float64)


def test_add_round(size, mode, rng):
    cdef:
        int i, old_round, status
        double[:] sample1 = random_double(size, rng)
        double[:] sample2 = random_double(size, rng)
        double res, std

    nfail = 0
    msg = []
    for i in range(size):
        old_round = fegetround()
        if old_round < 0:
            raise RuntimeError("Couldn't get rounding mode")
        try:
            if mode == 'up':
                res = add_round_up(sample1[i], sample2[i])
                status = fesetround(FE_UPWARD)
            elif mode == 'down':
                res = add_round_down(sample1[i], sample2[i])
                status = fesetround(FE_DOWNWARD)
            else:
                raise ValueError("Invalid rounding mode")
            if status != 0:
                raise RuntimeError("Failed to set rounding mode")
            std = sample1[i] + sample2[i]
        finally:
            fesetround(old_round)
        if isnan(res) and isnan(std):
            continue
        if res != std:
            nfail += 1
            msg.append("{:.21g} + {:.21g} = {:.21g} != {:.21g}"
                       .format(sample1[i], sample2[i], std, res))

    if nfail:
        s = "{}/{} failures with mode {}.".format(nfail, size, mode)
        msg = [s] + msg
        assert_(False, "\n".join(msg))


# Python wrappers for a few of the "double-double" C functions defined
# in cephes/dd_*.  The wrappers are not part of the public API; they are
# for use in scipy.special unit tests only.

def _dd_exp(double xhi, double xlo):
    cdef double2 x = dd_create(xhi, xlo)
    cdef double2 y = dd_exp(&x)
    return y.hi, y.lo


def _dd_log(double xhi, double xlo):
    cdef double2 x = dd_create(xhi, xlo)
    cdef double2 y = dd_log(&x)
    return y.hi, y.lo
