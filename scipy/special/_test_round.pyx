import numpy as np
from numpy.random import random_integers
from numpy.testing import assert_

cdef extern from "_round.h":
    double add_round_up(double, double) nogil
    double add_round_down(double, double) nogil
    int fesetround(int) nogil
    int fegetround() nogil
    int FE_UPWARD
    int FE_DOWNWARD

cdef extern from "numpy/npy_math.h":
    int npy_isnan(double) nogil


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


def test_add_round(size, mode):
    cdef:
        int i, old_round, status
        double[:] sample1 = random_integers(low=np.iinfo(np.int64).min,
                                            high=np.iinfo(np.int64).max,
                                            size=size).view(np.float64)
        double[:] sample2 = random_integers(low=np.iinfo(np.int64).min,
                                            high=np.iinfo(np.int64).max,
                                            size=size).view(np.float64)
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
        if npy_isnan(res) and npy_isnan(std):
            continue
        if res != std:
            nfail += 1
            msg.append("{:.21g} + {:.21g} = {:.21g} != {:.21g}"
                       .format(sample1[i], sample2[i], std, res))

    if nfail:
        s = "{}/{} failures with mode {}.".format(nfail, size, mode)
        msg = [s] + msg
        assert_(False, "\n".join(msg))
