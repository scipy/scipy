from numpy import random
from numpy.testing import assert_

cdef extern from "_round.h":
    double add_round_up(double, double) nogil
    double add_round_down(double, double) nogil
    int fesetround(int) nogil
    int fegetround() nogil
    int FE_UPWARD
    int FE_DOWNWARD


def have_fenv():
    return True if fegetround() == 0 else False


def test_add_round(low, high, size, mode):
    cdef:
        int i, old_round
        double[:] sample1 = random.uniform(low=low, high=high, size=size)
        double[:] sample2 = random.uniform(low=low, high=high, size=size)
        double res, std

    nfail = 0
    msg = []
    for i in range(size):
        old_round = fegetround()
        try:
            if mode == 'up':
                res = add_round_up(sample1[i], sample2[i])
                fesetround(FE_UPWARD)
            elif mode == 'down':
                res = add_round_down(sample1[i], sample2[i])
                fesetround(FE_DOWNWARD)
            else:
                raise ValueError("Invalid rounding mode")
            std = sample1[i] + sample2[i]
        finally:
            fesetround(old_round)
        if res != std:
            nfail += 1
            msg.append("{:.21g} + {:.21g} = {:.21g} != {:.21g}"
                       .format(sample1[i], sample2[i], std, res))
    if nfail:
        s = "{}/{} failures with mode {}.".format(nfail, size, mode)
        msg = [s] + msg
        assert_(False, "\n".join(msg))
