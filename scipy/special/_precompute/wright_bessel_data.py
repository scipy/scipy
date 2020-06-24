"""Compute a grid of values for Wright's generalized Bessel function
and save the values to data files for use in tests. Using mpmath directly in
tests would take too long.

This takes about 10 minutes to run on a 2.7 GHz i7 Macbook Pro.
"""
from functools import lru_cache
import os
from time import time

import numpy as np
from scipy.special._mptestutils import mpf2float

try:
    import mpmath as mp
except ImportError:
    pass

# exp_inf: smallest value x for which exp(x) == inf
exp_inf = 709.78271289338403


# 64 Byte per value
@lru_cache(maxsize=100_000)
def rgamma_cached(x, dps):
    with mp.workdps(dps):
        return mp.rgamma(x)


def mp_wright_bessel(a, b, x, dps=50, maxterms=1000):
    """Compute Wright' generalized Bessel function as Series with mpmath.
    """
    with mp.workdps(dps):
        a, b, x = mp.mpf(a), mp.mpf(b), mp.mpf(x)
        res = mp.nsum(lambda k: x**k / mp.fac(k)
                      * rgamma_cached(a * k + b, dps=dps),
                      [0, mp.inf],
                      tol=dps, method='s', steps=[maxterms]
                      )
        return mpf2float(res)


def main():
    t0 = time()
    print(__doc__)
    pwd = os.path.dirname(__file__)
    eps = np.finfo(float).eps * 100

    a_range = np.array([eps,
                        1e-4 * (1 - eps), 1e-4, 1e-4 * (1 + eps),
                        1e-3 * (1 - eps), 1e-3, 1e-3 * (1 + eps),
                        1 * (1 - eps), 1, 1 * (1 + eps),
                        2, 4.999, 5, 10])
    b_range = np.array([0, eps, 1e-10, 1e-5, 0.1, 1, 2, 10, 100])
    x_range = np.array([0, eps, 1 - eps, 1, 1 + eps,
                        1.5,
                        2 - eps, 2, 2 + eps,
                        10 * (1 - eps), 10, 10 * (1 + eps),
                        100 * (1 - eps), 100, 100 * (1 + eps),
                        500,
                        exp_inf,
                        1e3, 1e5, 1e10, 1e20])

    a_range, b_range, x_range = np.meshgrid(a_range, b_range, x_range)

    # filter out some values, especially too large x
    bool_filter = ~((a_range < 5e-3) & (x_range >= exp_inf))
    bool_filter = bool_filter & ~((a_range < 0.2) & (x_range > exp_inf))
    bool_filter = bool_filter & ~((a_range < 0.5) & (x_range > 1e3))
    bool_filter = bool_filter & ~((a_range < 0.56) & (x_range > 5e3))
    bool_filter = bool_filter & ~((a_range < 1) & (x_range > 1e4))
    bool_filter = bool_filter & ~((a_range < 1.4) & (x_range > 1e5))
    bool_filter = bool_filter & ~((a_range < 1.8) & (x_range > 1e6))
    bool_filter = bool_filter & ~((a_range < 2.2) & (x_range > 1e7))
    bool_filter = bool_filter & ~((a_range < 2.5) & (x_range > 1e8))
    bool_filter = bool_filter & ~((a_range < 2.9) & (x_range > 1e9))
    bool_filter = bool_filter & ~((a_range < 3.3) & (x_range > 1e10))
    bool_filter = bool_filter & ~((a_range < 3.7) & (x_range > 1e11))
    bool_filter = bool_filter & ~((a_range < 4) & (x_range > 1e12))
    bool_filter = bool_filter & ~((a_range < 4.4) & (x_range > 1e13))
    bool_filter = bool_filter & ~((a_range < 4.7) & (x_range > 1e14))
    bool_filter = bool_filter & ~((a_range < 5.1) & (x_range > 1e15))
    bool_filter = bool_filter & ~((a_range < 5.4) & (x_range > 1e16))
    bool_filter = bool_filter & ~((a_range < 5.8) & (x_range > 1e17))
    bool_filter = bool_filter & ~((a_range < 6.2) & (x_range > 1e18))
    bool_filter = bool_filter & ~((a_range < 6.2) & (x_range > 1e18))
    bool_filter = bool_filter & ~((a_range < 6.5) & (x_range > 1e19))
    bool_filter = bool_filter & ~((a_range < 6.9) & (x_range > 1e20))

    # filter out known values that do not meet the required numerical accuracy
    # see test test_wright_data_grid_failures
    bool_filter = bool_filter & ~((np.abs(a_range - 1) <= 1e-3)
                                  & (b_range == 10) & (x_range >= 500))
    bool_filter = bool_filter & ~((a_range == 1) & (b_range == 100)
                                  & (x_range >= 100000))
    bool_filter = bool_filter & ~((a_range >= 1) & (a_range <= 2)
                                  & (b_range == 100) & (x_range >= 100000))

    # and flatten
    a_range = a_range[bool_filter].flatten()
    b_range = b_range[bool_filter].flatten()
    x_range = x_range[bool_filter].flatten()

    dataset = []
    print(f"Computing {x_range.size} single points.")
    for i in range(x_range.size):
        a = a_range.flatten()[i]
        b = b_range.flatten()[i]
        x = x_range.flatten()[i]
        # take care of difficult corner cases
        maxterms = 1000
        if a < 1e-6 and x >= exp_inf/10:
            maxterms = 2000
        f = mp_wright_bessel(a, b, x, maxterms=maxterms)
        dataset.append((a, b, x, f))
    dataset = np.array(dataset)

    filename = os.path.join(pwd, '..', 'tests', 'data', 'local',
                            'wright_bessel.txt')
    np.savetxt(filename, dataset)

    print("{:.1f} minutes elapsed".format((time() - t0)/60))


if __name__ == "__main__":
    main()
