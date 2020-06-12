"""Compute a grid of values for Wright's generalized Bessel function
and save the values to data files for use in tests. Using mpmath directly in
tests would take too long.

This takes about 30 minutes to run on a 2.7 GHz Macbook Pro.
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
def rgamma_cached(x):
    return mp.rgamma(x)


def mp_wright_bessel(a, b, x, dps=50, maxterms=1000):
    """Compute Wright' generalized Bessel function as Series with mpmath.
    """
    with mp.workdps(dps):
        a, b, x = mp.mpf(a), mp.mpf(b), mp.mpf(x)
        res = mp.nsum(lambda k: x**k / mp.fac(k) * rgamma_cached(a * k + b),
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
                        1e-3 - eps, 1e-3, 1e-3 + eps,
                        1e-4 - eps, 1e-4, 1e-4 + eps])
    b_range = np.array([0, eps, 1e-10, 1e-5, 0.1, 1, 2, 10, 100])
    x_range = np.array([0, eps, 1 - eps, 1, 1 + eps,
                        1.5,
                        2 - eps, 2, 2 + eps,
                        10 * (1 - eps), 10, 10 * (1 + eps),
                        100 * (1 - eps), 100, 100 * (1 + eps),
                        600,
                        709.78271289338402,  # largest value np(x) is not inf
                        1e3, 1e5, 1e10, 1e20])

    a_range, b_range, x_range = np.meshgrid(a_range, b_range, x_range)
    # filter out some values
    b_filter = ~((a_range < 1e-6) & (x_range >= exp_inf))
    b_filter = b_filter | ((a_range < 1e-4) & (x_range >= 700))

    # and flatten
    a_range = a_range[b_filter].flatten()
    b_range = b_range[b_filter].flatten()
    x_range = x_range[b_filter].flatten()

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
