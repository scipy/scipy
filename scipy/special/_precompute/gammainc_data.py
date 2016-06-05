"""
Compute gammainc for large arguments and parameters and save the
values in a data file for use in tests. We can't just compare to
mpmath's gammainc in test_mpmath.TestSystematic because it would take
too long.

Note that mpmath's gammainc is computed using hypercomb, but since it
doesn't allow the user to increase the maximum number of terms used in
the series it doesn't converge for many arguments. To get around this
we copy the mpmath implementation but use more terms.

This takes about 14 minutes to run on a 2.3 GHz Macbook Pro with 4GB
ram.

Sources:
[1] Fredrik Johansson and others. mpmath: a Python library for
    arbitrary-precision floating-point arithmetic (version 0.19),
    December 2013. http://mpmath.org/.

"""
from __future__ import division, print_function, absolute_import

import os
from time import time
import numpy as np
from numpy import pi

from scipy.special._mptestutils import mpf2float

try:
    import mpmath as mp
except ImportError:
    try:
        import sympy.mpmath as mp
    except ImportError:
        pass


def gammainc(a, x, dps=50, maxterms=10**8):
    """
    Compute gammainc exactly like mpmath does but allow for more
    summands in hypercomb. See

    mpmath/functions/expintegrals.py#L134
    
    in the mpmath github repository.

    """
    with mp.workdps(dps):
        z, a, b = mp.mpf(a), mp.mpf(x), mp.mpf(x)
        G = [z]
        negb = mp.fneg(b, exact=True)

        def h(z):
            T1 = [mp.exp(negb), b, z], [1, z, -1], [], G, [1], [1+z], b
            return (T1,)

        res = mp.hypercomb(h, [z], maxterms=maxterms)
        return mpf2float(res)


def main():
    # It would be nice to have data for larger values, but either this
    # requires prohibitively large precision (dps > 800) or mpmath has
    # a bug. For example, gammainc(1e20, 1e20, dps=800) returns a
    # value around 0.03, while the true value should be close to 0.5
    # (DLMF 8.12.15).
    rmax = 14
    t0 = time()
    print(__doc__)
    # Region where 0.6 <= x/a <= 1. The transition to the asymptotic
    # series begins at x/a = 0.7.
    r = np.logspace(4, rmax, 30)
    theta = np.logspace(np.log10(pi/4), np.log10(np.arctan(0.6)), 30)
    r, theta = np.meshgrid(r, theta)
    a, x = r*np.cos(theta), r*np.sin(theta)
    a, x = a.flatten(), x.flatten()
    dataset = []
    for i, (a0, x0) in enumerate(zip(a, x)):
        dataset.append((a0, x0, gammainc(a0, x0)))
    dataset = np.array(dataset)

    fn = os.path.join(os.path.dirname(__file__), '..', 'tests',
                      'data', 'local', 'gammainc.txt')
    np.savetxt(fn, dataset)
    print("{} minutes elapsed".format((time() - t0)/60))


if __name__ == "__main__":
    main()
