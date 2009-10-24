"""
Test Scipy functions versus mpmath, if available.

"""
import numpy as np
from numpy.testing import *
from test_data import Data
import scipy.special as sc

try:
    import mpmath
    def mpmath_check(func): return func
except ImportError:
    try:
        import sympy.mpmath as mpmath
        def mpmath_check(func): return func
    except ImportError:
        def mpmath_check(func):
            return dec.skipif(True, "mpmath library not found")(func)

#------------------------------------------------------------------------------

@mpmath_check
def test_expi_complex():
    dataset = []
    for r in np.logspace(-99, 2, 10):
        for p in np.linspace(0, 2*np.pi, 30):
            z = r*np.exp(1j*p)
            dataset.append((z, mpmath.ei(z)))
    dataset = np.array(dataset, dtype=np.complex_)

    Data(sc.expi, dataset, 0, 1).check()
