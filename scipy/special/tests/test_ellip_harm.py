#
# Tests for the Ellipsoidal Harmonic Function,
# Distributed under the same license as SciPy itself.
#

from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_equal
from scipy.special import ellip_harm
from numpy import array, sqrt

from scipy.special._testutils import FuncData

def E01(h2, k2, s):
    return 1

def E11(h2, k2, s):
    return s

def E12(h2, k2, s):
    return sqrt(abs(s*s - h2))

def E13(h2, k2, s):
    return sqrt(abs(s*s - k2))

def E21(h2, k2, s):
    return s*s - 1/3*((h2 + k2) + sqrt(abs((h2 + k2)*(h2 + k2)-3*h2*k2)))

def E22(h2, k2, s):
    return s*s - 1/3*((h2 + k2) - sqrt(abs((h2 + k2)*(h2 + k2)-3*h2*k2)))

def E23(h2, k2, s):
    return s * sqrt(abs(s*s - h2))

def E24(h2, k2, s):
    return s * sqrt(abs(s*s - k2))

def E25(h2, k2, s):
    return sqrt(abs((s*s - h2)*(s*s - k2)))

def E31(h2, k2, s):
    return s*s*s - (s/5)*(2*(h2 + k2) + sqrt(4*(h2 + k2)*(h2 + k2) - 15*h2*k2))

def E32(h2, k2, s):
    return s*s*s - (s/5)*(2*(h2 + k2) - sqrt(4*(h2 + k2)*(h2 + k2) - 15*h2*k2))

def E33(h2, k2, s):
    return sqrt(abs(s*s - h2))*(s*s - 1/5*((h2 + 2*k2) + sqrt(abs((h2 + 2*k2)*(h2 + 2*k2) - 5*h2*k2))))

def E34(h2, k2, s):
    return sqrt(abs(s*s - h2))*(s*s - 1/5*((h2 + 2*k2) - sqrt(abs((h2 + 2*k2)*(h2 + 2*k2) - 5*h2*k2))))

def E35(h2, k2, s):
    return sqrt(abs(s*s - k2))*(s*s - 1/5*((2*h2 + k2) + sqrt(abs((2*h2 + k2)*(2*h2 + k2) - 5*h2*k2))))

def E36(h2, k2, s):
    return sqrt(abs(s*s - k2))*(s*s - 1/5*((2*h2 + k2) - sqrt(abs((2*h2 + k2)*(2*h2 + k2) - 5*h2*k2))))

def E37(h2, k2, s):
    return s * sqrt(abs((s*s - h2)*(s*s - k2)))

def test_values():
   
    assert_equal(ellip_harm(5,8,1,2,2.5,1,1), ellip_harm(5,8,1,2,2.5))

    data = [
        (5,8,0,1,2.5,1,1, E01(5, 8, 2.5)),
        (5,8,1,1,2.5,1,1, E11(5, 8, 2.5)),
        (5,8,1,2,2.5,1,1, E12(5, 8, 2.5)),
        (5,8,1,3,2.5,1,1, E13(5, 8, 2.5)),
        (5,8,2,1,2.5,1,1, E21(5, 8, 2.5)),
        (5,8,2,2,2.5,1,1, E22(5, 8, 2.5)),
        (5,8,2,3,2.5,1,1, E23(5, 8, 2.5)),
        (5,8,2,4,2.5,1,1, E24(5, 8, 2.5)),
        (5,8,2,5,2.5,1,1, E25(5, 8, 2.5)),
        (5,8,3,1,2.5,1,1, E31(5, 8, 2.5)),
        (5,8,3,2,2.5,1,1, E32(5, 8, 2.5)),
        (5,8,3,3,2.5,1,1, E33(5, 8, 2.5)),
        (5,8,3,4,2.5,1,1, E34(5, 8, 2.5)),
        (5,8,3,5,2.5,1,1, E35(5, 8, 2.5)),
        (5,8,3,6,2.5,1,1, E36(5, 8, 2.5)),
        (5,8,3,7,2.5,1,1, E37(5, 8, 2.5)),
        (6.25,36,0,1,7.3,1,1, E01(6.25, 36, 7.3)),
        (6.25,36,1,1,7.3,1,1, E11(6.25, 36, 7.3)),
        (6.25,36,1,2,7.3,1,1, E12(6.25, 36, 7.3)),
        (6.25,36,1,3,7.3,1,1, E13(6.25, 36, 7.3)),
        (6.25,36,2,1,7.3,1,1, E21(6.25, 36, 7.3)),
        (6.25,36,2,2,7.3,1,1, E22(6.25, 36, 7.3)),
        (6.25,36,2,3,7.3,1,1, E23(6.25, 36, 7.3)),
        (6.25,36,2,4,7.3,1,1, E24(6.25, 36, 7.3)),
        (6.25,36,2,5,7.3,1,1, E25(6.25, 36, 7.3)),
        (6.25,36,3,1,7.3,1,1, E31(6.25, 36, 7.3)),
        (6.25,36,3,2,7.3,1,1, E32(6.25, 36, 7.3)),
        (6.25,36,3,3,7.3,1,1, E33(6.25, 36, 7.3)),
        (6.25,36,3,4,7.3,1,1, E34(6.25, 36, 7.3)),
        (6.25,36,3,5,7.3,1,1, E35(6.25, 36, 7.3)),
        (6.25,36,3,6,7.3,1,1, E36(6.25, 36, 7.3)),
        (6.25,36,3,7,7.3,1,1, E37(6.25, 36, 7.3)),
        ]

    data = array(data, dtype=float)

    def w(a, b, c, d, e, f, g):
        return ellip_harm(a, b, c, d, e, f, g)
    olderr = np.seterr(all='ignore')
    try:
        FuncData(w, data, (0,1,2,3,4,5,6), 7, rtol=1e-10, atol=1e-13).check()  
    finally:
        np.seterr(**olderr)

