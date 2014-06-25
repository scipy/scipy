#
# Tests for the Ellipsoidal Harmonic Function,
# Distributed under the same license as SciPy itself.
#

from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_equal, assert_almost_equal, assert_allclose
from scipy.special import ellip_harm, ellip_harm_2, ellip_normal
from scipy.integrate import quad
from numpy import array, sqrt, pi
from scipy.special._testutils import FuncData
        
def test_ellip_norm():

    def G01(h2, k2):
        return 4*pi

    def G11(h2, k2):
        return 4*pi*h2*k2/3

    def G12(h2, k2):
        return 4*pi*h2*(k2 - h2)/3

    def G13(h2, k2):
        return 4*pi*k2*(k2 - h2)/3

    def G22(h2, k2):
        res = 2*(h2**4 + k2**4) - 4*h2*k2*(h2**2 + k2**2) + 6*h2**2*k2**2 + sqrt(h2**2 + k2**2 - h2*k2)*(-2*(h2**3 + k2**3) + 3*h2*k2*(h2 + k2)) 
        return 16*pi/405*res

    def G21(h2, k2):
        res = 2*(h2**4 + k2**4) - 4*h2*k2*(h2**2 + k2**2) + 6*h2**2*k2**2 + sqrt(h2**2 + k2**2 - h2*k2)*(2*(h2**3 + k2**3) - 3*h2*k2*(h2 + k2)) 
        return 16*pi/405*res

    def G23(h2, k2):
        return 4*pi*h2**2*k2*(k2 - h2)/15

    def G24(h2, k2):
        return 4*pi*h2*k2**2*(k2 - h2)/15

    def G25(h2, k2):
        return 4*pi*h2*k2*(k2 - h2)**2/15

    def G32(h2, k2):
        res = 16*(h2**4 + k2**4) - 36*h2*k2*(h2**2 + k2**2) + 46*h2**2*k2**2 + sqrt(4*(h2**2 + k2**2) - 7*h2*k2)*(-8*(h2**3 + k2**3) + 11*h2*k2*(h2 + k2))
        return 16*pi/13125*k2*h2*res

    def G31(h2, k2):
        res = 16*(h2**4 + k2**4) - 36*h2*k2*(h2**2 + k2**2) + 46*h2**2*k2**2 + sqrt(4*(h2**2 + k2**2) - 7*h2*k2)*(8*(h2**3 + k2**3) - 11*h2*k2*(h2 + k2))
        return 16*pi/13125*h2*k2*res

    def G34(h2, k2):
        res = 6*h2**4 + 16*k2**4 - 12*h2**3*k2 - 28*h2*k2**3 + 34*h2**2*k2**2 + sqrt(h2**2 + 4*k2**2 - h2*k2)*(-6*h2**3 - 8*k2**3 + 9*h2**2*k2 + 13*h2*k2**2)
        return 16*pi/13125*h2*(k2 - h2)*res

    def G33(h2, k2):
        res = 6*h2**4 + 16*k2**4 - 12*h2**3*k2 - 28*h2*k2**3 + 34*h2**2*k2**2 + sqrt(h2**2 + 4*k2**2 - h2*k2)*(6*h2**3 + 8*k2**3 - 9*h2**2*k2 - 13*h2*k2**2)
        return 16*pi/13125*h2*(k2 - h2)*res

    def G36(h2, k2):
        res = 16*h2**4 + 6*k2**4 - 28*h2**3*k2 - 12*h2*k2**3 + 34*h2**2*k2**2 + sqrt(4*h2**2 + k2**2 - h2*k2)*(-8*h2**3 - 6*k2**3 + 13*h2**2*k2 + 9*h2*k2**2)
        return 16*pi/13125*k2*(k2 - h2)*res

    def G35(h2, k2):
        res = 16*h2**4 + 6*k2**4 - 28*h2**3*k2 - 12*h2*k2**3 + 34*h2**2*k2**2 + sqrt(4*h2**2 + k2**2 - h2*k2)*(8*h2**3 + 6*k2**3 - 13*h2**2*k2 - 9*h2*k2**2)
        return 16*pi/13125*k2*(k2 - h2)*res

    def G37(h2, k2):
        return 4*pi*h2**2*k2**2*(k2 - h2)**2/105

    assert_almost_equal(ellip_normal(5,8,0,1), G01(5, 8))
    assert_almost_equal(ellip_normal(5,8,1,1), G11(5, 8))
    assert_almost_equal(ellip_normal(5,8,1,2), G12(5, 8))
    assert_almost_equal(ellip_normal(5,8,1,3), G13(5, 8))
    assert_almost_equal(ellip_normal(5,8,2,1), G21(5, 8))
    assert_almost_equal(ellip_normal(5,8,2,2), G22(5, 8))
    assert_almost_equal(ellip_normal(5,8,2,3), G23(5, 8))
    assert_almost_equal(ellip_normal(5,8,2,4), G24(5, 8))
    assert_almost_equal(ellip_normal(5,8,2,5), G25(5, 8))
    assert_almost_equal(ellip_normal(5,8,3,1), G31(5, 8))
    assert_almost_equal(ellip_normal(5,8,3,2), G32(5, 8))
    assert_almost_equal(ellip_normal(5,8,3,3), G33(5, 8))
    assert_almost_equal(ellip_normal(5,8,3,4), G34(5, 8))
    assert_almost_equal(ellip_normal(5,8,3,5), G35(5, 8))
    assert_almost_equal(ellip_normal(5,8,3,6), G36(5, 8))
    assert_almost_equal(ellip_normal(5,8,3,7), G37(5, 8))
    assert_almost_equal(ellip_normal(10.14,18.432,0,1), G01(10.14, 18.432))
    assert_almost_equal(ellip_normal(10.14,18.432,1,1), G11(10.14, 18.432))
    assert_almost_equal(ellip_normal(10.14,18.432,1,2), G12(10.14, 18.432))
    assert_almost_equal(ellip_normal(10.14,18.432,1,3), G13(10.14, 18.432))
    assert_almost_equal(ellip_normal(10.14,18.432,2,1), G21(10.14, 18.432))
    assert_almost_equal(ellip_normal(10.14,18.432,2,2), G22(10.14, 18.432))
    assert_almost_equal(ellip_normal(10.14,18.432,2,3), G23(10.14, 18.432))
    assert_almost_equal(ellip_normal(10.14,18.432,2,4), G24(10.14, 18.432))
    assert_almost_equal(ellip_normal(10.14,18.432,2,5), G25(10.14, 18.432))
    assert_almost_equal(ellip_normal(10.14,18.432,3,1), G31(10.14, 18.432))
    assert_almost_equal(ellip_normal(10.14,18.432,3,2), G32(10.14, 18.432))
    assert_almost_equal(ellip_normal(10.14,18.432,3,3), G33(10.14, 18.432))
    assert_almost_equal(ellip_normal(10.14,18.432,3,4), G34(10.14, 18.432))
    assert_almost_equal(ellip_normal(10.14,18.432,3,5), G35(10.14, 18.432))
    assert_almost_equal(ellip_normal(10.14,18.432,3,6), G36(10.14, 18.432))
    assert_almost_equal(ellip_normal(10.14,18.432,3,7), G37(10.14, 18.432))

def test_ellip_harm_2():

    def I1(h2, k2, s):
        res = ellip_harm_2(h2, k2, 1, 1, s)/(3 * ellip_harm(h2, k2, 1, 1, s)) + ellip_harm_2(h2, k2, 1, 2, s)/(3 * ellip_harm(h2, k2, 1, 2, s)) + ellip_harm_2(h2, k2, 1, 3, s)/(3 * ellip_harm(h2, k2, 1, 3, s))
        return res

    assert_almost_equal(I1(5, 8, 10), 1/(10*sqrt((100-5)*(100-8))))
    assert_almost_equal(ellip_harm_2(5, 8, 2, 1, 10), 0.00108056853382)
    assert_almost_equal(ellip_harm_2(5, 8, 2, 2, 10), 0.00105820513809)
    assert_almost_equal(ellip_harm_2(5, 8, 2, 3, 10), 0.00106058384743)
    assert_almost_equal(ellip_harm_2(5, 8, 2, 4, 10), 0.00106774492306)
    assert_almost_equal(ellip_harm_2(5, 8, 2, 5, 10), 0.00107976356454)

def test_ellip_harm():

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

