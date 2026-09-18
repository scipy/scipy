import numpy as np
from numpy import pi, log, sqrt
from numpy.testing import assert_, assert_allclose, assert_equal

from scipy.special._testutils import FuncData
import scipy.special as sc

# Euler-Mascheroni constant
euler = 0.57721566490153286


def test_consistency():
    # Make sure the implementation of digamma for real arguments
    # agrees with the implementation of digamma for complex arguments.

    # It's all poles after -1e16
    x = np.r_[-np.logspace(15, -30, 200), np.logspace(-30, 300, 200)]
    dataset = np.vstack((x + 0j, sc.digamma(x))).T
    FuncData(sc.digamma, dataset, 0, 1, rtol=5e-14, nan_ok=True).check()


def test_special_values():
    # Test special values from Gauss's digamma theorem. See
    #
    # https://en.wikipedia.org/wiki/Digamma_function

    dataset = [
        (1, -euler),
        (0.5, -2*log(2) - euler),
        (1/3, -pi/(2*sqrt(3)) - 3*log(3)/2 - euler),
        (1/4, -pi/2 - 3*log(2) - euler),
        (1/6, -pi*sqrt(3)/2 - 2*log(2) - 3*log(3)/2 - euler),
        (1/8,
         -pi/2 - 4*log(2) - (pi + log(2 + sqrt(2)) - log(2 - sqrt(2)))/sqrt(2) - euler)
    ]

    dataset = np.asarray(dataset)
    FuncData(sc.digamma, dataset, 0, 1, rtol=1e-14).check()


def test_nonfinite():
    pts = [0.0, -0.0, np.inf]
    std = [-np.inf, np.inf, np.inf]
    assert_equal(sc.digamma(pts), std)
    assert_(all(np.isnan(sc.digamma([-np.inf, -1]))))


def test_digammainv_roundtrip():
    x = np.logspace(-20, 20, 500)
    assert_allclose(sc.digammainv(sc.digamma(x)), x, rtol=2e-14)


def test_digammainv_special_values():
    # Reference values for digamma were computed with mpmath.
    dataset = np.array([
        (1e-100, -1e100),
        (1e-50, -1e50),
        (1e-20, -1e20),
        (1e-10, -10000000000.577215),
        (1e-3, -1000.5755719318103),
        (1e-1, -10.423754940411076),
        (1/4, -4.2274535333762655),
        (1/3, -3.1320337800208065),
        (1/2, -1.9635100260214235),
        (1.0, -0.5772156649015329),
        (1.4616321449683623, 0.0),
        (10, 2.251752589066721),
        (1e5, 11.512920464961896),
        (1e50, 115.12925464970229),
        (1e200, 460.51701859880916),
    ])
    assert_allclose(sc.digammainv(dataset[:, 1]), dataset[:, 0], rtol=1e-13)


def test_digammainv_nonfinite():
    res = sc.digammainv([np.nan, np.inf, -np.inf])
    assert np.isnan(res[0])
    assert_equal(res[1:], [np.inf, 0.0])


def test_digammainv_ufunc():
    y = np.array([sc.digamma(0.5), 0.0, sc.digamma(10.0)])
    out = np.empty_like(y)
    res = sc.digammainv(y, out=out)
    assert res is out
    assert_allclose(out, [0.5, 1.4616321449683623, 10.0], rtol=1e-14)

    y32 = y.astype(np.float32)
    res32 = sc.digammainv(y32)
    assert_equal(res32.dtype, np.float32)
    assert_allclose(res32, [0.5, 1.4616321, 10.0], rtol=2e-6)
