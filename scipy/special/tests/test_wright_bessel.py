# Reference MPMATH implementation:
#
# import mpmath
# from mpmath import nsum
#
# def Wright_Series_MPMATH(a, b, z, dps=50, method='r+s+e', steps=[1000]):
#    """Compute Wright' generalized Bessel function as Series.
#
#    This uses mpmath for arbitrary precision.
#    """
#    with mpmath.workdps(dps):
#        res = nsum(lambda k: z**k/mpmath.fac(k) * mpmath.rgamma(a*k+b),
#                          [0, mpmath.inf],
#                          tol=dps, method=method, steps=steps
#                          )
#
#    return res

import pytest
import numpy as np
from numpy.testing import assert_equal, assert_allclose

import scipy.special as sc
from scipy.special import rgamma
from scipy.special.wright_bessel import wright_bessel


@pytest.mark.parametrize('a', [0, 1e-6, 0.1, 0.5, 1, 10])
@pytest.mark.parametrize('b', [0, 1e-6, 0.1, 0.5, 1, 10])
def test_wright_bessel_zero(a, b):
    """Test at x = 0."""
    assert_equal(wright_bessel(a, b, 0.), rgamma(b))


@pytest.mark.parametrize('b', [0, 1e-6, 0.1, 0.5, 1, 10])
@pytest.mark.parametrize('x', [0, 1e-6, 0.1, 0.5, 1])
def test_wright_bessel_iv(b, x):
    """Test relation of wright_bessel and modified bessel function iv.

    iv(z) = (1/2*z)**v * Phi(1, v+1; 1/4*z**2).
    See https://dlmf.nist.gov/10.46.E2
    """
    if x != 0:
        v = b - 1
        wb = wright_bessel(1, v + 1, x**2 / 4.)
        # Note: iv(v, x) has precision of less than 1e-12 for some cases
        # e.g v=1-1e-6 and x=1e-06)
        assert_allclose(np.power(x / 2., v) * wb,
                        sc.iv(v, x),
                        rtol=1e-11, atol=1e-11)


@pytest.mark.parametrize('a', [0, 1e-6, 0.1, 0.5, 1, 10])
@pytest.mark.parametrize('b', [1, 1 + 1e-3, 2, 5, 10])
@pytest.mark.parametrize('x', [0, 1e-6, 0.1, 0.5, 1, 5, 10, 100])
def test_wright_functional(a, b, x):
    """Test functional relation of wright_bessel.

    Phi(a, b-1, z) = a*z*Phi(a, b+a, z) + (b-1)*Phi(a, b, z)

    Note that d/dx Phi(a, b, x) = Phi(a, b-1, x)
    See Eq. (22) of
    B. Stankovic, On the Function of E. M. Wright,
    Publ. de l' Institut Mathematique, Beograd,
    Nouvelle S`er. 10 (1970), 113-124.
    """
    assert_allclose(wright_bessel(a, b - 1, x),
                    a*x*wright_bessel(a, b + a, x)
                    + (b - 1) * wright_bessel(a, b, x),
                    rtol=1e-8, atol=1e-8)


@pytest.mark.xfail
@pytest.mark.parametrize(
    'a, b, x, phi, accuracy',
    [[0.9999999999999778, 10, 500, 454308.62922200584, 1e-11],
     [0.9999999999999778, 10, 709.7827128933841, 521034418.0650804, 1e-11],
     [0.9999999999999778, 10, 1000, 2446819552259.9204, 1e-11],
     [1, 10, 500, 454308.6292213986, 1e-11],
     [1, 10, 709.7827128933841, 521034418.064186, 1e-11],
     [1, 10, 1000, 2446819552254.583, 1e-11],
     [1.0000000000000222, 10, 500, 454308.62922079145, 1e-11],
     [1.0000000000000222, 10, 709.7827128933841, 521034418.0632915, 1e-11],
     [1.0000000000000222, 10, 1000, 2446819552249.245, 1e-11],
     [1, 100, 100000, 1.0269334596230763e+22, 1e-11],
     [1.0000000000000222, 100, 100000, 1.0269334595866202e+22, 1e-11],
     [2, 100, 100000, 5.475089686699177e-153, 1e-11]])
def test_wright_data_grid_failures(a, b, x, phi, accuracy):
    """Test cases of test_data that do not reach relative accuracy of 1e-11"""
    assert wright_bessel(a, b, x) == pytest.approx(phi, rel=accuracy)
