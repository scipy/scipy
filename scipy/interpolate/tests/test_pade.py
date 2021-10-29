from numpy.testing import (assert_array_equal, assert_array_almost_equal)
from scipy.interpolate import pade


def test_pade_trivial():
    nump, denomp = pade([1.0], 0)
    assert_array_equal(nump.coef, [1.0])
    assert_array_equal(denomp.coef, [1.0])

    nump, denomp = pade([1.0], 0, 0)
    assert_array_equal(nump.coef, [1.0])
    assert_array_equal(denomp.coef, [1.0])


def test_pade_4term_exp():
    # First four Taylor coefficients of exp(x).
    an = [1.0, 1.0, 0.5, 1.0/6]

    nump, denomp = pade(an, 0)
    assert_array_almost_equal(nump.coef, an)
    assert_array_almost_equal(denomp.coef, [1.0])

    nump, denomp = pade(an, 1)
    assert_array_almost_equal(nump.coef, [1.0, 2.0/3, 1.0/6])
    assert_array_almost_equal(denomp.coef, [1.0, -1.0/3])

    nump, denomp = pade(an, 2)
    assert_array_almost_equal(nump.coef, [1.0, 1.0/3])
    assert_array_almost_equal(denomp.coef, [1.0, -2.0/3, 1.0/6])

    nump, denomp = pade(an, 3)
    assert_array_almost_equal(nump.coef, [1.0])
    assert_array_almost_equal(denomp.coef, [1.0, -1.0, 0.5, -1.0/6])

    # Testing inclusion of optional parameter
    nump, denomp = pade(an, 0, 3)
    assert_array_almost_equal(nump.coef, [1.0, 1.0, 0.5, 1.0/6])
    assert_array_almost_equal(denomp.coef, [1.0])

    nump, denomp = pade(an, 1, 2)
    assert_array_almost_equal(nump.coef, [1.0, 2.0/3, 1.0/6])
    assert_array_almost_equal(denomp.coef, [1.0, -1.0/3])

    nump, denomp = pade(an, 2, 1)
    assert_array_almost_equal(nump.coef, [1.0, 1.0/3])
    assert_array_almost_equal(denomp.coef, [1.0, -2.0/3, 1.0/6])

    nump, denomp = pade(an, 3, 0)
    assert_array_almost_equal(nump.coef, [1.0])
    assert_array_almost_equal(denomp.coef, [1.0, -1.0, 0.5, -1.0/6])

    # Testing reducing array.
    nump, denomp = pade(an, 0, 2)
    assert_array_almost_equal(nump.coef, [1.0, 1.0, 0.5])
    assert_array_almost_equal(denomp.coef, [1.0])

    nump, denomp = pade(an, 1, 1)
    assert_array_almost_equal(nump.coef, [1.0, 1.0/2])
    assert_array_almost_equal(denomp.coef, [1.0, -1.0/2])

    nump, denomp = pade(an, 2, 0)
    assert_array_almost_equal(nump.coef, [1.0])
    assert_array_almost_equal(denomp.coef, [1.0, -1.0, 1.0/2])


def test_pade_ints():
    # Simple test sequences (one of ints, one of floats).
    an_int = [1, 2, 3, 4]
    an_flt = [1.0, 2.0, 3.0, 4.0]

    # Make sure integer arrays give the same result as float arrays with same values.
    for i in range(0, len(an_int)):
        for j in range(0, len(an_int) - i):

            # Create float and int pade approximation for given order.
            nump_int, denomp_int = pade(an_int, i, j)
            nump_flt, denomp_flt = pade(an_flt, i, j)

            # Check that they are the same.
            assert_array_equal(nump_int.coef, nump_flt.coef)
            assert_array_equal(denomp_int.coef, denomp_flt.coef)


def test_pade_complex():
    # Test sequence with known solutions - see page 6 of 10.1109/PESGM.2012.6344759.
    # Variable x is parameter - these tests will work with any complex number.
    x = 0.2 + 0.6j
    xc = x.conjugate()
    an = [1.0, x, -x*xc, xc*(x**2) + x*(xc**2),
          -(x**3)*xc - 3*(x*xc)**2 - x*(xc**3)]

    nump, denomp = pade(an, 1, 1)
    assert_array_almost_equal(nump.coef, [1.0, x + xc])
    assert_array_almost_equal(denomp.coef, [1.0, xc])

    nump, denomp = pade(an, 1, 2)
    assert_array_almost_equal(nump.coef, [1.0, 2*x + xc, x**2])
    assert_array_almost_equal(denomp.coef, [1.0, x + xc])

    nump, denomp = pade(an, 2, 2)
    assert_array_almost_equal(nump.coef,
                              [1.0, 2*(x + xc), x**2 + x*xc + xc**2])
    assert_array_almost_equal(denomp.coef, [1.0, x + 2*xc, xc**2])
