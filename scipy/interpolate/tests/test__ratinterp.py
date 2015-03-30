from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_allclose, assert_array_less, assert_equal

from scipy._lib.six import xrange
from scipy.interpolate._ratinterp import floater_hormann, lsq_rational, _sort_and_average_duplicates
from scipy.special import factorial


def test_floater_hormann_weights():
    # Check that the vectorized computation for the weights is done
    # correctly

    def weights(x, d):
        # Non-vectorized computation
        n = len(x)
        u = np.zeros((n,), dtype=float)
        for i in xrange(n):
            for k in xrange(max(i-d, 0), min(i+1, n-d)):
                c = 1.0
                for p in xrange(k, k+d+1):
                    if p != i:
                        c /= abs(x[i] - x[p])
                u[i] += c
        u *= (-1)**(np.arange(n) - d)
        return u

    np.random.seed(1234)
    x = np.random.rand(31)
    y = np.random.rand(31)
    x.sort()

    for d in xrange(0, x.size+1):
        # Check weights
        u0 = weights(x, min(d, x.size-1))
        ip = floater_hormann(x, y, d=d)
        assert_allclose(ip.u, u0, rtol=1e-13,
                        err_msg="d=%d" % (d,))

        # Check interpolation
        assert_allclose(ip(x), y)


def test_floater_hormann_weights_uniform():
    # Check against explicit results on an uniform grid
    x = np.arange(11)

    def scale(d):
        return (-1)**(np.arange(len(x)) + d) * factorial(d)

    ip = floater_hormann(x, 0*x, d=0)
    assert_allclose(ip.u*scale(0), [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

    ip = floater_hormann(x, 0*x, d=1)
    assert_allclose(ip.u*scale(1), [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1])

    ip = floater_hormann(x, 0*x, d=2)
    assert_allclose(ip.u*scale(2), [1, 3, 4, 4, 4, 4, 4, 4, 4, 3, 1])

    ip = floater_hormann(x, 0*x, d=3)
    assert_allclose(ip.u*scale(3), [1, 4, 7, 8, 8, 8, 8, 8, 7, 4, 1])

    ip = floater_hormann(x, 0*x, d=4)
    assert_allclose(ip.u*scale(4), [1, 5, 11, 15, 16, 16, 16, 15, 11, 5, 1])


def test_lsq_rational_order():
    # Check that the order conditions for the polynomials apply
    np.random.seed(1234)
    x = np.random.rand(11)
    y = np.random.rand(11)

    for m in range(5, 11):
        ip = lsq_rational(x, y, m=m)
        denom = ip.get_denominator()
        numer = ip.get_numerator()

        # Check that the numerator polynomial is (approximatively) of order m
        assert_array_less(abs(numer.coeffs[:-(m+1)]), 1e-14 + 1e-6*abs(numer.coeffs[-(m+1):]).max(), repr(m))

        # Check that the denominator polynomial is (approximatively) of order 10-m
        k = 10 - m
        assert_array_less(abs(denom.coeffs[:-(k+1)]), 1e-14 + 1e-6*abs(denom.coeffs[-(k+1):]).max(), repr(m))

        # Check interpolation
        assert_allclose(ip(x), y)


def test_lsq_interpolation():
    np.random.seed(1234)

    f = lambda x: np.sin(x) / (1 + x**2)

    x = np.random.rand(31) + 1j*np.random.rand(31)
    y = f(x)

    ip = lsq_rational(x, y)
    xx = np.random.rand(1001) + 1j*np.random.rand(1001)

    # Interpolation of analytic functions converges fast, even for a
    # fairly small number of points.
    assert_allclose(ip(xx), f(xx), rtol=1e-5)


def test_floater_hormann_interpolation():
    np.random.seed(1234)

    x = np.linspace(0, 1, 51)
    y = 1 / (1 + x**2)
    h = x[1] - x[0]

    xx = np.random.rand(1001)
    yy = 1 / (1 + xx**2)

    # Check convergence order
    for d in range(0, 10):
        ip = floater_hormann(x, y, d=d)
        tol = 10*h**(d+1)
        msg = repr((d, abs(ip(xx)/yy - 1).max(), tol))
        assert_allclose(ip(xx), yy, atol=1e-10, rtol=tol, err_msg=msg)


def test_sort_and_average_duplicates():
    x, y = _sort_and_average_duplicates([1, 0, 2, 3, 4], [2, 1, 3, 4, 5])
    assert_equal(x, [0, 1, 2, 3, 4])
    assert_equal(y, [1, 2, 3, 4, 5])

    x, y = _sort_and_average_duplicates([2, 1, 1, 3, 4], [3, 1, 2, 4, 5])
    assert_equal(x, [1, 2, 3, 4])
    assert_allclose(y, [1.5, 3, 4, 5])

    x, y = _sort_and_average_duplicates([2, 1, 1, 1, 4], [3, 1, 2, 4, 5])
    assert_equal(x, [1, 2, 4])
    assert_allclose(y, [(1+2+4)/3, 3, 5])

    x, y = _sort_and_average_duplicates([4, 1, 1, 1, 4], [3, 1, 2, 4, 5])
    assert_equal(x, [1, 4])
    assert_allclose(y, [(1+2+4)/3, (3+5)/2])

    x, y = _sort_and_average_duplicates([1, 1, 1, 1, 1], [3, 1, 2, 4, 5])
    assert_equal(x, [1])
    assert_allclose(y, [(1+2+3+4+5)/5])

    x, y = _sort_and_average_duplicates([2, 1, 1, 1, 4], [[3, 1, 2, 4, 5]]*4)
    assert_equal(x, [1, 2, 4])
    assert_allclose(y, [[(1+2+4)/3, 3, 5]]*4)


def test_polynomials():
    def func(x):
        return (1 + 5*x + 2*x**3) / (3j + x + x**2)

    xx = np.linspace(0, 1, 7)
    ip = lsq_rational(xx, func(xx))

    # Check power basis form (low order so it still works)
    numer = ip.get_numerator()
    denom = ip.get_denominator()

    numer_coef = numer.coef / numer.coef[-1]
    denom_coef = denom.coef / numer.coef[-1]
    assert_allclose(numer_coef, [0, 0, 0, 2, 0, 5, 1], atol=1e-6)
    assert_allclose(denom_coef, [0, 0, 0, 0, 1, 1, 3j], atol=1e-6)

    xx = np.linspace(0, 1, 200)
    assert_allclose(numer(xx)/denom(xx), ip(xx), rtol=1e-8)

    def root_cmp(a, b, radius):
        a = a[abs(a) < radius]
        a.sort()
        b = b[abs(b) < radius]
        b.sort()
        assert_allclose(a, b)

    root_cmp(ip.numerator_roots, numer.roots, 1e3)
    root_cmp(ip.denominator_roots, denom.roots, 1e3)

    # Check barycentric form
    numer_bary = ip.get_numerator('barycentric')
    denom_bary = ip.get_denominator('barycentric')
    assert_allclose(numer_bary(xx), numer(xx))
    assert_allclose(denom_bary(xx), denom(xx))

    # Check barycentric form (high order should still work)
    xx = np.linspace(0, 1, 100)
    yy = np.cos(2*np.pi*xx)
    ip = floater_hormann(xx, yy)
    xx = np.linspace(0, 1, 5000)
    numer = ip.get_numerator('barycentric')
    denom = ip.get_denominator('barycentric')
    assert_allclose(numer(xx) / denom(xx), ip(xx))
