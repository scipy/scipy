import pytest
import numpy as np
from numpy.testing import assert_equal, assert_array_equal, assert_allclose, suppress_warnings

from scipy.special._ufuncs import sinpi, cospi, tanpi, cotpi


def test_integer_real_part():
    x = np.arange(-100, 101)
    y = np.hstack((-np.linspace(310, -30, 10), np.linspace(-30, 310, 10)))
    x, y = np.meshgrid(x, y)
    z = x + 1j*y
    # In the following we should be *exactly* right
    res = sinpi(z)
    assert_equal(res.real, 0.0)
    res = cospi(z)
    assert_equal(res.imag, 0.0)


def test_half_integer_real_part():
    x = np.arange(-100, 101) + 0.5
    y = np.hstack((-np.linspace(310, -30, 10), np.linspace(-30, 310, 10)))
    x, y = np.meshgrid(x, y)
    z = x + 1j*y
    # In the following we should be *exactly* right
    res = sinpi(z)
    assert_equal(res.imag, 0.0)
    res = cospi(z)
    assert_equal(res.real, 0.0)


@pytest.mark.skip("Temporary skip while gh-19526 is being resolved")
def test_intermediate_overlow():
    # Make sure we avoid overflow in situations where cosh/sinh would
    # overflow but the product with sin/cos would not
    sinpi_pts = [complex(1 + 1e-14, 227),
                 complex(1e-35, 250),
                 complex(1e-301, 445)]
    # Data generated with mpmath
    sinpi_std = [complex(-8.113438309924894e+295, -np.inf),
                 complex(1.9507801934611995e+306, np.inf),
                 complex(2.205958493464539e+306, np.inf)]
    with suppress_warnings() as sup:
        sup.filter(RuntimeWarning, "invalid value encountered in multiply")
        for p, std in zip(sinpi_pts, sinpi_std):
            res = sinpi(p)
            assert_allclose(res.real, std.real)
            assert_allclose(res.imag, std.imag)

    # Test for cosine, less interesting because cos(0) = 1.
    p = complex(0.5 + 1e-14, 227)
    std = complex(-8.113438309924894e+295, -np.inf)
    with suppress_warnings() as sup:
        sup.filter(RuntimeWarning, "invalid value encountered in multiply")
        res = cospi(p)
        assert_allclose(res.real, std.real)
        assert_allclose(res.imag, std.imag)

def test_sinpi():
    x = np.arange(-16, 16, 1 / 16)

    expected = np.sin(np.pi * x)
    actual = sinpi(x)

    assert_allclose(expected, actual, atol=1e-10)

def test_cospi():
    x = np.arange(-16, 16, 1 / 16)

    expected = np.cos(np.pi * x)
    actual = cospi(x)

    assert_allclose(expected, actual, atol=1e-10)

def test_tanpi():
    x = np.arange(-16, 16, 1 / 16)
    x = x[x % 0.5 != 0]

    expected = np.tan(np.pi * x)
    actual = tanpi(x)

    assert_allclose(expected, actual, atol=1e-10)

def test_cotpi():
    x = np.arange(-16, 16, 1 / 16)
    x = x[x % 0.5 != 0]

    expected = 1 / np.tan(np.pi * x)
    actual = cotpi(x)

    assert_allclose(expected, actual, atol=1e-10)

def test_tanpi_complex():
    x = np.arange(-4, 4, 1 / 16)
    y = np.array([-32, -16, -2, -1, -0.5, 0.5, 1, 2, 16, 32]) 

    x, y = np.meshgrid(x, y)

    z = x + 1j*y

    expected = np.tan(np.pi * z)
    actual = tanpi(z)

    assert_allclose(expected.real, actual.real, rtol=1e-7, atol=1e-7)
    assert_allclose(expected.imag, actual.imag, rtol=1e-7, atol=1e-7)

def test_tanpi_complex_near_pole():
    eps_list = [
        1/4, 1/8, 1/16,
        1/32, 1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048, 
        1/4096, 1/8192, 1/16384, 1/32768, 1/65536
    ]

    for eps in eps_list:
        for x in [-4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5]:
            z = np.array([
                complex(eps, eps), complex(eps, 0), complex(eps, -eps),
                complex(0, eps), complex(0, -eps),
                complex(-eps, eps), complex(-eps, 0), complex(-eps, -eps)
            ])

            actual = tanpi(z + x)
            expected = -1.0 / np.tan(z * np.pi)

            assert_allclose(expected.real, actual.real, rtol=1e-13, atol=1e-13)
            assert_allclose(expected.imag, actual.imag, rtol=1e-13, atol=1e-13)

def test_cotpi_complex():
    x = np.arange(-4, 4, 1 / 16)
    y = np.array([-32, -16, -2, -1, -0.5, 0.5, 1, 2, 16, 32]) 

    x, y = np.meshgrid(x, y)

    z = x + 1j*y

    expected = 1 / np.tan(np.pi * z)
    actual = cotpi(z)

    assert_allclose(expected.real, actual.real, rtol=1e-7, atol=1e-7)
    assert_allclose(expected.imag, actual.imag, rtol=1e-7, atol=1e-7)

def test_cotpi_complex_near_pole():
    eps_list = [
        1/4, 1/8, 1/16,
        1/32, 1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048, 
        1/4096, 1/8192, 1/16384, 1/32768, 1/65536
    ]

    for eps in eps_list:
        for x in [-4, -3, -2, -1, 0, 1, 2, 3, 4]:
            z = np.array([
                complex(eps, eps), complex(eps, 0), complex(eps, -eps),
                complex(0, eps), complex(0, -eps),
                complex(-eps, eps), complex(-eps, 0), complex(-eps, -eps)
            ])

            actual = cotpi(z + x)
            expected = 1.0 / np.tan(z * np.pi)

            assert_allclose(expected.real, actual.real, rtol=1e-13, atol=1e-13)
            assert_allclose(expected.imag, actual.imag, rtol=1e-13, atol=1e-13)

def test_sinpi_zero_sign_ieee754_2008():
    y = sinpi([-3.0, -2.0, -1.0, -0.0, +0.0, 1.0, 2.0, 3.0])
    assert_array_equal(y, np.zeros(8))
    assert_array_equal(np.signbit(y), [True, True, True, True, False, False, False, False])

def test_cospi_zero_sign_ieee754_2008():
    y = cospi([-2.5, -1.5, -0.5, 0.5, 1.5, 2.5])
    assert_array_equal(y, np.zeros(6))
    assert_array_equal(np.signbit(y), [False, False, False, False, False, False])

def test_tanpi_zero_sign_ieee754_2019():
    y = tanpi([-3.0, -2.0, -1.0, -0.0, +0.0, 1.0, 2.0, 3.0])
    assert_array_equal(y, np.zeros(8))
    assert_array_equal(np.signbit(y), [False, True, False, True, False, True, False, True])

def test_cotpi_zero_sign():
    y = cotpi([-3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5])
    assert_array_equal(y, np.zeros(8))
    assert_array_equal(np.signbit(y), [False, True, False, True, False, True, False, True])

def test_sinpi_zero_sign_ieee754_2008_complex():
    y = sinpi(np.array([-3.0, -2.0, -1.0, -0.0, +0.0, 1.0, 2.0, 3.0], dtype=np.complex128))
    assert_array_equal(y.real, np.zeros(8))
    assert_array_equal(np.signbit(y.real), [True, True, True, True, False, False, False, False])

def test_cospi_zero_sign_ieee754_2008_complex():
    y = cospi(np.array([-2.5, -1.5, -0.5, 0.5, 1.5, 2.5], dtype=np.complex128))
    assert_array_equal(y.real, np.zeros(6))
    assert_array_equal(np.signbit(y.real), [False, False, False, False, False, False])

def test_tanpi_zero_sign_ieee754_2019_complex():
    y = tanpi(np.array([-3.0, -2.0, -1.0, -0.0, +0.0, 1.0, 2.0, 3.0], dtype=np.complex128))
    assert_array_equal(y.real, np.zeros(8))
    assert_array_equal(np.signbit(y.real), [False, True, False, True, False, True, False, True])

def test_cotpi_zero_sign_complex():
    y = cotpi(np.array([-3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5], dtype=np.complex128))
    assert_array_equal(y.real, np.zeros(8))
    assert_array_equal(np.signbit(y.real), [False, True, False, True, False, True, False, True])
