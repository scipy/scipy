#
# Tests of spherical Bessel functions.
#
import warnings

import numpy as np
from numpy.testing import assert_allclose
import pytest
from numpy import sin, cos, sinh, cosh, exp, inf, nan, r_, pi

from scipy.special import spherical_jn, spherical_yn, spherical_in, spherical_kn
from scipy.integrate import quad

from scipy._lib._array_api import make_xp_test_case, xp_assert_close


@make_xp_test_case(spherical_jn)
class TestSphericalJn:
    def test_spherical_jn_exact(self, xp):
        # https://dlmf.nist.gov/10.49.E3
        # Note: exact expression is numerically stable only for small
        # n or z >> n.
        x_np = np.asarray([0.12, 1.23, 12.34, 123.45, 1234.5])
        x = xp.asarray(x_np)
        desired = xp.asarray((-1/x_np + 3/x_np**3)*sin(x_np) - 3/x_np**2*cos(x_np))
        xp_assert_close(spherical_jn(2, x), desired)

    def test_spherical_jn_recurrence_complex(self, xp):
        # https://dlmf.nist.gov/10.51.E1
        n = xp.asarray([1, 2, 3, 7, 12])
        x = xp.asarray(1.1 + 1.5j, dtype=xp.complex128)
        xp_assert_close(spherical_jn(n - 1, x) + spherical_jn(n + 1, x),
                        xp.astype(2*n + 1, xp.float64)/x*spherical_jn(n, x))

    def test_spherical_jn_recurrence_real(self, xp):
        # https://dlmf.nist.gov/10.51.E1
        n = xp.asarray([1, 2, 3, 7, 12])
        x = xp.asarray(0.12, dtype=xp.float64)
        xp_assert_close(spherical_jn(n - 1, x) + spherical_jn(n + 1,x),
                        xp.astype(2*n + 1, xp.float64)/x*spherical_jn(n, x))

    def test_spherical_jn_inf_real(self, xp):
        # https://dlmf.nist.gov/10.52.E3
        n = xp.asarray(6)
        x = xp.asarray([-inf, inf], dtype=xp.float64)
        xp_assert_close(spherical_jn(n, x), xp.asarray([0, 0], dtype=xp.float64))

    def test_spherical_jn_inf_complex(self):
        # https://dlmf.nist.gov/10.52.E3
        n = 7
        x = np.array([-inf + 0j, inf + 0j, inf*(1+1j)])
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", "invalid value encountered in multiply", RuntimeWarning)
            assert_allclose(spherical_jn(n, x), np.array([0, 0, inf*(1+1j)]))

    def test_spherical_jn_large_arg_1(self, xp):
        # https://github.com/scipy/scipy/issues/2165
        # Reference value computed using mpmath, via
        # besselj(n + mpf(1)/2, z)*sqrt(pi/(2*z))
        xp_assert_close(
            spherical_jn(2, xp.asarray(3350.507, dtype=xp.float64)),
            xp.asarray(-0.00029846226538040747, dtype=xp.float64),
            check_0d=False,
        )

    def test_spherical_jn_large_arg_2(self, xp):
        # https://github.com/scipy/scipy/issues/1641
        # Reference value computed using mpmath, via
        # besselj(n + mpf(1)/2, z)*sqrt(pi/(2*z))
        xp_assert_close(
            spherical_jn(2, xp.asarray(10000, dtype=xp.float64)),
            xp.asarray(3.0590002633029811e-05, dtype=xp.float64),
            check_0d=False,
        )

    def test_spherical_jn_at_zero(self, xp):
        # https://dlmf.nist.gov/10.52.E1
        # But note that n = 0 is a special case: j0 = sin(x)/x -> 1
        n = xp.asarray([0, 1, 2, 5, 10, 100])
        x = xp.asarray(0, dtype=xp.float64)
        xp_assert_close(
            spherical_jn(n, x), xp.asarray([1, 0, 0, 0, 0, 0], dtype=xp.float64)
        )


@make_xp_test_case(spherical_yn)
class TestSphericalYn:
    def test_spherical_yn_exact(self, xp):
        # https://dlmf.nist.gov/10.49.E5
        # Note: exact expression is numerically stable only for small
        # n or z >> n.
        x_np = np.asarray([0.12, 1.23, 12.34, 123.45, 1234.5])
        x = xp.asarray(x_np)
        desired = xp.asarray((1/x_np - 3/x_np**3)*cos(x_np) - 3/x_np**2*sin(x_np))
        xp_assert_close(spherical_yn(2, x), desired)

    def test_spherical_yn_recurrence_real(self, xp):
        # https://dlmf.nist.gov/10.51.E1
        n = xp.asarray([1, 2, 3, 7, 12])
        x = xp.asarray(0.12, dtype=xp.float64)
        xp_assert_close(spherical_yn(n - 1, x) + spherical_yn(n + 1, x),
                        xp.astype(2*n + 1, xp.float64)/x*spherical_yn(n, x))

    def test_spherical_yn_recurrence_complex(self, xp):
        # https://dlmf.nist.gov/10.51.E1
        n = xp.asarray([1, 2, 3, 7, 12])
        x = xp.asarray(1.1 + 1.5j, dtype=xp.complex128)
        xp_assert_close(spherical_yn(n - 1, x) + spherical_yn(n + 1, x),
                        xp.astype(2*n + 1, xp.float64)/x*spherical_yn(n, x))

    def test_spherical_yn_inf_real(self, xp):
        # https://dlmf.nist.gov/10.52.E3
        n = xp.asarray(6)
        x = xp.asarray([-inf, inf], dtype=xp.float64)
        xp_assert_close(spherical_yn(n, x),
                        xp.asarray([0, 0], dtype=xp.float64))

    def test_spherical_yn_inf_complex(self):
        # https://dlmf.nist.gov/10.52.E3
        n = 7
        x = np.array([-inf + 0j, inf + 0j, inf*(1+1j)])
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", "invalid value encountered in multiply", RuntimeWarning)
            assert_allclose(spherical_yn(n, x), np.array([0, 0, inf*(1+1j)]))

    def test_spherical_yn_at_zero(self, xp):
        # https://dlmf.nist.gov/10.52.E2
        n = xp.asarray([0, 1, 2, 5, 10, 100])
        x = xp.asarray(0, dtype=xp.float64)
        xp_assert_close(spherical_yn(n, x),
                        xp.full(n.shape, -inf, dtype=xp.float64))

    def test_spherical_yn_at_zero_complex(self, xp):
        # Consistently with numpy:
        # >>> -np.cos(0)/0
        # -inf
        # >>> -np.cos(0+0j)/(0+0j)
        # (-inf + nan*j)
        n = xp.asarray([0, 1, 2, 5, 10, 100])
        x = xp.asarray(0 + 0j, dtype=xp.complex128)
        xp_assert_close(
            spherical_yn(n, x), xp.full(n.shape, xp.nan, dtype=xp.complex128)
        )


@make_xp_test_case(spherical_jn, spherical_yn)
class TestSphericalJnYnCrossProduct:
    def test_spherical_jn_yn_cross_product_1(self, xp):
        # https://dlmf.nist.gov/10.50.E3
        n = xp.asarray([1, 5, 8])
        x = xp.asarray([0.1, 1, 10], dtype=xp.float64)
        left = (spherical_jn(n + 1, x) * spherical_yn(n, x) -
                spherical_jn(n, x) * spherical_yn(n + 1, x))
        right = 1/x**2
        xp_assert_close(left, right)

    def test_spherical_jn_yn_cross_product_2(self, xp):
        # https://dlmf.nist.gov/10.50.E3
        n = xp.asarray([1, 5, 8])
        x = xp.asarray([0.1, 1, 10], dtype=xp.float64)
        left = (spherical_jn(n + 2, x) * spherical_yn(n, x) -
                spherical_jn(n, x) * spherical_yn(n + 2, x))
        right = xp.astype(2*n + 3, xp.float64)/x**3
        xp_assert_close(left, right)


import numpy as np
import pytest
from scipy._lib._array_api import make_xp_test_case, xp_assert_close
from scipy.special import spherical_kn


@make_xp_test_case(spherical_kn)
class TestSphericalKn:
    def test_spherical_kn_exact(self, xp):
        # https://dlmf.nist.gov/10.49.E13
        x_np = np.array([0.12, 1.23, 12.34, 123.45])
        x = xp.asarray(x_np)
        expected = np.pi / 2 * np.exp(-x_np) * (
            1 / x_np + 3 / x_np**2 + 3 / x_np**3
        )
        xp_assert_close(spherical_kn(2, x), xp.asarray(expected))

    def test_spherical_kn_recurrence_real(self, xp):
        # https://dlmf.nist.gov/10.51.E4
        n_np = np.array([1, 2, 3, 7, 12])
        n = xp.asarray(n_np)
        x_np = 0.12
        x = xp.asarray(x_np)
        sign_lower = xp.asarray((-1.0) ** (n_np - 1))
        sign_upper = xp.asarray((-1.0) ** (n_np + 1))
        factor = xp.asarray((-1.0) ** n_np * (2 * n_np + 1) / x_np)
        xp_assert_close(
            sign_lower * spherical_kn(n - 1, x)
            - sign_upper * spherical_kn(n + 1, x),
            factor * spherical_kn(n, x),
        )

    def test_spherical_kn_recurrence_complex(self, xp):
        # https://dlmf.nist.gov/10.51.E4
        n_np = np.array([1, 2, 3, 7, 12])
        n = xp.asarray(n_np)
        x_np = 1.1 + 1.5j
        x = xp.asarray(x_np)
        sign_lower = xp.asarray((-1.0) ** (n_np - 1))
        sign_upper = xp.asarray((-1.0) ** (n_np + 1))
        factor = xp.asarray((-1.0) ** n_np * (2 * n_np + 1) / x_np)
        xp_assert_close(
            sign_lower * spherical_kn(n - 1, x)
            - sign_upper * spherical_kn(n + 1, x),
            factor * spherical_kn(n, x),
        )

    def test_spherical_kn_inf_real(self, xp):
        # https://dlmf.nist.gov/10.52.E6
        n = 5
        x = xp.asarray(np.array([-np.inf, np.inf]))
        expected = xp.asarray(np.array([-np.inf, 0]))
        xp_assert_close(spherical_kn(n, x), expected)

    def test_spherical_kn_inf_complex(self, xp):
        # https://dlmf.nist.gov/10.52.E6
        # The behavior at complex infinity depends on the sign of the real
        # part: if Re(z) >= 0, then the limit is 0; if Re(z) < 0, then it's
        # z*inf. This distinction cannot be captured, so we return nan.
        n = 7
        x = xp.asarray(np.array([-np.inf + 0j, np.inf + 0j, np.inf * (1 + 1j)]))
        expected = xp.asarray(np.array([-np.inf, 0, np.nan], dtype=np.complex128))
        xp_assert_close(spherical_kn(n, x), expected)

    def test_spherical_kn_at_zero(self, xp):
        # https://dlmf.nist.gov/10.52.E2
        n_np = np.array([0, 1, 2, 5, 10, 100])
        n = xp.asarray(n_np)
        x = xp.asarray(0)
        expected = xp.asarray(np.full(n_np.shape, np.inf))
        xp_assert_close(spherical_kn(n, x), expected)

    def test_spherical_kn_at_zero_complex(self, xp):
        # https://dlmf.nist.gov/10.52.E2
        n_np = np.array([0, 1, 2, 5, 10, 100])
        n = xp.asarray(n_np)
        x = xp.asarray(0 + 0j)
        expected = xp.asarray(np.full(n_np.shape, np.nan, dtype=np.complex128))
        xp_assert_close(spherical_kn(n, x), expected)


@make_xp_test_case(spherical_in)
class TestSphericalIn:
    def test_spherical_in_exact(self, xp):
        # https://dlmf.nist.gov/10.49.E9
        x_np = np.asarray([0.12, 1.23, 12.34, 123.45])
        x = xp.asarray(x_np)
        desired = xp.asarray(
            (1/x_np + 3/x_np**3)*sinh(x_np) - 3/x_np**2*cosh(x_np)
        )
        xp_assert_close(spherical_in(2, x), desired)

    def test_spherical_in_recurrence_real(self, xp):
        # https://dlmf.nist.gov/10.51.E4
        n = xp.asarray([1, 2, 3, 7, 12])
        x = xp.asarray(0.12, dtype=xp.float64)
        xp_assert_close(spherical_in(n - 1, x) - spherical_in(n + 1, x),
                        xp.astype(2*n + 1, xp.float64)/x*spherical_in(n, x))

    def test_spherical_in_recurrence_complex(self, xp):
        # https://dlmf.nist.gov/10.51.E1
        n = xp.asarray([1, 2, 3, 7, 12])
        x = xp.asarray(1.1 + 1.5j, dtype=xp.complex128)
        xp_assert_close(spherical_in(n - 1, x) - spherical_in(n + 1, x),
                        xp.astype(2*n + 1, xp.float64)/x*spherical_in(n, x))

    def test_spherical_in_inf_real(self, xp):
        # https://dlmf.nist.gov/10.52.E3
        n = xp.asarray(5)
        x = xp.asarray([-inf, inf], dtype=xp.float64)
        xp_assert_close(spherical_in(n, x),
                        xp.asarray([-inf, inf], dtype=xp.float64))

    def test_spherical_in_inf_complex(self, xp):
        # https://dlmf.nist.gov/10.52.E5
        # Ideally, i1n(n, 1j*inf) = 0 and i1n(n, (1+1j)*inf) = (1+1j)*inf, but
        # this appears impossible to achieve because C99 regards any complex
        # value with at least one infinite part as a complex infinity, so
        # 1j*inf cannot be distinguished from (1+1j)*inf. Therefore, nan is
        # the correct return value.
        n = xp.asarray(7)
        x = xp.asarray([-inf + 0j, inf + 0j, inf*(1+1j)],
                       dtype=xp.complex128)
        xp_assert_close(spherical_in(n, x),
                        xp.asarray([-inf, inf, nan], dtype=xp.complex128))

    def test_spherical_in_at_zero(self, xp):
        # https://dlmf.nist.gov/10.52.E1
        # But note that n = 0 is a special case: i0 = sinh(x)/x -> 1
        n = xp.asarray([0, 1, 2, 5, 10, 100])
        x = xp.asarray(0, dtype=xp.float64)
        xp_assert_close(
            spherical_in(n, x),
            xp.asarray([1, 0, 0, 0, 0, 0], dtype=xp.float64),
        )


class SphericalDerivativesTestCase:
    def fundamental_theorem(self, n, a, b):
        integral, tolerance = quad(lambda z: self.df(n, z), a, b)
        assert_allclose(integral,
                        self.f(n, b) - self.f(n, a),
                        atol=tolerance)

    @pytest.mark.slow
    def test_fundamental_theorem_0(self):
        self.fundamental_theorem(0, 3.0, 15.0)

    @pytest.mark.slow
    def test_fundamental_theorem_7(self):
        self.fundamental_theorem(7, 0.5, 1.2)


class TestSphericalJnDerivatives(SphericalDerivativesTestCase):
    def f(self, n, z):
        return spherical_jn(n, z)

    def df(self, n, z):
        return spherical_jn(n, z, derivative=True)

    def test_spherical_jn_d_zero(self):
        n = np.array([0, 1, 2, 3, 7, 15])
        assert_allclose(spherical_jn(n, 0, derivative=True),
                        np.array([0, 1/3, 0, 0, 0, 0]))


class TestSphericalYnDerivatives(SphericalDerivativesTestCase):
    def f(self, n, z):
        return spherical_yn(n, z)

    def df(self, n, z):
        return spherical_yn(n, z, derivative=True)


class TestSphericalInDerivatives(SphericalDerivativesTestCase):
    def f(self, n, z):
        return spherical_in(n, z)

    def df(self, n, z):
        return spherical_in(n, z, derivative=True)

    def test_spherical_in_d_zero(self):
        n = np.array([0, 1, 2, 3, 7, 15])
        spherical_in(n, 0, derivative=False)
        assert_allclose(spherical_in(n, 0, derivative=True),
                        np.array([0, 1/3, 0, 0, 0, 0]))


class TestSphericalKnDerivatives(SphericalDerivativesTestCase):
    def f(self, n, z):
        return spherical_kn(n, z)

    def df(self, n, z):
        return spherical_kn(n, z, derivative=True)


class TestSphericalOld:
    # These are tests from the TestSpherical class of test_basic.py,
    # rewritten to use spherical_* instead of sph_* but otherwise unchanged.

    def test_sph_in(self):
        # This test reproduces test_basic.TestSpherical.test_sph_in.
        i1n = np.empty((2,2))
        x = 0.2

        i1n[0][0] = spherical_in(0, x)
        i1n[0][1] = spherical_in(1, x)
        i1n[1][0] = spherical_in(0, x, derivative=True)
        i1n[1][1] = spherical_in(1, x, derivative=True)

        inp0 = (i1n[0][1])
        inp1 = (i1n[0][0] - 2.0/0.2 * i1n[0][1])
        assert_allclose(i1n[0], np.array([1.0066800127054699381,
                                          0.066933714568029540839]),
                        atol=1.5e-12, rtol=0.0)
        assert_allclose(i1n[1], [inp0, inp1], atol=1.5e-12, rtol=0)

    def test_sph_in_kn_order0(self):
        x = 1.
        sph_i0 = np.empty((2,))
        sph_i0[0] = spherical_in(0, x)
        sph_i0[1] = spherical_in(0, x, derivative=True)
        sph_i0_expected = np.array([np.sinh(x)/x,
                                    np.cosh(x)/x-np.sinh(x)/x**2])
        assert_allclose(r_[sph_i0], sph_i0_expected, atol=1.5e-7, rtol=0)

        sph_k0 = np.empty((2,))
        sph_k0[0] = spherical_kn(0, x)
        sph_k0[1] = spherical_kn(0, x, derivative=True)
        sph_k0_expected = np.array([0.5*pi*exp(-x)/x,
                                    -0.5*pi*exp(-x)*(1/x+1/x**2)])
        assert_allclose(r_[sph_k0], sph_k0_expected, atol=1.5e-7, rtol=0)

    def test_sph_jn(self):
        s1 = np.empty((2,3))
        x = 0.2

        s1[0][0] = spherical_jn(0, x)
        s1[0][1] = spherical_jn(1, x)
        s1[0][2] = spherical_jn(2, x)
        s1[1][0] = spherical_jn(0, x, derivative=True)
        s1[1][1] = spherical_jn(1, x, derivative=True)
        s1[1][2] = spherical_jn(2, x, derivative=True)

        s10 = -s1[0][1]
        s11 = s1[0][0]-2.0/0.2*s1[0][1]
        s12 = s1[0][1]-3.0/0.2*s1[0][2]
        assert_allclose(s1[0], [0.99334665397530607731,
                                0.066400380670322230863,
                                0.0026590560795273856680],
                        atol=1.5e-12, rtol=0)
        assert_allclose(s1[1], [s10, s11, s12], atol=1.5e-12, rtol=0)

    def test_sph_kn(self):
        kn = np.empty((2,3))
        x = 0.2

        kn[0][0] = spherical_kn(0, x)
        kn[0][1] = spherical_kn(1, x)
        kn[0][2] = spherical_kn(2, x)
        kn[1][0] = spherical_kn(0, x, derivative=True)
        kn[1][1] = spherical_kn(1, x, derivative=True)
        kn[1][2] = spherical_kn(2, x, derivative=True)

        kn0 = -kn[0][1]
        kn1 = -kn[0][0]-2.0/0.2*kn[0][1]
        kn2 = -kn[0][1]-3.0/0.2*kn[0][2]
        assert_allclose(kn[0], [6.4302962978445670140,
                                38.581777787067402086,
                                585.15696310385559829],
                        atol=1.5e-12, rtol=0)
        assert_allclose(kn[1], [kn0, kn1, kn2], atol=1.5e-9, rtol=0)

    def test_sph_yn(self):
        sy1 = spherical_yn(2, 0.2)
        sy2 = spherical_yn(0, 0.2)
        # previous values in the system
        assert_allclose(sy1, -377.52483, atol=1.5e-5, rtol=0)
        assert_allclose(sy2, -4.9003329, atol=1.5e-5, rtol=0)
        sphpy = (spherical_yn(0, 0.2) - 2*spherical_yn(2, 0.2))/3
        sy3 = spherical_yn(1, 0.2, derivative=True)
        # compare correct derivative val. (correct =-system val).
        assert_allclose(sy3, sphpy, atol=1.5e-4, rtol=0)


@pytest.mark.parametrize('derivative', [False, True])
@pytest.mark.parametrize('fun', [spherical_jn, spherical_in,
                                 spherical_yn, spherical_kn])
def test_negative_real_gh14582(derivative, fun):
    # gh-14582 reported that the spherical Bessel functions did not work
    # with negative real argument `z`. Check that this is resolved.
    rng = np.random.default_rng(3598435982345987234)
    size = 25
    n = rng.integers(0, 10, size=size)
    z = rng.standard_normal(size=size)
    res = fun(n, z, derivative=derivative)
    ref = fun(n, z+0j, derivative=derivative)
    assert_allclose(res, ref.real)
