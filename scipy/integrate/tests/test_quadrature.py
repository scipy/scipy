# mypy: disable-error-code="attr-defined"
import pytest
import numpy as np
from numpy import cos, sin, pi
from numpy.testing import (assert_equal, assert_almost_equal, assert_allclose,
                           assert_, suppress_warnings)

from scipy.integrate import (quadrature, romberg, romb, newton_cotes,
                             cumulative_trapezoid, cumtrapz, trapz, trapezoid,
                             quad, simpson, simps, fixed_quad, AccuracyWarning,
                             qmc_quad)
from scipy.integrate._tanhsinh import _tanhsinh
from scipy import stats, special as sc


class TestFixedQuad:
    def test_scalar(self):
        n = 4
        expected = 1/(2*n)
        got, _ = fixed_quad(lambda x: x**(2*n - 1), 0, 1, n=n)
        # quadrature exact for this input
        assert_allclose(got, expected, rtol=1e-12)

    def test_vector(self):
        n = 4
        p = np.arange(1, 2*n)
        expected = 1/(p + 1)
        got, _ = fixed_quad(lambda x: x**p[:, None], 0, 1, n=n)
        assert_allclose(got, expected, rtol=1e-12)


class TestQuadrature:
    def quad(self, x, a, b, args):
        raise NotImplementedError

    def test_quadrature(self):
        # Typical function with two extra arguments:
        def myfunc(x, n, z):       # Bessel function integrand
            return cos(n*x-z*sin(x))/pi
        val, err = quadrature(myfunc, 0, pi, (2, 1.8))
        table_val = 0.30614353532540296487
        assert_almost_equal(val, table_val, decimal=7)

    def test_quadrature_rtol(self):
        def myfunc(x, n, z):       # Bessel function integrand
            return 1e90 * cos(n*x-z*sin(x))/pi
        val, err = quadrature(myfunc, 0, pi, (2, 1.8), rtol=1e-10)
        table_val = 1e90 * 0.30614353532540296487
        assert_allclose(val, table_val, rtol=1e-10)

    def test_quadrature_miniter(self):
        # Typical function with two extra arguments:
        def myfunc(x, n, z):       # Bessel function integrand
            return cos(n*x-z*sin(x))/pi
        table_val = 0.30614353532540296487
        for miniter in [5, 52]:
            val, err = quadrature(myfunc, 0, pi, (2, 1.8), miniter=miniter)
            assert_almost_equal(val, table_val, decimal=7)
            assert_(err < 1.0)

    def test_quadrature_single_args(self):
        def myfunc(x, n):
            return 1e90 * cos(n*x-1.8*sin(x))/pi
        val, err = quadrature(myfunc, 0, pi, args=2, rtol=1e-10)
        table_val = 1e90 * 0.30614353532540296487
        assert_allclose(val, table_val, rtol=1e-10)

    def test_romberg(self):
        # Typical function with two extra arguments:
        def myfunc(x, n, z):       # Bessel function integrand
            return cos(n*x-z*sin(x))/pi
        val = romberg(myfunc, 0, pi, args=(2, 1.8))
        table_val = 0.30614353532540296487
        assert_almost_equal(val, table_val, decimal=7)

    def test_romberg_rtol(self):
        # Typical function with two extra arguments:
        def myfunc(x, n, z):       # Bessel function integrand
            return 1e19*cos(n*x-z*sin(x))/pi
        val = romberg(myfunc, 0, pi, args=(2, 1.8), rtol=1e-10)
        table_val = 1e19*0.30614353532540296487
        assert_allclose(val, table_val, rtol=1e-10)

    def test_romb(self):
        assert_equal(romb(np.arange(17)), 128)

    def test_romb_gh_3731(self):
        # Check that romb makes maximal use of data points
        x = np.arange(2**4+1)
        y = np.cos(0.2*x)
        val = romb(y)
        val2, err = quad(lambda x: np.cos(0.2*x), x.min(), x.max())
        assert_allclose(val, val2, rtol=1e-8, atol=0)

        # should be equal to romb with 2**k+1 samples
        with suppress_warnings() as sup:
            sup.filter(AccuracyWarning, "divmax .4. exceeded")
            val3 = romberg(lambda x: np.cos(0.2*x), x.min(), x.max(), divmax=4)
        assert_allclose(val, val3, rtol=1e-12, atol=0)

    def test_non_dtype(self):
        # Check that we work fine with functions returning float
        import math
        valmath = romberg(math.sin, 0, 1)
        expected_val = 0.45969769413185085
        assert_almost_equal(valmath, expected_val, decimal=7)

    def test_newton_cotes(self):
        """Test the first few degrees, for evenly spaced points."""
        n = 1
        wts, errcoff = newton_cotes(n, 1)
        assert_equal(wts, n*np.array([0.5, 0.5]))
        assert_almost_equal(errcoff, -n**3/12.0)

        n = 2
        wts, errcoff = newton_cotes(n, 1)
        assert_almost_equal(wts, n*np.array([1.0, 4.0, 1.0])/6.0)
        assert_almost_equal(errcoff, -n**5/2880.0)

        n = 3
        wts, errcoff = newton_cotes(n, 1)
        assert_almost_equal(wts, n*np.array([1.0, 3.0, 3.0, 1.0])/8.0)
        assert_almost_equal(errcoff, -n**5/6480.0)

        n = 4
        wts, errcoff = newton_cotes(n, 1)
        assert_almost_equal(wts, n*np.array([7.0, 32.0, 12.0, 32.0, 7.0])/90.0)
        assert_almost_equal(errcoff, -n**7/1935360.0)

    def test_newton_cotes2(self):
        """Test newton_cotes with points that are not evenly spaced."""

        x = np.array([0.0, 1.5, 2.0])
        y = x**2
        wts, errcoff = newton_cotes(x)
        exact_integral = 8.0/3
        numeric_integral = np.dot(wts, y)
        assert_almost_equal(numeric_integral, exact_integral)

        x = np.array([0.0, 1.4, 2.1, 3.0])
        y = x**2
        wts, errcoff = newton_cotes(x)
        exact_integral = 9.0
        numeric_integral = np.dot(wts, y)
        assert_almost_equal(numeric_integral, exact_integral)

    # ignore the DeprecationWarning emitted by the even kwd
    @pytest.mark.filterwarnings('ignore::DeprecationWarning')
    def test_simpson(self):
        y = np.arange(17)
        assert_equal(simpson(y), 128)
        assert_equal(simpson(y, dx=0.5), 64)
        assert_equal(simpson(y, x=np.linspace(0, 4, 17)), 32)

        y = np.arange(4)
        x = 2**y
        assert_equal(simpson(y, x=x, even='avg'), 13.875)
        assert_equal(simpson(y, x=x, even='first'), 13.75)
        assert_equal(simpson(y, x=x, even='last'), 14)

        # `even='simpson'`
        # integral should be exactly 21
        x = np.linspace(1, 4, 4)
        def f(x):
            return x**2

        assert_allclose(simpson(f(x), x=x, even='simpson'), 21.0)
        assert_allclose(simpson(f(x), x=x, even='avg'), 21 + 1/6)

        # integral should be exactly 114
        x = np.linspace(1, 7, 4)
        assert_allclose(simpson(f(x), dx=2.0, even='simpson'), 114)
        assert_allclose(simpson(f(x), dx=2.0, even='avg'), 115 + 1/3)

        # `even='simpson'`, test multi-axis behaviour
        a = np.arange(16).reshape(4, 4)
        x = np.arange(64.).reshape(4, 4, 4)
        y = f(x)
        for i in range(3):
            r = simpson(y, x=x, even='simpson', axis=i)
            it = np.nditer(a, flags=['multi_index'])
            for _ in it:
                idx = list(it.multi_index)
                idx.insert(i, slice(None))
                integral = x[tuple(idx)][-1]**3 / 3 - x[tuple(idx)][0]**3 / 3
                assert_allclose(r[it.multi_index], integral)

        # test when integration axis only has two points
        x = np.arange(16).reshape(8, 2)
        y = f(x)
        for even in ['simpson', 'avg', 'first', 'last']:
            r = simpson(y, x=x, even=even, axis=-1)

            integral = 0.5 * (y[:, 1] + y[:, 0]) * (x[:, 1] - x[:, 0])
            assert_allclose(r, integral)

        # odd points, test multi-axis behaviour
        a = np.arange(25).reshape(5, 5)
        x = np.arange(125).reshape(5, 5, 5)
        y = f(x)
        for i in range(3):
            r = simpson(y, x=x, axis=i)
            it = np.nditer(a, flags=['multi_index'])
            for _ in it:
                idx = list(it.multi_index)
                idx.insert(i, slice(None))
                integral = x[tuple(idx)][-1]**3 / 3 - x[tuple(idx)][0]**3 / 3
                assert_allclose(r[it.multi_index], integral)

        # Tests for checking base case
        x = np.array([3])
        y = np.power(x, 2)
        assert_allclose(simpson(y, x=x, axis=0), 0.0)
        assert_allclose(simpson(y, x=x, axis=-1), 0.0)

        x = np.array([3, 3, 3, 3])
        y = np.power(x, 2)
        assert_allclose(simpson(y, x=x, axis=0), 0.0)
        assert_allclose(simpson(y, x=x, axis=-1), 0.0)

        x = np.array([[1, 2, 4, 8], [1, 2, 4, 8], [1, 2, 4, 8]])
        y = np.power(x, 2)
        zero_axis = [0.0, 0.0, 0.0, 0.0]
        default_axis = [170 + 1/3] * 3   # 8**3 / 3 - 1/3
        assert_allclose(simpson(y, x=x, axis=0), zero_axis)
        # the following should be exact for even='simpson'
        assert_allclose(simpson(y, x=x, axis=-1), default_axis)

        x = np.array([[1, 2, 4, 8], [1, 2, 4, 8], [1, 8, 16, 32]])
        y = np.power(x, 2)
        zero_axis = [0.0, 136.0, 1088.0, 8704.0]
        default_axis = [170 + 1/3, 170 + 1/3, 32**3 / 3 - 1/3]
        assert_allclose(simpson(y, x=x, axis=0), zero_axis)
        assert_allclose(simpson(y, x=x, axis=-1), default_axis)

    def test_simpson_even_is_deprecated(self):
        x = np.linspace(0, 3, 4)
        y = x**2
        with pytest.deprecated_call():
            simpson(y, x=x, even='first')

    @pytest.mark.parametrize('droplast', [False, True])
    def test_simpson_2d_integer_no_x(self, droplast):
        # The inputs are 2d integer arrays.  The results should be
        # identical to the results when the inputs are floating point.
        y = np.array([[2, 2, 4, 4, 8, 8, -4, 5],
                      [4, 4, 2, -4, 10, 22, -2, 10]])
        if droplast:
            y = y[:, :-1]
        result = simpson(y, axis=-1)
        expected = simpson(np.array(y, dtype=np.float64), axis=-1)
        assert_equal(result, expected)

    def test_simps(self):
        # Basic coverage test for the alias
        y = np.arange(5)
        x = 2**y
        assert_allclose(
            simpson(y, x=x, dx=0.5),
            simps(y, x=x, dx=0.5)
        )


class TestCumulative_trapezoid:
    def test_1d(self):
        x = np.linspace(-2, 2, num=5)
        y = x
        y_int = cumulative_trapezoid(y, x, initial=0)
        y_expected = [0., -1.5, -2., -1.5, 0.]
        assert_allclose(y_int, y_expected)

        y_int = cumulative_trapezoid(y, x, initial=None)
        assert_allclose(y_int, y_expected[1:])

    def test_y_nd_x_nd(self):
        x = np.arange(3 * 2 * 4).reshape(3, 2, 4)
        y = x
        y_int = cumulative_trapezoid(y, x, initial=0)
        y_expected = np.array([[[0., 0.5, 2., 4.5],
                                [0., 4.5, 10., 16.5]],
                               [[0., 8.5, 18., 28.5],
                                [0., 12.5, 26., 40.5]],
                               [[0., 16.5, 34., 52.5],
                                [0., 20.5, 42., 64.5]]])

        assert_allclose(y_int, y_expected)

        # Try with all axes
        shapes = [(2, 2, 4), (3, 1, 4), (3, 2, 3)]
        for axis, shape in zip([0, 1, 2], shapes):
            y_int = cumulative_trapezoid(y, x, initial=3.45, axis=axis)
            assert_equal(y_int.shape, (3, 2, 4))
            y_int = cumulative_trapezoid(y, x, initial=None, axis=axis)
            assert_equal(y_int.shape, shape)

    def test_y_nd_x_1d(self):
        y = np.arange(3 * 2 * 4).reshape(3, 2, 4)
        x = np.arange(4)**2
        # Try with all axes
        ys_expected = (
            np.array([[[4., 5., 6., 7.],
                       [8., 9., 10., 11.]],
                      [[40., 44., 48., 52.],
                       [56., 60., 64., 68.]]]),
            np.array([[[2., 3., 4., 5.]],
                      [[10., 11., 12., 13.]],
                      [[18., 19., 20., 21.]]]),
            np.array([[[0.5, 5., 17.5],
                       [4.5, 21., 53.5]],
                      [[8.5, 37., 89.5],
                       [12.5, 53., 125.5]],
                      [[16.5, 69., 161.5],
                       [20.5, 85., 197.5]]]))

        for axis, y_expected in zip([0, 1, 2], ys_expected):
            y_int = cumulative_trapezoid(y, x=x[:y.shape[axis]], axis=axis,
                                         initial=None)
            assert_allclose(y_int, y_expected)

    def test_x_none(self):
        y = np.linspace(-2, 2, num=5)

        y_int = cumulative_trapezoid(y)
        y_expected = [-1.5, -2., -1.5, 0.]
        assert_allclose(y_int, y_expected)

        y_int = cumulative_trapezoid(y, initial=1.23)
        y_expected = [1.23, -1.5, -2., -1.5, 0.]
        assert_allclose(y_int, y_expected)

        y_int = cumulative_trapezoid(y, dx=3)
        y_expected = [-4.5, -6., -4.5, 0.]
        assert_allclose(y_int, y_expected)

        y_int = cumulative_trapezoid(y, dx=3, initial=1.23)
        y_expected = [1.23, -4.5, -6., -4.5, 0.]
        assert_allclose(y_int, y_expected)

    def test_cumtrapz(self):
        # Basic coverage test for the alias
        x = np.arange(3 * 2 * 4).reshape(3, 2, 4)
        y = x
        assert_allclose(cumulative_trapezoid(y, x, dx=0.5, axis=0, initial=0),
                        cumtrapz(y, x, dx=0.5, axis=0, initial=0),
                        rtol=1e-14)


class TestTrapezoid:
    """This function is tested in NumPy more extensive, just do some
    basic due diligence here."""
    def test_trapezoid(self):
        y = np.arange(17)
        assert_equal(trapezoid(y), 128)
        assert_equal(trapezoid(y, dx=0.5), 64)
        assert_equal(trapezoid(y, x=np.linspace(0, 4, 17)), 32)

        y = np.arange(4)
        x = 2**y
        assert_equal(trapezoid(y, x=x, dx=0.1), 13.5)

    def test_trapz(self):
        # Basic coverage test for the alias
        y = np.arange(4)
        x = 2**y
        assert_equal(trapezoid(y, x=x, dx=0.5, axis=0),
                     trapz(y, x=x, dx=0.5, axis=0))


class TestQMCQuad:
    def test_input_validation(self):
        message = "`func` must be callable."
        with pytest.raises(TypeError, match=message):
            qmc_quad("a duck", [0, 0], [1, 1])

        message = "`func` must evaluate the integrand at points..."
        with pytest.raises(ValueError, match=message):
            qmc_quad(lambda: 1, [0, 0], [1, 1])

        def func(x):
            assert x.ndim == 1
            return np.sum(x)
        message = "Exception encountered when attempting vectorized call..."
        with pytest.warns(UserWarning, match=message):
            qmc_quad(func, [0, 0], [1, 1])

        message = "`n_points` must be an integer."
        with pytest.raises(TypeError, match=message):
            qmc_quad(lambda x: 1, [0, 0], [1, 1], n_points=1024.5)

        message = "`n_estimates` must be an integer."
        with pytest.raises(TypeError, match=message):
            qmc_quad(lambda x: 1, [0, 0], [1, 1], n_estimates=8.5)

        message = "`qrng` must be an instance of scipy.stats.qmc.QMCEngine."
        with pytest.raises(TypeError, match=message):
            qmc_quad(lambda x: 1, [0, 0], [1, 1], qrng="a duck")

        message = "`qrng` must be initialized with dimensionality equal to "
        with pytest.raises(ValueError, match=message):
            qmc_quad(lambda x: 1, [0, 0], [1, 1], qrng=stats.qmc.Sobol(1))

        message = r"`log` must be boolean \(`True` or `False`\)."
        with pytest.raises(TypeError, match=message):
            qmc_quad(lambda x: 1, [0, 0], [1, 1], log=10)

    def basic_test(self, n_points=2**8, n_estimates=8, signs=np.ones(2)):

        ndim = 2
        mean = np.zeros(ndim)
        cov = np.eye(ndim)

        def func(x):
            return stats.multivariate_normal.pdf(x.T, mean, cov)

        rng = np.random.default_rng(2879434385674690281)
        qrng = stats.qmc.Sobol(ndim, seed=rng)
        a = np.zeros(ndim)
        b = np.ones(ndim) * signs
        res = qmc_quad(func, a, b, n_points=n_points,
                       n_estimates=n_estimates, qrng=qrng)
        ref = stats.multivariate_normal.cdf(b, mean, cov, lower_limit=a)
        atol = sc.stdtrit(n_estimates-1, 0.995) * res.standard_error  # 99% CI
        assert_allclose(res.integral, ref, atol=atol)
        assert np.prod(signs)*res.integral > 0

        rng = np.random.default_rng(2879434385674690281)
        qrng = stats.qmc.Sobol(ndim, seed=rng)
        logres = qmc_quad(lambda *args: np.log(func(*args)), a, b,
                          n_points=n_points, n_estimates=n_estimates,
                          log=True, qrng=qrng)
        assert_allclose(np.exp(logres.integral), res.integral, rtol=1e-14)
        assert np.imag(logres.integral) == (np.pi if np.prod(signs) < 0 else 0)
        assert_allclose(np.exp(logres.standard_error),
                        res.standard_error, rtol=1e-14, atol=1e-16)

    @pytest.mark.parametrize("n_points", [2**8, 2**12])
    @pytest.mark.parametrize("n_estimates", [8, 16])
    def test_basic(self, n_points, n_estimates):
        self.basic_test(n_points, n_estimates)

    @pytest.mark.parametrize("signs", [[1, 1], [-1, -1], [-1, 1], [1, -1]])
    def test_sign(self, signs):
        self.basic_test(signs=signs)

    @pytest.mark.parametrize("log", [False, True])
    def test_zero(self, log):
        message = "A lower limit was equal to an upper limit, so"
        with pytest.warns(UserWarning, match=message):
            res = qmc_quad(lambda x: 1, [0, 0], [0, 1], log=log)
        assert res.integral == (-np.inf if log else 0)
        assert res.standard_error == 0

    def test_flexible_input(self):
        # check that qrng is not required
        # also checks that for 1d problems, a and b can be scalars
        def func(x):
            return stats.norm.pdf(x, scale=2)

        res = qmc_quad(func, 0, 1)
        ref = stats.norm.cdf(1, scale=2) - stats.norm.cdf(0, scale=2)
        assert_allclose(res.integral, ref, 1e-2)


class TestTanhSinh:
    # Test problems from [1] Section 6
    def f1(self, t):
        return t * np.log(1 + t)

    f1.ref = 0.25
    f1.b = 1

    def f2(self, t):
        return t ** 2 * np.arctan(t)

    f2.ref = (np.pi - 2 + 2 * np.log(2)) / 12
    f2.b = 1

    def f3(self, t):
        return np.exp(t) * np.cos(t)

    f3.ref = (np.exp(np.pi / 2) - 1) / 2
    f3.b = np.pi / 2

    def f4(self, t):
        a = np.sqrt(2 + t ** 2)
        return np.arctan(a) / ((1 + t ** 2) * a)

    f4.ref = 5 * np.pi ** 2 / 96
    f4.b = 1

    def f5(self, t):
        return np.sqrt(t) * np.log(t)

    f5.ref = -4 / 9
    f5.b = 1

    def f6(self, t):
        return np.sqrt(1 - t ** 2)

    f6.ref = np.pi / 4
    f6.b = 1

    def f7(self, t):
        return np.sqrt(t) / np.sqrt(1 - t ** 2)

    f7.ref = 2 * np.sqrt(np.pi) * sc.gamma(3 / 4) / sc.gamma(1 / 4)
    f7.b = 1

    def f8(self, t):
        return np.log(t) ** 2

    f8.ref = 2
    f8.b = 1

    def f9(self, t):
        return np.log(np.cos(t))

    f9.ref = -np.pi * np.log(2) / 2
    f9.b = np.pi / 2

    def f10(self, t):
        return np.sqrt(np.tan(t))

    f10.ref = np.pi * np.sqrt(2) / 2
    f10.b = np.pi / 2

    def f11(self, t):
        return 1 / (1 + t ** 2)

    f11.ref = np.pi / 2
    f11.b = np.inf

    def f12(self, t):
        return np.exp(-t) / np.sqrt(t)

    f12.ref = np.sqrt(np.pi)
    f12.b = np.inf

    def f13(self, t):
        return np.exp(-t ** 2 / 2)

    f13.ref = np.sqrt(np.pi / 2)
    f13.b = np.inf

    def f14(self, t):
        return np.exp(-t) * np.cos(t)

    f14.ref = 0.5
    f14.b = np.inf

    def f15(self, t):
        return np.sin(t) / t

    f15.ref = np.pi / 2
    f15.b = np.inf

    def error(self, res, ref, log=False):
        err = abs(res - ref)

        if not log:
            return err

        with np.errstate(divide='ignore'):
            return np.log10(err)

    def test_input_validation(self):
        f = self.f1

        message = '`f` must be callable.'
        with pytest.raises(ValueError, match=message):
            _tanhsinh(42, 0, f.b)

        message = '...must be True or False.'
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, log=2)

        message = '...must be reals.'
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 1+1j, f.b)
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, np.nan)
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, atol='ekki')
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, rtol=pytest)
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, minweight=object())

        message = '...must be positive and finite.'
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, minweight=np.inf)

        message = '...must be non-negative and finite.'
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, rtol=-1)
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, atol=np.inf)

        message = '...may not be positive infinity.'
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, rtol=np.inf, log=True)
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, atol=np.inf, log=True)

        message = '...must be integers.'
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, maxlevel=object())
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, maxfun=1+1j)
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, minlevel="migratory coconut")

        message = '...must be non-negative.'
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, maxlevel=-1)
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, maxfun=-1)
        with pytest.raises(ValueError, match=message):
            _tanhsinh(f, 0, f.b, minlevel=-1)

    @pytest.mark.parametrize("limits, val", [
        [(0, np.inf), 0.5],  # b infinite
        [(-np.inf, 0), 0.5],  # a infinite
        [(-np.inf, np.inf), 1],  # a and b infinite
        [(np.inf, -np.inf), -1],  # flipped limits
        [(1, -1), stats.norm.cdf(-1) -  stats.norm.cdf(1)],  # flipped limits
    ])
    def test_integral_transforms(self, limits, val):
        # Check that the integral transforms are behaving for both log and
        # normal integration
        dist = stats.norm()

        res = _tanhsinh(dist.pdf, *limits)
        assert_allclose(res.integral, val)

        logres = _tanhsinh(dist.logpdf, *limits, log=True)
        assert_allclose(np.exp(logres.integral), val)
        # Transformation should not make the result complex unnecessarily
        assert (np.isreal(logres.integral) if val > 0
                else np.iscomplex(logres.integral))

        assert_allclose(np.exp(logres.error), res.error, atol=1e-16)

    # 15 skipped intentionally; it's very difficult numerically
    @pytest.mark.parametrize('f_number', range(1, 15))
    def test_basic(self, f_number):
        f = getattr(self, f"f{f_number}")
        res = _tanhsinh(f, 0, f.b)
        assert_allclose(res.integral, f.ref)

    def test_convergence(self):
        # demonstrate that number of accurate digits doubles each iteration
        f = self.f1
        last_logerr = 0
        for i in range(4):
            res = _tanhsinh(f, 0, f.b, minlevel=0, maxlevel=i)
            logerr = self.error(res.integral, f.ref, log=True)
            assert (logerr < last_logerr * 2 or logerr < -15.5)
            last_logerr = logerr

    def test_feval(self):
        # Test function evaluation count
        dist = stats.norm()
        def f(x):
            f.feval += len(x)
            return dist.pdf(x)

        f.feval = 0
        res = _tanhsinh(f, 0, np.inf)
        assert_allclose(res.integral, 0.5)
        assert res.feval == f.feval

        f.feval = 0
        res = _tanhsinh(f, -np.inf, np.inf)
        assert_allclose(res.integral, 1)
        assert res.feval == f.feval

    def test_options_and_result_attributes(self):
        # demonstrate that options are behaving as advertised and status
        # messages are as intended
        def f(x):
            f.calls += 1
            f.feval += len(x)
            return self.f2(x)
        f.ref = self.f2.ref
        f.b = self.f2.b
        default_rtol = 1e-12
        default_atol = f.ref * default_rtol  # effective default absolute tol

        # Test default options
        f.feval, f.calls = 0, 0
        ref = _tanhsinh(f, 0, f.b)
        assert self.error(ref.integral, f.ref) < ref.error < default_atol
        assert ref.feval == f.feval
        ref.calls = f.calls  # reference number of function calls
        assert ref.success is True
        assert ref.status == 0
        assert ref.message.startswith("The algorithm completed successfully")

        # Test `maxlevel` equal to required number of function evaluations
        # We should get all the same results
        f.feval, f.calls = 0, 0
        maxlevel = ref.calls + 1  # default is first call goes up to level 2
        res = _tanhsinh(f, 0, f.b, maxlevel=maxlevel)
        assert res == ref

        # Now reduce the maximum level. We won't meet tolerances.
        f.feval, f.calls = 0, 0
        maxlevel -= 1
        assert maxlevel >= 2  # can't compare errors otherwise
        res = _tanhsinh(f, 0, f.b, maxlevel=maxlevel)
        assert self.error(res.integral, f.ref) < res.error > default_atol
        assert res.feval == f.feval < ref.feval
        assert f.calls == ref.calls - 1
        assert res.success is False
        assert res.status == 1
        assert res.message.endswith("maximum level to be exceeded.")

        # Test `maxfun` equal to required number of function evaluations
        # We should get all the same results
        f.feval, f.calls = 0, 0
        maxfun = ref.feval
        res = _tanhsinh(f, 0, f.b, maxfun = maxfun)
        assert res == ref

        # Now reduce `maxfun`. We won't meet tolerances.
        f.feval, f.calls = 0, 0
        maxfun -= 1
        res = _tanhsinh(f, 0, f.b, maxfun=maxfun)
        assert self.error(res.integral, f.ref) < res.error > default_atol
        assert res.feval == f.feval < ref.feval
        assert f.calls == ref.calls - 1
        assert res.success is False
        assert res.status == 2
        assert res.message.endswith("evaluation limit to be exceeded.")

        # Take this result to be the new reference
        ref = res
        ref.calls = f.calls

        # Test `atol`
        f.feval, f.calls = 0, 0
        # With this tolerance, we should get the exact same result as ref
        atol = np.nextafter(ref.error, np.inf)
        res = _tanhsinh(f, 0, f.b, atol=atol)
        assert res.integral == ref.integral
        assert res.error == ref.error
        assert res.feval == f.feval == ref.feval
        assert f.calls == ref.calls
        # Except the result is considered to be successful
        assert res.success is True
        assert res.status == 0
        assert res.message.startswith("The algorithm completed successfully")

        f.feval, f.calls = 0, 0
        # With a tighter tolerance, we should get a more accurate result
        atol = np.nextafter(ref.error, -np.inf)
        res = _tanhsinh(f, 0, f.b, atol=atol)
        assert self.error(res.integral, f.ref) < res.error < atol
        assert res.feval == f.feval > ref.feval
        assert f.calls > ref.calls
        assert res.success is True
        assert res.status == 0
        assert res.message.startswith("The algorithm completed successfully")

        # Test `rtol`
        f.feval, f.calls = 0, 0
        # With this tolerance, we should get the exact same result as ref
        rtol = np.nextafter(ref.error/ref.integral, np.inf)
        res = _tanhsinh(f, 0, f.b, rtol=rtol)
        assert res.integral == ref.integral
        assert res.error == ref.error
        assert res.feval == f.feval == ref.feval
        assert f.calls == ref.calls
        # Except the result is considered to be successful
        assert res.success is True
        assert res.status == 0
        assert res.message.startswith("The algorithm completed successfully")

        f.feval, f.calls = 0, 0
        # With a tighter tolerance, we should get a more accurate result
        rtol = np.nextafter(ref.error/ref.integral, -np.inf)
        res = _tanhsinh(f, 0, f.b, rtol=rtol)
        assert self.error(res.integral, f.ref)/f.ref < res.error/res.integral < rtol
        assert res.feval == f.feval > ref.feval
        assert f.calls > ref.calls
        assert res.success is True
        assert res.status == 0
        assert res.message.startswith("The algorithm completed successfully")

        # Test `minweight` and status 3
        f = self.f11
        # The transformed integrand produces a NaN when evaluated close to zero
        # The weight is very small here, so with normal options, the NaN is
        # replaced with zero. However, if we set `minweight` too small, we get
        # status 3.
        res = _tanhsinh(f, 0, f.b, minweight=1e-300)
        assert res.success is False
        assert res.status == 3
        assert res.message.startswith("An invalid value")
        assert np.isnan(res.integral)
        assert np.isnan(res.error)

    @pytest.mark.parametrize('rtol', [1e-4, 1e-14])
    def test_log(self, rtol):
        # Test equivalence of log-integration and regular integration
        dist = stats.norm()

        # Problem with positive integrand/real log-integrand)
        res = _tanhsinh(dist.logpdf, -1, 2, log=True, rtol=np.log(rtol))
        ref = _tanhsinh(dist.pdf, -1, 2, rtol=rtol)
        assert_allclose(np.exp(res.integral), ref.integral)
        assert_allclose(np.exp(res.error), ref.error)
        assert res.feval == ref.feval

        # Problem with real integrand/complex log-integrand
        def f(x):
            return -dist.logpdf(x)*dist.pdf(x)

        def logf(x):
            return np.log(dist.logpdf(x) + 0j) + dist.logpdf(x) + np.pi * 1j

        res = _tanhsinh(logf, -np.inf, np.inf, log=True)
        ref = _tanhsinh(f, -np.inf, np.inf)
        assert_allclose(np.exp(res.integral), ref.integral)
        assert_allclose(np.exp(res.error), ref.error)
        assert res.feval == ref.feval

    def test_complex(self):
        # Test case with finite limits
        def f(x):
            return np.exp(1j * x)

        a, b = 0, np.pi/4
        res = _tanhsinh(f, a, b)
        ref = np.sqrt(2)/2 + (1-np.sqrt(2)/2)*1j
        assert_allclose(res.integral, ref)

        # Test case involving a few transformations
        dist1 = stats.norm(scale=1)
        dist2 = stats.norm(scale=2)
        def f(x):
            return dist1.pdf(x) + 1j*dist2.pdf(x)
        res = _tanhsinh(f, np.inf, -np.inf)
        assert_allclose(res.integral, -(1+1j))

    @pytest.mark.parametrize("maxlevel", range(4))
    def test_minlevel(self, maxlevel):
        # Verify that minlevel does not change the values at which the
        # integrand is evaluated or the integral/error estimates, only the
        # number of function calls
        def f(x):
            f.calls += 1
            f.feval += len(x)
            f.x = np.concatenate((f.x, x))
            return self.f2(x)
        f.feval, f.calls, f.x = 0, 0, np.array([])

        ref = _tanhsinh(f, 0, self.f2.b, minlevel=0, maxlevel=maxlevel)
        ref_x = np.sort(f.x)

        for minlevel in range(0, maxlevel + 1):
            f.feval, f.calls, f.x = 0, 0, np.array([])
            options = dict(minlevel=minlevel, maxlevel=maxlevel)
            res = _tanhsinh(f, 0, self.f2.b, **options)
            # Should be very close; all that has changed is the order of values
            assert_allclose(res.integral, ref.integral, rtol=4e-16)
            # Difference in absolute errors << magnitude of integral
            assert_allclose(res.error, ref.error, atol=4e-16 * ref.integral)
            # Note: previous line is not satisfied for all problems. See note
            # corresponding with `d4` in `_tanhsinh._estimate_error`.
            assert res.feval == f.feval == len(f.x)
            assert f.calls == maxlevel - minlevel + 1
            assert res.status == ref.status
            assert_equal(ref_x, np.sort(f.x))
