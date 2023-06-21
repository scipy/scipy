import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_equal, assert_array_less

from scipy import stats, optimize
from scipy.optimize import _chandrupatla as zeros
import scipy.optimize._tstutils


class TestChandrupatla:

    def f(self, q, dist, p):
        return dist.cdf(q) - p

    @pytest.mark.parametrize('p', [0.6, np.linspace(-0.05, 1.05, 10)])
    def test_basic(self, p):
        # Invert distribution CDF and compare against distrtibution `ppf`
        dist = stats.norm()
        res = zeros._chandrupatla(self.f, -5, 5, args=(dist, p))
        ref = dist.ppf(p)
        np.testing.assert_allclose(res.x, ref)
        assert res.x.shape == ref.shape

    @pytest.mark.parametrize('shape', [tuple(), (12,), (3, 4), (3, 2, 2)])
    def test_vectorization(self, shape):
        # Test for correct functionality, output shapes, and dtypes for various
        # input shapes.
        p = np.linspace(-0.05, 1.05, 12).reshape(shape) if shape else 0.6
        dist = stats.norm()
        args = (dist, p)

        @np.vectorize
        def chandrupatla_single(p):
            return zeros._chandrupatla(self.f, -5, 5, args=(dist, p))

        def f(*args, **kwargs):
            f.f_evals += 1
            return self.f(*args, **kwargs)
        f.f_evals = 0

        res = zeros._chandrupatla(f, -5, 5, args=args)
        refs = chandrupatla_single(p).ravel()

        ref_x = [ref.x for ref in refs]
        assert_allclose(res.x.ravel(), ref_x)
        assert_equal(res.x.shape, shape)

        ref_fun = [ref.fun for ref in refs]
        assert_allclose(res.fun.ravel(), ref_fun)
        assert_equal(res.fun.shape, shape)
        assert_equal(res.fun, self.f(res.x, *args))

        ref_success = [ref.success for ref in refs]
        assert_equal(res.success.ravel(), ref_success)
        assert_equal(res.success.shape, shape)
        assert np.issubdtype(res.success.dtype, np.bool_)

        ref_flag = [ref.status for ref in refs]
        assert_equal(res.status.ravel(), ref_flag)
        assert_equal(res.status.shape, shape)
        assert np.issubdtype(res.status.dtype, np.integer)

        ref_nfev = [ref.nfev for ref in refs]
        assert_equal(res.nfev.ravel(), ref_nfev)
        assert_equal(np.max(res.nfev), f.f_evals)
        assert_equal(res.nfev.shape, res.fun.shape)
        assert np.issubdtype(res.nfev.dtype, np.integer)

        ref_nit = [ref.nit for ref in refs]
        assert_equal(res.nit.ravel(), ref_nit)
        assert_equal(np.max(res.nit), f.f_evals-2)
        assert_equal(res.nit.shape, res.fun.shape)
        assert np.issubdtype(res.nit.dtype, np.integer)

        ref_xl = [ref.xl for ref in refs]
        assert_allclose(res.xl.ravel(), ref_xl)
        assert_equal(res.xl.shape, shape)

        ref_xr = [ref.xr for ref in refs]
        assert_allclose(res.xr.ravel(), ref_xr)
        assert_equal(res.xr.shape, shape)

        assert_array_less(res.xl, res.xr)
        finite = np.isfinite(res.x)
        assert np.all((res.x[finite] == res.xl[finite])
                      | (res.x[finite] == res.xr[finite]))

        ref_fl = [ref.fl for ref in refs]
        assert_allclose(res.fl.ravel(), ref_fl)
        assert_equal(res.fl.shape, shape)
        assert_allclose(res.fl, self.f(res.xl, *args))

        ref_fr = [ref.fr for ref in refs]
        assert_allclose(res.fr.ravel(), ref_fr)
        assert_equal(res.fr.shape, shape)
        assert_allclose(res.fr, self.f(res.xr, *args))

        assert np.all(np.abs(res.fun[finite]) ==
                      np.minimum(np.abs(res.fl[finite]),
                                 np.abs(res.fr[finite])))


    def test_flags(self):
        # Test cases that should produce different status flags; show that all
        # can be produced simultaneously.
        def f(x):
            return [x[0] - 2.5, x[1] - 10, (x[2]-0.1)**3, np.nan]

        res = zeros._chandrupatla(f, [0] * 4, [np.pi] * 4, maxiter=2)

        ref_flags = np.array([zeros._ECONVERGED, zeros._ESIGNERR,
                              zeros._ECONVERR, zeros._EVALUEERR])
        assert_equal(res.status, ref_flags)

    def test_convergence(self):
        # Test that the convergence tolerances behave as expected
        rng = np.random.default_rng(2585255913088665241)
        p = rng.random(size=3)
        dist = stats.norm()
        bracket = (-5, 5)
        args = (dist, p)
        kwargs0 = dict(args=args, xatol=0, xrtol=0, fatol=0, frtol=0)

        kwargs = kwargs0.copy()
        kwargs['xatol'] = 1e-3
        res1 = zeros._chandrupatla(self.f, *bracket, **kwargs)
        assert_array_less(res1.xr - res1.xl, 1e-3)
        kwargs['xatol'] = 1e-6
        res2 = zeros._chandrupatla(self.f, *bracket, **kwargs)
        assert_array_less(res2.xr - res2.xl, 1e-6)
        assert_array_less(res2.xr - res2.xl, res1.xr - res1.xl)

        kwargs = kwargs0.copy()
        kwargs['xrtol'] = 1e-3
        res1 = zeros._chandrupatla(self.f, *bracket, **kwargs)
        assert_array_less(res1.xr - res1.xl, 1e-3 * np.abs(res1.x))
        kwargs['xrtol'] = 1e-6
        res2 = zeros._chandrupatla(self.f, *bracket, **kwargs)
        assert_array_less(res2.xr - res2.xl, 1e-6 * np.abs(res2.x))
        assert_array_less(res2.xr - res2.xl, res1.xr - res1.xl)

        kwargs = kwargs0.copy()
        kwargs['fatol'] = 1e-3
        res1 = zeros._chandrupatla(self.f, *bracket, **kwargs)
        assert_array_less(np.abs(res1.fun), 1e-3)
        kwargs['fatol'] = 1e-6
        res2 = zeros._chandrupatla(self.f, *bracket, **kwargs)
        assert_array_less(np.abs(res2.fun), 1e-6)
        assert_array_less(np.abs(res2.fun), np.abs(res1.fun))

        kwargs = kwargs0.copy()
        kwargs['frtol'] = 1e-3
        x1, x2 = bracket
        f0 = np.minimum(abs(self.f(x1, *args)), abs(self.f(x2, *args)))
        res1 = zeros._chandrupatla(self.f, *bracket, **kwargs)
        assert_array_less(np.abs(res1.fun), 1e-3*f0)
        kwargs['frtol'] = 1e-6
        res2 = zeros._chandrupatla(self.f, *bracket, **kwargs)
        assert_array_less(np.abs(res2.fun), 1e-6*f0)
        assert_array_less(np.abs(res2.fun), np.abs(res1.fun))

    def test_maxiter_callback(self):
        # Test behavior of `maxiter` parameter and `callback` interface
        p = 0.612814
        dist = stats.norm()
        bracket = (-5, 5)
        maxiter = 5

        res = zeros._chandrupatla(self.f, *bracket, args=(dist, p),
                                  maxiter=maxiter)
        assert not np.any(res.success)
        assert np.all(res.nfev == maxiter+2)
        assert np.all(res.nit == maxiter)

        def callback(res):
            callback.iter += 1
            callback.res = res
            assert hasattr(res, 'x')
            if callback.iter == 0:
                # callback is called once with initial bracket
                assert res.xl, res.xr == bracket
            assert res.status == zeros._EINPROGRESS
            if callback.iter == maxiter:
                raise StopIteration
        callback.iter = -1  # callback called once before first iteration
        callback.res = None

        res2 = zeros._chandrupatla(self.f, *bracket, args=(dist, p),
                                   callback=callback)

        # terminating with callback is identical to terminating due to maxiter
        # (except for `status`)
        for key in res.keys():
            if key == 'status':
                assert res[key] == zeros._ECONVERR
                assert callback.res[key] == zeros._EINPROGRESS
                assert res2[key] == zeros._ECALLBACK
            else:
                assert res2[key] == callback.res[key] == res[key]

    @pytest.mark.parametrize('case', optimize._tstutils._CHANDRUPATLA_TESTS)
    def test_nit_expected(self, case):
        # Test that `_chandrupatla` implements Chandrupatla's algorithm:
        # in all 40 test cases, the number of iterations performed
        # matches the number reported in the original paper.
        f, bracket, root, nfeval, id = case
        # Chandrupatla's criterion is equivalent to
        # abs(x2-x1) < 4*abs(xmin)*xrtol + xatol, but we use the more standard
        # abs(x2-x1) < abs(xmin)*xrtol + xatol. Therefore, set xrtol to 4x
        # that used by Chandrupatla in tests.
        res = zeros._chandrupatla(f, *bracket, xrtol=4e-10, xatol=1e-5)
        assert_allclose(res.fun, f(root), rtol=1e-8, atol=2e-3)
        assert_equal(res.nfev, nfeval)

    @pytest.mark.parametrize("dtype", (np.float16, np.float32, np.float64))
    def test_dtype(self, dtype):
        # Test that dtypes are preserved

        root = 0.622
        def f(x):
            return ((x - root) ** 3).astype(dtype)

        res = zeros._chandrupatla(f, dtype(-3), dtype(5), xatol=1e-3)
        assert res.x.dtype == dtype
        assert_allclose(res.x, root, atol=1e-3)

    def test_input_validation(self):
        # Test input validation for appropriate error messages

        message = '`func` must be callable.'
        with pytest.raises(ValueError, match=message):
            zeros._chandrupatla(None, -4, 4)

        message = 'Bracket and function output must be real numbers.'
        with pytest.raises(ValueError, match=message):
            zeros._chandrupatla(lambda x: x, -4+1j, 4)

        message = "shape mismatch: objects cannot be broadcast"
        # raised by `np.broadcast, but the traceback is readable IMO
        with pytest.raises(ValueError, match=message):
            zeros._chandrupatla(lambda x: x, [-2, -3], [3, 4, 5])
        with pytest.raises(ValueError, match=message):
            zeros._chandrupatla(lambda x: [x[0], x[1], x[1]], [-3, -3], [5, 5])

        message = 'Tolerances must be non-negative scalars.'
        with pytest.raises(ValueError, match=message):
            zeros._chandrupatla(lambda x: x, -4, 4, xatol=-1)
        with pytest.raises(ValueError, match=message):
            zeros._chandrupatla(lambda x: x, -4, 4, xrtol=None)
        with pytest.raises(ValueError, match=message):
            zeros._chandrupatla(lambda x: x, -4, 4, fatol='ekki')
        with pytest.raises(ValueError, match=message):
            zeros._chandrupatla(lambda x: x, -4, 4, frtol=None)

        message = '`maxiter` must be a non-negative integer.'
        with pytest.raises(ValueError, match=message):
            zeros._chandrupatla(lambda x: x, -4, 4, maxiter=1.5)
        with pytest.raises(ValueError, match=message):
            zeros._chandrupatla(lambda x: x, -4, 4, maxiter=-1)

        message = '`callback` must be callable.'
        with pytest.raises(ValueError, match=message):
            zeros._chandrupatla(lambda x: x, -4, 4, callback='shrubbery')

    def test_special_cases(self):
        # Test edge cases and other special cases

        # Test that integers are not passed to `f`
        # (otherwise this would overflow)
        def f(x):
            # assert np.issubdtype(x.dtype, np.floating)
            return x ** 99 - 1

        res = zeros._chandrupatla(f, -7, 5)
        assert res.success
        assert_allclose(res.x, 1)

        # Test that if both ends of bracket equal root, algorithm reports
        # convergence
        def f(x):
            return x**2 - 1

        res = zeros._chandrupatla(f, 1, 1)
        # assert res.success
        assert_equal(res.x, 1)

        def f(x):
            return 1/x

        with np.errstate(invalid='ignore'):
            res = zeros._chandrupatla(f, np.inf, np.inf)
        assert res.success
        assert_equal(res.x, np.inf)

        # Test maxiter = 0. Should do nothing to bracket.
        def f(x):
            return x**3 - 1

        bracket = (-3, 5)
        res = zeros._chandrupatla(f, *bracket, maxiter=0)
        assert res.xl, res.xr == bracket
        assert res.nit == 0
        assert res.nfev == 2
        assert res.status == -2
        assert res.x == -3  # best so far

        # Test maxiter = 1
        res = zeros._chandrupatla(f, *bracket, maxiter=1)
        assert res.success
        assert res.status == 0
        assert res.nit == 1
        assert res.nfev == 3
        assert_allclose(res.x, 1)

        # Test scalar `args` (not in tuple)
        def f(x, c):
            return c*x - 1

        res = zeros._chandrupatla(f, -1, 1, args=3)
        assert_allclose(res.x, 1/3)

        # # TODO: Test zero tolerance
        # # ~~What's going on here - why are iterations repeated?~~
        # # tl goes to zero when xatol=xrtol=0. When function is nearly linear,
        # # this causes convergence issues.
        # def f(x):
        #     return np.cos(x)
        #
        # res = zeros._chandrupatla(f, 0, np.pi, xatol=0, xrtol=0)
        # assert res.nit < 100
        # xp = np.nextafter(res.x, np.inf)
        # xm = np.nextafter(res.x, -np.inf)
        # assert np.abs(res.fun) < np.abs(f(xp))
        # assert np.abs(res.fun) < np.abs(f(xm))
