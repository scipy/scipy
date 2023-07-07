import pytest

from functools import lru_cache

from numpy.testing import (assert_warns, assert_,
                           assert_allclose,
                           assert_equal,
                           assert_array_equal,
                           assert_array_less,
                           suppress_warnings)
import numpy as np
from numpy import finfo, power, nan, isclose, sqrt, exp, sin, cos

from scipy import stats, optimize
from scipy.optimize import (_zeros_py as zeros, newton, root_scalar,
                            OptimizeResult)

from scipy._lib._util import getfullargspec_no_self as _getfullargspec

# Import testing parameters
from scipy.optimize._tstutils import get_tests, functions as tstutils_functions

TOL = 4*np.finfo(float).eps  # tolerance

_FLOAT_EPS = finfo(float).eps

bracket_methods = [zeros.bisect, zeros.ridder, zeros.brentq, zeros.brenth,
                   zeros.toms748]
gradient_methods = [zeros.newton]
all_methods = bracket_methods + gradient_methods  # noqa

# A few test functions used frequently:
# # A simple quadratic, (x-1)^2 - 1
def f1(x):
    return x ** 2 - 2 * x - 1


def f1_1(x):
    return 2 * x - 2


def f1_2(x):
    return 2.0 + 0 * x


def f1_and_p_and_pp(x):
    return f1(x), f1_1(x), f1_2(x)


# Simple transcendental function
def f2(x):
    return exp(x) - cos(x)


def f2_1(x):
    return exp(x) + sin(x)


def f2_2(x):
    return exp(x) + cos(x)


# lru cached function
@lru_cache
def f_lrucached(x):
    return x


class TestScalarRootFinders:
    # Basic tests for all scalar root finders

    xtol = 4 * np.finfo(float).eps
    rtol = 4 * np.finfo(float).eps

    def _run_one_test(self, tc, method, sig_args_keys=None,
                      sig_kwargs_keys=None, **kwargs):
        method_args = []
        for k in sig_args_keys or []:
            if k not in tc:
                # If a,b not present use x0, x1. Similarly for f and func
                k = {'a': 'x0', 'b': 'x1', 'func': 'f'}.get(k, k)
            method_args.append(tc[k])

        method_kwargs = dict(**kwargs)
        method_kwargs.update({'full_output': True, 'disp': False})
        for k in sig_kwargs_keys or []:
            method_kwargs[k] = tc[k]

        root = tc.get('root')
        func_args = tc.get('args', ())

        try:
            r, rr = method(*method_args, args=func_args, **method_kwargs)
            return root, rr, tc
        except Exception:
            return root, zeros.RootResults(nan, -1, -1, zeros._EVALUEERR), tc

    def run_tests(self, tests, method, name, known_fail=None, **kwargs):
        r"""Run test-cases using the specified method and the supplied signature.

        Extract the arguments for the method call from the test case
        dictionary using the supplied keys for the method's signature."""
        # The methods have one of two base signatures:
        # (f, a, b, **kwargs)  # newton
        # (func, x0, **kwargs)  # bisect/brentq/...
        sig = _getfullargspec(method)  # FullArgSpec with args, varargs, varkw, defaults, ...
        assert_(not sig.kwonlyargs)
        nDefaults = len(sig.defaults)
        nRequired = len(sig.args) - nDefaults
        sig_args_keys = sig.args[:nRequired]
        sig_kwargs_keys = []
        if name in ['secant', 'newton', 'halley']:
            if name in ['newton', 'halley']:
                sig_kwargs_keys.append('fprime')
                if name in ['halley']:
                    sig_kwargs_keys.append('fprime2')
            kwargs['tol'] = self.xtol
        else:
            kwargs['xtol'] = self.xtol
            kwargs['rtol'] = self.rtol

        results = [list(self._run_one_test(
            tc, method, sig_args_keys=sig_args_keys,
            sig_kwargs_keys=sig_kwargs_keys, **kwargs)) for tc in tests]
        # results= [[true root, full output, tc], ...]

        known_fail = known_fail or []
        notcvgd = [elt for elt in results if not elt[1].converged]
        notcvgd = [elt for elt in notcvgd if elt[-1]['ID'] not in known_fail]
        notcvged_IDS = [elt[-1]['ID'] for elt in notcvgd]
        assert_equal([len(notcvged_IDS), notcvged_IDS], [0, []])

        # The usable xtol and rtol depend on the test
        tols = {'xtol': self.xtol, 'rtol': self.rtol}
        tols.update(**kwargs)
        rtol = tols['rtol']
        atol = tols.get('tol', tols['xtol'])

        cvgd = [elt for elt in results if elt[1].converged]
        approx = [elt[1].root for elt in cvgd]
        correct = [elt[0] for elt in cvgd]
        # See if the root matches the reference value
        notclose = [[a] + elt for a, c, elt in zip(approx, correct, cvgd) if
                    not isclose(a, c, rtol=rtol, atol=atol)
                    and elt[-1]['ID'] not in known_fail]
        # If not, evaluate the function and see if is 0 at the purported root
        fvs = [tc['f'](aroot, *tc.get('args', tuple()))
               for aroot, c, fullout, tc in notclose]
        notclose = [[fv] + elt for fv, elt in zip(fvs, notclose) if fv != 0]
        assert_equal([notclose, len(notclose)], [[], 0])

    def run_collection(self, collection, method, name, smoothness=None,
                       known_fail=None, **kwargs):
        r"""Run a collection of tests using the specified method.

        The name is used to determine some optional arguments."""
        tests = get_tests(collection, smoothness=smoothness)
        self.run_tests(tests, method, name, known_fail=known_fail, **kwargs)


class TestBracketMethods(TestScalarRootFinders):
    @pytest.mark.parametrize('method', bracket_methods)
    @pytest.mark.parametrize('function', tstutils_functions)
    def test_basic_root_scalar(self, method, function):
        # Tests bracketing root finders called via `root_scalar` on a small
        # set of simple problems, each of which has a root at `x=1`. Checks for
        # converged status and that the root was found.
        a, b = .5, sqrt(3)

        r = root_scalar(function, method=method.__name__, bracket=[a, b], x0=a,
                        xtol=self.xtol, rtol=self.rtol)
        assert r.converged
        assert_allclose(r.root, 1.0, atol=self.xtol, rtol=self.rtol)

    @pytest.mark.parametrize('method', bracket_methods)
    @pytest.mark.parametrize('function', tstutils_functions)
    def test_basic_individual(self, method, function):
        # Tests individual bracketing root finders on a small set of simple
        # problems, each of which has a root at `x=1`. Checks for converged
        # status and that the root was found.
        a, b = .5, sqrt(3)
        root, r = method(function, a, b, xtol=self.xtol, rtol=self.rtol,
                         full_output=True)

        assert r.converged
        assert_allclose(root, 1.0, atol=self.xtol, rtol=self.rtol)

    @pytest.mark.parametrize('method', bracket_methods)
    def test_aps_collection(self, method):
        self.run_collection('aps', method, method.__name__, smoothness=1)

    @pytest.mark.parametrize('method', [zeros.bisect, zeros.ridder,
                                        zeros.toms748])
    def test_chandrupatla_collection(self, method):
        known_fail = {'fun7.4'} if method == zeros.ridder else {}
        self.run_collection('chandrupatla', method, method.__name__,
                            known_fail=known_fail)

    @pytest.mark.parametrize('method', bracket_methods)
    def test_lru_cached_individual(self, method):
        # check that https://github.com/scipy/scipy/issues/10846 is fixed
        # (`root_scalar` failed when passed a function that was `@lru_cache`d)
        a, b = -1, 1
        root, r = method(f_lrucached, a, b, full_output=True)
        assert r.converged
        assert_allclose(root, 0)


class TestChandrupatla(TestScalarRootFinders):

    def f(self, q, p):
        return stats.norm.cdf(q) - p

    @pytest.mark.parametrize('p', [0.6, np.linspace(-0.05, 1.05, 10)])
    def test_basic(self, p):
        # Invert distribution CDF and compare against distrtibution `ppf`
        res = zeros._chandrupatla(self.f, -5, 5, args=(p,))
        ref = stats.norm().ppf(p)
        np.testing.assert_allclose(res.x, ref)
        assert res.x.shape == ref.shape

    @pytest.mark.parametrize('shape', [tuple(), (12,), (3, 4), (3, 2, 2)])
    def test_vectorization(self, shape):
        # Test for correct functionality, output shapes, and dtypes for various
        # input shapes.
        p = np.linspace(-0.05, 1.05, 12).reshape(shape) if shape else 0.6
        args = (p,)

        @np.vectorize
        def chandrupatla_single(p):
            return zeros._chandrupatla(self.f, -5, 5, args=(p,))

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
        def f(xs, js):
            funcs = [lambda x: x - 2.5,
                     lambda x: x - 10,
                     lambda x: (x - 0.1)**3,
                     lambda x: np.nan]

            return [funcs[j](x) for x, j in zip(xs, js)]

        args = (np.arange(4, dtype=np.int64),)
        res = zeros._chandrupatla(f, [0]*4, [np.pi]*4, args=args, maxiter=2)

        ref_flags = np.array([zeros._ECONVERGED, zeros._ESIGNERR,
                              zeros._ECONVERR, zeros._EVALUEERR])
        assert_equal(res.status, ref_flags)

    def test_convergence(self):
        # Test that the convergence tolerances behave as expected
        rng = np.random.default_rng(2585255913088665241)
        p = rng.random(size=3)
        bracket = (-5, 5)
        args = (p,)
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
        bracket = (-5, 5)
        maxiter = 5

        def f(q, p):
            res = stats.norm().cdf(q) - p
            f.x = q
            f.fun = res
            return res
        f.x = None
        f.fun = None

        res = zeros._chandrupatla(f, *bracket, args=(p,),
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
            else:
                # Ensure that attributes are updating each iteration
                assert f.x[0] in {res.xl, res.xr}
                assert f.fun[0] in {res.fl, res.fr}
            assert res.status == zeros._EINPROGRESS
            if callback.iter == maxiter:
                raise StopIteration
        callback.iter = -1  # callback called once before first iteration
        callback.res = None

        res2 = zeros._chandrupatla(f, *bracket, args=(p,),
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

        message = 'Abscissae and function output must be real numbers.'
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
        assert res.success
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


class TestNewton(TestScalarRootFinders):
    def test_newton_collections(self):
        known_fail = ['aps.13.00']
        known_fail += ['aps.12.05', 'aps.12.17']  # fails under Windows Py27
        for collection in ['aps', 'complex']:
            self.run_collection(collection, zeros.newton, 'newton',
                                smoothness=2, known_fail=known_fail)

    def test_halley_collections(self):
        known_fail = ['aps.12.06', 'aps.12.07', 'aps.12.08', 'aps.12.09',
                      'aps.12.10', 'aps.12.11', 'aps.12.12', 'aps.12.13',
                      'aps.12.14', 'aps.12.15', 'aps.12.16', 'aps.12.17',
                      'aps.12.18', 'aps.13.00']
        for collection in ['aps', 'complex']:
            self.run_collection(collection, zeros.newton, 'halley',
                                smoothness=2, known_fail=known_fail)

    def test_newton(self):
        for f, f_1, f_2 in [(f1, f1_1, f1_2), (f2, f2_1, f2_2)]:
            x = zeros.newton(f, 3, tol=1e-6)
            assert_allclose(f(x), 0, atol=1e-6)
            x = zeros.newton(f, 3, x1=5, tol=1e-6)  # secant, x0 and x1
            assert_allclose(f(x), 0, atol=1e-6)
            x = zeros.newton(f, 3, fprime=f_1, tol=1e-6)   # newton
            assert_allclose(f(x), 0, atol=1e-6)
            x = zeros.newton(f, 3, fprime=f_1, fprime2=f_2, tol=1e-6)  # halley
            assert_allclose(f(x), 0, atol=1e-6)

    def test_newton_by_name(self):
        r"""Invoke newton through root_scalar()"""
        for f, f_1, f_2 in [(f1, f1_1, f1_2), (f2, f2_1, f2_2)]:
            r = root_scalar(f, method='newton', x0=3, fprime=f_1, xtol=1e-6)
            assert_allclose(f(r.root), 0, atol=1e-6)
        for f, f_1, f_2 in [(f1, f1_1, f1_2), (f2, f2_1, f2_2)]:
            r = root_scalar(f, method='newton', x0=3, xtol=1e-6)  # without f'
            assert_allclose(f(r.root), 0, atol=1e-6)

    def test_secant_by_name(self):
        r"""Invoke secant through root_scalar()"""
        for f, f_1, f_2 in [(f1, f1_1, f1_2), (f2, f2_1, f2_2)]:
            r = root_scalar(f, method='secant', x0=3, x1=2, xtol=1e-6)
            assert_allclose(f(r.root), 0, atol=1e-6)
            r = root_scalar(f, method='secant', x0=3, x1=5, xtol=1e-6)
            assert_allclose(f(r.root), 0, atol=1e-6)
        for f, f_1, f_2 in [(f1, f1_1, f1_2), (f2, f2_1, f2_2)]:
            r = root_scalar(f, method='secant', x0=3, xtol=1e-6)  # without x1
            assert_allclose(f(r.root), 0, atol=1e-6)

    def test_halley_by_name(self):
        r"""Invoke halley through root_scalar()"""
        for f, f_1, f_2 in [(f1, f1_1, f1_2), (f2, f2_1, f2_2)]:
            r = root_scalar(f, method='halley', x0=3,
                            fprime=f_1, fprime2=f_2, xtol=1e-6)
            assert_allclose(f(r.root), 0, atol=1e-6)

    def test_root_scalar_fail(self):
        message = 'fprime2 must be specified for halley'
        with pytest.raises(ValueError, match=message):
            root_scalar(f1, method='halley', fprime=f1_1, x0=3, xtol=1e-6)  # no fprime2
        message = 'fprime must be specified for halley'
        with pytest.raises(ValueError, match=message):
            root_scalar(f1, method='halley', fprime2=f1_2, x0=3, xtol=1e-6)  # no fprime

    def test_array_newton(self):
        """test newton with array"""

        def f1(x, *a):
            b = a[0] + x * a[3]
            return a[1] - a[2] * (np.exp(b / a[5]) - 1.0) - b / a[4] - x

        def f1_1(x, *a):
            b = a[3] / a[5]
            return -a[2] * np.exp(a[0] / a[5] + x * b) * b - a[3] / a[4] - 1

        def f1_2(x, *a):
            b = a[3] / a[5]
            return -a[2] * np.exp(a[0] / a[5] + x * b) * b**2

        a0 = np.array([
            5.32725221, 5.48673747, 5.49539973,
            5.36387202, 4.80237316, 1.43764452,
            5.23063958, 5.46094772, 5.50512718,
            5.42046290
        ])
        a1 = (np.sin(range(10)) + 1.0) * 7.0
        args = (a0, a1, 1e-09, 0.004, 10, 0.27456)
        x0 = [7.0] * 10
        x = zeros.newton(f1, x0, f1_1, args)
        x_expected = (
            6.17264965, 11.7702805, 12.2219954,
            7.11017681, 1.18151293, 0.143707955,
            4.31928228, 10.5419107, 12.7552490,
            8.91225749
        )
        assert_allclose(x, x_expected)
        # test halley's
        x = zeros.newton(f1, x0, f1_1, args, fprime2=f1_2)
        assert_allclose(x, x_expected)
        # test secant
        x = zeros.newton(f1, x0, args=args)
        assert_allclose(x, x_expected)

    def test_array_newton_complex(self):
        def f(x):
            return x + 1+1j

        def fprime(x):
            return 1.0

        t = np.full(4, 1j)
        x = zeros.newton(f, t, fprime=fprime)
        assert_allclose(f(x), 0.)

        # should work even if x0 is not complex
        t = np.ones(4)
        x = zeros.newton(f, t, fprime=fprime)
        assert_allclose(f(x), 0.)

        x = zeros.newton(f, t)
        assert_allclose(f(x), 0.)

    def test_array_secant_active_zero_der(self):
        """test secant doesn't continue to iterate zero derivatives"""
        x = zeros.newton(lambda x, *a: x*x - a[0], x0=[4.123, 5],
                         args=[np.array([17, 25])])
        assert_allclose(x, (4.123105625617661, 5.0))

    def test_array_newton_integers(self):
        # test secant with float
        x = zeros.newton(lambda y, z: z - y ** 2, [4.0] * 2,
                         args=([15.0, 17.0],))
        assert_allclose(x, (3.872983346207417, 4.123105625617661))
        # test integer becomes float
        x = zeros.newton(lambda y, z: z - y ** 2, [4] * 2, args=([15, 17],))
        assert_allclose(x, (3.872983346207417, 4.123105625617661))

    def test_array_newton_zero_der_failures(self):
        # test derivative zero warning
        assert_warns(RuntimeWarning, zeros.newton,
                     lambda y: y**2 - 2, [0., 0.], lambda y: 2 * y)
        # test failures and zero_der
        with pytest.warns(RuntimeWarning):
            results = zeros.newton(lambda y: y**2 - 2, [0., 0.],
                                   lambda y: 2*y, full_output=True)
            assert_allclose(results.root, 0)
            assert results.zero_der.all()
            assert not results.converged.any()

    def test_newton_combined(self):
        def f1(x):
            return x ** 2 - 2 * x - 1
        def f1_1(x):
            return 2 * x - 2
        def f1_2(x):
            return 2.0 + 0 * x

        def f1_and_p_and_pp(x):
            return x**2 - 2*x-1, 2*x-2, 2.0

        sol0 = root_scalar(f1, method='newton', x0=3, fprime=f1_1)
        sol = root_scalar(f1_and_p_and_pp, method='newton', x0=3, fprime=True)
        assert_allclose(sol0.root, sol.root, atol=1e-8)
        assert_equal(2*sol.function_calls, sol0.function_calls)

        sol0 = root_scalar(f1, method='halley', x0=3, fprime=f1_1, fprime2=f1_2)
        sol = root_scalar(f1_and_p_and_pp, method='halley', x0=3, fprime2=True)
        assert_allclose(sol0.root, sol.root, atol=1e-8)
        assert_equal(3*sol.function_calls, sol0.function_calls)

    def test_newton_full_output(self):
        # Test the full_output capability, both when converging and not.
        # Use simple polynomials, to avoid hitting platform dependencies
        # (e.g., exp & trig) in number of iterations

        x0 = 3
        expected_counts = [(6, 7), (5, 10), (3, 9)]

        for derivs in range(3):
            kwargs = {'tol': 1e-6, 'full_output': True, }
            for k, v in [['fprime', f1_1], ['fprime2', f1_2]][:derivs]:
                kwargs[k] = v

            x, r = zeros.newton(f1, x0, disp=False, **kwargs)
            assert_(r.converged)
            assert_equal(x, r.root)
            assert_equal((r.iterations, r.function_calls), expected_counts[derivs])
            if derivs == 0:
                assert r.function_calls <= r.iterations + 1
            else:
                assert_equal(r.function_calls, (derivs + 1) * r.iterations)

            # Now repeat, allowing one fewer iteration to force convergence failure
            iters = r.iterations - 1
            x, r = zeros.newton(f1, x0, maxiter=iters, disp=False, **kwargs)
            assert_(not r.converged)
            assert_equal(x, r.root)
            assert_equal(r.iterations, iters)

            if derivs == 1:
                # Check that the correct Exception is raised and
                # validate the start of the message.
                with pytest.raises(
                    RuntimeError,
                    match='Failed to converge after %d iterations, value is .*' % (iters)):
                    x, r = zeros.newton(f1, x0, maxiter=iters, disp=True, **kwargs)

    def test_deriv_zero_warning(self):
        def func(x):
            return x ** 2 - 2.0
        def dfunc(x):
            return 2 * x
        assert_warns(RuntimeWarning, zeros.newton, func, 0.0, dfunc, disp=False)
        with pytest.raises(RuntimeError, match='Derivative was zero'):
            zeros.newton(func, 0.0, dfunc)

    def test_newton_does_not_modify_x0(self):
        # https://github.com/scipy/scipy/issues/9964
        x0 = np.array([0.1, 3])
        x0_copy = x0.copy()  # Copy to test for equality.
        newton(np.sin, x0, np.cos)
        assert_array_equal(x0, x0_copy)

    def test_gh17570_defaults(self):
        # Previously, when fprime was not specified, root_scalar would default
        # to secant. When x1 was not specified, secant failed.
        # Check that without fprime, the default is secant if x1 is specified
        # and newton otherwise.
        res_newton_default = root_scalar(f1, method='newton', x0=3, xtol=1e-6)
        res_secant_default = root_scalar(f1, method='secant', x0=3, x1=2,
                                         xtol=1e-6)
        # `newton` uses the secant method when `x1` and `x2` are specified
        res_secant = newton(f1, x0=3, x1=2, tol=1e-6, full_output=True)[1]

        # all three found a root
        assert_allclose(f1(res_newton_default.root), 0, atol=1e-6)
        assert res_newton_default.root.shape == tuple()
        assert_allclose(f1(res_secant_default.root), 0, atol=1e-6)
        assert res_secant_default.root.shape == tuple()
        assert_allclose(f1(res_secant.root), 0, atol=1e-6)
        assert res_secant.root.shape == tuple()

        # Defaults are correct
        assert (res_secant_default.root
                == res_secant.root
                != res_newton_default.iterations)
        assert (res_secant_default.iterations
                == res_secant_default.function_calls - 1  # true for secant
                == res_secant.iterations
                != res_newton_default.iterations
                == res_newton_default.function_calls/2)  # newton 2-point diff


def test_gh_5555():
    root = 0.1

    def f(x):
        return x - root

    methods = [zeros.bisect, zeros.ridder]
    xtol = rtol = TOL
    for method in methods:
        res = method(f, -1e8, 1e7, xtol=xtol, rtol=rtol)
        assert_allclose(root, res, atol=xtol, rtol=rtol,
                        err_msg='method %s' % method.__name__)


def test_gh_5557():
    # Show that without the changes in 5557 brentq and brenth might
    # only achieve a tolerance of 2*(xtol + rtol*|res|).

    # f linearly interpolates (0, -0.1), (0.5, -0.1), and (1,
    # 0.4). The important parts are that |f(0)| < |f(1)| (so that
    # brent takes 0 as the initial guess), |f(0)| < atol (so that
    # brent accepts 0 as the root), and that the exact root of f lies
    # more than atol away from 0 (so that brent doesn't achieve the
    # desired tolerance).
    def f(x):
        if x < 0.5:
            return -0.1
        else:
            return x - 0.6

    atol = 0.51
    rtol = 4 * _FLOAT_EPS
    methods = [zeros.brentq, zeros.brenth]
    for method in methods:
        res = method(f, 0, 1, xtol=atol, rtol=rtol)
        assert_allclose(0.6, res, atol=atol, rtol=rtol)


def test_brent_underflow_in_root_bracketing():
    # Tetsing if an interval [a,b] brackets a zero of a function
    # by checking f(a)*f(b) < 0 is not reliable when the product
    # underflows/overflows. (reported in issue# 13737)

    underflow_scenario = (-450.0, -350.0, -400.0)
    overflow_scenario = (350.0, 450.0, 400.0)

    for a, b, root in [underflow_scenario, overflow_scenario]:
        c = np.exp(root)
        for method in [zeros.brenth, zeros.brentq]:
            res = method(lambda x: np.exp(x)-c, a, b)
            assert_allclose(root, res)


class TestRootResults:
    r = zeros.RootResults(root=1.0, iterations=44, function_calls=46, flag=0)

    def test_repr(self):
        expected_repr = ("      converged: True\n           flag: converged"
                         "\n function_calls: 46\n     iterations: 44\n"
                         "           root: 1.0")
        assert_equal(repr(self.r), expected_repr)

    def test_type(self):
        assert isinstance(self.r, OptimizeResult)


def test_complex_halley():
    """Test Halley's works with complex roots"""
    def f(x, *a):
        return a[0] * x**2 + a[1] * x + a[2]

    def f_1(x, *a):
        return 2 * a[0] * x + a[1]

    def f_2(x, *a):
        retval = 2 * a[0]
        try:
            size = len(x)
        except TypeError:
            return retval
        else:
            return [retval] * size

    z = complex(1.0, 2.0)
    coeffs = (2.0, 3.0, 4.0)
    y = zeros.newton(f, z, args=coeffs, fprime=f_1, fprime2=f_2, tol=1e-6)
    # (-0.75000000000000078+1.1989578808281789j)
    assert_allclose(f(y, *coeffs), 0, atol=1e-6)
    z = [z] * 10
    coeffs = (2.0, 3.0, 4.0)
    y = zeros.newton(f, z, args=coeffs, fprime=f_1, fprime2=f_2, tol=1e-6)
    assert_allclose(f(y, *coeffs), 0, atol=1e-6)


def test_zero_der_nz_dp():
    """Test secant method with a non-zero dp, but an infinite newton step"""
    # pick a symmetrical functions and choose a point on the side that with dx
    # makes a secant that is a flat line with zero slope, EG: f = (x - 100)**2,
    # which has a root at x = 100 and is symmetrical around the line x = 100
    # we have to pick a really big number so that it is consistently true
    # now find a point on each side so that the secant has a zero slope
    dx = np.finfo(float).eps ** 0.33
    # 100 - p0 = p1 - 100 = p0 * (1 + dx) + dx - 100
    # -> 200 = p0 * (2 + dx) + dx
    p0 = (200.0 - dx) / (2.0 + dx)
    with suppress_warnings() as sup:
        sup.filter(RuntimeWarning, "RMS of")
        x = zeros.newton(lambda y: (y - 100.0)**2, x0=[p0] * 10)
    assert_allclose(x, [100] * 10)
    # test scalar cases too
    p0 = (2.0 - 1e-4) / (2.0 + 1e-4)
    with suppress_warnings() as sup:
        sup.filter(RuntimeWarning, "Tolerance of")
        x = zeros.newton(lambda y: (y - 1.0) ** 2, x0=p0, disp=False)
    assert_allclose(x, 1)
    with pytest.raises(RuntimeError, match='Tolerance of'):
        x = zeros.newton(lambda y: (y - 1.0) ** 2, x0=p0, disp=True)
    p0 = (-2.0 + 1e-4) / (2.0 + 1e-4)
    with suppress_warnings() as sup:
        sup.filter(RuntimeWarning, "Tolerance of")
        x = zeros.newton(lambda y: (y + 1.0) ** 2, x0=p0, disp=False)
    assert_allclose(x, -1)
    with pytest.raises(RuntimeError, match='Tolerance of'):
        x = zeros.newton(lambda y: (y + 1.0) ** 2, x0=p0, disp=True)


def test_array_newton_failures():
    """Test that array newton fails as expected"""
    # p = 0.68  # [MPa]
    # dp = -0.068 * 1e6  # [Pa]
    # T = 323  # [K]
    diameter = 0.10  # [m]
    # L = 100  # [m]
    roughness = 0.00015  # [m]
    rho = 988.1  # [kg/m**3]
    mu = 5.4790e-04  # [Pa*s]
    u = 2.488  # [m/s]
    reynolds_number = rho * u * diameter / mu  # Reynolds number

    def colebrook_eqn(darcy_friction, re, dia):
        return (1 / np.sqrt(darcy_friction) +
                2 * np.log10(roughness / 3.7 / dia +
                             2.51 / re / np.sqrt(darcy_friction)))

    # only some failures
    with pytest.warns(RuntimeWarning):
        result = zeros.newton(
            colebrook_eqn, x0=[0.01, 0.2, 0.02223, 0.3], maxiter=2,
            args=[reynolds_number, diameter], full_output=True
        )
        assert not result.converged.all()
    # they all fail
    with pytest.raises(RuntimeError):
        result = zeros.newton(
            colebrook_eqn, x0=[0.01] * 2, maxiter=2,
            args=[reynolds_number, diameter], full_output=True
        )


# this test should **not** raise a RuntimeWarning
def test_gh8904_zeroder_at_root_fails():
    """Test that Newton or Halley don't warn if zero derivative at root"""

    # a function that has a zero derivative at it's root
    def f_zeroder_root(x):
        return x**3 - x**2

    # should work with secant
    r = zeros.newton(f_zeroder_root, x0=0)
    assert_allclose(r, 0, atol=zeros._xtol, rtol=zeros._rtol)
    # test again with array
    r = zeros.newton(f_zeroder_root, x0=[0]*10)
    assert_allclose(r, 0, atol=zeros._xtol, rtol=zeros._rtol)

    # 1st derivative
    def fder(x):
        return 3 * x**2 - 2 * x

    # 2nd derivative
    def fder2(x):
        return 6*x - 2

    # should work with newton and halley
    r = zeros.newton(f_zeroder_root, x0=0, fprime=fder)
    assert_allclose(r, 0, atol=zeros._xtol, rtol=zeros._rtol)
    r = zeros.newton(f_zeroder_root, x0=0, fprime=fder,
                     fprime2=fder2)
    assert_allclose(r, 0, atol=zeros._xtol, rtol=zeros._rtol)
    # test again with array
    r = zeros.newton(f_zeroder_root, x0=[0]*10, fprime=fder)
    assert_allclose(r, 0, atol=zeros._xtol, rtol=zeros._rtol)
    r = zeros.newton(f_zeroder_root, x0=[0]*10, fprime=fder,
                     fprime2=fder2)
    assert_allclose(r, 0, atol=zeros._xtol, rtol=zeros._rtol)

    # also test that if a root is found we do not raise RuntimeWarning even if
    # the derivative is zero, EG: at x = 0.5, then fval = -0.125 and
    # fder = -0.25 so the next guess is 0.5 - (-0.125/-0.5) = 0 which is the
    # root, but if the solver continued with that guess, then it will calculate
    # a zero derivative, so it should return the root w/o RuntimeWarning
    r = zeros.newton(f_zeroder_root, x0=0.5, fprime=fder)
    assert_allclose(r, 0, atol=zeros._xtol, rtol=zeros._rtol)
    # test again with array
    r = zeros.newton(f_zeroder_root, x0=[0.5]*10, fprime=fder)
    assert_allclose(r, 0, atol=zeros._xtol, rtol=zeros._rtol)
    # doesn't apply to halley


def test_gh_8881():
    r"""Test that Halley's method realizes that the 2nd order adjustment
    is too big and drops off to the 1st order adjustment."""
    n = 9

    def f(x):
        return power(x, 1.0/n) - power(n, 1.0/n)

    def fp(x):
        return power(x, (1.0-n)/n)/n

    def fpp(x):
        return power(x, (1.0-2*n)/n) * (1.0/n) * (1.0-n)/n

    x0 = 0.1
    # The root is at x=9.
    # The function has positive slope, x0 < root.
    # Newton succeeds in 8 iterations
    rt, r = newton(f, x0, fprime=fp, full_output=True)
    assert r.converged
    # Before the Issue 8881/PR 8882, halley would send x in the wrong direction.
    # Check that it now succeeds.
    rt, r = newton(f, x0, fprime=fp, fprime2=fpp, full_output=True)
    assert r.converged


def test_gh_9608_preserve_array_shape():
    """
    Test that shape is preserved for array inputs even if fprime or fprime2 is
    scalar
    """
    def f(x):
        return x**2

    def fp(x):
        return 2 * x

    def fpp(x):
        return 2

    x0 = np.array([-2], dtype=np.float32)
    rt, r = newton(f, x0, fprime=fp, fprime2=fpp, full_output=True)
    assert r.converged

    x0_array = np.array([-2, -3], dtype=np.float32)
    # This next invocation should fail
    with pytest.raises(IndexError):
        result = zeros.newton(
            f, x0_array, fprime=fp, fprime2=fpp, full_output=True
        )

    def fpp_array(x):
        return np.full(np.shape(x), 2, dtype=np.float32)

    result = zeros.newton(
        f, x0_array, fprime=fp, fprime2=fpp_array, full_output=True
    )
    assert result.converged.all()


@pytest.mark.parametrize(
    "maximum_iterations,flag_expected",
    [(10, zeros.CONVERR), (100, zeros.CONVERGED)])
def test_gh9254_flag_if_maxiter_exceeded(maximum_iterations, flag_expected):
    """
    Test that if the maximum iterations is exceeded that the flag is not
    converged.
    """
    result = zeros.brentq(
        lambda x: ((1.2*x - 2.3)*x + 3.4)*x - 4.5,
        -30, 30, (), 1e-6, 1e-6, maximum_iterations,
        full_output=True, disp=False)
    assert result[1].flag == flag_expected
    if flag_expected == zeros.CONVERR:
        # didn't converge because exceeded maximum iterations
        assert result[1].iterations == maximum_iterations
    elif flag_expected == zeros.CONVERGED:
        # converged before maximum iterations
        assert result[1].iterations < maximum_iterations


def test_gh9551_raise_error_if_disp_true():
    """Test that if disp is true then zero derivative raises RuntimeError"""

    def f(x):
        return x*x + 1

    def f_p(x):
        return 2*x

    assert_warns(RuntimeWarning, zeros.newton, f, 1.0, f_p, disp=False)
    with pytest.raises(
            RuntimeError,
            match=r'^Derivative was zero\. Failed to converge after \d+ iterations, value is [+-]?\d*\.\d+\.$'):
        zeros.newton(f, 1.0, f_p)
    root = zeros.newton(f, complex(10.0, 10.0), f_p)
    assert_allclose(root, complex(0.0, 1.0))


@pytest.mark.parametrize('solver_name',
                         ['brentq', 'brenth', 'bisect', 'ridder', 'toms748'])
def test_gh3089_8394(solver_name):
    # gh-3089 and gh-8394 reported that bracketing solvers returned incorrect
    # results when they encountered NaNs. Check that this is resolved.
    def f(x):
        return np.nan

    solver = getattr(zeros, solver_name)
    with pytest.raises(ValueError, match="The function value at x..."):
        solver(f, 0, 1)


@pytest.mark.parametrize('method',
                         ['brentq', 'brenth', 'bisect', 'ridder', 'toms748'])
def test_gh18171(method):
    # gh-3089 and gh-8394 reported that bracketing solvers returned incorrect
    # results when they encountered NaNs. Check that `root_scalar` returns
    # normally but indicates that convergence was unsuccessful. See gh-18171.
    def f(x):
        f._count += 1
        return np.nan
    f._count = 0

    res = root_scalar(f, bracket=(0, 1), method=method)
    assert res.converged is False
    assert res.flag.startswith("The function value at x")
    assert res.function_calls == f._count
    assert str(res.root) in res.flag


@pytest.mark.parametrize('solver_name',
                         ['brentq', 'brenth', 'bisect', 'ridder', 'toms748'])
@pytest.mark.parametrize('rs_interface', [True, False])
def test_function_calls(solver_name, rs_interface):
    # There do not appear to be checks that the bracketing solvers report the
    # correct number of function evaluations. Check that this is the case.
    solver = ((lambda f, a, b, **kwargs: root_scalar(f, bracket=(a, b)))
              if rs_interface else getattr(zeros, solver_name))

    def f(x):
        f.calls += 1
        return x**2 - 1
    f.calls = 0

    res = solver(f, 0, 10, full_output=True)

    if rs_interface:
        assert res.function_calls == f.calls
    else:
        assert res[1].function_calls == f.calls


def test_gh_14486_converged_false():
    """Test that zero slope with secant method results in a converged=False"""
    def lhs(x):
        return x * np.exp(-x*x) - 0.07

    with pytest.warns(RuntimeWarning, match='Tolerance of'):
        res = root_scalar(lhs, method='secant', x0=-0.15, x1=1.0)
    assert not res.converged
    assert res.flag == 'convergence error'

    with pytest.warns(RuntimeWarning, match='Tolerance of'):
        res = newton(lhs, x0=-0.15, x1=1.0, disp=False, full_output=True)[1]
    assert not res.converged
    assert res.flag == 'convergence error'


@pytest.mark.parametrize('solver_name',
                         ['brentq', 'brenth', 'bisect', 'ridder', 'toms748'])
@pytest.mark.parametrize('rs_interface', [True, False])
def test_gh5584(solver_name, rs_interface):
    # gh-5584 reported that an underflow can cause sign checks in the algorithm
    # to fail. Check that this is resolved.
    solver = ((lambda f, a, b, **kwargs: root_scalar(f, bracket=(a, b)))
              if rs_interface else getattr(zeros, solver_name))

    def f(x):
        return 1e-200*x

    # Report failure when signs are the same
    with pytest.raises(ValueError, match='...must have different signs'):
        solver(f, -0.5, -0.4, full_output=True)

    # Solve successfully when signs are different
    res = solver(f, -0.5, 0.4, full_output=True)
    res = res if rs_interface else res[1]
    assert res.converged
    assert_allclose(res.root, 0, atol=1e-8)

    # Solve successfully when one side is negative zero
    res = solver(f, -0.5, float('-0.0'), full_output=True)
    res = res if rs_interface else res[1]
    assert res.converged
    assert_allclose(res.root, 0, atol=1e-8)


def test_gh13407():
    # gh-13407 reported that the message produced by `scipy.optimize.toms748`
    # when `rtol < eps` is incorrect, and also that toms748 is unusual in
    # accepting `rtol` as low as eps while other solvers raise at 4*eps. Check
    # that the error message has been corrected and that `rtol=eps` can produce
    # a lower function value than `rtol=4*eps`.
    def f(x):
        return x**3 - 2*x - 5

    xtol = 1e-300
    eps = np.finfo(float).eps
    x1 = zeros.toms748(f, 1e-10, 1e10, xtol=xtol, rtol=1*eps)
    f1 = f(x1)
    x4 = zeros.toms748(f, 1e-10, 1e10, xtol=xtol, rtol=4*eps)
    f4 = f(x4)
    assert f1 < f4

    # using old-style syntax to get exactly the same message
    message = fr"rtol too small \({eps/2:g} < {eps:g}\)"
    with pytest.raises(ValueError, match=message):
        zeros.toms748(f, 1e-10, 1e10, xtol=xtol, rtol=eps/2)


def test_newton_complex_gh10103():
    # gh-10103 reported a problem when `newton` is pass a Python complex x0,
    # no `fprime` (secant method), and no `x1` (`x1` must be constructed).
    # Check that this is resolved.
    def f(z):
        return z - 1
    res = newton(f, 1+1j)
    assert_allclose(res, 1, atol=1e-12)

    res = root_scalar(f, x0=1+1j, x1=2+1.5j, method='secant')
    assert_allclose(res.root, 1, atol=1e-12)


@pytest.mark.parametrize('method', all_methods)
def test_maxiter_int_check_gh10236(method):
    # gh-10236 reported that the error message when `maxiter` is not an integer
    # was difficult to interpret. Check that this was resolved (by gh-10907).
    message = "'float' object cannot be interpreted as an integer"
    with pytest.raises(TypeError, match=message):
        method(f1, 0.0, 1.0, maxiter=72.45)
