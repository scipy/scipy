import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from scipy.optimize import biteopt, rosen, Bounds, OptimizeResult


class TestSolver:
    def test_rosenbrock(self):
        bounds = Bounds(lb=[-5.0, -5.0], ub=[5.0, 5.0])
        res = biteopt(rosen, bounds, rng=1234)

        assert_allclose(res.x, [1.0, 1.0], rtol=1e-3, atol=1e-3)
        assert_allclose(res.fun, 0.0, atol=1e-6)

    def test_styblinski_tang(self):
        def styblinski_tang(pos):
            x, y = pos
            return 0.5 * (x**4 - 16 * x**2 + 5 * x
                          + y**4 - 16 * y**2 + 5 * y)

        bounds = [(-4.0, 4.0), (-4.0, 4.0)]
        res = biteopt(styblinski_tang, bounds, rng=1234)

        assert_allclose(res.x, [-2.903534, -2.903534], rtol=1e-3, atol=1e-3)
        assert_allclose(res.fun, -39.16617 * 2, rtol=1e-4)

    def test_args_are_passed(self):
        # Extra fixed parameters supplied via ``args`` must reach the
        # objective. Here they shift the location of the minimum.
        def shifted(x, center):
            return np.sum((x - center) ** 2)

        center = np.array([1.5, -2.0])
        bounds = [(-5.0, 5.0), (-5.0, 5.0)]
        res = biteopt(shifted, bounds, args=(center,), rng=0)

        assert_allclose(res.x, center, rtol=1e-3, atol=1e-3)
        assert_allclose(res.fun, 0.0, atol=1e-6)

    def test_single_precision_objective(self):
        # An objective returning a single-precision (float32) value must be
        # handled correctly: the wrapper promotes it to a Python float.
        def rosen_f32(x):
            return np.float32(rosen(x))

        bounds = [(-5.0, 5.0), (-5.0, 5.0)]
        res = biteopt(rosen_f32, bounds, rng=1234)

        assert_allclose(res.x, [1.0, 1.0], rtol=1e-2, atol=1e-2)
        assert np.isfinite(res.fun)

    def test_nfev_counts_objective_calls_correctly(self):
        bounds = [(-5.0, 5.0), (-5.0, 5.0)]

        class CountingObjective:
            def __init__(self):
                self.nfev = 0

            def __call__(self, x):
                self.nfev += 1
                return rosen(x)

        objective = CountingObjective()
        res = biteopt(objective, bounds, maxfun=200, rng=0)
        assert res.nfev == objective.nfev

    def test_biteopt_stays_within_bounds(self):
        lower_bound = 2
        upper_bound = 3
        def fun(x):
            assert (x >= lower_bound).all()
            assert (x <= upper_bound).all()
            return np.square(x).sum()
        bounds = [(lower_bound, upper_bound)] * 3
        res = biteopt(fun, bounds, rng=0)
        assert np.all(res.x >= lower_bound)
        assert np.all(res.x <= upper_bound)
        assert_allclose(res.x, lower_bound, rtol=1e-3, atol=1e-3)

    @pytest.mark.parametrize("depth", [1, 2, 4, 9])
    def test_depth(self, depth):
        bounds = [(-5.0, 5.0), (-5.0, 5.0)]
        res = biteopt(rosen, bounds, depth=depth, rng=0, maxfun=10000)
        assert_allclose(res.x, [1.0, 1.0], rtol=1e-3, atol=1e-3)

    @pytest.mark.parametrize("depth", [1, 2, 4, 9])
    def test_maxfun_respected_across_depths(self, depth):
        # biteopt internally scales its iteration count by sqrt(depth);
        # the wrapper must pre-divide so nfev never exceeds the budget.
        maxfun = 500
        bounds = [(-5.0, 5.0), (-5.0, 5.0)]
        res = biteopt(rosen, bounds, maxfun=maxfun, depth=depth, rng=0)

        assert res.nfev <= maxfun

    def test_nan_objective_does_not_crash(self):
        # A NaN return is sanitized to a large sentinel internally, so the
        # optimization must complete without raising.
        res = biteopt(lambda x: np.nan, [(-5.0, 5.0), (-5.0, 5.0)], rng=0)
        assert isinstance(res, OptimizeResult)

    def test_inf_objective_does_not_crash(self):
        # An objective returning +inf everywhere must not crash the optimizer;
        # the run must still complete and return a result.
        res = biteopt(lambda x: np.inf, [(-5.0, 5.0), (-5.0, 5.0)], rng=0)
        assert isinstance(res, OptimizeResult)

    def test_f_min_stops_early(self):
        # With a generous f_min threshold the optimizer should stop well before
        # exhausting its full iteration budget, using far fewer evaluations
        # than an unconstrained run.
        target = 1.0
        bounds = [(-5.0, 5.0), (-5.0, 5.0)]
        res_full = biteopt(rosen, bounds, rng=0)
        res_early = biteopt(rosen, bounds, f_min=target, rng=0)

        assert res_early.fun <= target
        assert res_early.nfev < res_full.nfev
        assert res_early.success

    def test_f_min_unreachable_reports_failure(self):
        # An unreachable f_min (below the global minimum of 0) can never be met,
        # so the run is reported as unsuccessful.
        bounds = [(-5.0, 5.0), (-5.0, 5.0)]
        res = biteopt(rosen, bounds, f_min=-1.0, rng=0)

        assert res.fun > -1.0
        assert res.success is False

    @pytest.mark.parametrize("exc_type", [ValueError, KeyError, OverflowError])
    def test_objective_exception_propagates(self, exc_type):
        # An exception raised inside the objective must propagate out of the
        # biteopt call unchanged, preserving the exact exception type.
        def objective(x):
            raise exc_type("Error")

        bounds = [(-5.0, 5.0), (-5.0, 5.0)]
        with pytest.raises(exc_type, match="Error"):
            biteopt(objective, bounds, rng=0)


class TestInputValidation:
    def setup_method(self):
        self.bounds = [(-5.0, 5.0), (-5.0, 5.0)]

    @pytest.mark.parametrize("not_callable", [42, "func", None, [1, 2, 3]])
    def test_func_not_callable(self, not_callable):
        with pytest.raises(TypeError, match="func must be callable"):
            biteopt(not_callable, self.bounds)

    def test_bounds_invalid_type(self):
        with pytest.raises(ValueError, match="bounds must be a sequence"):
            biteopt(rosen, 5)

    def test_bounds_inconsistent(self):
        with pytest.raises(ValueError, match="Bounds are not consistent"):
            biteopt(rosen, [(5.0, -5.0), (-5.0, 5.0)])

    def test_bounds_inf(self):
        with pytest.raises(ValueError, match="Bounds must not be inf"):
            biteopt(rosen, [(-np.inf, 5.0), (-5.0, 5.0)])

    def test_bounds_empty(self):
        with pytest.raises(ValueError, match="at least one"):
            biteopt(rosen, [])

    @pytest.mark.parametrize("maxfun, exc_type, msg", [
        (0, ValueError, "must be an integer not less than 1"),
        (-1, ValueError, "must be an integer not less than 1"),
        (1.5, TypeError, "must be an integer"),
    ])
    def test_invalid_maxfun(self, maxfun, exc_type, msg):
        with pytest.raises(exc_type, match=msg):
            biteopt(rosen, self.bounds, maxfun=maxfun)

    @pytest.mark.parametrize("kwargs, exc_type, msg", [
        ({"depth": 0}, ValueError, "must be an integer not less than 1"),
        ({"depth": -1}, ValueError, "must be an integer not less than 1"),
        ({"depth": 37}, ValueError, r"depth must be an integer in \[1, 36\]"),
        ({"depth": 1.5}, TypeError, "must be an integer"),
    ])
    def test_invalid_depth(self, kwargs, exc_type, msg):
        with pytest.raises(exc_type, match=msg):
            biteopt(rosen, self.bounds, **kwargs)

    def test_f_min_not_float(self):
        with pytest.raises(ValueError, match="float"):
            biteopt(rosen, self.bounds, f_min="not a float")

    def test_func_returns_non_scalar(self):
        # The objective function must return a scalar; if it returns an array
        # or other non-scalar that cannot be converted to a single value, 
        # it must raise a ValueError.
        def non_scalar(x):
            return np.array([1.0, 2.0])

        with pytest.raises(ValueError, match="must return a scalar value"):
            biteopt(non_scalar, self.bounds, rng=0)


class TestRNG:
    def setup_method(self):
        self.bounds = [(-5.0, 5.0), (-5.0, 5.0)]

    def test_reproducible(self):
        # Two generators seeded identically must yield bit-for-bit identical
        # results.
        res1 = biteopt(rosen, self.bounds, rng=np.random.default_rng(42))
        res2 = biteopt(rosen, self.bounds, rng=np.random.default_rng(42))
        assert_array_equal(res1.x, res2.x)
        assert res1.fun == res2.fun
        assert res1.nfev == res2.nfev

    def test_generator_equivalent_to_seed(self):
        # Passing a Generator must be accepted and must be equivalent to
        # passing the integer seed used to construct that Generator.
        res_seed = biteopt(rosen, self.bounds, rng=7)
        res_gen = biteopt(rosen, self.bounds, rng=np.random.default_rng(7))
        assert_array_equal(res_seed.x, res_gen.x)
        assert res_seed.fun == res_gen.fun

    def test_different_seeds_differ(self):
        # Generators seeded differently drive different PRNG trajectories, so
        # the returned minimizers should not be identical. Use a small budget
        # so the runs do not both collapse onto the exact same float optimum.
        res1 = biteopt(rosen, self.bounds, maxfun=50,
                       rng=np.random.default_rng(1))
        res2 = biteopt(rosen, self.bounds, maxfun=50,
                       rng=np.random.default_rng(2))
        assert not np.array_equal(res1.x, res2.x)

    def test_rng_none_runs(self):
        # An unseeded run must still complete successfully.
        res = biteopt(rosen, self.bounds, rng=None)
        assert res.success is True

    def test_accepts_bit_generator(self):
        # default_rng also accepts a bare BitGenerator.
        res = biteopt(rosen, self.bounds, rng=np.random.PCG64(0))
        assert_allclose(res.x, [1.0, 1.0], rtol=1e-3, atol=1e-3)
