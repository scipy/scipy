import math
import scipy
import itertools

import pytest

from scipy._lib._array_api import array_namespace, xp_assert_close, xp_size, np_compat
from scipy.conftest import array_api_compatible

from scipy.integrate import cubature

from scipy.integrate._rules import (
    Rule, FixedRule,
    NestedFixedRule,
    GaussLegendreQuadrature, GaussKronrodQuadrature,
    GenzMalikCubature,
)

pytestmark = [pytest.mark.usefixtures("skip_xp_backends"),]
skip_xp_backends = pytest.mark.skip_xp_backends

# The integrands ``genz_malik_1980_*`` come from the paper:
#   A.C. Genz, A.A. Malik, Remarks on algorithm 006: An adaptive algorithm for
#   numerical integration over an N-dimensional rectangular region, Journal of
#   Computational and Applied Mathematics, Volume 6, Issue 4, 1980, Pages 295-302,
#   ISSN 0377-0427, https://doi.org/10.1016/0771-050X(80)90039-X.


def basic_1d_integrand(x, n, xp):
    x_reshaped = xp.reshape(x, (-1, 1, 1))
    n_reshaped = xp.reshape(n, (1, -1, 1))

    return x_reshaped**n_reshaped


def basic_1d_integrand_exact(n, xp):
    # Exact only for integration over interval [0, 2].
    return xp.reshape(2**(n+1)/(n+1), (-1, 1))


def basic_nd_integrand(x, n, xp):
    return xp.reshape(xp.sum(x, axis=-1), (-1, 1))**xp.reshape(n, (1, -1))


def basic_nd_integrand_exact(n, xp):
    # Exact only for integration over interval [0, 2].
    return (-2**(3+n) + 4**(2+n))/((1+n)*(2+n))


def genz_malik_1980_f_1(x, r, alphas, xp):
    r"""
    .. math:: f_1(\mathbf x) = \cos\left(2\pi r + \sum^n_{i = 1}\alpha_i x_i\right)

    .. code-block:: mathematica

        genzMalik1980f1[x_List, r_, alphas_List] := Cos[2*Pi*r + Total[x*alphas]]
    """

    npoints, ndim = x.shape[0], x.shape[-1]

    alphas_reshaped = alphas[None, ...]
    x_reshaped = xp.reshape(x, (npoints, *([1]*(len(alphas.shape) - 1)), ndim))

    return xp.cos(2*math.pi*r + xp.sum(alphas_reshaped * x_reshaped, axis=-1))


def genz_malik_1980_f_1_exact(a, b, r, alphas, xp):
    ndim = xp_size(a)
    a = xp.reshape(a, (*([1]*(len(alphas.shape) - 1)), ndim))
    b = xp.reshape(b, (*([1]*(len(alphas.shape) - 1)), ndim))

    return (
        (-2)**ndim
        * 1/xp.prod(alphas, axis=-1)
        * xp.cos(2*math.pi*r + xp.sum(alphas * (a+b) * 0.5, axis=-1))
        * xp.prod(xp.sin(alphas * (a-b)/2), axis=-1)
    )


def genz_malik_1980_f_1_random_args(rng, shape, xp):
    r = xp.asarray(rng.random(shape[:-1]))
    alphas = xp.asarray(rng.random(shape))

    difficulty = 9
    normalisation_factors = xp.sum(alphas, axis=-1)[..., None]
    alphas = difficulty * alphas / normalisation_factors

    return (r, alphas)


def genz_malik_1980_f_2(x, alphas, betas, xp):
    r"""
    .. math:: f_2(\mathbf x) = \prod^n_{i = 1} (\alpha_i^2 + (x_i - \beta_i)^2)^{-1}

    .. code-block:: mathematica

        genzMalik1980f2[x_List, alphas_List, betas_List] :=
            1/Times @@ ((alphas^2 + (x - betas)^2))
    """
    npoints, ndim = x.shape[0], x.shape[-1]

    alphas_reshaped = alphas[None, ...]
    betas_reshaped = betas[None, ...]

    x_reshaped = xp.reshape(x, (npoints, *([1]*(len(alphas.shape) - 1)), ndim))

    return 1/xp.prod(alphas_reshaped**2 + (x_reshaped-betas_reshaped)**2, axis=-1)


def genz_malik_1980_f_2_exact(a, b, alphas, betas, xp):
    ndim = xp_size(a)
    a = xp.reshape(a, (*([1]*(len(alphas.shape) - 1)), ndim))
    b = xp.reshape(b, (*([1]*(len(alphas.shape) - 1)), ndim))

    # `xp` is the unwrapped namespace, so `.atan` won't work for `xp = np` and np<2.
    xp_test = array_namespace(a)

    return (
        (-1)**ndim * 1/xp.prod(alphas, axis=-1)
        * xp.prod(
            xp_test.atan((a - betas)/alphas) - xp_test.atan((b - betas)/alphas),
            axis=-1,
        )
    )


def genz_malik_1980_f_2_random_args(rng, shape, xp):
    ndim = shape[-1]
    alphas = xp.asarray(rng.random(shape))
    betas = xp.asarray(rng.random(shape))

    difficulty = 25.0
    products = xp.prod(alphas**xp.asarray(-2.0), axis=-1)
    normalisation_factors = (products**xp.asarray(1 / (2*ndim)))[..., None]
    alphas = alphas * normalisation_factors * math.pow(difficulty, 1 / (2*ndim))

    # Adjust alphas from distribution used in Genz and Malik 1980 since denominator
    # is very small for high dimensions.
    alphas *= 10

    return alphas, betas


def genz_malik_1980_f_3(x, alphas, xp):
    r"""
    .. math:: f_3(\mathbf x) = \exp\left(\sum^n_{i = 1} \alpha_i x_i\right)

    .. code-block:: mathematica

        genzMalik1980f3[x_List, alphas_List] := Exp[Dot[x, alphas]]
    """

    npoints, ndim = x.shape[0], x.shape[-1]

    alphas_reshaped = alphas[None, ...]
    x_reshaped = xp.reshape(x, (npoints, *([1]*(len(alphas.shape) - 1)), ndim))

    return xp.exp(xp.sum(alphas_reshaped * x_reshaped, axis=-1))


def genz_malik_1980_f_3_exact(a, b, alphas, xp):
    ndim = xp_size(a)
    a = xp.reshape(a, (*([1]*(len(alphas.shape) - 1)), ndim))
    b = xp.reshape(b, (*([1]*(len(alphas.shape) - 1)), ndim))

    return (
        (-1)**ndim * 1/xp.prod(alphas, axis=-1)
        * xp.prod(xp.exp(alphas * a) - xp.exp(alphas * b), axis=-1)
    )


def genz_malik_1980_f_3_random_args(rng, shape, xp):
    alphas = xp.asarray(rng.random(shape))
    normalisation_factors = xp.sum(alphas, axis=-1)[..., None]
    difficulty = 12.0
    alphas = difficulty * alphas / normalisation_factors

    return (alphas,)


def genz_malik_1980_f_4(x, alphas, xp):
    r"""
    .. math:: f_4(\mathbf x) = \left(1 + \sum^n_{i = 1} \alpha_i x_i\right)^{-n-1}

    .. code-block:: mathematica
        genzMalik1980f4[x_List, alphas_List] :=
            (1 + Dot[x, alphas])^(-Length[alphas] - 1)
    """

    npoints, ndim = x.shape[0], x.shape[-1]

    alphas_reshaped = alphas[None, ...]
    x_reshaped = xp.reshape(x, (npoints, *([1]*(len(alphas.shape) - 1)), ndim))

    return (1 + xp.sum(alphas_reshaped * x_reshaped, axis=-1))**(-ndim-1)


def genz_malik_1980_f_4_exact(a, b, alphas, xp):
    ndim = xp_size(a)

    def F(x):
        x_reshaped = xp.reshape(x, (*([1]*(len(alphas.shape) - 1)), ndim))

        return (
            (-1)**ndim/xp.prod(alphas, axis=-1)
            / math.factorial(ndim)
            / (1 + xp.sum(alphas * x_reshaped, axis=-1))
        )

    return _eval_indefinite_integral(F, a, b, xp)


def _eval_indefinite_integral(F, a, b, xp):
    """
    Calculates a definite integral from points `a` to `b` by summing up over the corners
    of the corresponding hyperrectangle.
    """

    ndim = xp_size(a)
    points = xp.stack([a, b], axis=0)

    out = 0
    for ind in itertools.product(range(2), repeat=ndim):
        selected_points = xp.asarray([points[i, j] for i, j in zip(ind, range(ndim))])
        out += pow(-1, sum(ind) + ndim) * F(selected_points)

    return out


def genz_malik_1980_f_4_random_args(rng, shape, xp):
    ndim = shape[-1]

    alphas = xp.asarray(rng.random(shape))
    normalisation_factors = xp.sum(alphas, axis=-1)[..., None]
    difficulty = 14.0
    alphas = (difficulty / ndim) * alphas / normalisation_factors

    return (alphas,)


def genz_malik_1980_f_5(x, alphas, betas, xp):
    r"""
    .. math::

        f_5(\mathbf x) = \exp\left(-\sum^n_{i = 1} \alpha^2_i (x_i - \beta_i)^2\right)

    .. code-block:: mathematica

        genzMalik1980f5[x_List, alphas_List, betas_List] :=
            Exp[-Total[alphas^2 * (x - betas)^2]]
    """

    npoints, ndim = x.shape[0], x.shape[-1]

    alphas_reshaped = alphas[None, ...]
    betas_reshaped = betas[None, ...]

    x_reshaped = xp.reshape(x, (npoints, *([1]*(len(alphas.shape) - 1)), ndim))

    return xp.exp(
        -xp.sum(alphas_reshaped**2 * (x_reshaped - betas_reshaped)**2, axis=-1)
    )


def genz_malik_1980_f_5_exact(a, b, alphas, betas, xp):
    ndim = xp_size(a)
    a = xp.reshape(a, (*([1]*(len(alphas.shape) - 1)), ndim))
    b = xp.reshape(b, (*([1]*(len(alphas.shape) - 1)), ndim))

    return (
        (1/2)**ndim
        * 1/xp.prod(alphas, axis=-1)
        * (math.pi**(ndim/2))
        * xp.prod(
            scipy.special.erf(alphas * (betas - a))
            + scipy.special.erf(alphas * (b - betas)),
            axis=-1,
        )
    )


def genz_malik_1980_f_5_random_args(rng, shape, xp):
    alphas = xp.asarray(rng.random(shape))
    betas = xp.asarray(rng.random(shape))

    difficulty = 21.0
    normalisation_factors = xp.sqrt(xp.sum(alphas**xp.asarray(2.0), axis=-1))[..., None]
    alphas = alphas / normalisation_factors * math.sqrt(difficulty)

    return alphas, betas


@array_api_compatible
class TestCubature:
    """
    Tests related to the interface of `cubature`.
    """

    @pytest.mark.parametrize("rule_str", [
        "gauss-kronrod",
        "genz-malik",
        "gk21",
        "gk15",
    ])
    def test_pass_str(self, rule_str, xp):
        n = xp.arange(5, dtype=xp.float64)
        a = xp.asarray([0, 0], dtype=xp.float64)
        b = xp.asarray([2, 2], dtype=xp.float64)

        res = cubature(basic_nd_integrand, a, b, rule=rule_str, args=(n, xp))

        xp_assert_close(
            res.estimate,
            basic_nd_integrand_exact(n, xp),
            rtol=1e-1,
            atol=0,
        )

    @skip_xp_backends(np_only=True,
                      reason='array-likes only supported for NumPy backend')
    def test_pass_array_like_not_array(self, xp):
        n = np_compat.arange(5, dtype=np_compat.float64)
        a = [0]
        b = [2]

        res = cubature(
            basic_1d_integrand,
            a,
            b,
            args=(n, xp)
        )

        xp_assert_close(
            res.estimate,
            basic_1d_integrand_exact(n, xp),
            rtol=1e-1,
            atol=0,
        )

    def test_stops_after_max_subdivisions(self, xp):
        a = xp.asarray([0])
        b = xp.asarray([1])
        rule = BadErrorRule()

        res = cubature(
            basic_1d_integrand,  # Any function would suffice
            a,
            b,
            rule=rule,
            max_subdivisions=10,
            args=(xp.arange(5, dtype=xp.float64), xp),
        )

        assert res.subdivisions == 10
        assert res.status == "not_converged"

    def test_a_and_b_must_be_1d(self, xp):
        a = xp.asarray([[0]], dtype=xp.float64)
        b = xp.asarray([[1]], dtype=xp.float64)

        with pytest.raises(Exception, match="`a` and `b` must be 1D arrays"):
            cubature(basic_1d_integrand, a, b, args=(xp,))

    def test_a_and_b_must_be_nonempty(self, xp):
        a = xp.asarray([])
        b = xp.asarray([])

        with pytest.raises(Exception, match="`a` and `b` must be nonempty"):
            cubature(basic_1d_integrand, a, b, args=(xp,))

    def test_zero_width_limits(self, xp):
        n = xp.arange(5, dtype=xp.float64)

        a = xp.asarray([0], dtype=xp.float64)
        b = xp.asarray([0], dtype=xp.float64)

        res = cubature(
            basic_1d_integrand,
            a,
            b,
            args=(n, xp),
        )

        xp_assert_close(
            res.estimate,
            xp.asarray([[0], [0], [0], [0], [0]], dtype=xp.float64),
            rtol=1e-1,
            atol=0,
        )


@pytest.mark.parametrize("rtol", [1e-4])
@pytest.mark.parametrize("atol", [1e-5])
@pytest.mark.parametrize("rule", [
    "gk15",
    "gk21",
    "genz-malik",
])
@array_api_compatible
class TestCubatureProblems:
    """
    Tests that `cubature` gives the correct answer.
    """

    @pytest.mark.parametrize("problem", [
        # -- f1 --
        (
            # Function to integrate, like `f(x, *args)`
            genz_malik_1980_f_1,

            # Exact solution, like `exact(a, b, *args)`
            genz_malik_1980_f_1_exact,

            # Coordinates of `a`
            [0],

            # Coordinates of `b`
            [10],

            # Arguments to pass to `f` and `exact`
            (
                1/4,
                [5],
            )
        ),
        (
            genz_malik_1980_f_1,
            genz_malik_1980_f_1_exact,
            [0, 0],
            [1, 1],
            (
                1/4,
                [2, 4],
            ),
        ),
        (
            genz_malik_1980_f_1,
            genz_malik_1980_f_1_exact,
            [0, 0],
            [5, 5],
            (
                1/2,
                [2, 4],
            )
        ),
        (
            genz_malik_1980_f_1,
            genz_malik_1980_f_1_exact,
            [0, 0, 0],
            [5, 5, 5],
            (
                1/2,
                [1, 1, 1],
            )
        ),

        # -- f2 --
        (
            genz_malik_1980_f_2,
            genz_malik_1980_f_2_exact,
            [-1],
            [1],
            (
                [5],
                [4],
            )
        ),
        (
            genz_malik_1980_f_2,
            genz_malik_1980_f_2_exact,

            [0, 0],
            [10, 50],
            (
                [-3, 3],
                [-2, 2],
            ),
        ),
        (
            genz_malik_1980_f_2,
            genz_malik_1980_f_2_exact,
            [0, 0, 0],
            [1, 1, 1],
            (
                [1, 1, 1],
                [1, 1, 1],
            )
        ),
        (
            genz_malik_1980_f_2,
            genz_malik_1980_f_2_exact,
            [0, 0, 0],
            [1, 1, 1],
            (
                [2, 3, 4],
                [2, 3, 4],
            )
        ),
        (
            genz_malik_1980_f_2,
            genz_malik_1980_f_2_exact,
            [-1, -1, -1],
            [1, 1, 1],
            (
                [1, 1, 1],
                [2, 2, 2],
            )
        ),
        (
            genz_malik_1980_f_2,
            genz_malik_1980_f_2_exact,
            [-1, -1, -1, -1],
            [1, 1, 1, 1],
            (
                [1, 1, 1, 1],
                [1, 1, 1, 1],
            )
        ),

        # -- f3 --
        (
            genz_malik_1980_f_3,
            genz_malik_1980_f_3_exact,
            [-1],
            [1],
            (
                [1/2],
            ),
        ),
        (
            genz_malik_1980_f_3,
            genz_malik_1980_f_3_exact,
            [0, -1],
            [1, 1],
            (
                [5, 5],
            ),
        ),
        (
            genz_malik_1980_f_3,
            genz_malik_1980_f_3_exact,
            [-1, -1, -1],
            [1, 1, 1],
            (
                [1, 1, 1],
            ),
        ),

        # -- f4 --
        (
            genz_malik_1980_f_4,
            genz_malik_1980_f_4_exact,
            [0],
            [2],
            (
                [1],
            ),
        ),
        (
            genz_malik_1980_f_4,
            genz_malik_1980_f_4_exact,
            [0, 0],
            [2, 1],
            ([1, 1],),
        ),
        (
            genz_malik_1980_f_4,
            genz_malik_1980_f_4_exact,
            [0, 0, 0],
            [1, 1, 1],
            ([1, 1, 1],),
        ),

        # -- f5 --
        (
            genz_malik_1980_f_5,
            genz_malik_1980_f_5_exact,
            [-1],
            [1],
            (
                [-2],
                [2],
            ),
        ),
        (
            genz_malik_1980_f_5,
            genz_malik_1980_f_5_exact,
            [-1, -1],
            [1, 1],
            (
                [2, 3],
                [4, 5],
            ),
        ),
        (
            genz_malik_1980_f_5,
            genz_malik_1980_f_5_exact,
            [-1, -1],
            [1, 1],
            (
                [-1, 1],
                [0, 0],
            ),
        ),
        (
            genz_malik_1980_f_5,
            genz_malik_1980_f_5_exact,
            [-1, -1, -1],
            [1, 1, 1],
            (
                [1, 1, 1],
                [1, 1, 1],
            ),
        ),
    ])
    def test_scalar_output(self, problem, rule, rtol, atol, xp):
        f, exact, a, b, args = problem

        a = xp.asarray(a, dtype=xp.float64)
        b = xp.asarray(b, dtype=xp.float64)
        args = tuple(xp.asarray(arg, dtype=xp.float64) for arg in args)

        ndim = xp_size(a)

        if rule == "genz-malik" and ndim < 2:
            pytest.skip("Genz-Malik cubature does not support 1D integrals")

        res = cubature(
            f,
            a,
            b,
            rule=rule,
            rtol=rtol,
            atol=atol,
            args=(*args, xp),
        )

        assert res.status == "converged"

        est = res.estimate
        exact_sol = exact(a, b, *args, xp)

        xp_assert_close(
            est,
            exact_sol,
            rtol=rtol,
            atol=atol,
            err_msg=f"estimate_error={res.error}, subdivisions={res.subdivisions}",
        )

    @pytest.mark.parametrize("problem", [
        (
            # Function to integrate, like `f(x, *args)`
            genz_malik_1980_f_1,

            # Exact solution, like `exact(a, b, *args)`
            genz_malik_1980_f_1_exact,

            # Function that generates random args of a certain shape.
            genz_malik_1980_f_1_random_args,
        ),
        (
            genz_malik_1980_f_2,
            genz_malik_1980_f_2_exact,
            genz_malik_1980_f_2_random_args,
        ),
        (
            genz_malik_1980_f_3,
            genz_malik_1980_f_3_exact,
            genz_malik_1980_f_3_random_args
        ),
        (
            genz_malik_1980_f_4,
            genz_malik_1980_f_4_exact,
            genz_malik_1980_f_4_random_args
        ),
        (
            genz_malik_1980_f_5,
            genz_malik_1980_f_5_exact,
            genz_malik_1980_f_5_random_args,
        ),
    ])
    @pytest.mark.parametrize("shape", [
        (2,),
        (3,),
        (4,),
        (1, 2),
        (1, 3),
        (1, 4),
        (3, 2),
        (3, 4, 2),
        (2, 1, 3),
    ])
    def test_array_output(self, problem, rule, shape, rtol, atol, xp):
        rng = np_compat.random.default_rng(1)
        ndim = shape[-1]

        if rule == "genz-malik" and ndim < 2:
            pytest.skip("Genz-Malik cubature does not support 1D integrals")

        if rule == "genz-malik" and ndim >= 5:
            pytest.mark.slow("Gauss-Kronrod is slow in >= 5 dim")

        f, exact, random_args = problem
        args = random_args(rng, shape, xp)

        a = xp.asarray([0] * ndim, dtype=xp.float64)
        b = xp.asarray([1] * ndim, dtype=xp.float64)

        res = cubature(
            f,
            a,
            b,
            rule=rule,
            rtol=rtol,
            atol=atol,
            args=(*args, xp),
        )

        est = res.estimate
        exact_sol = exact(a, b, *args, xp)

        xp_assert_close(
            est,
            exact_sol,
            rtol=rtol,
            atol=atol,
            err_msg=f"estimate_error={res.error}, subdivisions={res.subdivisions}",
        )

        err_msg = (f"estimate_error={res.error}, "
                   f"subdivisions= {res.subdivisions}, "
                   f"true_error={xp.abs(res.estimate - exact_sol)}")
        assert res.status == "converged", err_msg

        assert res.estimate.shape == shape[:-1]


@array_api_compatible
class TestRules:
    """
    Tests related to the general Rule interface (currently private).
    """

    @pytest.mark.parametrize("problem", [
        (
            # 2D problem, 1D rule
            [0, 0],
            [1, 1],
            GaussKronrodQuadrature,
            (21,),
        ),
        (
            # 1D problem, 2D rule
            [0],
            [1],
            GenzMalikCubature,
            (2,),
        )
    ])
    def test_incompatible_dimension_raises_error(self, problem, xp):
        a, b, quadrature, quadrature_args = problem
        rule = quadrature(*quadrature_args, xp=xp)

        a = xp.asarray(a, dtype=xp.float64)
        b = xp.asarray(b, dtype=xp.float64)

        with pytest.raises(Exception, match="incompatible dimension"):
            rule.estimate(basic_1d_integrand, a, b, args=(xp,))

    def test_estimate_with_base_classes_raise_error(self, xp):
        a = xp.asarray([0])
        b = xp.asarray([1])

        for base_class in [Rule(), FixedRule()]:
            with pytest.raises(Exception):
                base_class.estimate(basic_1d_integrand, a, b, args=(xp,))


@array_api_compatible
class TestRulesQuadrature:
    """
    Tests underlying quadrature rules (ndim == 1).
    """

    @pytest.mark.parametrize(("rule", "rule_args"), [
        (GaussLegendreQuadrature, (3,)),
        (GaussLegendreQuadrature, (5,)),
        (GaussLegendreQuadrature, (10,)),
        (GaussKronrodQuadrature, (15,)),
        (GaussKronrodQuadrature, (21,)),
    ])
    def test_base_1d_quadratures_simple(self, rule, rule_args, xp):
        quadrature = rule(*rule_args, xp=xp)

        n = xp.arange(5, dtype=xp.float64)

        def f(x):
            x_reshaped = xp.reshape(x, (-1, 1, 1))
            n_reshaped = xp.reshape(n, (1, -1, 1))

            return x_reshaped**n_reshaped

        a = xp.asarray([0])
        b = xp.asarray([2])

        exact = xp.reshape(2**(n+1)/(n+1), (-1, 1))
        estimate = quadrature.estimate(f, a, b)

        xp_assert_close(
            estimate,
            exact,
            rtol=1e-1,
            atol=0,
        )

    @pytest.mark.parametrize(("rule_pair", "rule_pair_args"), [
        ((GaussLegendreQuadrature, GaussLegendreQuadrature), (10, 5)),
    ])
    def test_base_1d_quadratures_error_from_difference(self, rule_pair, rule_pair_args,
                                                       xp):
        n = xp.arange(5, dtype=xp.float64)
        a = xp.asarray([0], dtype=xp.float64)
        b = xp.asarray([2], dtype=xp.float64)

        higher = rule_pair[0](rule_pair_args[0], xp=xp)
        lower = rule_pair[1](rule_pair_args[1], xp=xp)

        rule = NestedFixedRule(higher, lower)
        res = cubature(
            basic_1d_integrand,
            a, b,
            rule=rule,
            rtol=1e-1,
            args=(n, xp),
        )

        xp_assert_close(
            res.estimate,
            basic_1d_integrand_exact(n, xp),
            rtol=1e-1,
            atol=0,
        )

    @pytest.mark.parametrize("quadrature", [
        GaussLegendreQuadrature
    ])
    def test_one_point_fixed_quad_impossible(self, quadrature, xp):
        with pytest.raises(Exception):
            quadrature(1, xp=xp)


@array_api_compatible
class TestRulesCubature:
    """
    Tests underlying cubature rules (ndim >= 2).
    """

    @pytest.mark.parametrize("ndim", range(2, 11))
    def test_genz_malik_func_evaluations(self, ndim, xp):
        """
        Tests that the number of function evaluations required for Genz-Malik cubature
        matches the number in Genz and Malik 1980.
        """

        nodes, _ = GenzMalikCubature(ndim, xp=xp).nodes_and_weights

        assert nodes.shape[0] == (2**ndim) + 2*ndim**2 + 2*ndim + 1

    def test_genz_malik_1d_raises_error(self, xp):
        with pytest.raises(Exception, match="only defined for ndim >= 2"):
            GenzMalikCubature(1, xp=xp)


class BadErrorRule(Rule):
    """
    A rule with fake high error so that cubature will keep on subdividing.
    """

    def estimate(self, f, a, b, args=()):
        xp = array_namespace(a, b)
        underlying = GaussLegendreQuadrature(10, xp=xp)

        return underlying.estimate(f, a, b, args)

    def estimate_error(self, f, a, b, args=()):
        xp = array_namespace(a, b)
        return xp.asarray(1e6, dtype=xp.float64)
