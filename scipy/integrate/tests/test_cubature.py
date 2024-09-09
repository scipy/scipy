import math
import scipy
import itertools

import pytest
import numpy as np

from numpy.testing import assert_allclose

from scipy.integrate import cubature

from scipy.integrate._rules import (
    Rule, FixedRule,
    NestedFixedRule,
    GaussLegendreQuadrature, GaussKronrodQuadrature,
    GenzMalikCubature,
)

# The integrands ``genz_malik_1980_*`` come from the paper:
#   A.C. Genz, A.A. Malik, Remarks on algorithm 006: An adaptive algorithm for
#   numerical integration over an N-dimensional rectangular region, Journal of
#   Computational and Applied Mathematics, Volume 6, Issue 4, 1980, Pages 295-302,
#   ISSN 0377-0427, https://doi.org/10.1016/0771-050X(80)90039-X.


def basic_1d_integrand(x, n):
    x_reshaped = x.reshape(-1, 1, 1)
    n_reshaped = n.reshape(1, -1, 1)

    return np.power(x_reshaped, n_reshaped)


def basic_1d_integrand_exact(n):
    return (2**(n+1)/(n+1)).reshape(-1, 1)


def basic_nd_integrand(x, n):
    return np.power(np.sum(x, axis=-1).reshape(-1, 1), n.reshape(1, -1))


def basic_nd_integrand_exact(n):
    return (-2**(3+n) + 4**(2+n))/((1+n)*(2+n))


def genz_malik_1980_f_1(x, r, alphas):
    r"""
    .. math:: f_1(\mathbf x) = \cos\left(2\pi r + \sum^n_{i = 1}\alpha_i x_i\right)

    .. code-block:: mathematica

        genzMalik1980f1[x_List, r_, alphas_List] := Cos[2*Pi*r + Total[x*alphas]]
    """

    npoints, ndim = x.shape[0], x.shape[-1]

    alphas_reshaped = alphas[np.newaxis, :]
    x_reshaped = x.reshape(npoints, *([1]*(len(alphas.shape) - 1)), ndim)

    return np.cos(2*np.pi*r + np.sum(alphas_reshaped * x_reshaped, axis=-1))


def genz_malik_1980_f_1_exact(a, b, r, alphas):
    ndim = len(a)
    a = a.reshape(*([1]*(len(alphas.shape) - 1)), ndim)
    b = b.reshape(*([1]*(len(alphas.shape) - 1)), ndim)

    return (-2)**ndim * 1/np.prod(alphas, axis=-1) \
        * np.cos(2*np.pi*r + np.sum(alphas * (a+b)/2, axis=-1)) \
        * np.prod(np.sin(alphas * (a-b)/2), axis=-1)


def genz_malik_1980_f_1_random_args(shape):
    r = np.random.rand(*shape[:-1])
    alphas = np.random.rand(*shape)

    difficulty = 9
    normalisation_factors = np.expand_dims(np.sum(alphas, axis=-1), -1)
    alphas = difficulty * alphas / normalisation_factors

    return (r, alphas)


def genz_malik_1980_f_2(x, alphas, betas):
    r"""
    .. math:: f_2(\mathbf x) = \prod^n_{i = 1} (\alpha_i^2 + (x_i - \beta_i)^2)^{-1}

    .. code-block:: mathematica

        genzMalik1980f2[x_List, alphas_List, betas_List] :=
            1/Times @@ ((alphas^2 + (x - betas)^2))
    """
    npoints, ndim = x.shape[0], x.shape[-1]

    alphas_reshaped = alphas[np.newaxis, :]
    betas_reshaped = betas[np.newaxis, :]

    x_reshaped = x.reshape(npoints, *([1]*(len(alphas.shape) - 1)), ndim)

    return 1/np.prod(alphas_reshaped**2 + (x_reshaped-betas_reshaped)**2, axis=-1)


def genz_malik_1980_f_2_exact(a, b, alphas, betas):
    ndim = len(a)
    a = a.reshape(*([1]*(len(alphas.shape) - 1)), ndim)
    b = b.reshape(*([1]*(len(alphas.shape) - 1)), ndim)

    return (-1)**ndim * 1/np.prod(alphas, axis=-1) \
        * np.prod(
            np.arctan((a - betas)/alphas) - np.arctan((b - betas)/alphas),
            axis=-1
        )


def genz_malik_1980_f_2_random_args(shape):
    ndim = shape[-1]
    alphas = np.random.rand(*shape)
    betas = np.random.rand(*shape)

    difficulty = 25
    products = np.prod(np.power(alphas, -2), axis=-1)
    normalisation_factors = np.expand_dims(np.power(products, 1 / (2*ndim)), axis=-1)
    alphas = alphas \
        * normalisation_factors \
        / np.power(difficulty, 1 / (2*ndim))

    assert_allclose(np.prod(np.power(alphas, -2), axis=-1), difficulty)

    # Adjust alphas from distribution used in Genz and Malik 1980 since denominator
    # is very small for high dimensions.
    alphas *= 10

    return alphas, betas


def genz_malik_1980_f_3(x, alphas):
    r"""
    .. math:: f_3(\mathbf x) = \exp\left(\sum^n_{i = 1} \alpha_i x_i\right)

    .. code-block:: mathematica

        genzMalik1980f3[x_List, alphas_List] := Exp[Dot[x, alphas]]
    """

    npoints, ndim = x.shape[0], x.shape[-1]

    alphas_reshaped = alphas[np.newaxis, :]
    x_reshaped = x.reshape(npoints, *([1]*(len(alphas.shape) - 1)), ndim)

    return np.exp(np.sum(alphas_reshaped * x_reshaped, axis=-1))


def genz_malik_1980_f_3_exact(a, b, alphas):
    ndim = len(a)
    a = a.reshape(*([1]*(len(alphas.shape) - 1)), ndim)
    b = b.reshape(*([1]*(len(alphas.shape) - 1)), ndim)

    return (-1)**ndim * 1/np.prod(alphas, axis=-1) \
        * np.prod(np.exp(alphas * a) - np.exp(alphas * b), axis=-1)


def genz_malik_1980_f_3_random_args(shape):
    alphas = np.random.rand(*shape)
    normalisation_factors = np.expand_dims(np.sum(alphas, axis=-1), -1)
    difficulty = 12
    alphas = difficulty * alphas / normalisation_factors

    return (alphas,)


def genz_malik_1980_f_4(x, alphas):
    r"""
    .. math:: f_4(\mathbf x) = \left(1 + \sum^n_{i = 1} \alpha_i x_i\right)^{-n-1}

    .. code-block:: mathematica
        genzMalik1980f4[x_List, alphas_List] :=
            (1 + Dot[x, alphas])^(-Length[alphas] - 1)
    """

    npoints, ndim = x.shape[0], x.shape[-1]

    alphas_reshaped = alphas[np.newaxis, :]
    x_reshaped = x.reshape(npoints, *([1]*(len(alphas.shape) - 1)), ndim)

    return ((1 + np.sum(alphas_reshaped * x_reshaped, axis=-1))**(-ndim-1))


def genz_malik_1980_f_4_exact(a, b, alphas):
    ndim = len(a)

    def F(x):
        x_reshaped = x.reshape(*([1]*(len(alphas.shape) - 1)), ndim)

        return (-1)**ndim/np.prod(alphas, axis=-1) / (
            math.factorial(ndim) * (1 + np.sum(alphas * x_reshaped, axis=-1))
        )

    return _eval_indefinite_integral(F, a, b)


def _eval_indefinite_integral(F, a, b):
    """
    Calculates a definite integral from points `a` to `b` by summing up over the corners
    of the corresponding hyperrectangle.
    """

    ndim = len(a)
    points = np.stack([a, b], axis=0)

    out = 0
    for ind in itertools.product(range(2), repeat=ndim):
        out += pow(-1, sum(ind) + ndim) * F(points[ind, tuple(range(ndim))])

    return out


def genz_malik_1980_f_4_random_args(shape):
    ndim = shape[-1]

    alphas = np.random.rand(*shape)
    normalisation_factors = np.expand_dims(np.sum(alphas, axis=-1), -1)
    difficulty = 14
    alphas = (difficulty/ndim) * alphas / normalisation_factors

    return (alphas,)


def genz_malik_1980_f_5(x, alphas, betas):
    r"""
    .. math::

        f_5(\mathbf x) = \exp\left(-\sum^n_{i = 1} \alpha^2_i (x_i - \beta_i)^2\right)

    .. code-block:: mathematica

        genzMalik1980f5[x_List, alphas_List, betas_List] :=
            Exp[-Total[alphas^2 * (x - betas)^2]]
    """

    npoints, ndim = x.shape[0], x.shape[-1]

    alphas_reshaped = alphas[np.newaxis, :]
    betas_reshaped = betas[np.newaxis, :]

    x_reshaped = x.reshape(npoints, *([1]*(len(alphas.shape) - 1)), ndim)

    return np.exp(
        -np.sum(alphas_reshaped**2 * (x_reshaped - betas_reshaped)**2, axis=-1)
    )


def genz_malik_1980_f_5_exact(a, b, alphas, betas):
    ndim = len(a)
    a = a.reshape(*([1]*(len(alphas.shape) - 1)), ndim)
    b = b.reshape(*([1]*(len(alphas.shape) - 1)), ndim)

    return (1/2)**ndim * 1/np.prod(alphas, axis=-1) \
        * np.power(np.pi, ndim/2) \
        * np.prod(
            scipy.special.erf(alphas * (betas - a))
            + scipy.special.erf(alphas * (b - betas)),
            axis=-1
        )


def genz_malik_1980_f_5_random_args(shape):
    alphas = np.random.rand(*shape)
    betas = np.random.rand(*shape)

    difficulty = 21
    normalisation_factors = np.expand_dims(
        np.sqrt(np.sum(np.power(alphas, 2), axis=-1)), -1
    )
    alphas = alphas / normalisation_factors * np.sqrt(difficulty)

    return alphas, betas


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
    def test_pass_str(self, rule_str):
        n = np.arange(5)
        a = np.array([0, 0])
        b = np.array([2, 2])

        res = cubature(basic_nd_integrand, a, b, rule_str, args=(n,))

        assert_allclose(
            res.estimate,
            basic_nd_integrand_exact(n),
            rtol=1e-1,
            atol=0,
        )

    def test_pass_list_not_array(self):
        n = np.arange(5)

        a = [0]
        b = [2]

        res = cubature(
            basic_1d_integrand,
            a,
            b,
            args=(n,)
        )

        assert_allclose(
            res.estimate,
            basic_1d_integrand_exact(n),
            rtol=1e-1,
            atol=0,
        )

    def test_stops_after_max_subdivisions(self):
        a = np.array([0])
        b = np.array([1])
        rule = BadErrorRule()

        res = cubature(
            basic_1d_integrand,  # Any function would suffice
            a,
            b,
            rule,
            max_subdivisions=10,
            args=(np.arange(5),),
        )

        assert res.subdivisions == 10
        assert res.status == "not_converged"

    def test_a_and_b_must_be_1d(self):
        a = np.array([[0]])
        b = np.array([[1]])

        with pytest.raises(Exception, match="`a` and `b` must be 1D arrays"):
            cubature(basic_1d_integrand, a, b)


@pytest.mark.parametrize("rtol", [1e-4])
@pytest.mark.parametrize("atol", [1e-5])
@pytest.mark.parametrize("rule", [
    "gk15",
    "gk21",
    "genz-malik",
])
class TestCubatureProblems:
    """
    Tests that `cubature` gives the correct answer.
    """

    problems_scalar_output = [
        # -- f1 --
        (
            # Function to integrate, like `f(x, *args)`
            genz_malik_1980_f_1,

            # Exact solution, like `exact(a, b, *args)`
            genz_malik_1980_f_1_exact,

            # Coordinates of `a`
            np.array([0]),

            # Coordinates of `b`
            np.array([10]),

            # Arguments to pass to `f` and `exact`
            (
                np.array([1/4]),
                np.array([5]),
            )
        ),
        (
            genz_malik_1980_f_1,
            genz_malik_1980_f_1_exact,
            np.array([0, 0]),
            np.array([1, 1]),
            (
                np.array([1/4]),
                np.array([2, 4]),
            ),
        ),
        (
            genz_malik_1980_f_1,
            genz_malik_1980_f_1_exact,
            np.array([0, 0]),
            np.array([5, 5]),
            (
                np.array([1/2]),
                np.array([2, 4]),
            )
        ),
        (
            genz_malik_1980_f_1,
            genz_malik_1980_f_1_exact,
            np.array([0, 0, 0]),
            np.array([10, 10, 10]),
            (
                np.array([1/2]),
                np.array([1, 1, 1]),
            )
        ),

        # -- f2 --
        (
            genz_malik_1980_f_2,
            genz_malik_1980_f_2_exact,
            np.array([-1]),
            np.array([1]),
            (
                np.array([5]),
                np.array([4]),
            )
        ),
        (
            genz_malik_1980_f_2,
            genz_malik_1980_f_2_exact,

            np.array([10, 50]),
            np.array([10, 50]),
            (
                np.array([-3, 3]),
                np.array([-2, 2])
            ),
        ),
        (
            genz_malik_1980_f_2,
            genz_malik_1980_f_2_exact,
            np.array([0, 0, 0]),
            np.array([1, 1, 1]),
            (
                np.array([1, 1, 1]),
                np.array([1, 1, 1]),
            )
        ),
        (
            genz_malik_1980_f_2,
            genz_malik_1980_f_2_exact,
            np.array([0, 0, 0]),
            np.array([1, 1, 1]),
            (
                np.array([2, 3, 4]),
                np.array([2, 3, 4]),
            )
        ),
        (
            genz_malik_1980_f_2,
            genz_malik_1980_f_2_exact,
            np.array([-1, -1, -1]),
            np.array([1, 1, 1]),
            (
                np.array([1, 1, 1]),
                np.array([2, 2, 2]),
            )
        ),
        (
            genz_malik_1980_f_2,
            genz_malik_1980_f_2_exact,
            np.array([-1, -1, -1, -1]),
            np.array([1, 1, 1, 1]),
            (
                np.array([1, 1, 1, 1]),
                np.array([1, 1, 1, 1]),
            )
        ),

        # -- f3 --
        (
            genz_malik_1980_f_3,
            genz_malik_1980_f_3_exact,
            np.array([-1]),
            np.array([1]),
            (
                np.array([1/2]),
            ),
        ),
        (
            genz_malik_1980_f_3,
            genz_malik_1980_f_3_exact,
            np.array([0, -1]),
            np.array([1, 1]),
            (
                np.array([5, 5]),
            ),
        ),
        (
            genz_malik_1980_f_3,
            genz_malik_1980_f_3_exact,
            np.array([-1, -1, -1]),
            np.array([1, 1, 1]),
            (
                np.array([1, 1, 1]),
            ),
        ),

        # -- f4 --
        (
            genz_malik_1980_f_4,
            genz_malik_1980_f_4_exact,
            np.array([0]),
            np.array([2]),
            (
                np.array([1]),
            ),
        ),
        (
            genz_malik_1980_f_4,
            genz_malik_1980_f_4_exact,
            np.array([0, 0]),
            np.array([2, 1]),
            (np.array([1, 1]),),
        ),
        (
            genz_malik_1980_f_4,
            genz_malik_1980_f_4_exact,
            np.array([0, 0, 0]),
            np.array([1, 1, 1]),
            (np.array([1, 1, 1]),),
        ),

        # -- f5 --
        (
            genz_malik_1980_f_5,
            genz_malik_1980_f_5_exact,
            np.array([-1]),
            np.array([1]),
            (
                np.array([-2]),
                np.array([2]),
            ),
        ),
        (
            genz_malik_1980_f_5,
            genz_malik_1980_f_5_exact,
            np.array([-1, -1]),
            np.array([1, 1]),
            (
                np.array([2, 3]),
                np.array([4, 5]),
            ),
        ),
        (
            genz_malik_1980_f_5,
            genz_malik_1980_f_5_exact,
            np.array([-1, -1]),
            np.array([1, 1]),
            (
                np.array([-1, 1]),
                np.array([0, 0]),
            ),
        ),
        (
            genz_malik_1980_f_5,
            genz_malik_1980_f_5_exact,
            np.array([-1, -1, -1]),
            np.array([1, 1, 1]),
            (
                np.array([1, 1, 1]),
                np.array([1, 1, 1]),
            ),
        ),
    ]

    problem_array_output = [
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
    ]

    @pytest.mark.parametrize("problem", problems_scalar_output)
    def test_scalar_output(self, problem, rule, rtol, atol):
        f, exact, a, b, args = problem

        ndim = len(a)

        if rule == "genz-malik" and ndim < 2:
            pytest.skip("Genz-Malik cubature does not support 1D integrals")

        res = cubature(
            f,
            a,
            b,
            rule,
            rtol,
            atol,
            args=args,
        )

        assert res.status == "converged"

        assert_allclose(
            res.estimate,
            exact(a, b, *args),
            rtol=rtol,
            atol=atol,
            err_msg=f"estimate_error={res.error}, subdivisions={res.subdivisions}"
        )

    @pytest.mark.parametrize("problem", problem_array_output)
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
    def test_array_output(self, problem, rule, shape, rtol, atol):
        np.random.seed(1)
        ndim = shape[-1]

        if rule == "genz-malik" and ndim < 2:
            pytest.skip("Genz-Malik cubature does not support 1D integrals")

        if rule == "genz-malik" and ndim >= 5:
            pytest.mark.slow("Gauss-Kronrod is slow in >= 5 dim")

        f, exact, random_args = problem
        args = random_args(shape)

        a = np.array([0] * ndim)
        b = np.array([1] * ndim)

        res = cubature(
            f,
            a,
            b,
            rule,
            rtol,
            atol,
            args=args,
        )

        assert_allclose(
            res.estimate,
            exact(a, b, *args),
            rtol=rtol,
            atol=atol,
            err_msg=f"estimate_error={res.error}, subdivisions={res.subdivisions}"
        )

        err_msg = (f"estimate_error={res.error}, "
                   f"subdivisions= {res.subdivisions}, "
                   f"true_error={np.abs(res.estimate - exact(a, b, *args))}")
        assert res.status == "converged", err_msg

        assert res.estimate.shape == shape[:-1]


class TestRules:
    """
    Tests related to the general Rule interface (currently private).
    """

    @pytest.mark.parametrize("problem", [
        (
            # 2D problem, 1D rule
            np.array([0, 0]),
            np.array([1, 1]),
            GaussKronrodQuadrature(21),
        ),
        (
            # 2D problem, 1D rule
            np.array([0, 0]),
            np.array([1, 1]),
            NestedFixedRule(GaussKronrodQuadrature(21), GaussKronrodQuadrature(21)),
        ),
        (
            # 1D problem, 2D rule
            np.array([0]),
            np.array([1]),
            GenzMalikCubature(2),
        ),
    ])
    def test_incompatible_dimension_raises_error(self, problem):
        a, b, rule = problem

        with pytest.raises(Exception, match="incompatible dimension"):
            rule.estimate(basic_1d_integrand, a, b)

    def test_estimate_with_base_classes_raise_error(self):
        a = np.array([0])
        b = np.array([1])

        for base_class in [Rule(), FixedRule()]:
            with pytest.raises(Exception):
                base_class.estimate(basic_1d_integrand, a, b)


class TestRulesQuadrature:
    """
    Tests underlying quadrature rules (ndim == 1).
    """

    @pytest.mark.parametrize("quadrature", [
        GaussLegendreQuadrature(3),
        GaussLegendreQuadrature(5),
        GaussLegendreQuadrature(10),
        GaussKronrodQuadrature(15),
        GaussKronrodQuadrature(21),
    ])
    def test_base_1d_quadratures_simple(self, quadrature):
        n = np.arange(5)

        def f(x):
            x_reshaped = x.reshape(-1, 1, 1)
            n_reshaped = n.reshape(1, -1, 1)

            return np.power(x_reshaped, n_reshaped)

        a = np.array([0])
        b = np.array([2])

        exact = (2**(n+1)/(n+1)).reshape(-1, 1)
        estimate = quadrature.estimate(f, a, b)

        assert_allclose(
            estimate,
            exact,
            rtol=1e-1,
            atol=0,
        )

    @pytest.mark.parametrize("quadrature_pair", [
        (GaussLegendreQuadrature(10), GaussLegendreQuadrature(5))
    ])
    def test_base_1d_quadratures_error_from_difference(self, quadrature_pair):
        n = np.arange(5)
        a = np.array([0])
        b = np.array([2])

        rule = NestedFixedRule(
            higher=quadrature_pair[0],
            lower=quadrature_pair[1]
        )

        res = cubature(basic_1d_integrand, a, b, rule, rtol=1e-1, args=(np.arange(5),))

        assert_allclose(
            res.estimate,
            basic_1d_integrand_exact(n),
            rtol=1e-1,
            atol=0,
        )

    @pytest.mark.parametrize("quadrature", [
        GaussLegendreQuadrature
    ])
    def test_one_point_fixed_quad_impossible(self, quadrature):
        with pytest.raises(Exception):
            quadrature(1)


class TestRulesCubature:
    """
    Tests underlying cubature rules (ndim >= 2).
    """

    @pytest.mark.parametrize("ndim", range(2, 11))
    def test_genz_malik_func_evaluations(self, ndim):
        """
        Tests that the number of function evaluations required for Genz-Malik cubature
        matches the number in Genz and Malik 1980.
        """

        nodes, _ = GenzMalikCubature(ndim).nodes_and_weights

        assert nodes.shape[0] == (2**ndim) + 2*ndim**2 + 2*ndim + 1

    def test_genz_malik_1d_raises_error(self):
        with pytest.raises(Exception, match="only defined for ndim >= 2"):
            GenzMalikCubature(1)


class BadErrorRule(Rule):
    """
    A rule with fake high error so that cubature will keep on subdividing.
    """

    underlying = GaussLegendreQuadrature(10)

    def estimate(self, f, a, b, args=()):
        return self.underlying.estimate(f, a, b, args)

    def estimate_error(self, f, a, b, args=()):
        return 1e6
