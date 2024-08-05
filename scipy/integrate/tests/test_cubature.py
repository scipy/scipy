import math
import scipy
import itertools

import pytest
import numpy as np

from numpy.testing import assert_allclose

from scipy.integrate._rules import (
    Rule, FixedRule, ProductNestedFixed,
    NestedFixedRule, NestedRule,
    NewtonCotesQuad, GaussLegendreQuad, GaussKronrodQuad, GenzMalikCub
)

from scipy.integrate._cubature import cub


def genz_malik_1980_f_1(x, r, alphas):
    r"""
    ``f_1`` from Genz and Malik 1980.

    Notes
    -----

    .. math:: f_1(\mathbf x) = \cos\left(2\pi r + \sum^n_{i = 1}\alpha_i x_i\right)

    .. code-block:: mathematica

        genzMalik1980f1[x_List, r_, alphas_List] := Cos[2*Pi*r + Total[x*alphas]]

    References
    ----------
    [1] A.C. Genz, A.A. Malik, Remarks on algorithm 006: An adaptive algorithm for
        numerical integration over an N-dimensional rectangular region, Journal of
        Computational and Applied Mathematics, Volume 6, Issue 4, 1980, Pages 295-302,
        ISSN 0377-0427, https://doi.org/10.1016/0771-050X(80)90039-X.
    """

    ndim = x.shape[-1]
    num_eval_points = x.shape[0]

    alphas_reshaped = alphas[np.newaxis, :]
    x_reshaped = x.reshape(num_eval_points, *([1]*(len(alphas.shape) - 1)), ndim)

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
    `f_2` from Genz and Malik 1980.

    Notes
    -----

    .. math:: f_2(\mathbf x) = \prod^n_{i = 1} (\alpha_i^2 + (x_i - \beta_i)^2)^{-1}

    .. code-block:: mathematica

        genzMalik1980f2[x_List, alphas_List, betas_List] :=
            1/Times @@ ((alphas^2 + (x - betas)^2))

    References
    ----------
    [1] A.C. Genz, A.A. Malik, Remarks on algorithm 006: An adaptive algorithm for
        numerical integration over an N-dimensional rectangular region, Journal of
        Computational and Applied Mathematics, Volume 6, Issue 4, 1980, Pages 295-302,
        ISSN 0377-0427, https://doi.org/10.1016/0771-050X(80)90039-X.
    """
    ndim = x.shape[-1]
    num_eval_points = x.shape[0]

    alphas_reshaped = alphas[np.newaxis, :]
    betas_reshaped = betas[np.newaxis, :]

    x_reshaped = x.reshape(num_eval_points, *([1]*(len(alphas.shape) - 1)), ndim)

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
    `f_3` from Genz and Malik 1980.

    Notes
    -----

    .. math:: f_3(\mathbf x) = \exp\left(\sum^n_{i = 1} \alpha_i x_i\right)

    .. code-block:: mathematica

        genzMalik1980f3[x_List, alphas_List] := Exp[Dot[x, alphas]]

    References
    ----------
    [1] A.C. Genz, A.A. Malik, Remarks on algorithm 006: An adaptive algorithm for
        numerical integration over an N-dimensional rectangular region, Journal of
        Computational and Applied Mathematics, Volume 6, Issue 4, 1980, Pages 295-302,
        ISSN 0377-0427, https://doi.org/10.1016/0771-050X(80)90039-X.
    """

    ndim = x.shape[-1]
    num_eval_points = x.shape[0]

    alphas_reshaped = alphas[np.newaxis, :]
    x_reshaped = x.reshape(num_eval_points, *([1]*(len(alphas.shape) - 1)), ndim)

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
    `f_4` from Genz and Malik 1980.

    Notes
    -----

    .. math:: f_4(\mathbf x) = \left(1 + \sum^n_{i = 1} \alpha_i x_i\right)^{-n-1}

    .. code-block:: mathematica
        genzMalik1980f4[x_List, alphas_List] :=
            (1 + Dot[x, alphas])^(-Length[alphas] - 1)

    References
    ----------
    [1] A.C. Genz, A.A. Malik, Remarks on algorithm 006: An adaptive algorithm for
        numerical integration over an N-dimensional rectangular region, Journal of
        Computational and Applied Mathematics, Volume 6, Issue 4, 1980, Pages 295-302,
        ISSN 0377-0427, https://doi.org/10.1016/0771-050X(80)90039-X.
    """

    ndim = x.shape[-1]
    num_eval_points = x.shape[0]

    alphas_reshaped = alphas[np.newaxis, :]
    x_reshaped = x.reshape(num_eval_points, *([1]*(len(alphas.shape) - 1)), ndim)

    return ((1 + np.sum(alphas_reshaped * x_reshaped, axis=-1))**(-ndim-1))


def genz_malik_1980_f_4_exact(a, b, alphas):
    ndim = len(a)

    def F(x):
        x_reshaped = x.reshape(*([1]*(len(alphas.shape) - 1)), ndim)

        return (-1)**ndim/np.prod(alphas, axis=-1) / (
            math.factorial(ndim) * (1 + np.sum(alphas * x_reshaped, axis=-1))
        )

    return _eval_indefinite_integral(F, a, b)


def genz_malik_1980_f_4_random_args(shape):
    ndim = shape[-1]

    alphas = np.random.rand(*shape)
    normalisation_factors = np.expand_dims(np.sum(alphas, axis=-1), -1)
    difficulty = 14
    alphas = (difficulty/ndim) * alphas / normalisation_factors

    return (alphas,)


def genz_malik_1980_f_5(x, alphas, betas):
    r"""
    `f_5` from Genz and Malik 1980.

    Notes
    -----

    .. math::

        f_5(\mathbf x) = \exp\left(-\sum^n_{i = 1} \alpha^2_i (x_i - \beta_i)^2\right)

    .. code-block:: mathematica

        genzMalik1980f5[x_List, alphas_List, betas_List] :=
            Exp[-Total[alphas^2 * (x - betas)^2]]

    References
    ----------
    [1] A.C. Genz, A.A. Malik, Remarks on algorithm 006: An adaptive algorithm for
        numerical integration over an N-dimensional rectangular region, Journal of
        Computational and Applied Mathematics, Volume 6, Issue 4, 1980, Pages 295-302,
        ISSN 0377-0427, https://doi.org/10.1016/0771-050X(80)90039-X.
    """

    ndim = x.shape[-1]
    num_eval_points = x.shape[0]

    alphas_reshaped = alphas[np.newaxis, :]
    betas_reshaped = betas[np.newaxis, :]

    x_reshaped = x.reshape(num_eval_points, *([1]*(len(alphas.shape) - 1)), ndim)

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
        np.array([10, 10]),
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


@pytest.mark.parametrize("problem", problems_scalar_output)
@pytest.mark.parametrize("quadrature", [
    GaussKronrodQuad(15),
    GaussKronrodQuad(21),
    GenzMalikCub
])
@pytest.mark.parametrize("rtol", [1e-4])
@pytest.mark.parametrize("atol", [1e-5])
def test_cub_scalar_output(problem, quadrature, rtol, atol):
    f, exact, a, b, args = problem

    ndim = len(a)

    if quadrature is GenzMalikCub and ndim < 2:
        pytest.skip("Genz-Malik cubature does not support 1D integrals")

    if isinstance(quadrature, GaussKronrodQuad):
        rule = ProductNestedFixed([quadrature] * ndim)
    elif quadrature is GenzMalikCub and ndim >= 2:
        rule = GenzMalikCub(ndim)
    else:
        raise "Unknown quadrature rule specified"

    res = cub(
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


problems_array_output = [
    (
        # Function to integrate, like `f(x, *args)`
        genz_malik_1980_f_1,

        # Exact solution, like `exact(a, b, *args)`
        genz_malik_1980_f_1_exact,

        # Function that generates random args of a certain shape, like `random(shape)`.
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


@pytest.mark.parametrize("problem", problems_array_output)
@pytest.mark.parametrize("quadrature", [
    GaussKronrodQuad(15),
    GaussKronrodQuad(21),
    GenzMalikCub
])
@pytest.mark.parametrize("shape", [
    (2,),
    (3,),
    (4,),
    (1, 2),
    (1, 3),
    (1, 4),
    (1, 5),
    (3, 2),
    (3, 4, 2),
    (2, 1, 3),
])
@pytest.mark.parametrize("rtol", [1e-3])
@pytest.mark.parametrize("atol", [1e-4])
def test_cub_array_output(problem, quadrature, shape, rtol, atol):
    np.random.seed(1)
    ndim = shape[-1]

    if quadrature is GenzMalikCub and ndim < 2:
        pytest.skip("Genz-Malik cubature does not support 1D integrals")

    if isinstance(quadrature, GaussKronrodQuad):
        rule = ProductNestedFixed([quadrature] * ndim)
    elif quadrature is GenzMalikCub and ndim >= 2:
        rule = GenzMalikCub(ndim)
    else:
        raise "Unknown quadrature rule specified"

    f, exact, random_args = problem
    args = random_args(shape)

    a = np.array([0] * ndim)
    b = np.array([1] * ndim)

    res = cub(
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

    assert res.status == "converged", f"estimate_error={res.error}, subdivisions=\
{res.subdivisions}, true_error={np.abs(res.estimate - exact(a, b, *args))}"
    assert res.estimate.shape == shape[:-1]


@pytest.mark.parametrize("ndim", range(2, 11))
def test_genz_malik_func_evaluations(ndim):
    """
    Tests that the number of function evaluations required for Genz-Malik cubature
    matches the number in Genz and Malik 1980.
    """

    nodes, _ = GenzMalikCub(ndim).nodes_and_weights

    assert nodes.shape[0] == (2**ndim) + 2*ndim**2 + 2*ndim + 1


@pytest.mark.parametrize("quadrature", [
    NewtonCotesQuad(3),
    NewtonCotesQuad(5),
    NewtonCotesQuad(10),
    NewtonCotesQuad(3, open=True),
    NewtonCotesQuad(5, open=True),
    NewtonCotesQuad(10, open=True),
    GaussLegendreQuad(3),
    GaussLegendreQuad(5),
    GaussLegendreQuad(10),
    GaussKronrodQuad(15),
    GaussKronrodQuad(21),
])
def test_base_1d_quadratures_simple(quadrature):
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


@pytest.mark.parametrize("quadrature_str", [
    "gk21",
    "gk15",
    "trapezoid"
])
def test_can_pass_str_to_cub(quadrature_str):
    n = np.arange(5)

    def f(x):
        x_reshaped = x.reshape(-1, 1, 1)
        n_reshaped = n.reshape(1, -1, 1)

        return np.power(x_reshaped, n_reshaped)

    a = np.array([0])
    b = np.array([2])

    exact = (2**(n+1)/(n+1)).reshape(-1, 1)
    res = cub(f, a, b, quadrature_str)

    assert_allclose(
        res.estimate,
        exact,
        rtol=1e-1,
        atol=0,
    )


def test_can_pass_list_to_cub():
    n = np.arange(5)

    def f(x):
        x_reshaped = x.reshape(-1, 1, 1)
        n_reshaped = n.reshape(1, -1, 1)

        return np.power(x_reshaped, n_reshaped)

    a = [0]
    b = [2]

    exact = (2**(n+1)/(n+1)).reshape(-1, 1)
    res = cub(f, a, b)

    assert_allclose(
        res.estimate,
        exact,
        rtol=1e-1,
        atol=0,
    )


@pytest.mark.parametrize("quadrature_pair", [
    (NewtonCotesQuad(10), NewtonCotesQuad(5)),
    (GaussLegendreQuad(10), GaussLegendreQuad(5))
])
def test_base_1d_quadratures_error_from_difference(quadrature_pair):
    n = np.arange(5)

    def f(x):
        x_reshaped = x.reshape(-1, 1, 1)
        n_reshaped = n.reshape(1, -1, 1)

        return np.power(x_reshaped, n_reshaped)

    a = np.array([0])
    b = np.array([2])

    exact = (2**(n+1)/(n+1)).reshape(-1, 1)

    rule = NestedFixedRule(
        higher=quadrature_pair[0],
        lower=quadrature_pair[1]
    )

    res = cub(f, a, b, rule, rtol=1e-1)

    assert_allclose(
        res.estimate,
        exact,
        rtol=1e-1,
        atol=0,
    )


@pytest.mark.parametrize("rule", [
    ProductNestedFixed([
        NestedFixedRule(NewtonCotesQuad(10), NewtonCotesQuad(8)),
        NestedFixedRule(NewtonCotesQuad(10), NewtonCotesQuad(8)),
    ]),
    ProductNestedFixed([
        NestedFixedRule(GaussLegendreQuad(10), NewtonCotesQuad(5)),
        NestedFixedRule(GaussLegendreQuad(10), GaussLegendreQuad(5)),
    ]),
    ProductNestedFixed([
        GaussKronrodQuad(21),
        GaussKronrodQuad(21),
    ]),
    GenzMalikCub(2),
])
def test_cub_with_kwargs(rule):
    np.random.seed(1)

    f = genz_malik_1980_f_1
    r, alphas = genz_malik_1980_f_1_random_args((3, 2))

    a = np.array([0, 0])
    b = np.array([1, 1])

    res = cub(f, a, b, rule, kwargs={
        "r": r,
        "alphas": alphas
    })
    exact = genz_malik_1980_f_1_exact(a, b, r, alphas)

    assert_allclose(
        res.estimate,
        exact,
        rtol=1e-3,
        atol=1e-4,
    )


def test_stops_after_max_subdivisions():
    # Define a cubature rule with fake high error so that cub will keep on subdividing

    class BadError(Rule):
        underlying = GaussLegendreQuad(10)

        def estimate(self, f, a, b, args=(), kwargs=dict()):
            return self.underlying.estimate(f, a, b, args, kwargs)

        def estimate_error(self, f, a, b, args=(), kwargs=dict()):
            return 1e6

    def f(x):
        return x

    a = np.array([0])
    b = np.array([1])
    rule = BadError()

    res = cub(
        f,
        a,
        b,
        rule,
        max_subdivisions=10,
    )

    assert res.subdivisions == 10
    assert not res.success
    assert res.status == "not_converged"


@pytest.mark.parametrize("quadrature", [
    NewtonCotesQuad,
    GaussLegendreQuad
])
def test_one_point_fixed_quad_impossible(quadrature):
    with pytest.raises(Exception):
        quadrature(1)


def test_estimate_with_base_classes_raise_error():
    def f(x):
        return x

    a = np.array([0])
    b = np.array([1])

    for base_class in [Rule(), FixedRule()]:
        with pytest.raises(Exception):
            base_class.estimate(f, a, b)


def test_genz_malik_1d_raises_error():
    with pytest.raises(Exception, match="only defined for ndim >= 2"):
        GenzMalikCub(1)


@pytest.mark.parametrize("problem", [
    (
        # 2D problem, 1D rule
        np.array([0, 0]),
        np.array([1, 1]),
        GaussKronrodQuad(21),
    ),
    (
        # 2D problem, 1D rule
        np.array([0, 0]),
        np.array([1, 1]),
        NestedRule(NewtonCotesQuad(10), NewtonCotesQuad(5)),
    ),
    (
        # 1D problem, 2D rule
        np.array([0]),
        np.array([1]),
        GenzMalikCub(2),
    ),
])
def test_incompatible_dimension_raises_error(problem):
    def f(x):
        return x

    a, b, rule = problem

    with pytest.raises(Exception, match="incompatible dimension"):
        cub(f, a, b, rule)


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
