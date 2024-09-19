import math
import scipy
import itertools

import pytest

import numpy as np  # TODO: remove

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
    return xp.reshape(2**(n+1)/(n+1), (-1, 1))


def basic_nd_integrand(x, n, xp):
    return xp.reshape(xp.sum(x, axis=-1), (-1, 1))**xp.reshape(n, (1, -1))


def basic_nd_integrand_exact(n, xp):
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


def f_gaussian(x, alphas):
    r"""
    .. math::

        f(\mathbf x) = \exp\left(-\sum^n_{i = 1} (\alpha_i x_i)^2 \right)
    """
    npoints, ndim = x.shape[0], x.shape[-1]
    alphas_reshaped = alphas[np.newaxis, :]
    x_reshaped = x.reshape(npoints, *([1]*(len(alphas.shape) - 1)), ndim)

    return np.exp(
        -np.sum((alphas_reshaped * x_reshaped)**2, axis=-1)
    )


def f_gaussian_exact(a, b, alphas):
    # Exact only when `a` and `b` are one of:
    #   (-oo, oo), or
    #   (0, oo), or
    #   (-oo, 0)
    # `alphas` can be arbitrary.

    ndim = len(a)
    double_infinite_count = 0
    semi_infinite_count = 0

    for i in range(len(a)):
        if np.isinf(a[i]) and np.isinf(b[i]):   # doubly-infinite
            double_infinite_count += 1
        elif np.isinf(a[i]) != np.isinf(b[i]):  # exclusive or, so semi-infinite
            semi_infinite_count += 1

    return (math.sqrt(np.pi) ** ndim) / (
        2**semi_infinite_count * np.prod(alphas, axis=-1)
    )


def f_gaussian_random_args(shape):
    np.random.seed(12)

    alphas = np.random.rand(*shape)

    # If alphas are very close to 0 this makes the problem very difficult due to large
    # values of ``f``.
    alphas *= 100

    return alphas


def f_modified_gaussian(x_arr, n):
    r"""
    .. math::

        f(x, y, z, w) = x^n \sqrt{y} \exp(-y-z^2-w^2)
    """
    x, y, z, w = x_arr[:, 0], x_arr[:, 1], x_arr[:, 2], x_arr[:, 3]
    res = (x ** n[:, np.newaxis]) * np.sqrt(y) * np.exp(-y-z**2-w**2)

    return res.T


def f_modified_gaussian_exact(a, b, n):
    # Exact only for the limits
    #   a = (0, 0, -oo, -oo)
    #   b = (1, oo, oo, oo)
    # but defined here as a function to match the format of the other integrands.
    return 1/(2 + 2*n) * math.pi ** (3/2)


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

        res = cubature(basic_nd_integrand, a, b, rule_str, args=(n, xp))

        xp_assert_close(
            res.estimate,
            basic_nd_integrand_exact(n, xp),
            rtol=1e-1,
            atol=0,
        )

    @skip_xp_backends(np_only=True,
                      reasons=['array-likes only supported for NumPy backend'])
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
            rule,
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

    def test_limits_other_way_around(self, xp):
        n = xp.arange(5)

        a = xp.asarray([2])
        b = xp.asarray([0])

        res = cubature(
            basic_1d_integrand,
            a,
            b,
            args=(n, xp),
        )

        np.testing.assert_allclose(
            res.estimate,
            -basic_1d_integrand_exact(n, xp),
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

    problems_scalar_output = [
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
            [10, 10, 10],
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

            [10, 50],
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
            rule,
            rtol,
            atol,
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
            rule,
            rtol,
            atol,
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

    @pytest.mark.parametrize("problem", [
        (
            # Function to integrate
            lambda x: x,

            # Exact value
            50,

            # Coordinates of `a`
            np.array([0]),

            # Coordinates of `b`
            np.array([10]),

            # Points by which to split up the initial region
            None,
        ),
        (
            lambda x: np.sin(x)/x,
            scipy.special.sici(1)[0] + scipy.special.sici(2)[0],
            np.array([-1]),
            np.array([2]),
            [np.array([0])],
        ),
        (
            lambda x: np.ones(len(x)),
            1,
            np.array([0, 0, 0]),
            np.array([1, 1, 1]),
            [
                np.array([0.5, 0.5, 0.5]),
            ],
        ),
        (
            lambda x: np.ones(len(x)),
            1,
            np.array([0, 0, 0]),
            np.array([1, 1, 1]),
            [
                np.array([0.25, 0.25, 0.25]),
                np.array([0.5, 0.5, 0.5]),
            ],
        ),
        (
            lambda x: np.ones(len(x)),
            1,
            np.array([0, 0, 0]),
            np.array([1, 1, 1]),
            [
                np.array([0.1, 0.25, 0.5]),
                np.array([0.25, 0.25, 0.25]),
                np.array([0.5, 0.5, 0.5]),
            ],
        )
    ])
    def test_break_points(self, problem, rule, rtol, atol):
        f, exact, a, b, points = problem
        ndim = a.size

        if rule == "genz-malik" and ndim < 2:
            pytest.skip("Genz-Malik cubature does not support 1D integrals")

        if rule == "genz-malik" and ndim >= 5:
            pytest.mark.slow("Gauss-Kronrod is slow in >= 5 dim")

        res = cubature(
            f,
            a,
            b,
            rule,
            rtol,
            atol,
            points=points,
        )

        np.testing.assert_allclose(
            res.estimate,
            exact,
            rtol=rtol,
            atol=atol,
            err_msg=f"estimate_error={res.error}, subdivisions={res.subdivisions}"
        )

        err_msg = (f"estimate_error={res.error}, "
                   f"subdivisions= {res.subdivisions}, "
                   f"true_error={np.abs(res.estimate - exact)}")
        assert res.status == "converged", err_msg

    @pytest.mark.parametrize("problem", [
        (
            # Function to integrate
            f_gaussian,

            # Exact solution
            f_gaussian_exact,

            # Arguments passed to f
            f_gaussian_random_args((1, 1)),

            # Limits, have to match the shape of the arguments
            np.array([-np.inf]),  # a
            np.array([np.inf]),   # b
        ),
        (
            f_gaussian,
            f_gaussian_exact,
            f_gaussian_random_args((1, 1),),
            np.array([0]),
            np.array([np.inf]),
        ),
        (
            f_gaussian,
            f_gaussian_exact,
            (f_gaussian_random_args((1, 1)),),
            np.array([-np.inf]),
            np.array([0]),
        ),
        (
            f_gaussian,
            f_gaussian_exact,
            (f_gaussian_random_args((1, 4)),),
            np.array([0, 0, -np.inf, -np.inf]),
            np.array([np.inf, np.inf, np.inf, np.inf]),
        ),
        (
            f_gaussian,
            f_gaussian_exact,
            (f_gaussian_random_args((1, 4)),),
            np.array([-np.inf, -np.inf, -np.inf, -np.inf]),
            np.array([0, 0, np.inf, np.inf]),
        ),

        # This particular problem can be slow
        pytest.param(
            (
                # f(x, y, z, w) = x^n * sqrt(y) * exp(-y-z**2-w**2) for n in [0,1,2,3]
                f_modified_gaussian,

                # This exact solution is for the below limits, not in general
                f_modified_gaussian_exact,

                (np.arange(4),),
                np.array([0, 0, -np.inf, -np.inf]),
                np.array([1, np.inf, np.inf, np.inf])
            ),

            marks=pytest.mark.slow,
        )
    ])
    def test_infinite_limits(self, problem, rule, rtol, atol):
        f, exact, args, a, b = problem

        ndim = a.size

        if rule == "genz-malik" and ndim < 2:
            pytest.skip("Genz-Malik cubature does not support 1D integrals")

        if rule == "genz-malik" and ndim >= 4:
            pytest.mark.slow("Gauss-Kronrod is slow in >= 5 dim")

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

        np.testing.assert_allclose(
            res.estimate,
            exact(a, b, *args),
            rtol=rtol,
            atol=atol,
            err_msg=f"error_estimate={res.error}, subdivisions={res.subdivisions}"
        )

    @pytest.mark.parametrize(
        "problem",
        [
            (
                # Integrate
                #   f(x_1, x_2, x_3) = x_1
                # with limits
                #   a = [0, 0, 0]
                #   b = [1, 1, 1]
                # f
                lambda x: x[:, 0],
                # exact solution
                0.5,
                # args
                (),
                # a & b
                np.array([0, 0]),
                np.array([1, 1]),
                # region, here just constants
                [
                    lambda x: (
                        np.zeros(x.shape[0]).reshape(-1, 1),
                        np.ones(x.shape[0]).reshape(-1, 1),
                    ),
                ],
            ),
            (
                # Integrate
                #   f(x_1, x_2, x_3) = 1
                # with limits
                #   a = [0, 0, 0]
                #   b = [1, 1, x_1 * x_2]
                lambda x: np.ones(x.shape[0]),
                0.25,
                (),
                np.array([0, 0]),
                np.array([1, 1]),
                [
                    lambda x: (
                        np.zeros(x.shape[0]).reshape(-1, 1),
                        (x[:, 0] * x[:, 1]).reshape(-1, 1),
                    ),
                ],
            ),
            (
                # Integrate
                #   f(x_1, ..., x_n) = 1
                # with limits
                #   a = [-1, -sqrt(1 - x_1**2)]
                #   b = [1, sqrt(1 - x_1**2)]
                lambda x: np.ones(x.shape[0]),
                np.pi,
                (),
                np.array([-1]),
                np.array([1]),
                [
                    lambda x: (
                        -np.sqrt(1 - x[:, 0] ** 2).reshape(-1, 1),
                        np.sqrt(1 - x[:, 0] ** 2).reshape(-1, 1),
                    ),
                ],
            ),
            (
                # Integrate
                #   f(x_1, x_2) = 1
                # with limits
                #   a = [-inf, -exp(-x_1**2)]
                #   b = [inf, exp(x_1**2)]
                lambda x: np.ones(x.shape[0]),
                2 * math.sqrt(np.pi),
                (),
                np.array([-np.inf]),
                np.array([np.inf]),
                [
                    lambda x: (
                        -np.exp(-x[:, 0] ** 2).reshape(-1, 1),
                        np.exp(-x[:, 0] ** 2).reshape(-1, 1),
                    ),
                ],
            ),
            (
                # Integrate
                #   f(x_1, x_2) = n
                # with limits
                #   a = [-inf, -exp(-x_1**2)]
                #   b = [inf, exp(x_1**2)]
                # for n = range(5)
                lambda x, n: np.tile(n, (x.shape[0], 1)),
                np.arange(5) * 2 * math.sqrt(np.pi),
                (np.arange(5),),
                np.array([-np.inf]),
                np.array([np.inf]),
                [
                    lambda x: (
                        -np.exp(-x[:, 0] ** 2).reshape(-1, 1),
                        np.exp(-x[:, 0] ** 2).reshape(-1, 1),
                    ),
                ],
            ),
            (
                # Integrate
                #   f(x, y, z) = n
                # with limits
                #   a = [-1, -sqrt(1-x^2), -sqrt(1-x^2-y^2)],
                #   b = [1, sqrt(1-x^2), sqrt(1-x^2-y^2)]
                # for n = range(5)
                # i.e. integrate a unit ball with densities 0, 1, 2, 3, 4
                lambda x, n: np.tile(n, (x.shape[0], 1)),
                4 / 3 * math.pi * np.arange(5),
                (np.arange(5),),
                np.array([-1]),
                np.array([1]),
                [
                    lambda x: (
                        -np.sqrt(1 - x[:, 0] ** 2).reshape(-1, 1),
                        np.sqrt(1 - x[:, 0] ** 2).reshape(-1, 1),
                    ),
                    lambda x: (
                        -np.sqrt(1 - x[:, 0] ** 2 - x[:, 1] ** 2).reshape(-1, 1),
                        np.sqrt(1 - x[:, 0] ** 2 - x[:, 1] ** 2).reshape(-1, 1),
                    ),
                ],
            ),
            (
                # Integrate
                #   f(x, y, z, w) = n
                # with limits
                #   a = [-1, -1, -sqrt(1-x^2), -sqrt(1-x^2)],
                #   b = [1, 1, sqrt(1-x^2), sqrt(1-x^2)]
                # for n = range(5)
                lambda x, n: np.tile(n, (x.shape[0], 1)),
                32 / 3 * np.arange(5),
                (np.arange(5),),
                np.array([-1, -1]),
                np.array([1, 1]),
                [
                    lambda x: (
                        -np.repeat(np.sqrt(1 - x[:, 0] ** 2), 2).reshape(-1, 2),
                        np.repeat(np.sqrt(1 - x[:, 0] ** 2), 2).reshape(-1, 2),
                    ),
                ],
            ),
            (
                # Integrate
                #   f(x, y, z, w) = n
                # with limits
                #   a = [0, 0, 0, 0],
                #   b = [1, x, x, x]
                # for n = range(5)
                lambda x, n: np.tile(n, (x.shape[0], 1)),
                1 / 4 * np.arange(5),
                (np.arange(5),),
                np.array([0]),
                np.array([1]),
                [
                    lambda x: (
                        np.zeros((x.shape[0], 3)),
                        np.repeat(x[:, 0], 3).reshape(-1, 3),
                    ),
                ],
            ),
            (
                # Integrate
                #   f(x, y, z, w, p) = n
                # with limits
                #   a = [0, 0, 0, 0, 0],
                #   b = [1, x, x, xy, xy]
                # for n = range(5)
                lambda x, n: np.tile(n, (x.shape[0], 1)),
                1 / 21 * np.arange(5),
                (np.arange(5),),
                np.array([0]),
                np.array([1]),
                [
                    lambda x: (
                        np.zeros((x.shape[0], 2)),
                        np.repeat(x[:, 0], 2).reshape(-1, 2),
                    ),
                    lambda x: (
                        np.zeros((x.shape[0], 2)),
                        np.repeat(x[:, 0] * x[:, 1], 2).reshape(-1, 2),
                    ),
                ],
            ),
            (
                # Integrate
                #   f(x, y) = n
                # with limits
                #   a = [0, -sqrt(1-x^2)],
                #   b = [1, sqrt(1-x^2)]
                # for n = range(5)
                lambda x, n: np.tile(n, (x.shape[0], 1)),
                math.pi * np.arange(5),
                (np.arange(5),),
                np.array([]),  # Here the initial limits are empty
                np.array([]),
                [
                    lambda x: (
                        -np.ones((x.shape[0], 1)),
                        np.ones((x.shape[0], 1)),
                    ),
                    lambda x: (
                        -np.sqrt(1 - x[:, 0] ** 2).reshape(-1, 1),
                        np.sqrt(1 - x[:, 0] ** 2).reshape(-1, 1),
                    ),
                ],
            ),
            pytest.param(
                (
                    # Integrate
                    #   f(x, y, z, w) = n
                    # with limits
                    #   a = [-1, -sqrt(1-x^2), -sqrt(1-x^2-y^2), -sqrt(1-x^2-y^2-z^2)],
                    #   b = [1, sqrt(1-x^2), sqrt(1-x^2-y^2), sqrt(1-x^2-y^2-z^2)]
                    # for n = range(5)
                    # i.e. integration over a hypersphere
                    lambda x, n: np.tile(n, (x.shape[0], 1)),
                    math.pi**2 / 2 * np.arange(5),
                    (np.arange(5),),
                    np.array([-1]),
                    np.array([1]),
                    [
                        lambda x: (
                            -np.sqrt(1 - x[:, 0] ** 2).reshape(-1, 1),
                            np.sqrt(1 - x[:, 0] ** 2).reshape(-1, 1),
                        ),
                        lambda x: (
                            -np.sqrt(1 - x[:, 0] ** 2 - x[:, 1] ** 2).reshape(-1, 1),
                            np.sqrt(1 - x[:, 0] ** 2 - x[:, 1] ** 2).reshape(-1, 1),
                        ),
                        lambda x: (
                            -np.sqrt(
                                1 - x[:, 0] ** 2 - x[:, 1] ** 2 - x[:, 2] ** 2
                            ).reshape(-1, 1),
                            np.sqrt(
                                1 - x[:, 0] ** 2 - x[:, 1] ** 2 - x[:, 2] ** 2
                            ).reshape(-1, 1),
                        ),
                    ],
                ),
                marks=pytest.mark.slow,
            ),
        ],
    )
    def test_cub_func_limits(self, problem, rule, rtol, atol):
        f, exact, args, a_outer, b_outer, region = problem

        res = cubature(
            f,
            a_outer,
            b_outer,
            rule,
            rtol,
            atol,
            args=args,
            region=region,
        )

        assert res.status == "converged"

        np.testing.assert_allclose(
            res.estimate,
            exact,
            rtol=rtol,
            atol=atol,
            err_msg=f"error_estimate={res.error}, subdivisions={res.subdivisions}"
        )


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
            rule,
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
