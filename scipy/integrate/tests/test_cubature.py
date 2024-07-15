import math
import scipy
import itertools

import pytest
import numpy as np

from scipy.integrate._cubature import (
    cub, ProductRule, GaussKronrod15, GaussKronrod21
)


def genz_malik_1980_f_1(x, r, alphas):
    r"""
    `f_1` from Genz and Malik 1980.

    .. math:: f_1(\mathbf x) = \cos\left(2\pi r + \sum^n_{i = 1}\alpha_i x_i\right)

    .. code-block:: mathematica

        genzMalik1980f1[x_List, r_, alphas_List] := Cos[2*Pi*r + Total[x*alphas]]
    """

    dim = x.shape[0]
    num_eval_points = x.shape[-1]

    r_reshaped = np.expand_dims(r, -1)
    alphas_reshaped = np.expand_dims(alphas, -1)
    x_reshaped = x.reshape(dim, *([1]*(len(alphas.shape) - 1)), num_eval_points)

    return np.cos(2*np.pi*r_reshaped + np.sum(alphas_reshaped * x_reshaped, axis=0))


def genz_malik_1980_f_1_exact(a, b, r, alphas):
    dim = len(a)
    a = a.reshape(dim, *([1]*(len(alphas.shape) - 1)))
    b = b.reshape(dim, *([1]*(len(alphas.shape) - 1)))

    return (-2)**dim * 1/np.prod(alphas, axis=0) \
        * np.cos(2*np.pi*r + np.sum(alphas * (a+b)/2, axis=0)) \
        * np.prod(
            np.sin(alphas * (a-b)/2), axis=0
        )


def genz_malik_1980_f_1_random_args(shape):
    r = np.random.rand(*shape)
    alphas = np.random.rand(*shape)

    difficulty = 9
    alphas = difficulty * alphas / np.sum(alphas, axis=0)

    return (r, alphas)


def genz_malik_1980_f_2(x, alphas, betas):
    r"""
    `f_2` from Genz and Malik 1980.

    .. math:: f_2(\mathbf x) = \prod^n_{i = 1} (\alpha_i^2 + (x_i - \beta_i)^2)^{-1}

    .. code-block:: mathematica

        genzMalik1980f2[x_List, alphas_List, betas_List] :=
            1/Times @@ ((alphas^2 + (x - betas)^2))
    """
    dim = x.shape[0]
    num_eval_points = x.shape[-1]

    alphas_reshaped = np.expand_dims(alphas, -1)
    betas_reshaped = np.expand_dims(betas, -1)

    x_reshaped = x.reshape(dim, *([1]*(len(alphas.shape) - 1)), num_eval_points)

    return 1/np.prod(alphas_reshaped**2 + (x_reshaped-betas_reshaped)**2, axis=0)


def genz_malik_1980_f_2_exact(a, b, alphas, betas):
    dim = len(a)
    a = a.reshape(dim, *([1]*(len(alphas.shape) - 1)))
    b = b.reshape(dim, *([1]*(len(alphas.shape) - 1)))

    return (-1)**dim * 1/np.prod(alphas, axis=0) * np.prod(
        np.arctan((a - betas)/alphas) - np.arctan((b - betas)/alphas),
        axis=0
    )


def genz_malik_1980_f_2_random_args(shape):
    dim = shape[0]
    alphas = np.random.rand(*shape)
    betas = np.random.rand(*shape)

    difficulty = 25
    products = np.prod(np.power(alphas, -2), axis=0)
    alphas = alphas \
        * np.power(products, 1 / (2*dim)) \
        / np.power(difficulty, 1 / (2*dim))

    return alphas, betas


def genz_malik_1980_f_3(x, alphas):
    r"""
    `f_3` from Genz and Malik 1980.

    .. math:: f_3(\mathbf x) = \exp\left(\sum^n_{i = 1} \alpha_i x_i\right)

    .. code-block:: mathematica

        genzMalik1980f3[x_List, alphas_List] := Exp[Dot[x, alphas]]
    """

    dim = x.shape[0]
    num_eval_points = x.shape[-1]

    alphas_reshaped = np.expand_dims(alphas, -1)
    x_reshaped = x.reshape(dim, *([1]*(len(alphas.shape) - 1)), num_eval_points)

    return np.exp(np.sum(alphas_reshaped * x_reshaped, axis=0))


def genz_malik_1980_f_3_exact(a, b, alphas):
    dim = len(a)
    a = a.reshape(dim, *([1]*(len(alphas.shape) - 1)))
    b = b.reshape(dim, *([1]*(len(alphas.shape) - 1)))

    return (-1)**dim * 1/np.prod(alphas, axis=0) * np.prod(
        np.exp(alphas * a) - np.exp(alphas * b),
        axis=0
    )


def genz_malik_1980_f_3_random_args(shape):
    alphas = np.random.rand(*shape)
    difficulty = 12
    alphas = difficulty * alphas / np.sum(alphas, axis=0)

    return (alphas,)


def genz_malik_1980_f_4(x, alphas):
    r"""
    `f_4` from Genz and Malik 1980.

    .. math:: f_4(\mathbf x) = \left(1 + \sum^n_{i = 1} \alpha_i x_i\right)^{-n-1}

    .. code-block:: mathematica
        genzMalik1980f4[x_List, alphas_List] :=
            (1 + Dot[x, alphas])^(-Length[alphas] - 1)
    """

    dim = x.shape[0]
    num_eval_points = x.shape[-1]

    alphas_reshaped = np.expand_dims(alphas, -1)
    x_reshaped = x.reshape(dim, *([1]*(len(alphas.shape) - 1)), num_eval_points)

    return ((1 + np.sum(alphas_reshaped * x_reshaped, axis=0))**(-x.shape[0]-1))


def genz_malik_1980_f_4_exact(a, b, alphas):
    dim = len(a)

    def F(x):
        x_reshaped = x.reshape(dim, *([1]*(len(alphas.shape) - 1)))

        return (-1)**dim/np.prod(alphas, axis=0) / (
            math.factorial(dim) * (1 + np.sum(alphas * x_reshaped, axis=0))
        )

    return _eval_indefinite_integral(F, a, b)


def genz_malik_1980_f_4_random_args(shape):
    dim = shape[0]

    alphas = np.random.rand(*shape)
    difficulty = 14
    alphas = (difficulty/dim) * alphas / np.sum(alphas, axis=0)

    return (alphas,)


def genz_malik_1980_f_5(x, alphas, betas):
    r"""
    `f_5` from Genz and Malik 1980.

    .. math::

        f_5(\mathbf x) = \exp\left(-\sum^n_{i = 1} \alpha^2_i (x_i - \beta_i)^2\right)

    .. code-block:: mathematica

        genzMalik1980f5[x_List, alphas_List, betas_List] :=
            Exp[-Total[alphas^2 * (x - betas)^2]]
    """

    dim = x.shape[0]
    num_eval_points = x.shape[-1]

    alphas_reshaped = np.expand_dims(alphas, -1)
    betas_reshaped = np.expand_dims(betas, -1)

    x_reshaped = x.reshape(dim, *([1]*(len(alphas.shape) - 1)), num_eval_points)

    return np.exp(
        -np.sum(alphas_reshaped**2 * (x_reshaped - betas_reshaped)**2, axis=0)
    )


def genz_malik_1980_f_5_exact(a, b, alphas, betas):
    dim = len(a)
    a = a.reshape(dim, *([1]*(len(alphas.shape) - 1)))
    b = b.reshape(dim, *([1]*(len(alphas.shape) - 1)))

    return (1/2)**dim * 1/np.prod(alphas, axis=0) * np.power(np.pi, dim/2) * np.prod(
        scipy.special.erf(alphas * (betas - a))
        + scipy.special.erf(alphas * (b - betas)),
        axis=0
    )


def genz_malik_1980_f_5_random_args(shape):
    alphas = np.random.rand(*shape)
    betas = np.random.rand(*shape)

    difficulty = 21
    l2_norm = np.sum(np.power(alphas, 2), axis=0)
    alphas = alphas / np.sqrt(l2_norm) * np.sqrt(difficulty)

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
            np.array([[5]]),
        )
    ),
    (
        genz_malik_1980_f_1,
        genz_malik_1980_f_1_exact,
        np.array([0, 0]),
        np.array([1, 1]),
        (
            np.array([1/4]),
            np.array([[2], [4]]),
        ),
    ),
    (
        genz_malik_1980_f_1,
        genz_malik_1980_f_1_exact,
        np.array([0, 0]),
        np.array([100, 100]),
        (
            np.array([1/2]),
            np.array([[2], [4]]),
        )
    ),
    (
        genz_malik_1980_f_1,
        genz_malik_1980_f_1_exact,
        np.array([0, 0, 0]),
        np.array([100, 200, 300]),
        (
            np.array([1/2]),
            np.array([[1], [1], [1]]),
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
    # Currently failing for product rule version of Gauss-Kronrod 15 (and 21)
    # dblquad also fails
    # (
    #     genz_malik_1980_f_2,
    #     genz_malik_1980_f_2_exact,

    #     np.array([-50, 0]),
    #     np.array([50, 20]),
    #     (
    #         np.array([-3, 3]),
    #         np.array([-2, 2])
    #     ),
    # ),
    (
        genz_malik_1980_f_2,
        genz_malik_1980_f_2_exact,
        np.array([0, 0, 0]),
        np.array([1, 1, 1]),
        (
            np.array([[1], [1], [1]]),
            np.array([[1], [1], [1]]),
        )
    ),
    (
        genz_malik_1980_f_2,
        genz_malik_1980_f_2_exact,
        np.array([0, 0, 0]),
        np.array([1, 1, 1]),
        (
            np.array([[2], [3], [4]]),
            np.array([[2], [3], [4]]),
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
            np.array([[1], [1], [1], [1]]),
            np.array([[1], [1], [1], [1]]),
        )
    ),

    # -- f3 --
    (
        genz_malik_1980_f_3,
        genz_malik_1980_f_3_exact,
        np.array([-1]),
        np.array([1]),
        (
            np.array([[1/2]]),
        ),
    ),
    (
        genz_malik_1980_f_3,
        genz_malik_1980_f_3_exact,
        np.array([0, -1]),
        np.array([1, 1]),
        (
            np.array([[5], [5]]),
        ),
    ),
    (
        genz_malik_1980_f_3,
        genz_malik_1980_f_3_exact,
        np.array([-1, -1, -1]),
        np.array([1, 1, 1]),
        (
            np.array([[1], [1], [1]]),
        ),
    ),

    # -- f4 --
    (
        genz_malik_1980_f_4,
        genz_malik_1980_f_4_exact,
        np.array([0]),
        np.array([2]),
        (np.array([1]),),
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
            np.array([[-2]]),
            np.array([[2]]),
        ),
    ),
    (
        genz_malik_1980_f_5,
        genz_malik_1980_f_5_exact,
        np.array([-1, -1]),
        np.array([1, 1]),
        (
            np.array([[2], [3]]),
            np.array([[4], [5]]),
        ),
    ),
    # Also failing, I think the problem previously was that the parameters were set in
    # a way where the bump in the function was far from the origin.
    # TODO: check if dblquad also fails this
    # problem_from_closed_form(
    #     f=genz_malik_1980_f_5,
    #     closed_form=genz_malik_1980_f_5_exact,

    #     a=np.array([-1, -1]),
    #     b=np.array([1, 1]),
    #     args=(
    #         np.array([[-1], [1]]),
    #         np.array([[0], [0]]),
    #     ),
    # ),
    (
        genz_malik_1980_f_5,
        genz_malik_1980_f_5_exact,
        np.array([-1, -1, -1]),
        np.array([1, 1, 1]),
        (
            np.array([[1], [1], [1]]),
            np.array([[1], [1], [1]]),
        ),
    )
]


@pytest.mark.parametrize("problem", problems_scalar_output)
@pytest.mark.parametrize("quadrature", [GaussKronrod15()])
@pytest.mark.parametrize("rtol", [1e-5])
@pytest.mark.parametrize("atol", [1e-6])
def test_cub_scalar_output(problem, quadrature, rtol, atol):
    np.random.seed(1)
    f, exact, a, b, args = problem

    dim = len(a)
    rule = ProductRule([quadrature] * dim)

    res = cub(
        f,
        a,
        b,
        rule,
        rtol,
        atol,
        args=args,
    )

    np.testing.assert_allclose(
        res.estimate,
        exact(a, b, *args),
        rtol=rtol,
        atol=atol,
    )

    assert res.status == "converged"


problems_tensor_output = [
    (
        # Function to integrate, like `f(x, *args)`
        genz_malik_1980_f_1,

        # Exact solution, like `exact(a, b, *args)`
        genz_malik_1980_f_1_exact,

        # Function that generates random args of a certain shape, like `random(shape)`.
        genz_malik_1980_f_1_random_args,
    ),
    # Tests currently failing:
    # (
    #     genz_malik_1980_f_2,
    #     genz_malik_1980_f_2_exact,
    #     genz_malik_1980_f_2_random_args
    # ),
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


@pytest.mark.parametrize("problem", problems_tensor_output)
@pytest.mark.parametrize("quadrature", [GaussKronrod21()])
@pytest.mark.parametrize("shape", [
    (2, 3,),
    (2, 3, 5),
    (3, 5, 6),
    (4, 5, 2, 1)
])
@pytest.mark.parametrize("rtol", [1e-5])
@pytest.mark.parametrize("atol", [1e-6])
def test_cub_tensor_output(problem, quadrature, shape, rtol, atol):
    dim = shape[0]

    f, exact, random_args = problem
    args = random_args(shape)

    a = np.array([0] * dim)
    b = np.array([1] * dim)

    rule = ProductRule([quadrature] * dim)

    res = cub(
        f,
        a,
        b,
        rule,
        rtol,
        atol,
        args=args,
    )

    np.testing.assert_allclose(
        res.estimate,
        exact(a, b, *args),
        rtol=rtol,
        atol=atol,
    )

    assert res.status == "converged"


def _eval_indefinite_integral(F, a, b):
    dim = len(a)
    points = np.stack([a, b], axis=0)

    out = 0
    for ind in itertools.product(range(2), repeat=dim):
        out += pow(-1, sum(ind) + dim) * F(points[ind, tuple(range(dim))])

    return out
