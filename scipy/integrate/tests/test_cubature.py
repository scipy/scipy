import math

import pytest
import scipy
import numpy as np

from scipy.integrate._cubature import (
    cub, ProductRule, GaussKronrod15, GaussKronrod21, _cartesian_product
)


class Problem:
    """
    Represents a general problem to be solved by quadrature algorithm.
    """
    def __init__(self, dim, f, a, b, args, exact):
        self.dim = dim
        self.f = f
        self.a = a
        self.b = b
        self.args = args
        self.exact = exact


def assert_integral_estimate_close(problem, res, rtol, atol):
    np.testing.assert_allclose(
        res.estimate,
        problem.exact,
        rtol=rtol,
        atol=atol,
        verbose=True,
        err_msg=f"""failed to find approx. integral of {problem.f.__name__} with \
rtol={rtol}, atol={atol} and points \
a={problem.a}
b={problem.b}""",
    )


def problem_from_closed_form(dim, f, closed_form, a, b, args):
    """
    Creates a new Problem given a function `f` and the corresponding `closed_form`.
    """
    return Problem(
        dim=dim,
        f=f,
        a=a,
        b=b,
        args=args,
        exact=closed_form(dim, a, b, *args)
    )


# TODO: I think this actually calculates the transpose of f
def genz_malik_1980_f_1(x, r, alphas):
    r"""
    `f_1` from Genz and Malik 1980.

    .. math:: f_1(\mathbf x) = \cos\left(2\pi r + \sum^n_{i = 1}\alpha_i x_i\right)

    .. code-block:: mathematica

        genzMalik1980f1[x_List, r_, alphas_List] := Cos[2*Pi*r + Total[x*alphas]]
    """

    return np.cos(np.expand_dims(2*np.pi*r.T, -1) + (alphas.T @ x))


# TODO: I think this actually calculates the transpose of the exact soln.
def genz_malik_1980_f_1_exact(dim, a, b, r, alphas):
    return (-2)**dim * 1/np.prod(alphas, axis=0).T \
        * np.cos(2*np.pi*r.T + alphas.T @ (a+b)/2) \
        * np.prod(np.sin(alphas.T * (a-b)/2), axis=-1)


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

    # Add extra dimension to parameters for eval_points
    alphas_reshaped = np.expand_dims(alphas, -1)
    betas_reshaped = np.expand_dims(betas, -1)

    x_reshaped = x.reshape(dim, *([1]*(len(alphas.shape) - 1)), num_eval_points)

    return 1/np.prod(alphas_reshaped**2 + (x_reshaped-betas_reshaped)**2, axis=0)


def genz_malik_1980_f_2_exact(dim, a, b, alphas, betas):
    # Expand a and b to fit the dimension of the parameters
    a = a.reshape(dim, *([1]*(len(alphas.shape) - 1)))
    b = b.reshape(dim, *([1]*(len(alphas.shape) - 1)))

    return (-1)**dim * 1/np.prod(alphas, axis=0) * np.prod(
        np.arctan((a - betas)/alphas) - np.arctan((b - betas)/alphas),
        axis=0
    )


def genz_malik_1980_f_3(x, alphas):
    r"""
    `f_3` from Genz and Malik 1980.

    .. math:: f_3(\mathbf x) = \exp\left(\sum^n_{i = 1} \alpha_i x_i\right)

    .. code-block:: mathematica

        genzMalik1980f3[x_List, alphas_List] := Exp[Dot[x, alphas]]
    """

    dim = x.shape[0]
    num_eval_points = x.shape[-1]

    # Add extra dimension to parameters for eval_points
    alphas_reshaped = np.expand_dims(alphas, -1)

    x_reshaped = x.reshape(dim, *([1]*(len(alphas.shape) - 1)), num_eval_points)

    return np.exp(np.sum(alphas_reshaped * x_reshaped, axis=0))


def genz_malik_1980_f_3_exact(dim, a, b, alphas):
    # Expand a and b to fit the dimension of the parameters
    a = a.reshape(dim, *([1]*(len(alphas.shape) - 1)))
    b = b.reshape(dim, *([1]*(len(alphas.shape) - 1)))

    return (-1)**dim * 1/np.prod(alphas, axis=0) * np.prod(
        np.exp(alphas * a) - np.exp(alphas * b),
        axis=0
    )


def genz_malik_1980_f_4(x, alphas):
    dim = x.shape[0]
    num_eval_points = x.shape[-1]
    alphas_reshaped = np.expand_dims(alphas, -1)
    x_reshaped = x.reshape(dim, *([1]*(len(alphas.shape) - 1)), num_eval_points)

    return ((1 + np.sum(alphas_reshaped * x_reshaped, axis=0))**(-x.shape[0]-1))


def genz_malik_1980_f_4_exact(dim, a, b, alphas):
    alphas_reshaped = np.expand_dims(alphas, -1)

    def F(x):
        x_reshaped = x.reshape(dim, *([1]*(len(alphas.shape) - 1)), 1)

        return 1/np.prod(alphas_reshaped, axis=0) / (
            math.factorial(dim) * (1 + np.sum(alphas_reshaped * x_reshaped, axis=0))
        )

    points = _cartesian_product(np.array([a, b]).T).T
    res = 0

    for i, point in enumerate(points):
        zeroes = bin(i).count('1')
        sign = (-1)**zeroes
        res += sign * F(point.reshape(-1, 1))

    return res


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

    # Add extra dimension to parameters for eval_points
    alphas_reshaped = np.expand_dims(alphas, -1)
    betas_reshaped = np.expand_dims(alphas, -1)

    x_reshaped = x.reshape(dim, *([1]*(len(alphas.shape) - 1)), num_eval_points)

    return np.exp(
        -np.sum(alphas_reshaped**2 * (x_reshaped - betas_reshaped)**2, axis=0)
    )


# TODO: not 100% sure this is the correct closed form
def genz_malik_1980_f_5_exact(dim, a, b, alphas, betas):
    # Expand a and b to fit the dimension of the parameters
    a = a.reshape(dim, *([1]*(len(alphas.shape) - 1)))
    b = b.reshape(dim, *([1]*(len(alphas.shape) - 1)))

    return (1/2)**dim * 1/np.prod(alphas, axis=0) * np.power(np.pi, dim/2) * np.prod(
        scipy.special.erf(alphas * (betas - a))
        + scipy.special.erf(alphas * (b - betas)),
        axis=0
    )


problems_scalar_output = [
    # Problems based on test integrals in Genz and Malik 1980:
    # -- f1 --

    problem_from_closed_form(
        dim=1,
        f=genz_malik_1980_f_1,
        closed_form=genz_malik_1980_f_1_exact,

        a=np.array([0]),
        b=np.array([10]),
        args=(
            np.array([1/4]),
            np.array([[5]]),
        ),
    ),
    problem_from_closed_form(
        dim=2,
        f=genz_malik_1980_f_1,
        closed_form=genz_malik_1980_f_1_exact,

        a=np.array([0, 0]),
        b=np.array([1, 1]),
        args=(
            np.array([1/4]),
            np.array([[2], [4]]),
        ),
    ),
    problem_from_closed_form(
        dim=2,
        f=genz_malik_1980_f_1,
        closed_form=genz_malik_1980_f_1_exact,

        a=np.array([0, 0]),
        b=np.array([100, 100]),
        args=(
            np.array([1/2]),
            np.array([[2], [4]]),
        )
    ),
    problem_from_closed_form(
        dim=3,
        f=genz_malik_1980_f_1,
        closed_form=genz_malik_1980_f_1_exact,

        a=np.array([0, 0, 0]),
        b=np.array([100, 200, 300]),
        args=(
            np.array([1/2]),
            np.array([[1], [1], [1]]),
        )
    ),

    # -- f2 --

    problem_from_closed_form(
        dim=1,
        f=genz_malik_1980_f_2,
        closed_form=genz_malik_1980_f_2_exact,

        a=np.array([-1]),
        b=np.array([1]),
        args=(
            np.array([5]),
            np.array([4]),
        )
    ),
    # Currently failing for product rule version of Gauss-Kronrod 15 (and 21)
    # dblquad also fails
    # problem_with_analytic_sol(
    #     dim=2,
    #     f=genz_malik_1980_f_2,
    #     exact_func=genz_malik_1980_f_2_exact,
    #     a=np.array([-50, 0]),
    #     b=np.array([50, 20]),
    #     alphas=np.array([-3, 3]),
    #     betas=np.array([-2, 2])
    # ),
    problem_from_closed_form(
        dim=3,
        f=genz_malik_1980_f_2,
        closed_form=genz_malik_1980_f_2_exact,

        a=np.array([0, 0, 0]),
        b=np.array([1, 1, 1]),
        args=(
            np.array([[1], [1], [1]]),
            np.array([[1], [1], [1]]),
        )
    ),
    problem_from_closed_form(
        dim=3,
        f=genz_malik_1980_f_2,
        closed_form=genz_malik_1980_f_2_exact,

        a=np.array([0, 0, 0]),
        b=np.array([1, 1, 1]),
        args=(
            np.array([[2], [3], [4]]),
            np.array([[2], [3], [4]]),
        )
    ),
    problem_from_closed_form(
        dim=3,
        f=genz_malik_1980_f_2,
        closed_form=genz_malik_1980_f_2_exact,

        a=np.array([-1, -1, -1]),
        b=np.array([1, 1, 1]),
        args=(
            np.array([1, 1, 1]),
            np.array([2, 2, 2]),
        )
    ),
    problem_from_closed_form(
        dim=4,
        f=genz_malik_1980_f_2,
        closed_form=genz_malik_1980_f_2_exact,

        a=np.array([-1, -1, -1, -1]),
        b=np.array([1, 1, 1, 1]),
        args=(
            np.array([[1], [1], [1], [1]]),
            np.array([[1], [1], [1], [1]]),
        )
    ),

    # # -- f_3 --

    problem_from_closed_form(
        dim=1,
        f=genz_malik_1980_f_3,
        closed_form=genz_malik_1980_f_3_exact,

        a=np.array([-1]),
        b=np.array([1]),
        args=(
            np.array([[1/2]]),
        ),
    ),
    problem_from_closed_form(
        dim=2,
        f=genz_malik_1980_f_3,
        closed_form=genz_malik_1980_f_3_exact,

        a=np.array([0, -1]),
        b=np.array([1, 1]),
        args=(
            np.array([[5], [5]]),
        ),
    ),
    problem_from_closed_form(
        dim=3,
        f=genz_malik_1980_f_3,
        closed_form=genz_malik_1980_f_3_exact,

        a=np.array([-1, -1, -1]),
        b=np.array([1, 1, 1]),
        args=(
            np.array([[1], [1], [1]]),
        ),
    ),

    # # -- f_4 --

    problem_from_closed_form(
        dim=1,
        f=genz_malik_1980_f_4,
        closed_form=genz_malik_1980_f_4_exact,
        a=np.array([0]),
        b=np.array([2]),
        args=(np.array([1]),),
    ),
    problem_from_closed_form(
        dim=2,
        f=genz_malik_1980_f_4,
        closed_form=genz_malik_1980_f_4_exact,
        a=np.array([0, 0]),
        b=np.array([2, 1]),
        args=(np.array([1, 1]),),
    ),
    problem_from_closed_form(
        dim=3,
        f=genz_malik_1980_f_4,
        closed_form=genz_malik_1980_f_4_exact,
        a=np.array([0, 0, 0]),
        b=np.array([1, 1, 1]),
        args=(np.array([1, 1, 1]),),
    ),

    # -- f_5 --

    problem_from_closed_form(
        dim=1,
        f=genz_malik_1980_f_5,
        closed_form=genz_malik_1980_f_5_exact,

        a=np.array([-1]),
        b=np.array([1]),
        args=(
            np.array([[-2]]),
            np.array([[2]]),
        ),
    ),
    problem_from_closed_form(
        dim=2,
        f=genz_malik_1980_f_5,
        closed_form=genz_malik_1980_f_5_exact,

        a=np.array([-1, -1]),
        b=np.array([1, 1]),
        args=(
            np.array([[2], [3]]),
            np.array([[4], [5]]),
        ),
    ),
    # Also failing, I think the problem previously was that the parameters were set in
    # a way where the bump in the function was far from the origin.
    # TODO: check if dblquad also fails this
    # problem_from_closed_form(
    #     dim=2,
    #     f=genz_malik_1980_f_5,
    #     closed_form=genz_malik_1980_f_5_exact,

    #     a=np.array([-1, -1]),
    #     b=np.array([1, 1]),
    #     args=(
    #         np.array([[-1], [1]]),
    #         np.array([[0], [0]]),
    #     ),
    # ),
    problem_from_closed_form(
        dim=3,
        f=genz_malik_1980_f_5,
        closed_form=genz_malik_1980_f_5_exact,

        a=np.array([-1, -1, -1]),
        b=np.array([1, 1, 1]),
        args=(
            np.array([[1], [1], [1]]),
            np.array([[1], [1], [1]]),
        ),
    )
]


@pytest.mark.parametrize("problem", problems_scalar_output)
@pytest.mark.parametrize("quadrature", [GaussKronrod15()])
@pytest.mark.parametrize("rtol", [1e-5])
@pytest.mark.parametrize("atol", [1e-8])
def test_cub_scalar_output(problem, quadrature, rtol, atol):
    rule = ProductRule([quadrature] * problem.dim)

    res = cub(
        problem.f,
        problem.a,
        problem.b,
        rule,
        rtol,
        atol,
        args=problem.args,
    )

    assert_integral_estimate_close(problem, res, rtol, atol)

    assert res.status == "converged"


@pytest.mark.skip()
@pytest.mark.parametrize("quadrature", [GaussKronrod15()])
@pytest.mark.parametrize("rtol", [1e-5])
@pytest.mark.parametrize("atol", [1e-8])
@pytest.mark.parametrize("dim", [1, 2, 3])
@pytest.mark.parametrize("args_shape", [(2,), (2, 3), (2, 3, 5), (2, 3, 5, 7)])
def test_genz_malik_1980_f_1_arbitrary_shape(quadrature, rtol, atol, dim, args_shape):
    np.random.seed(1)
    difficulty = 9
    rule = ProductRule([quadrature] * dim)

    a = np.array([0]*dim)
    b = np.array([1]*dim)
    alphas = np.random.rand(dim, *args_shape)
    r = np.random.rand(*args_shape)

    # Adjust the parameters as in the paper
    alphas = difficulty * alphas / np.sum(alphas, axis=0)

    assert np.allclose(np.sum(alphas, axis=0), difficulty)

    problem = problem_from_closed_form(
        dim=dim,
        f=genz_malik_1980_f_1,
        closed_form=genz_malik_1980_f_1_exact,
        a=a,
        b=b,
        args=(
            r,
            alphas,
        ),
    )

    res = cub(
        problem.f,
        problem.a,
        problem.b,
        rule,
        rtol,
        atol,
        args=problem.args,
    )

    assert_integral_estimate_close(problem, res, rtol, atol)


@pytest.mark.skip()
@pytest.mark.parametrize("quadrature", [GaussKronrod21()])
@pytest.mark.parametrize("rtol", [1e-5])
@pytest.mark.parametrize("atol", [1e-8])
@pytest.mark.parametrize("dim", [1, 2, 3])
@pytest.mark.parametrize("args_shape", [(2,), (2, 3), (2, 3, 5), (2, 3, 5, 7)])
def test_genz_malik_1980_f_2_arbitrary_shape(quadrature, rtol, atol, dim, args_shape):
    np.random.seed(1)
    difficulty = 25
    rule = ProductRule([quadrature] * dim)

    a = np.array([0]*dim)
    b = np.array([1]*dim)
    alphas = np.random.rand(dim, *args_shape)
    betas = np.random.rand(dim, *args_shape)

    # Adjust the parameters as in the paper

    products = np.prod(np.power(alphas, -2), axis=0)
    alphas = alphas \
        * np.power(products, 1 / (2*dim)) \
        / np.power(difficulty, 1 / (2*dim))

    assert np.allclose(np.prod(np.power(alphas, -2), axis=0), difficulty)

    problem = problem_from_closed_form(
        dim=dim,
        f=genz_malik_1980_f_2,
        closed_form=genz_malik_1980_f_2_exact,
        a=a,
        b=b,
        args=(
            alphas,
            betas,
        ),
    )

    res = cub(
        problem.f,
        problem.a,
        problem.b,
        rule,
        rtol,
        atol,
        args=problem.args,
    )

    assert_integral_estimate_close(problem, res, rtol, atol)


@pytest.mark.skip()
@pytest.mark.parametrize("quadrature", [GaussKronrod21()])
@pytest.mark.parametrize("rtol", [1e-5])
@pytest.mark.parametrize("atol", [1e-8])
@pytest.mark.parametrize("dim", [1, 2, 3])
@pytest.mark.parametrize("args_shape", [(2,), (2, 3), (2, 3, 5), (2, 3, 5, 7)])
def test_genz_malik_1980_f_3_arbitrary_shape(quadrature, rtol, atol, dim, args_shape):
    np.random.seed(1)
    difficulty = 12
    rule = ProductRule([quadrature] * dim)

    a = np.array([0]*dim)
    b = np.array([1]*dim)
    alphas = np.random.rand(dim, *args_shape)

    # Adjust the parameters as in the paper
    alphas = difficulty * alphas / np.sum(alphas, axis=0)

    assert np.allclose(np.sum(alphas, axis=0), difficulty)

    problem = problem_from_closed_form(
        dim=dim,
        f=genz_malik_1980_f_3,
        closed_form=genz_malik_1980_f_3_exact,
        a=a,
        b=b,
        args=(
            alphas,
        ),
    )

    res = cub(
        problem.f,
        problem.a,
        problem.b,
        rule,
        rtol,
        atol,
        args=problem.args,
    )

    assert_integral_estimate_close(problem, res, rtol, atol)


@pytest.mark.skip()
@pytest.mark.parametrize("quadrature", [GaussKronrod21()])
@pytest.mark.parametrize("rtol", [1e-5])
@pytest.mark.parametrize("atol", [1e-8])
@pytest.mark.parametrize("dim", [1, 2, 3])
@pytest.mark.parametrize("args_shape", [(2,), (2, 3), (2, 3, 5), (2, 3, 5, 7)])
def test_genz_malik_1980_f_5_arbitrary_shape(quadrature, rtol, atol, dim, args_shape):
    np.random.seed(1)
    difficulty = 21
    rule = ProductRule([quadrature] * dim)

    a = np.array([0]*dim)
    b = np.array([1]*dim)
    alphas = np.random.rand(dim, *args_shape)
    betas = np.random.rand(dim, *args_shape)

    # Adjust the parameters as in the paper
    l2_norm = np.sum(np.power(alphas, 2), axis=0)
    alphas = alphas / np.sqrt(l2_norm) * np.sqrt(difficulty)

    assert np.allclose(np.sum(np.power(alphas, 2), axis=0), difficulty)

    problem = problem_from_closed_form(
        dim=dim,
        f=genz_malik_1980_f_5,
        closed_form=genz_malik_1980_f_5_exact,
        a=a,
        b=b,
        args=(
            alphas,
            betas,
        ),
    )

    res = cub(
        problem.f,
        problem.a,
        problem.b,
        rule,
        rtol,
        atol,
        args=problem.args,
    )

    assert_integral_estimate_close(problem, res, rtol, atol)


# TODO: test that product rules are calculated properly
# TODO: test that inconsistent dimensions are reported
