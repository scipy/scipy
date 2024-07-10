import scipy
import math

import pytest
import numpy as np

# import scipy
from scipy.integrate._cubature import cub, ProductRule, GaussKronrod15


class Problem:
    """Represents a general problem to be solved by quadrature algorithm"""
    def __init__(self, dim, f, a, b, exact):
        self.dim = dim
        self.f = f
        self.a = a
        self.b = b
        self.exact = exact


def problem_from_closed_form(dim, parameterized_f, closed_form, a, b, **params):
    """
    Creates a new Problem given a function `parameterised_f(**params)` which will return
    a function `f` to be integrated specified in terms of some parameters, and also
    given a `closed_form` which returns the exact solution to the problem in terms of
    the bounds of integration `a` and `b`, and the parameters used to construct `f`.
    """

    return Problem(
        dim=dim,
        f=parameterized_f(**params),
        a=a,
        b=b,
        exact=closed_form(dim, a, b, **params)
    )


def genz_malik_1980_f_1(r, alphas):
    r"""
    `f_1` from Genz and Malik 1980.

    .. math:: f_1(\mathbf x) = \cos\left(2\pi r + \sum^n_{i = 1}\alpha_i x_i\right)

    .. code-block:: mathematica

        genzMalik1980f1[x_List, r_, alphas_List] := Cos[2*Pi*r + Total[x*alphas]]
    """
    def f(x):
        alphas_reshaped = alphas[:, np.newaxis]
        return np.cos(2*np.pi*r + np.sum(x * alphas_reshaped, axis=0))

    return f


def genz_malik_1980_f_1_exact(dim, a, b, r, alphas):
    return (-2)**dim * 1/math.prod(alphas) \
        * np.cos(2*np.pi*r + np.sum((a * alphas) + (b * alphas)) / 2) \
        * np.prod(np.sin((a * alphas - b * alphas)/2))


def genz_malik_1980_f_2(alphas, betas):
    r"""
    `f_2` from Genz and Malik 1980.

    .. math:: f_2(\mathbf x) = \prod^n_{i = 1} (\alpha_i^2 + (x_i - \beta_i)^2)^{-1}

    .. code-block:: mathematica

        genzMalik1980f2[x_List, alphas_List, betas_List] :=
            1/Times @@ ((alphas^2 + (x - betas)^2))
    """

    def f(x):
        alphas_reshaped = alphas[:, np.newaxis]
        betas_reshaped = betas[:, np.newaxis]
        return (1/np.prod(alphas_reshaped**2 + (x - betas_reshaped)**2, axis=0))

    return f


def genz_malik_1980_f_2_exact(dim, a, b, alphas, betas):
    return (-1)**dim * 1/math.prod(alphas) * np.prod(
        np.arctan((a - betas)/alphas) - np.arctan((b - betas)/alphas)
    )


def genz_malik_1980_f_3(alphas):
    r"""
    `f_3` from Genz and Malik 1980.

    .. math:: f_3(\mathbf x) = \exp\left(\sum^n_{i = 1} \alpha_i x_i\right)

    .. code-block:: mathematica

        genzMalik1980f3[x_List, alphas_List] := Exp[Dot[x, alphas]]
    """

    def f(x):
        alphas_reshaped = alphas[:, np.newaxis]
        return np.exp(np.sum(alphas_reshaped * x, axis=0))

    return f


# TODO: not 100% sure this is the correct closed form
def genz_malik_1980_f_3_exact(dim, a, b, alphas):
    return (-1)**dim * 1/math.prod(alphas) \
        * np.prod(np.exp(alphas * a) - np.exp(alphas * b))


def genz_malik_1980_f_4(alphas):
    r"""
    `f_4` from Genz and Malik 1980.

    .. math:: f_4(\mathbf x) = \left(1 + \sum^n_{i = 1} \alpha_i x_i\right)^{-n-1}

    .. code-block:: mathematica

        genzMalik1980f4[x_List, alphas_List] :=
            (1 + Dot[x, alphas])^(-Length[alphas] - 1)
    """

    def f(x):
        alphas_reshaped = alphas[:, np.newaxis]
        return ((1 + np.sum(alphas_reshaped * x, axis=0))**(-x.shape[0]-1))

    return f


def genz_malik_1980_f_5(alphas, betas):
    r"""
    `f_5` from Genz and Malik 1980.

    .. math::

        f_5(\mathbf x) = \exp\left(-\sum^n_{i = 1} \alpha^2_i (x_i - \beta_i)^2\right)

    .. code-block:: mathematica

        genzMalik1980f5[x_List, alphas_List, betas_List] :=
            Exp[-Total[alphas^2 * (x - betas)^2]]
    """

    def f(x):
        alphas_reshaped = alphas[:, np.newaxis]
        betas_reshaped = betas[:, np.newaxis]
        return np.exp(-np.sum((alphas_reshaped**2) * (x - betas_reshaped)**2, axis=0))

    return f


# TODO: not 100% sure this is the correct closed form
def genz_malik_1980_f_5_exact(dim, a, b, alphas, betas):
    return (1/2)**dim * 1/math.prod(alphas) * np.power(np.pi, dim/2) * np.prod(
        scipy.special.erf(alphas * (betas - a))
        + scipy.special.erf(alphas * (b - betas))
    )


problems = [
    # Problems based on test integrals in Genz and Malik 1980:
    # -- f1 --

    problem_from_closed_form(
        dim=1,
        parameterized_f=genz_malik_1980_f_1,
        closed_form=genz_malik_1980_f_1_exact,

        a=np.array([0]),
        b=np.array([10]),
        r=1/4,
        alphas=np.array([5])
    ),
    problem_from_closed_form(
        dim=2,
        parameterized_f=genz_malik_1980_f_1,
        closed_form=genz_malik_1980_f_1_exact,

        a=np.array([0, 0]),
        b=np.array([1, 1]),
        r=1/4,
        alphas=np.array([2, 4]),
    ),
    problem_from_closed_form(
        dim=2,
        parameterized_f=genz_malik_1980_f_1,
        closed_form=genz_malik_1980_f_1_exact,

        a=np.array([0, 0]),
        b=np.array([100, 100]),
        r=1/2,
        alphas=np.array([2, 4]),
    ),
    problem_from_closed_form(
        dim=3,
        parameterized_f=genz_malik_1980_f_1,
        closed_form=genz_malik_1980_f_1_exact,

        a=np.array([0, 0, 0]),
        b=np.array([100, 200, 300]),
        r=1/2,
        alphas=np.array([1, 1, 1]),
    ),

    # -- f2 --

    problem_from_closed_form(
        dim=1,
        parameterized_f=genz_malik_1980_f_2,
        closed_form=genz_malik_1980_f_2_exact,

        a=np.array([-1]),
        b=np.array([1]),
        alphas=np.array([5]),
        betas=np.array([4])
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
        parameterized_f=genz_malik_1980_f_2,
        closed_form=genz_malik_1980_f_2_exact,

        a=np.array([0, 0, 0]),
        b=np.array([1, 1, 1]),
        alphas=np.array([1, 1, 1]),
        betas=np.array([1, 1, 1]),
    ),
    problem_from_closed_form(
        dim=3,
        parameterized_f=genz_malik_1980_f_2,
        closed_form=genz_malik_1980_f_2_exact,

        a=np.array([0, 0, 0]),
        b=np.array([1, 1, 1]),
        alphas=np.array([2, 3, 4]),
        betas=np.array([2, 3, 4]),
    ),
    problem_from_closed_form(
        dim=3,
        parameterized_f=genz_malik_1980_f_2,
        closed_form=genz_malik_1980_f_2_exact,

        a=np.array([-1, -1, -1]),
        b=np.array([1, 1, 1]),
        alphas=np.array([1, 1, 1]),
        betas=np.array([2, 2, 2]),
    ),
    problem_from_closed_form(
        dim=4,
        parameterized_f=genz_malik_1980_f_2,
        closed_form=genz_malik_1980_f_2_exact,

        a=np.array([-1, -1, -1, -1]),
        b=np.array([1, 1, 1, 1]),
        alphas=np.array([1, 1, 1, 1]),
        betas=np.array([1, 1, 1, 1]),
    ),

    # -- f_3 --

    problem_from_closed_form(
        dim=1,
        parameterized_f=genz_malik_1980_f_3,
        closed_form=genz_malik_1980_f_3_exact,

        a=np.array([-1]),
        b=np.array([1]),
        alphas=np.array([1/2]),
    ),
    problem_from_closed_form(
        dim=2,
        parameterized_f=genz_malik_1980_f_3,
        closed_form=genz_malik_1980_f_3_exact,

        a=np.array([0, -1]),
        b=np.array([1, 1]),
        alphas=np.array([5, 5]),
    ),
    problem_from_closed_form(
        dim=3,
        parameterized_f=genz_malik_1980_f_3,
        closed_form=genz_malik_1980_f_3_exact,

        a=np.array([-1, -1, -1]),
        b=np.array([1, 1, 1]),
        alphas=np.array([1, 1, 1]),
    ),

    # -- f_4 --

    Problem(
        dim=1,
        f=genz_malik_1980_f_4(
            alphas=np.array([1])
        ),
        a=np.array([0]),
        b=np.array([2]),
        exact=2/3,
    ),
    Problem(
        dim=2,
        f=genz_malik_1980_f_4(
            alphas=np.array([1, 1]),
        ),
        a=np.array([0, 0]),
        b=np.array([2, 1]),
        exact=5/24,
    ),
    Problem(
        dim=3,
        f=genz_malik_1980_f_4(
            alphas=np.array([1, 1, 1]),
        ),
        a=np.array([0, 0, 0]),
        b=np.array([1, 1, 1]),
        exact=1/24,
    ),

    # -- f_5 --

    problem_from_closed_form(
        dim=1,
        parameterized_f=genz_malik_1980_f_5,
        closed_form=genz_malik_1980_f_5_exact,

        a=np.array([-1]),
        b=np.array([1]),
        alphas=np.array([-2]),
        betas=np.array([2]),
    ),
    problem_from_closed_form(
        dim=2,
        parameterized_f=genz_malik_1980_f_5,
        closed_form=genz_malik_1980_f_5_exact,

        a=np.array([-1, -1]),
        b=np.array([1, 1]),
        alphas=np.array([2, 3]),
        betas=np.array([4, 5]),
    ),
    problem_from_closed_form(
        dim=3,
        parameterized_f=genz_malik_1980_f_5,
        closed_form=genz_malik_1980_f_5_exact,

        a=np.array([-1, -1, -1]),
        b=np.array([1, 1, 1]),
        alphas=np.array([1, 1, 1]),
        betas=np.array([1, 1, 1]),
    )
]


@pytest.mark.parametrize("problem", problems)
@pytest.mark.parametrize("quadrature", [GaussKronrod15()])
@pytest.mark.parametrize("rtol", [1e-5])
@pytest.mark.parametrize("atol", [1e-8])
def test_cub_product_rules(problem, quadrature, rtol, atol):
    rule = ProductRule([quadrature] * problem.dim)

    res = cub(
        problem.f,
        problem.a,
        problem.b,
        rule,
        rtol,
        atol,
    )

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

    assert res.status == "converged"


# TODO: test that product rules are calculated properly
# TODO: test that inconsistent dimensions are reported
