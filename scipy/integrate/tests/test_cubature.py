import pytest

import numpy as np
import scipy

from scipy.integrate._cubature import cub, ProductRule, GaussKronrod21
import scipy.special


class Problem:
    """Represents a general problem to be solved by quadrature algorithm"""
    def __init__(self, dim, f, ranges, exact):
        self.dim = dim
        self.f = f
        self.ranges = ranges
        self.exact = exact


def genz_malik_1980_f_1(r, alphas):
    r"""
    `f_1` from Genz and Malik 1980.

    .. math:: f_1(\mathbf x) = \cos\left(2\pi r + \sum^n_{i = 1}\alpha_i x_i\right)

    .. code-block:: mathematica

        genzMalik1980f1[x_List, r_, alphas_List] := Cos[2*Pi*r + Total[x*alphas]]
    """
    def f(x):
        return np.cos(2*np.pi*r + np.sum(x * alphas, axis=1))[:, np.newaxis]

    return f


def genz_malik_1980_f_2(alphas, betas):
    r"""
    `f_2` from Genz and Malik 1980.

    .. math:: f_2(\mathbf x) = \prod^n_{i = 1} (\alpha_i^2 + (x_i - \beta_i)^2)^{-1}

    .. code-block:: mathematica

        genzMalik1980f2[x_List, alphas_List, betas_List] :=
            1/Times @@ ((alphas^2 + (x - betas)^2))
    """

    def f(x):
        return (1/np.prod(alphas**2 + (x - betas)**2, axis=1))[:, np.newaxis]

    return f


def genz_malik_1980_f_3(alphas):
    r"""
    `f_3` from Genz and Malik 1980.

    .. math:: f_3(\mathbf x) = \exp\left(\sum^n_{i = 1} \alpha_i x_i\right)

    .. code-block:: mathematica

        genzMalik1980f3[x_List, alphas_List] := Exp[Dot[x, alphas]]
    """

    def f(x):
        return np.exp(np.sum(alphas * x, axis=1))[:, np.newaxis]

    return f


def genz_malik_1980_f_4(alphas):
    r"""
    `f_4` from Genz and Malik 1980.

    .. math:: f_4(\mathbf x) = \left(1 + \sum^n_{i = 1} \alpha_i x_i\right)^{-n-1}

    .. code-block:: mathematica

        genzMalik1980f4[x_List, alphas_List] :=
            (1 + Dot[x, alphas])^(-Length[alphas] - 1)
    """

    def f(x):
        return ((1 + np.sum(alphas * x, axis=1))**(-x.shape[-1]-1))[:, np.newaxis]

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
        return np.exp(-np.sum((alphas**2) * (x - betas)**2, axis=1))[:, np.newaxis]

    return f


problems = [
    # Problems from Genz and Malik 1980:

    Problem(
        # Three dimensional problem
        dim=3,

        # Function to integrate
        f=genz_malik_1980_f_1(
            r=1/2,
            alphas=np.array([1, 1, 1]),
        ),

        # Ranges to integrate over
        ranges=np.array([[0, 1], [0, 1], [0, 1]]),

        # Exact answer to the problem
        exact=-8 * np.cos(3/2) * np.sin(1/2) ** 3,
    ),
    Problem(
        dim=3,
        f=genz_malik_1980_f_2(
            alphas=np.array([1, 1, 1]),
            betas=np.array([1, 1, 1]),
        ),
        ranges=np.array([[0, 1], [0, 1], [0, 1]]),
        exact=np.pi ** 3 / 64,
    ),
    Problem(
        dim=3,
        f=genz_malik_1980_f_2(
            alphas=np.array([2, 3, 4]),
            betas=np.array([2, 3, 4]),
        ),
        ranges=np.array([[0, 1], [0, 1], [0, 1]]),
        exact=(
            1/1536
            * (np.pi - 4*np.arctan(1/2))
            * (np.pi - 4*np.arctan(2/3))
            * (np.pi - 4*np.arctan(3/4))
        ),
    ),
    Problem(
        dim=3,
        f=genz_malik_1980_f_2(
            alphas=np.array([1, 1, 1]),
            betas=np.array([2, 2, 2]),
        ),
        ranges=np.array([[-1, 1], [-1, 1], [-1, 1]]),
        exact=-1/64 * (np.pi - 4*np.arctan(3))**3,
    ),
    Problem(
        dim=4,
        f=genz_malik_1980_f_2(
            alphas=np.array([1, 1, 1, 1]),
            betas=np.array([1, 1, 1, 1]),
        ),
        ranges=np.array([[-1, 1], [-1, 1], [-1, 1], [-1, 1]]),
        exact=np.arctan(2)**4,
    ),
    Problem(
        dim=3,
        f=genz_malik_1980_f_3(
            alphas=np.array([1, 1, 1]),
        ),
        ranges=np.array([[-1, 1], [-1, 1], [-1, 1]]),
        exact=(-1+np.e**2)**3 / np.e**3
    ),
    Problem(
        dim=3,
        f=genz_malik_1980_f_4(
            alphas=np.array([1, 1, 1]),
        ),
        ranges=np.array([[0, 1], [0, 1], [0, 1]]),
        exact=1/24,
    ),
    Problem(
        dim=3,
        f=genz_malik_1980_f_5(
            alphas=np.array([1, 1, 1]),
            betas=np.array([1, 1, 1])
        ),
        ranges=np.array([[-1, 1], [-1, 1], [-1, 1]]),
        exact=1/8 * np.pi**(3/2) * scipy.special.erf(2)**3,
    ),
]


@pytest.mark.parametrize("problem", problems)
@pytest.mark.parametrize("tol", [1e-3, 1e-5, 1e-7])
def test_cub(problem, tol):
    rule = ProductRule([GaussKronrod21()] * problem.dim)

    integral_estimate, _ = cub(
        problem.f,
        problem.ranges,
        rule,
        tol,
    )

    np.testing.assert_allclose(
        integral_estimate,
        problem.exact,
        rtol=0,
        atol=tol,
        verbose=True,
        err_msg=f"""failed to find approx. integral of {problem.f.__name__} with \
rtol=TODO, atol={tol} and ranges {np.array_str(problem.ranges, max_line_width=np.inf)} \
""",
    )


# TODO: more test integrals
# TODO: test that product rules are calculated properly
# TODO: test that inconsistent dimensions are reported
