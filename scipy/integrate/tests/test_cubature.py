import scipy

import pytest
import numpy as np

# import scipy
from scipy.integrate._cubature import cub, ProductRule, GaussKronrod15


class Problem:
    """Represents a general problem to be solved by quadrature algorithm"""
    def __init__(self, dim, f, a, b, exact, failing=False):
        self.dim = dim
        self.f = f
        self.a = a
        self.b = b
        self.exact = exact
        self.failing = failing


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


problems = [
    # Problems based on test integrals in Genz and Malik 1980:

    Problem(
        # Three dimensional problem
        dim=1,

        # Function to integrate
        f=genz_malik_1980_f_1(
            r=1/4,
            alphas=np.array([5]),
        ),

        # Coordinates of hypercube
        a=np.array([0]),
        b=np.array([10]),

        # Exact answer to the problem
        exact=-2/5 * np.sin(25)**2
    ),
    Problem(
        dim=2,
        f=genz_malik_1980_f_1(
            r=1/4,
            alphas=np.array([2, 4]),
        ),
        a=np.array([0, 0]),
        b=np.array([1, 1]),
        exact=-np.cos(1)*(1 + 2*np.cos(2))*np.sin(1)**3,
    ),
    Problem(
        dim=2,
        f=genz_malik_1980_f_1(
            r=1/2,
            alphas=np.array([2, 4]),
        ),
        a=np.array([0, 0]),
        b=np.array([100, 100]),
        exact=1/4 * np.sin(200) * (np.sin(200) - np.sin(400)),
    ),
    Problem(
        dim=3,
        f=genz_malik_1980_f_1(
            r=1/2,
            alphas=np.array([1, 1, 1]),
        ),
        a=np.array([0, 0, 0]),
        b=np.array([100, 200, 300]),
        exact=4*np.sin(100)*np.sin(150)*(np.sin(250) - np.sin(350))
    ),
    Problem(
        dim=1,
        f=genz_malik_1980_f_2(
            alphas=np.array([5]),
            betas=np.array([4]),
        ),
        a=np.array([-1]),
        b=np.array([1]),
        exact=1/20 * (np.pi - 4*np.arctan(3/5)),
    ),
    # Currently failing for product rule version of Gauss-Kronrod 15 (and 21)
    # dblquad also fails
    Problem(
        dim=2,
        f=genz_malik_1980_f_2(
            alphas=np.array([-3, 3]),
            betas=np.array([-2, 2])
        ),
        a=np.array([-50, 0]),
        b=np.array([50, 20]),
        exact=(
            1/9 * (np.arctan(2/3) + np.arctan(6)) * (np.arctan(16) + np.arctan(52/3))
        ),
        failing=True,
    ),
    Problem(
        dim=3,
        f=genz_malik_1980_f_2(
            alphas=np.array([1, 1, 1]),
            betas=np.array([1, 1, 1]),
        ),
        a=np.array([0, 0, 0]),
        b=np.array([1, 1, 1]),
        exact=np.pi ** 3 / 64,
    ),
    Problem(
        dim=3,
        f=genz_malik_1980_f_2(
            alphas=np.array([2, 3, 4]),
            betas=np.array([2, 3, 4]),
        ),
        a=np.array([0, 0, 0]),
        b=np.array([1, 1, 1]),
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
        a=np.array([-1, -1, -1]),
        b=np.array([1, 1, 1]),
        exact=-1/64 * (np.pi - 4*np.arctan(3))**3,
    ),
    Problem(
        dim=4,
        f=genz_malik_1980_f_2(
            alphas=np.array([1, 1, 1, 1]),
            betas=np.array([1, 1, 1, 1]),
        ),
        a=np.array([-1, -1, -1, -1]),
        b=np.array([1, 1, 1, 1]),
        exact=np.arctan(2)**4,
    ),
    Problem(
        dim=1,
        f=genz_malik_1980_f_3(
            alphas=np.array([1/2]),
        ),
        a=np.array([-1]),
        b=np.array([1]),
        exact=4*np.sinh(1/2),
    ),
    Problem(
        dim=2,
        f=genz_malik_1980_f_3(
            alphas=np.array([5, 5]),
        ),
        a=np.array([0, -1]),
        b=np.array([1, 1]),
        exact=(-1 + np.e**5)**2 * (1 + np.e**5)/(25*np.e**5)
    ),
    Problem(
        dim=3,
        f=genz_malik_1980_f_3(
            alphas=np.array([1, 1, 1]),
        ),
        a=np.array([-1, -1, -1]),
        b=np.array([1, 1, 1]),
        exact=(-1+np.e**2)**3 / np.e**3
    ),
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
    Problem(
        dim=1,
        f=genz_malik_1980_f_5(
            alphas=np.array([-2]),
            betas=np.array([2]),
        ),
        a=np.array([-1]),
        b=np.array([1]),
        exact=1/4 * np.sqrt(np.pi) * (-scipy.special.erf(2) + scipy.special.erf(6)),
    ),
    Problem(
        dim=2,
        f=genz_malik_1980_f_5(
            alphas=np.array([2, 3]),
            betas=np.array([4, 5]),
        ),
        a=np.array([-1, -1]),
        b=np.array([1, 1]),
        exact=1/24 * np.pi
                   * (scipy.special.erf(6) - scipy.special.erf(10))
                   * (scipy.special.erf(12) - scipy.special.erf(18)),
    ),
    Problem(
        dim=3,
        f=genz_malik_1980_f_5(
            alphas=np.array([1, 1, 1]),
            betas=np.array([1, 1, 1])
        ),
        a=np.array([-1, -1, -1]),
        b=np.array([1, 1, 1]),
        exact=1/8 * np.pi**(3/2) * scipy.special.erf(2)**3,
    ),
]


@pytest.mark.parametrize("problem", problems)
@pytest.mark.parametrize("quadrature", [GaussKronrod15()])
@pytest.mark.parametrize("rtol", [1e-5])
@pytest.mark.parametrize("atol", [1e-8])
def test_cub_product_rules(problem, quadrature, rtol, atol):
    rule = ProductRule([quadrature] * problem.dim)

    if problem.failing:
        pytest.xfail()
        return

    res = cub(
        problem.f,
        problem.a,
        problem.b,
        rule,
        rtol,
        atol,
    )

    # The estimate should be within the specified tolerance of the solution
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
