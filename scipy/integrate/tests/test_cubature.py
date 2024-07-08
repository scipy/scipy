import numpy as np
import scipy

from scipy.integrate._cubature import cub, ProductRule, GaussKronrod21
import scipy.special


class Problem:
    """Represents a general problem to be solved by a numerical integrator"""
    def __init__(self, dim, f, ranges, exact, tol, rule):
        self.dim = dim
        self.f = f
        self.ranges = ranges
        self.exact = exact
        self.tol = tol
        self.rule = rule

    def calculate(self):
        integral_estimate, error_estimate = cub(
            self.f,
            self.ranges,
            self.rule,
            self.tol,
        )

        self.integral_estimate = integral_estimate
        self.error_estimate = error_estimate

    def assert_within_tol(self):
        np.testing.assert_allclose(
            self.integral_estimate,
            self.exact,
            rtol=0,
            atol=self.tol,
            verbose=True,
            err_msg=f"""failed to find approx. integral of {self.f.__name__} with \
rtol=TODO, atol={self.tol} with ranges {self.ranges}"""
        )


def genz_malik_1980_f_1(r, alphas):
    def f(x):
        return np.cos(2*np.pi*r + np.sum(x * alphas, axis=1))[:, np.newaxis]

    return f


def test_genz_malik_1980_f_1():
    def genz_malik_1980_f_1(r, alphas):
        def f(x):
            return np.cos(2*np.pi*r + np.sum(x * alphas, axis=1))[:, np.newaxis]

        return f

    problem = Problem(
        # Three dimensional problem
        dim=3,

        # Test function in terms of these parameters
        f=genz_malik_1980_f_1(
            r=1/2,
            alphas=np.array([1, 1, 1]),
        ),

        # Ranges to integrate over
        ranges=np.array([[0, 1], [0, 1], [0, 1]]),

        # Exact answer to the problem
        exact=-8 * np.cos(3/2) * np.sin(1/2) ** 3,

        # Tolerance in which we wish to attempt solving the problem
        tol=1e-8,

        # Rule to use
        rule=ProductRule([GaussKronrod21(), GaussKronrod21(), GaussKronrod21()]),
    )
    problem.calculate()
    problem.assert_within_tol()


def test_genz_malik_1980_f_2():
    def genz_malik_1980_f_2(alphas, betas):
        def f(x):
            return (1/np.prod(alphas**2 + (x - betas)**2, axis=1))[:, np.newaxis]

        return f

    problem = Problem(
        # Three dimensional problem
        dim=3,

        # Test function in terms of these parameters
        f=genz_malik_1980_f_2(
            alphas=np.array([1, 1, 1]),
            betas=np.array([1, 1, 1]),
        ),

        # Ranges to integrate over
        ranges=np.array([[0, 1], [0, 1], [0, 1]]),

        # Exact answer to the problem
        exact=np.pi ** 3 / 64,

        # Tolerance in which we wish to attempt solving the problem
        tol=1e-8,

        # Rule to use
        rule=ProductRule([GaussKronrod21(), GaussKronrod21(), GaussKronrod21()]),
    )
    problem.calculate()
    problem.assert_within_tol()

    problem = Problem(
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
        tol=1e-8,
        rule=ProductRule([GaussKronrod21(), GaussKronrod21(), GaussKronrod21()]),
    )
    problem.calculate()
    problem.assert_within_tol()

    problem = Problem(
        dim=3,
        f=genz_malik_1980_f_2(
            alphas=np.array([1, 1, 1]),
            betas=np.array([2, 2, 2]),
        ),
        ranges=np.array([[-1, 1], [-1, 1], [-1, 1]]),
        exact=(
            -1/64 * (np.pi - 4*np.arctan(3))**3
        ),
        tol=1e-8,
        rule=ProductRule([GaussKronrod21(), GaussKronrod21(), GaussKronrod21()]),
    )
    problem.calculate()
    problem.assert_within_tol()

    problem = Problem(
        dim=4,
        f=genz_malik_1980_f_2(
            alphas=np.array([1, 1, 1, 1]),
            betas=np.array([1, 1, 1, 1]),
        ),
        ranges=np.array([[-1, 1], [-1, 1], [-1, 1], [-1, 1]]),
        exact=(
            np.arctan(2)**4
        ),
        tol=1e-8,
        rule=ProductRule([
            GaussKronrod21(),
            GaussKronrod21(),
            GaussKronrod21(),
            GaussKronrod21(),
        ]),
    )
    problem.calculate()
    problem.assert_within_tol()


def test_genz_malik_1980_f_3():
    def genz_malik_1980_f_3(alphas):
        def f(x):
            return np.exp(np.sum(alphas * x, axis=1))[:, np.newaxis]

        return f

    problem = Problem(
        dim=3,
        f=genz_malik_1980_f_3(
            alphas=np.array([1, 1, 1]),
        ),
        ranges=np.array([[-1, 1], [-1, 1], [-1, 1]]),
        exact=(
            (-1+np.e**2)**3 / np.e**3
        ),
        tol=1e-8,
        rule=ProductRule([
            GaussKronrod21(),
            GaussKronrod21(),
            GaussKronrod21(),
        ]),
    )
    problem.calculate()
    problem.assert_within_tol()


# TODO: more test integrals
# TODO: test that product rules are calculated properly
# TODO: test that inconsistent dimensions are reported


def test_genz_malik_1980_f_4():
    def genz_malik_1980_f_4(alphas):
        def f(x):
            return ((1 + np.sum(alphas * x, axis=1))**(-x.shape[-1]-1))[:, np.newaxis]

        return f

    problem = Problem(
        dim=3,
        f=genz_malik_1980_f_4(
            alphas=np.array([1, 1, 1]),
        ),
        ranges=np.array([[0, 1], [0, 1], [0, 1]]),
        exact=(
            1/24
        ),
        tol=1e-8,
        rule=ProductRule([
            GaussKronrod21(),
            GaussKronrod21(),
            GaussKronrod21(),
        ]),
    )
    problem.calculate()
    problem.assert_within_tol()


def test_genz_malik_1980_f_5():
    def genz_malik_1980_f_5(alphas, betas):
        def f(x):
            return np.exp(-np.sum((alphas**2) * (x - betas)**2, axis=1))[:, np.newaxis]

        return f

    problem = Problem(
        dim=3,
        f=genz_malik_1980_f_5(
            alphas=np.array([1, 1, 1]),
            betas=np.array([1, 1, 1])
        ),
        ranges=np.array([[-1, 1], [-1, 1], [-1, 1]]),
        exact=(
            1/8 * np.pi**(3/2) * scipy.special.erf(2)**3
        ),
        tol=1e-8,
        rule=ProductRule([
            GaussKronrod21(),
            GaussKronrod21(),
            GaussKronrod21(),
        ]),
    )
    problem.calculate()
    problem.assert_within_tol()
