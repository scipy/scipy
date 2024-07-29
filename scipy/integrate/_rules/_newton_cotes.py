import numpy as np

from functools import cached_property

from ._base import FixedCub


class NewtonCotesQuad(FixedCub):
    """
    Newton-Cotes cubature. Newton-Cotes rules consist of function evaluations at equally
    spaced nodes.

    Newton-Cotes has no error estimator. It is possible to give it an error estimator
    by using DualEstimateCubature to estimate the error as the difference between two
    Newton-Cotes rules. See Examples.

    Newton-Cotes is a 1D rule. To use it for multidimensional integrals, it will be
    necessary to take the product of multiple Newton-Cotes rules. See Examples.

    Parameters
    ----------
    npoints : int
        Number of equally spaced evaluation points.

    open : bool, default=False
        Whether this should be an open rule. Open rules do not include the endpoints
        as nodes.

    Attributes
    ----------
    npoints : int
        Number of equally spaced evaluation points.

    open : bool, default=False
        Whether this is an open rule. Open rules do not include the endpoints as nodes.

    Examples
    --------
    Evaluating a simple integral, first without error estimation and then with error
    estimation:

    >>> import numpy as np
    >>> from scipy.integrate._cubature import *
    >>> def f(x):
    ...     return np.sin(np.sqrt(x))
    >>> rule = NewtonCotes(10)
    >>> a, b = np.array([0]), np.array([np.pi * 2])
    >>> rule.estimate(f, a, b) # Very inaccuracte
     array([0.60003617])
    >>> rule.error_estimate(f, a, b) is None
     True
    >>> rule_with_err_est = DualEstimateCubature(NewtonCotes(10), NewtonCotes(8))
    >>> rule_with_err_est.error_estimate(f, a, b) # Error is high
     array([6.21045267])

    Evaluating a 2D integral, using the product of NewtonCotes with an error estimator:

    >>> import numpy as np
    >>> from scipy.integrate._cubature import *
    >>> def f(x):
    ...     # f(x) = cos(x_1) + cos(x_2)
    ...     return np.sum(np.cos(x), axis=0)
    >>> rule_with_err_est = DualEstimateCubature(NewtonCotes(10), NewtonCotes(8))
    >>> rule_2d = Product([rule_with_err_est, rule_with_err_est])
    >>> a, b = np.array([0, 0]), np.array([1, 1])
    >>> rule_2d.estimate(f, a, b) # True value 2*sin(1), approximately 1.6829
     np.float64(1.6829419696151342)
    >>> rule_2d.error_estimate(f, a, b)
     np.float64(6.823492881835591e-10)
    """

    def __init__(self, npoints, open=False):
        if npoints < 2:
            raise ValueError(
                "at least 2 points required for Newton-Cotes cubature"
            )

        self.npoints = npoints
        self.open = open

    @cached_property
    def rule(self):
        if self.open:
            h = 2/self.npoints
            nodes = np.linspace(-1 + h, 1 - h, num=self.npoints).reshape(1, -1)
        else:
            nodes = np.linspace(-1, 1, num=self.npoints).reshape(1, -1)

        weights = _newton_cotes_weights(nodes.reshape(-1))

        return nodes, weights


def _newton_cotes_weights(points):
    order = len(points) - 1

    a = np.vander(points, increasing=True)
    a = np.transpose(a)

    i = np.arange(order + 1)
    b = (1 - np.power(-1, i + 1)) / (2 * (i + 1))

    # TODO: figure out why this 2x is necessary
    return 2*np.linalg.solve(a, b)
