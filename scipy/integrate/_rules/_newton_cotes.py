import numpy as np

from functools import cached_property

from ._base import FixedRule


class NewtonCotesQuadrature(FixedRule):
    """
    Newton-Cotes quadrature.

    Newton-Cotes is a 1D rule. To use it for multidimensional integrals, it will be
    necessary to take the product of multiple Newton-Cotes rules using `ProductFixed`.
    See Examples.

    Parameters
    ----------
    npoints : int
        Number of nodes for the higher-order rule.

    Examples
    --------
    Evaluate a 1D integral. Note in this example that ``f`` returns an array, so the
    estimates will also be arrays, despite the fact that this is a 1D problem.

    >>> import numpy as np
    >>> from scipy.integrate import cubature
    >>> from scipy.integrate._rules import NewtonCotesQuadrature
    >>> def f(x):
    ...     return np.cos(x)
    >>> rule = NewtonCotesQuadrature(21) # Use 21-point GaussLegendre
    >>> a, b = np.array([0]), np.array([1])
    >>> rule.estimate(f, a, b) # True value sin(1), approximately 0.84147
     array([0.84147098])
    >>> rule.estimate_error(f, a, b)
     array([1.11022302e-16])

    Evaluate a 2D integral. Note that in this example ``f`` returns a float, so the
    estimates will also be floats.

    >>> import numpy as np
    >>> from scipy.integrate import cubature
    >>> from scipy.integrate._rules import ProductFixed, NewtonCotesQuadrature
    >>> def f(x):
    ...     # f(x) = cos(x_1) + cos(x_2)
    ...     return np.sum(np.cos(x), axis=-1)
    >>> rule = ProductFixed(
    ...     [NewtonCotesQuadrature(15), NewtonCotesQuadrature(15)]
    ... ) # Use 15-point Newton-Cotes
    >>> a, b = np.array([0, 0]), np.array([1, 1])
    >>> rule.estimate(f, a, b) # True value 2*sin(1), approximately 1.6829
     np.float64(1.682941969615793)
    >>> rule.estimate_error(f, a, b)
     np.float64(2.220446049250313e-16)
    """

    def __init__(self, npoints, open=False):
        if npoints < 2:
            raise ValueError(
                "at least 2 points required for Newton-Cotes cubature"
            )

        self.npoints = npoints
        self.open = open

    @cached_property
    def nodes_and_weights(self):
        if self.open:
            h = 2/self.npoints
            nodes = np.linspace(-1 + h, 1 - h, num=self.npoints)
        else:
            nodes = np.linspace(-1, 1, num=self.npoints)

        weights = _newton_cotes_weights(nodes)

        return nodes, weights


def _newton_cotes_weights(points):
    order = len(points) - 1

    a = np.vander(points, increasing=True)
    a = np.transpose(a)

    i = np.arange(order + 1)
    b = (1 - np.power(-1, i + 1)) / (i + 1)

    return np.linalg.solve(a, b)
