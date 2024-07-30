from functools import cached_property

from scipy.special import roots_legendre

from ._base import FixedCub


class GaussLegendreQuad(FixedCub):
    """
    Gauss-Legendre quadrature.

    Gauss-Legendre is a 1D rule. To use it for multidimensional integrals, it will be
    necessary to take the FixedProductCub of multiple Gauss-Legendre rules. See
    Examples.

    Parameters
    ----------
    npoints : int
        Number of nodes for the higher-order rule.

    Examples
    --------
    Evaluate a 1D integral. Note in this example that ``f`` returns an array, so the
    estimates will also be arrays, despite the fact that this is a 1D problem.

    >>> import numpy as np
    >>> from scipy.integrate._cubature import cub
    >>> from scipy.integrate._rules import GaussLegendreQuad
    >>> def f(x):
    ...     return np.cos(x)
    >>> rule = GaussLegendreQuad(21) # Use 21-point GaussLegendre
    >>> a, b = np.array([0]), np.array([1])
    >>> rule.estimate(f, a, b) # True value sin(1), approximately 0.84147
     array([0.84147098])
    >>> rule.error_estimate(f, a, b)
     array([1.11022302e-16])

    Evaluate a 2D integral. Note that in this example ``f`` returns a float, so the
    estimates will also be floats.

    >>> import numpy as np
    >>> from scipy.integrate._cubature import cub
    >>> from scipy.integrate._rules import FixedProductCub, GaussLegendreQuad
    >>> def f(x):
    ...     # f(x) = cos(x_1) + cos(x_2)
    ...     return np.sum(np.cos(x), axis=0)
    >>> rule = FixedProductCub(
    ...     [GaussLegendreQuad(15), GaussLegendreQuad(15)]
    ... ) # Use 15-point GaussKronrod
    >>> a, b = np.array([0, 0]), np.array([1, 1])
    >>> rule.estimate(f, a, b) # True value 2*sin(1), approximately 1.6829
     np.float64(1.682941969615793)
    >>> rule.error_estimate(f, a, b)
     np.float64(2.220446049250313e-16)
    """

    def __init__(self, npoints):
        if npoints < 2:
            raise ValueError(
                "At least 2 nodes required for Gauss-Legendre cubature"
            )

        self.npoints = npoints

    @cached_property
    def nodes_and_weights(self):
        return roots_legendre(self.npoints)
