from functools import cached_property

from scipy.special import roots_legendre

from ._base import FixedRule


class GaussLegendreQuadrature(FixedRule):
    """
    Gauss-Legendre quadrature.

    Parameters
    ----------
    npoints : int
        Number of nodes for the higher-order rule.

    Examples
    --------
    Evaluate a 1D integral. Note in this example that ``f`` returns an array, so the
    estimates will also be arrays.

    >>> import numpy as np
    >>> from scipy.integrate import cubature
    >>> from scipy.integrate._rules import GaussLegendreQuadrature
    >>> def f(x):
    ...     return np.cos(x)
    >>> rule = GaussLegendreQuadrature(21) # Use 21-point GaussLegendre
    >>> a, b = np.array([0]), np.array([1])
    >>> rule.estimate(f, a, b) # True value sin(1), approximately 0.84147
     array([0.84147098])
    >>> rule.estimate_error(f, a, b)
     array([1.11022302e-16])
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
