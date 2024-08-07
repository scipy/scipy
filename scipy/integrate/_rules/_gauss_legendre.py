from functools import cached_property

from scipy.special import roots_legendre

from ._base import FixedRule


class GaussLegendreQuad(FixedRule):
    """
    Gauss-Legendre quadrature.

    Gauss-Legendre is a 1D rule. To use it for multidimensional integrals, it will be
    necessary to take the product of multiple Gauss-Legendre rules using `ProductFixed`.
    See Examples.

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
    >>> from scipy.integrate._rules import GaussLegendreQuad
    >>> def f(x):
    ...     return np.cos(x)
    >>> rule = GaussLegendreQuad(21) # Use 21-point GaussLegendre
    >>> a, b = np.array([0]), np.array([1])
    >>> rule.estimate(f, a, b) # True value sin(1), approximately 0.84147
     array([0.84147098])
    >>> rule.estimate_error(f, a, b)
     array([1.11022302e-16])

    Evaluate a 2D integral. Note that in this example ``f`` returns a float, so the
    estimates will also be floats.

    >>> import numpy as np
    >>> from scipy.integrate import cubature
    >>> from scipy.integrate._rules import ProductFixed, GaussLegendreQuad
    >>> def f(x):
    ...     # f(x) = cos(x_1) + cos(x_2)
    ...     return np.sum(np.cos(x), axis=-1)
    >>> rule = ProductFixed(
    ...     [GaussLegendreQuad(15), GaussLegendreQuad(15)]
    ... )
    >>> a, b = np.array([0, 0]), np.array([1, 1])
    >>> rule.estimate(f, a, b) # True value 2*sin(1), approximately 1.6829
     np.float64(1.682941969615793)
    >>> rule.estimate_error(f, a, b)
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
