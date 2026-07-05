import numpy as np

from scipy._lib._array_api import xp_compat_namespace

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

    xp : array_namespace, optional
        The namespace for the node and weight arrays. Default is None, where NumPy is
        used.

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

    def __init__(self, npoints, xp=None):
        if npoints < 2:
            raise ValueError(
                "At least 2 nodes required for Gauss-Legendre cubature"
            )

        self.npoints = npoints

        self.xp = xp_compat_namespace(xp)

    @cached_property
    def nodes_and_weights(self):
        # The nodes and weights are constants, kept as NumPy (host) data so
        # their values exist independently of any array library's default
        # device; they are converted onto the namespace and device of the
        # integration limits at apply time (see `_cached_cast`, gh-22680).
        nodes, weights = roots_legendre(self.npoints)

        return (
            np.asarray(nodes, dtype=np.float64),
            np.asarray(weights, dtype=np.float64),
        )
