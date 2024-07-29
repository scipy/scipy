from functools import cached_property

from scipy.special import roots_legendre

from ._base import FixedCub


class GaussLegendreQuad(FixedCub):
    """
    Gauss-Legendre cubature.

    Gauss-Legendre has no error estimator. It is possible to give it an error estimator
    by using DualEstimateCubature to estimate the error as the difference between two
    Gauss-Legendre rules. See Examples.

    Gauss-Legendre is a 1D rule. To use it for multidimensional integrals, it will be
    necessary to take the product of multiple Gauss-Legendre rules. See Examples.

    Parameters
    ----------
    npoints : int
        Number of nodes.

    Attributes
    ----------
    npoints : int
        Number of nodes

    Examples
    --------
    Evaluating a simple integral, first without error estimation and then with error
    estimation:

    TODO
    """

    def __init__(self, npoints):
        if npoints < 2:
            raise ValueError(
                "At least 2 nodes required for Gauss-Legendre cubature"
            )

        self.npoints = npoints

    @cached_property
    def rule(self):
        nodes, weights = roots_legendre(self.npoints)

        # roots_legendre returns nodes as array of shape (self.npoints,), need to
        # reshape so the spacial dimension is first
        return nodes.reshape(1, -1), weights
