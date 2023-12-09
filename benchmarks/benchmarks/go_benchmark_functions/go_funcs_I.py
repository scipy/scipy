from numpy import sin, sum
from .go_benchmark import Benchmark


class Infinity(Benchmark):

    r"""
    Infinity objective function.

    This class defines the Infinity [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Infinity}}(x) = \sum_{i=1}^{n} x_i^{6} 
        \left [ \sin\left ( \frac{1}{x_i} \right ) + 2 \right ]


    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-1, 1]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x_i = 0` for
    :math:`i = 1, ..., n`

    .. [1] Gavana, A. Global Optimization Benchmarks and AMPGO retrieved 2015
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-1.0] * self.N, [1.0] * self.N))

        self.global_optimum = [[1e-16 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum(x ** 6.0 * (sin(1.0 / x) + 2.0))
