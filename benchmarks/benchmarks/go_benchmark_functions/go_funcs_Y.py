# -*- coding: utf-8 -*-
from numpy import abs, sum, cos, pi
from .go_benchmark import Benchmark


class YaoLiu04(Benchmark):

    r"""
    Yao-Liu 4 objective function.

    This class defines the Yao-Liu function 4 [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{YaoLiu04}}(x) = {max}_i \left\{ \left | x_i \right | ,
                                 1 \leq i \leq n \right\}

    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-10, 10]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x_i = 0` for
    :math:`i = 1, ..., n`

    .. [1] Yao X., Liu Y. (1997) Fast evolution strategies.
    In: Angeline P.J., Reynolds R.G., McDonnell J.R., Eberhart R. (eds)
    Evolutionary Programming VI. EP 1997.
    Lecture Notes in Computer Science, vol 1213. Springer, Berlin, Heidelberg

    .. [2] Mishra, S. Global Optimization by Differential Evolution and
    Particle Swarm Methods: Evaluation on Some Benchmark Functions.
    Munich Personal RePEc Archive, 2006, 1005

    TODO line 1201.  Gavana code and documentation differ.
    max(abs(x)) != abs(max(x))
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-10.0] * self.N, [10.0] * self.N))

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return abs(x).max()


class YaoLiu09(Benchmark):

    r"""
    Yao-Liu 9 objective function.

    This class defines the Yao-Liu [1]_ function 9 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{YaoLiu09}}(x) = \sum_{i=1}^n \left [ x_i^2
                                 - 10 \cos(2 \pi x_i ) + 10 \right ]

    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-5.12, 5.12]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x_i = 0` for
    :math:`i = 1, ..., n`

    .. [1] Yao X., Liu Y. (1997) Fast evolution strategies.
    In: Angeline P.J., Reynolds R.G., McDonnell J.R., Eberhart R. (eds)
    Evolutionary Programming VI. EP 1997.
    Lecture Notes in Computer Science, vol 1213. Springer, Berlin, Heidelberg

    .. [2] Gavana, A. Global Optimization Benchmarks and AMPGO retrieved 2015
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-5.12] * self.N, [5.12] * self.N))

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum(x ** 2.0 - 10.0 * cos(2 * pi * x) + 10)
