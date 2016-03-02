# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import

from numpy import asarray, atleast_2d, floor, arange, sin, sqrt, prod, sum, round
from .go_benchmark import Benchmark


class Katsuura(Benchmark):

    r"""
    Katsuura objective function.

    This class defines the Katsuura [1]_ global optimization problem. This is a
    multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Katsuura}}(x) = \prod_{i=0}^{n-1} \left [ 1 +
        (i+1) \sum_{k=1}^{d} \lfloor (2^k x_i) \rfloor 2^{-k} \right ]


    Where, in this exercise, :math:`d = 32`.

    Here, :math:`n` represents the number of dimensions and 
    :math:`x_i \in [0, 100]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 1` for :math:`x_i = 0` for
    :math:`i = 1, ..., n`.

    .. [1] Adorio, E. MVF - "Multivariate Test Functions Library in C for
    Unconstrained Global Optimization", 2005
    .. [2] Gavana, A. Global Optimization Benchmarks and AMPGO retrieved 2015

    TODO: Adorio has wrong global minimum.  Adorio uses round, Gavana docstring
    uses floor, but Gavana code uses round.  We'll use round...
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N, [100.0] * self.N)

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.custom_bounds = [(0, 1), (0, 1)]
        self.fglob = 1.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        d = 32
        k = atleast_2d(arange(1, d + 1)).T
        i = arange(0., self.N * 1.)
        inner = round(2 ** k * x) * (2. ** (-k))
        return prod(sum(inner, axis=0) * (i + 1) + 1)


class Keane(Benchmark):

    r"""
    Keane objective function.

    This class defines the Keane [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Keane}}(x) = \frac{\sin^2(x_1 - x_2)\sin^2(x_1 + x_2)}
        {\sqrt{x_1^2 + x_2^2}}


    with :math:`x_i \in [0, 10]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = 0.0` for 
    :math:`x = [7.85396153, 7.85396135]`.

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    TODO: Jamil #69, there is no way that the function can have a negative
    value.  Everything is squared.  I think that they have the wrong solution.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[7.85396153, 7.85396135]]
        self.custom_bounds = [(-1, 0.34), (-1, 0.34)]
        self.fglob = 0.

    def fun(self, x, *args):
        self.nfev += 1

        val = sin(x[0] - x[1]) ** 2 * sin(x[0] + x[1]) ** 2
        return val / sqrt(x[0] ** 2 + x[1] ** 2)


class Kowalik(Benchmark):

    r"""
    Kowalik objective function.

    This class defines the Kowalik [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Kowalik}}(x) = \sum_{i=0}^{10} \left [ a_i
        - \frac{x_1 (b_i^2 + b_i x_2)} {b_i^2 + b_i x_3 + x_4} \right ]^2

    Where:

    .. math::

        \begin{matrix}
        a = [4, 2, 1, 1/2, 1/4 1/8, 1/10, 1/12, 1/14, 1/16] \\
        b = [0.1957, 0.1947, 0.1735, 0.1600, 0.0844, 0.0627,
             0.0456, 0.0342, 0.0323, 0.0235, 0.0246]\\
        \end{matrix}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \in 
    [-5, 5]` for :math:`i = 1, ..., 4`.

    *Global optimum*: :math:`f(x) = 0.00030748610` for :math:`x = 
    [0.192833, 0.190836, 0.123117, 0.135766]`.

    ..[1] http://www.itl.nist.gov/div898/strd/nls/data/mgh09.shtml
    """

    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.0] * self.N, [5.0] * self.N)
        self.global_optimum = [[0.192833, 0.190836, 0.123117, 0.135766]]
        self.fglob = 0.00030748610

        self.a = asarray([4.0, 2.0, 1.0, 1 / 2.0, 1 / 4.0, 1 / 6.0, 1 / 8.0,
                          1 / 10.0, 1 / 12.0, 1 / 14.0, 1 / 16.0])
        self.b = asarray([0.1957, 0.1947, 0.1735, 0.1600, 0.0844, 0.0627,
                          0.0456, 0.0342, 0.0323, 0.0235, 0.0246])

    def fun(self, x, *args):
        self.nfev += 1

        vec = self.b - (x[0] * (self.a ** 2 + self.a * x[1])
                   / (self.a ** 2 + self.a * x[2] + x[3]))
        return sum(vec ** 2)
