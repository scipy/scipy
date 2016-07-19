# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import

from numpy import abs, sum, sign, arange
from .go_benchmark import Benchmark


class Zacharov(Benchmark):

    r"""
    Zacharov objective function.

    This class defines the Zacharov [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Zacharov}}(x) = \sum_{i=1}^{n} x_i^2 + \left ( \frac{1}{2}
                                 \sum_{i=1}^{n} i x_i \right )^2
                                 + \left ( \frac{1}{2} \sum_{i=1}^{n} i x_i 
                                 \right )^4

    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-5, 10]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x_i = 0` for
    :math:`i = 1, ..., n`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-5.0] * self.N, [10.0] * self.N))
        self.custom_bounds = ([-1, 1], [-1, 1])

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        u = sum(x ** 2)
        v = sum(arange(1, self.N + 1) * x)
        return u + (0.5 * v) ** 2 + (0.5 * v) ** 4


class ZeroSum(Benchmark):

    r"""
    ZeroSum objective function.

    This class defines the ZeroSum [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{ZeroSum}}(x) = \begin{cases}
                                0 & \textrm{if} \sum_{i=1}^n x_i = 0 \\
                                1 + \left(10000 \left |\sum_{i=1}^n x_i\right|
                                \right)^{0.5} & \textrm{otherwise}
                                \end{cases}

    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-10, 10]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0` where :math:`\sum_{i=1}^n x_i = 0`

    .. [1] Gavana, A. Global Optimization Benchmarks and AMPGO retrieved 2015
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-10.0] * self.N, [10.0] * self.N))

        self.global_optimum = [[]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        if abs(sum(x)) < 3e-16:
            return 0.0
        return 1.0 + (10000.0 * abs(sum(x))) ** 0.5


class Zettl(Benchmark):

    r"""
    Zettl objective function.

    This class defines the Zettl [1]_ global optimization problem. This is a
    multimodal minimization problem defined as follows:

    .. math::

       f_{\text{Zettl}}(x) = \frac{1}{4} x_{1} + \left(x_{1}^{2} - 2 x_{1}
                             + x_{2}^{2}\right)^{2}


    with :math:`x_i \in [-1, 5]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = -0.0037912` for :math:`x = [-0.029896, 0.0]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-5.0] * self.N, [10.0] * self.N))

        self.global_optimum = [[-0.02989597760285287, 0.0]]
        self.fglob = -0.003791237220468656

    def fun(self, x, *args):
        self.nfev += 1

        return (x[0] ** 2 + x[1] ** 2 - 2 * x[0]) ** 2 + 0.25 * x[0]


class Zimmerman(Benchmark):

    r"""
    Zimmerman objective function.

    This class defines the Zimmerman [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Zimmerman}}(x) = \max \left[Zh1(x), Zp(Zh2(x))
                                  \textrm{sgn}(Zh2(x)), Zp(Zh3(x))
                                  \textrm{sgn}(Zh3(x)),
                                  Zp(-x_1)\textrm{sgn}(x_1),
                                  Zp(-x_2)\textrm{sgn}(x_2) \right]


    Where, in this exercise:

    .. math::

        \begin{cases}
        Zh1(x) = 9 - x_1 - x_2 \\
        Zh2(x) = (x_1 - 3)^2 + (x_2 - 2)^2 \\
        Zh3(x) = x_1x_2 - 14 \\
        Zp(t) = 100(1 + t)
        \end{cases}


    Where :math:`x` is a vector and :math:`t` is a scalar.

    Here, :math:`x_i \in [0, 100]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x = [7, 2]`

    .. [1] Gavana, A. Global Optimization Benchmarks and AMPGO retrieved 2015

    TODO implementation from Gavana
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([0.0] * self.N, [100.0] * self.N))
        self.custom_bounds = ([0.0, 8.0], [0.0, 8.0])

        self.global_optimum = [[7.0, 2.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        Zh1 = lambda x: 9.0 - x[0] - x[1]
        Zh2 = lambda x: (x[0] - 3.0) ** 2.0 + (x[1] - 2.0) ** 2.0 - 16.0
        Zh3 = lambda x: x[0] * x[1] - 14.0
        Zp = lambda x: 100.0 * (1.0 + x)

        return max(Zh1(x),
                   Zp(Zh2(x)) * sign(Zh2(x)),
                   Zp(Zh3(x)) * sign(Zh3(x)),
                   Zp(-x[0]) * sign(x[0]),
                   Zp(-x[1]) * sign(x[1]))


class Zirilli(Benchmark):

    r"""
    Zettl objective function.

    This class defines the Zirilli [1]_ global optimization problem. This is a
    unimodal minimization problem defined as follows:

    .. math::

        f_{\text{Zirilli}}(x) = 0.25x_1^4 - 0.5x_1^2 + 0.1x_1 + 0.5x_2^2

    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-10, 10]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = -0.3523` for :math:`x = [-1.0465, 0]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-10.0] * self.N, [10.0] * self.N))
        self.custom_bounds = ([-2.0, 2.0], [-2.0, 2.0])

        self.global_optimum = [[-1.0465, 0.0]]
        self.fglob = -0.35238603

    def fun(self, x, *args):
        self.nfev += 1

        return 0.25 * x[0] ** 4 - 0.5 * x[0] ** 2 + 0.1 * x[0] + 0.5 * x[1] ** 2
