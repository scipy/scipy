# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import

from numpy import abs, sum, sin, cos, pi, sqrt
from .go_benchmark import Benchmark


class Ursem01(Benchmark):

    r"""
    Ursem 1 objective function.

    This class defines the Ursem 1 global optimization problem. This is a
    unimodal minimization problem defined as follows:

    .. math::

        f_{\text{Ursem01}}(x) = - \sin(2x_1 - 0.5 \pi) - 3 \cos(x_2) - 0.5 x_1

    with :math:`x_1 \in [-2.5, 3]` and :math:`x_2 \in [-2, 2]`.

    *Global optimum*: :math:`f(x) = -4.81681406371` for
    :math:`x = [1.69714, 0.0]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = [(-2.5, 3.0), (-2.0, 2.0)]

        self.global_optimum = [[1.69714, 0.0]]
        self.fglob = -4.81681406371

    def fun(self, x, *args):
        self.nfev += 1

        return (-sin(2 * x[0] - 0.5 * pi) - 3.0 * cos(x[1]) - 0.5 * x[0])


class Ursem03(Benchmark):

    r"""
    Ursem 3 objective function.

    This class defines the Ursem 3 global optimization problem. This is a
    multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Ursem03}}(x) = - \sin(2.2 \pi x_1 + 0.5 \pi) 
                                \frac{2 - \lvert x_1 \rvert}{2}
                                \frac{3 - \lvert x_1 \rvert}{2}
                                - \sin(2.2 \pi x_2 + 0.5 \pi)
                                \frac{2 - \lvert x_2 \rvert}{2}
                                \frac{3 - \lvert x_2 \rvert}{2}

    with :math:`x_1 \in [-2, 2]`, :math:`x_2 \in [-1.5, 1.5]`.

    *Global optimum*: :math:`f(x) = -3` for :math:`x = [0, 0]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = [(-2, 2), (-1.5, 1.5)]

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = -3.0

    def fun(self, x, *args):
        self.nfev += 1

        u = -(sin(2.2 * pi * x[0] + 0.5 * pi)
              * ((2.0 - abs(x[0])) / 2.0) * ((3.0 - abs(x[0])) / 2))
        v = -(sin(2.2 * pi * x[1] + 0.5 * pi)
              * ((2.0 - abs(x[1])) / 2) * ((3.0 - abs(x[1])) / 2))
        return u + v


class Ursem04(Benchmark):

    r"""
    Ursem 4 objective function.

    This class defines the Ursem 4 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Ursem04}}(x) = -3 \sin(0.5 \pi x_1 + 0.5 \pi)
                                \frac{2 - \sqrt{x_1^2 + x_2 ^ 2}}{4}

    with :math:`x_i \in [-2, 2]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = -1.5` for :math:`x = [0, 0]` for
    :math:`i = 1, 2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-2.0] * self.N, [2.0] * self.N)

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = -1.5

    def fun(self, x, *args):
        self.nfev += 1

        return (-3 * sin(0.5 * pi * x[0] + 0.5 * pi)
                * (2 - sqrt(x[0] ** 2 + x[1] ** 2)) / 4)


class UrsemWaves(Benchmark):

    r"""
    Ursem Waves objective function.

    This class defines the Ursem Waves global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{UrsemWaves}}(x) = -0.9x_1^2 + (x_2^2 - 4.5x_2^2)x_1x_2
                                   + 4.7 \cos \left[ 2x_1 - x_2^2(2 + x_1)
                                   \right ] \sin(2.5 \pi x_1)

    with :math:`x_1 \in [-0.9, 1.2]`, :math:`x_2 \in [-1.2, 1.2]`.

    *Global optimum*: :math:`f(x) = -8.5536` for :math:`x = [1.2, 1.2]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = [(-0.9, 1.2), (-1.2, 1.2)]

        self.global_optimum = [[1.2 for _ in range(self.N)]]
        self.fglob = -8.5536

    def fun(self, x, *args):
        self.nfev += 1

        u = -0.9 * x[0] ** 2
        v = (x[1] ** 2 - 4.5 * x[1] ** 2) * x[0] * x[1]
        w = 4.7 * cos(2 * x[0] - x[1] ** 2 * (2 + x[0])) * sin(2.5 * pi * x[0])
        return (u + v + w)
