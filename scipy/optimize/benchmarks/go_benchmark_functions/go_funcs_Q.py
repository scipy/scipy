# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import
import numpy as np
from numpy import (abs, arctan2, asarray, cos, exp, floor, log, log10,
                   arange, pi, prod, roll, sign, sin, sqrt, sum, where,
                   zeros, tan, tanh, dot)

from scipy.misc import factorial
from .go_benchmark import Benchmark


class Qing(Benchmark):

    """
    Qing objective function.

    This class defines the Qing global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Qing}}(\\mathbf{x}) = \\sum_{i=1}^{n} (x_i^2 - i)^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-500, 500]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = \\pm \\sqrt(i)` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-500.0] * self.N,
                           [500.0] * self.N)
        self.custom_bounds = [(-2, 2), (-2, 2)]
        self.global_optimum = [[sqrt(_) for _ in range(1, self.N + 1)]]
        self.fglob = 0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        i = arange(1, self.N + 1)
        return sum((x ** 2.0 - i) ** 2.0)


class Quadratic(Benchmark):

    """
    Quadratic objective function.

    This class defines the Quadratic global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Quadratic}}(\\mathbf{x}) = -3803.84 - 138.08x_1 - 232.92x_2 + 128.08x_1^2 + 203.64x_2^2 + 182.25x_1x_2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -3873.72418` for :math:`\\mathbf{x} = [0.19388, 0.48513]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.custom_bounds = [(0, 1), (0, 1)]
        self.global_optimum = [[0.19388, 0.48513]]
        self.fglob = -3873.72418
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return (-3803.84 - 138.08 * x[0] - 232.92 * x[1] + 128.08 * x[0] ** 2.0
                + 203.64 * x[1] ** 2.0 + 182.25 * x[0] * x[1])


class Quintic(Benchmark):

    """
    Quintic objective function.

    This class defines the Quintic global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Quintic}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\left|{x_{i}^{5} - 3 x_{i}^{4} + 4 x_{i}^{3} + 2 x_{i}^{2} - 10 x_{i} -4}\\right|


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = -1` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.custom_bounds = [(-2, 2), (-2, 2)]

        self.global_optimum = [[-1.0 for _ in range(self.N)]]
        self.fglob = 0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum(abs(x ** 5 - 3 * x ** 4 + 4 * x ** 3 + 2 * x ** 2
                       - 10 * x - 4))
