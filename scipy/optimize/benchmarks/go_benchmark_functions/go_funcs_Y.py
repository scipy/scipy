# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import
import numpy as np
from numpy import (abs, arctan2, asarray, cos, exp, floor, log, log10,
                   arange, pi, prod, roll, sign, sin, sqrt, sum, where,
                   zeros, tan, tanh, dot)

from scipy.misc import factorial
from .go_benchmark import Benchmark


class YaoLiu04(Benchmark):

    """
    Yao-Liu 4 objective function.

    This class defines the Yao-Liu function 4 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{YaoLiu04}}(\\mathbf{x}) = {max}_i \\left\{ \\left | x_i \\right | , 1 \\leq i \\leq n \\right\}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return abs(x.max())


class YaoLiu09(Benchmark):

    """
    Yao-Liu 9 objective function.

    This class defines the Yao-Liu function 9 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{YaoLiu09}}(\\mathbf{x}) = \\sum_{i=1}^n \\left [ x_i^2 - 10 \\cos(2 \\pi x_i ) + 10 \\right ]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5.12, 5.12]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.12] * self.N, [5.12] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum(x ** 2.0 - 10.0 * cos(2 * pi * x) + 10)
