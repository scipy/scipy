# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import
import numpy as np
from numpy import (abs, arctan2, asarray, cos, exp, floor, log, log10,
                   arange, pi, prod, roll, sign, sin, sqrt, sum, where,
                   zeros, tan, tanh, dot)

from scipy.misc import factorial
from .go_benchmark import Benchmark


class VenterSobiezcczanskiSobieski(Benchmark):

    """
    Venter Sobiezcczanski-Sobieski objective function.

    This class defines the Venter Sobiezcczanski-Sobieski global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{VenterSobiezcczanskiSobieski}}(\\mathbf{x}) = x_1^2 - 100 \\cos^2(x_1) - 100 \\cos(x_1^2/30) + x_2^2 - 100 \\cos^2(x_2) - 100 \\cos(x_2^2/30)

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-50, 50]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -400` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-50.0] * self.N, [50.0] * self.N)
        self.custom_bounds = ([-10, 10], [-10, 10])

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = -400

    def fun(self, x, *args):
        self.nfev += 1

        u = x[0] ** 2.0 - 100.0 * cos(x[0]) ** 2.0
        v = -100.0 * cos(x[0] ** 2.0 / 30.0) + x[1] ** 2.0
        w = - 100.0 * cos(x[1]) ** 2.0 - 100.0 * cos(x[1] ** 2.0 / 30.0)
        return u + v + w


class Vincent(Benchmark):

    """
    Vincent objective function.

    This class defines the Vincent global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Vincent}}(\\mathbf{x}) = - \\sum_{i=1}^{n} \\sin(10 \\log(x))

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0.25, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = -n` for :math:`x_i = 7.70628098` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.25] * self.N, [10.0] * self.N)

        self.global_optimum = [[7.70628098 for _ in range(self.N)]]
        self.fglob = -float(self.N)
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return -sum(sin(10.0 * log(x)))
