# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import
import numpy as np
from numpy import (abs, arctan2, asarray, cos, exp, floor, log, log10,
                   arange, pi, prod, roll, sign, sin, sqrt, sum, where,
                   zeros, tan, tanh, dot)

from scipy.misc import factorial
from .go_benchmark import Benchmark


class FreudensteinRoth(Benchmark):

    """
    FreudensteinRoth objective function.

    This class defines the Freudenstein & Roth global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{FreudensteinRoth}}(\\mathbf{x}) =  \\left\{x_1 - 13 + \\left[(5 - x_2)x_2 - 2 \\right] x_2 \\right\}^2 + \\left \{x_1 - 29 + \\left[(x_2 + 1)x_2 - 14 \\right] x_2 \\right\}^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [5, 4]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.custom_bounds = [(-3, 3), (-5, 5)]

        self.global_optimum = [[5.0, 4.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        f1 = (-13.0 + x[0] + ((5.0 - x[1]) * x[1] - 2.0) * x[1]) ** 2
        f2 = (-29.0 + x[0] + ((x[1] + 1.0) * x[1] - 14.0) * x[1]) ** 2

        return f1 + f2
