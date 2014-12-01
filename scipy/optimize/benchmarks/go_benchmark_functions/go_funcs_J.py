# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import
import numpy as np
from numpy import (abs, arctan2, asarray, cos, exp, floor, log, log10,
                   arange, pi, prod, roll, sign, sin, sqrt, sum, where,
                   zeros, tan, tanh, dot)

from scipy.misc import factorial
from .go_benchmark import Benchmark


class JennrichSampson(Benchmark):

    """
    Jennrich-Sampson objective function.

    This class defines the Jennrich-Sampson global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{JennrichSampson}}(\\mathbf{x}) = \\sum_{i=1}^{10} \\left [2 + 2i - (e^{ix_1} + e^{ix_2}) \\right ]^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 124.3621824` for :math:`\\mathbf{x} = [0.257825, 0.257825]`.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-1.0] * self.N, [1.0] * self.N)

        self.global_optimum = [[0.257825, 0.257825]]
        self.custom_bounds = [(-1, 0.34), (-1, 0.34)]
        self.fglob = 124.3621824

    def fun(self, x, *args):
        self.nfev += 1

        i = arange(1, 11)
        return sum((2 + 2 * i - (exp(i * x[0]) + exp(i * x[1]))) ** 2)


class Judge(Benchmark):

    """
    Judge objective function.

    This class defines the Judge global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Judge}}(\\mathbf{x}) = \\sum_{i=1}^{20} \\left [ \\left (x_1 + A_i x_2 + B x_2^2 \\right ) - C_i \\right ]^2

    Where, in this exercise:

    .. math::

        \\begin{cases} C = [4.284, 4.149, 3.877, 0.533, 2.211, 2.389, 2.145,  3.231, 1.998, 1.379, 2.106, 1.428, 1.011, 2.179, 2.858, 1.388, 1.651, 1.593, 1.046, 2.152] \\\\
        A = [0.286, 0.973, 0.384, 0.276, 0.973, 0.543, 0.957, 0.948, 0.543, 0.797, 0.936, 0.889, 0.006, 0.828, 0.399, 0.617, 0.939, 0.784, 0.072, 0.889] \\\\
        B = [0.645, 0.585, 0.310, 0.058, 0.455, 0.779, 0.259, 0.202, 0.028, 0.099, 0.142, 0.296, 0.175, 0.180, 0.842, 0.039, 0.103, 0.620, 0.158, 0.704] \\end{cases}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 16.0817307` for :math:`\\mathbf{x} = [0.86479, 1.2357]`.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[0.86479, 1.2357]]
        self.custom_bounds = [(-2.0, 2.0), (-2.0, 2.0)]
        self.fglob = 16.0817307
        self.c = asarray([4.284, 4.149, 3.877, 0.533, 2.211, 2.389, 2.145,
                          3.231, 1.998, 1.379, 2.106, 1.428, 1.011, 2.179,
                          2.858, 1.388, 1.651, 1.593, 1.046, 2.152])

        self.a = asarray([0.286, 0.973, 0.384, 0.276, 0.973, 0.543, 0.957,
                          0.948, 0.543, 0.797, 0.936, 0.889, 0.006, 0.828,
                          0.399, 0.617, 0.939, 0.784, 0.072, 0.889])

        self.b = asarray([0.645, 0.585, 0.310, 0.058, 0.455, 0.779, 0.259,
                          0.202, 0.028, 0.099, 0.142, 0.296, 0.175, 0.180,
                          0.842, 0.039, 0.103, 0.620, 0.158, 0.704])

    def fun(self, x, *args):
        self.nfev += 1

        return sum(((x[0] + x[1] * self.a + (x[1] ** 2.0) * self.b) - self.c)
                    ** 2.0)
