# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import
import numpy as np
from numpy import (abs, arctan2, asarray, cos, exp, floor, log, log10,
                   arange, pi, prod, roll, sign, sin, sqrt, sum, where,
                   zeros, tan, tanh, dot)

from scipy.misc import factorial
from .go_benchmark import Benchmark


class XinSheYang01(Benchmark):

    """
    Xin-She Yang 1 objective function.

    This class defines the Xin-She Yang 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{XinSheYang01}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\epsilon_i \\lvert x_i \\rvert^i


    The variable :math:`\\epsilon_i, (i=1,...,n)` is a random variable uniformly distributed in :math:`[0, 1]`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.0] * self.N, [5.0] * self.N)
        self.custom_bounds = ([-2, 2], [-2, 2])

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        i = arange(1.0, self.N + 1.0)
        return sum(np.random.random() * (abs(x) ** i))


class XinSheYang02(Benchmark):

    """
    Xin-She Yang 2 objective function.

    This class defines the Xin-She Yang 2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{XinSheYang02}}(\\mathbf{x}) = \\frac{\\sum_{i=1}^{n} \\lvert{x_{i}}\\rvert}{e^{\\sum_{i=1}^{n} \\sin\\left(x_{i}^{2.0}\\right)}}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-2\\pi, 2\\pi]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-2 * pi] * self.N,
                           [2 * pi] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum(abs(x)) * exp(-sum(sin(x ** 2.0)))


class XinSheYang03(Benchmark):

    """
    Xin-She Yang 3 objective function.

    This class defines the Xin-She Yang 3 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{XinSheYang03}}(\\mathbf{x}) = e^{-\\sum_{i=1}^{n} (x_i/\\beta)^{2m}} - 2e^{-\\sum_{i=1}^{n} x_i^2} \\prod_{i=1}^{n} \\cos^2(x_i)


    Where, in this exercise, :math:`\\beta = 15` and :math:`m = 3`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-20, 20]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = -1` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-20.0] * self.N, [20.0] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = -1.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        beta, m = 15.0, 5.0
        u = sum((x / beta) ** (2 * m))
        v = sum(x ** 2)
        w = prod(cos(x) ** 2)

        return exp(-u) - 2 * exp(-v) * w


class XinSheYang04(Benchmark):

    """
    Xin-She Yang 4 objective function.

    This class defines the Xin-She Yang 4 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{XinSheYang04}}(\\mathbf{x}) = \\left[ \\sum_{i=1}^{n} \\sin^2(x_i) - e^{-\\sum_{i=1}^{n} x_i^2} \\right ] e^{-\\sum_{i=1}^{n} \\sin^2 \\sqrt{ \\lvert x_i \\rvert }}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = -1` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = -1.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        u = sum(sin(x) ** 2)
        v = sum(x ** 2)
        w = sum(sin(sqrt(abs(x))) ** 2)
        return (u - exp(-v)) * exp(-w)


class Xor(Benchmark):

    """
    Xor objective function.

    This class defines the Xor global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Xor}}(\\mathbf{x}) = \\left[ 1 + \\exp \\left( - \\frac{x_7}{1 + \\exp(-x_1 - x_2 - x_5)} - \\frac{x_8}{1 + \\exp(-x_3 - x_4 - x_6)} - x_9 \\right ) \\right ]^{-2} \\\\
       + \\left [ 1 + \\exp \\left( -\\frac{x_7}{1 + \\exp(-x_5)} - \\frac{x_8}{1 + \\exp(-x_6)} - x_9 \\right ) \\right] ^{-2} \\\\
       + \\left [1 - \\left\\{1 + \\exp \\left(-\\frac{x_7}{1 + \\exp(-x_1 - x_5)} - \\frac{x_8}{1 + \\exp(-x_3 - x_6)} - x_9 \\right ) \\right\\}^{-1} \\right ]^2 \\\\
       + \\left [1 - \\left\\{1 + \\exp \\left(-\\frac{x_7}{1 + \\exp(-x_2 - x_5)} - \\frac{x_8}{1 + \\exp(-x_4 - x_6)} - x_9 \\right ) \\right\\}^{-1} \\right ]^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,...,9`.

    *Global optimum*: :math:`f(x_i) = 0.9597588` for :math:`\\mathbf{x} = [1, -1, 1, -1, -1, 1, 1, -1, 0.421134]`

    """

    def __init__(self, dimensions=9):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-1.0] * self.N, [1.0] * self.N)

        self.global_optimum = [[1.0, -1.0, 1.0,
                               -1.0, -1.0, 1.0, 1.0, -1.0, 0.421134]]
        self.fglob = 0.9597588

    def fun(self, x, *args):
        self.nfev += 1

        F11 = x[6] / (1.0 + exp(-x[0] - x[1] - x[4]))
        F12 = x[7] / (1.0 + exp(-x[2] - x[3] - x[5]))
        F1 = (1.0 + exp(-F11 - F12 - x[8])) ** (-2)
        F21 = x[6] / (1.0 + exp(-x[4]))
        F22 = x[7] / (1.0 + exp(-x[5]))
        F2 = (1.0 + exp(-F21 - F22 - x[8])) ** (-2)
        F31 = x[6] / (1.0 + exp(-x[0] - x[4]))
        F32 = x[7] / (1.0 + exp(-x[2] - x[5]))
        F3 = (1.0 - (1.0 + exp(-F31 - F32 - x[8])) ** (-1)) ** 2
        F41 = x[6] / (1.0 + exp(-x[1] - x[4]))
        F42 = x[7] / (1.0 + exp(-x[3] - x[5]))
        F4 = (1.0 - (1.0 + exp(-F41 - F42 - x[8])) ** (-1)) ** 2

        return F1 + F2 + F3 + F4
