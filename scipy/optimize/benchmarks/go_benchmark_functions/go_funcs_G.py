# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import
import numpy as np
from numpy import (abs, arctan2, asarray, cos, exp, floor, log, log10,
                   arange, pi, prod, roll, sign, sin, sqrt, sum, where,
                   zeros, tan, tanh, dot)

from scipy.misc import factorial
from .go_benchmark import Benchmark


class Gear(Benchmark):

    """
    Gear objective function.

    This class defines the Gear global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Gear}}(\\mathbf{x}) = \\left \\{ \\frac{1.0}{6.931} - \\frac{\\lfloor x_1\\rfloor \\lfloor x_2 \\rfloor } {\\lfloor x_3 \\rfloor \\lfloor x_4 \\rfloor } \\right\\}^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [12, 60]` for :math:`i=1,...,4`.

    *Global optimum*: :math:`f(x_i) = 2.7 \\cdot 10^{-12}` for :math:`\\mathbf{x} = [16, 19, 43, 49]`, where the various
    :math:`x_i` may be permuted.

    """

    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([12.0] * self.N, [60.0] * self.N)
        self.global_optimum = [[16, 19, 43, 49]]
        self.fglob = 2.7e-12

    def fun(self, x, *args):
        self.nfev += 1

        return (1. / 6.931
                - floor(x[0]) * floor(x[1]) / floor(x[2]) / floor(x[3])) ** 2


class Giunta(Benchmark):

    """
    Giunta objective function.

    This class defines the Giunta global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Giunta}}(\\mathbf{x}) = 0.6 + \\sum_{i=1}^{n} \\left[\\sin^{2}\\left(1 - \\frac{16}{15} x_i\\right) - \\frac{1}{50} \\sin\\left(4 - \\frac{64}{15} x_i\\right) - \\sin\\left(1 - \\frac{16}{15} x_i\\right)\\right]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0.06447042053690566` for :math:`\\mathbf{x} = [0.4673200277395354, 0.4673200169591304]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-1.0] * self.N, [1.0] * self.N)

        self.global_optimum = [[0.4673200277395354, 0.4673200169591304]]
        self.fglob = 0.06447042053690566

    def fun(self, x, *args):
        self.nfev += 1

        arg = 16 * x / 15.0 - 1
        return 0.6 + sum(sin(arg) + sin(arg) ** 2 + sin(4 * arg) / 50)


class GoldsteinPrice(Benchmark):

    """
    Goldstein-Price objective function.

    This class defines the Goldstein-Price global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{GoldsteinPrice}}(\\mathbf{x}) = \\left[ 1+(x_1+x_2+1)^2(19-14x_1+3x_1^2-14x_2+6x_1x_2+3x_2^2) \\right] \\left[ 30+(2x_1-3x_2)^2(18-32x_1+12x_1^2+48x_2-36x_1x_2+27x_2^2) \\right]


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-2, 2]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 3` for :math:`\\mathbf{x} = [0, -1]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-2.0] * self.N, [2.0] * self.N)

        self.global_optimum = [[0., -1.]]
        self.fglob = 3.0

    def fun(self, x, *args):
        self.nfev += 1

        a = 1 + (x[0] + x[1] + 1) ** 2 * \
            (19 - 14 * x[0] + 3 * x[0] ** 2 -
             14 * x[1] + 6 * x[0] * x[1] + 3 * x[1] ** 2)
        b = 30 + (2 * x[0] - 3 * x[1]) ** 2 * \
            (18 - 32 * x[0] + 12 * x[0] ** 2 +
             48 * x[1] - 36 * x[0] * x[1] + 27 * x[1] ** 2)
        return a * b


class Griewank(Benchmark):

    """
    Griewank objective function.

    This class defines the Griewank global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Griewank}}(\\mathbf{x}) = \\frac{1}{4000}\\sum_{i=1}^n x_i^2 - \\prod_{i=1}^n\\cos\\left(\\frac{x_i}{\\sqrt{i}}\\right) + 1

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-600, 600]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-600.0] * self.N,
                           [600.0] * self.N)
        self.custom_bounds = [(-50, 50), (-50, 50)]

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        i = arange(1., np.size(x))
        return sum(x ** 2 / 4000) - prod(cos(x / sqrt(i))) + 1


class Gulf(Benchmark):

    """
    Gulf objective function.

    This class defines the Gulf global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Gulf}}(\\mathbf{x}) = \\sum_{i=1}^99 \\left( e^{-\\frac{\\lvert y_i - x_2 \\rvert^{x_3}}{x_1}    }  - t_i \\right)

    Where, in this exercise:

    .. math::

       t_i = i/100 \\\\
       y_i = 25 + [-50 \\log(t_i)]^{2/3}


    Here, :math:`x_i \\in [0, 60]` for :math:`i=1,2,3`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [50, 25, 1.5]`

    """

    def __init__(self, dimensions=3):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N, [50.0] * self.N)

        self.global_optimum = [[50.0, 25.0, 1.5]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        m = 99.
        i = arange(1., m + 1)
        y = 25 + (-50 * log(i / 100.)) ** (2 / 3.)
        vec = (exp(-((abs(y - x[1])) ** x[2] / x[0])) - i / 100.)
        return sum(vec ** 2)
