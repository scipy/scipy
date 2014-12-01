# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import
import numpy as np
from numpy import (abs, arctan2, asarray, cos, exp, floor, log, log10,
                   arange, pi, prod, roll, sign, sin, sqrt, sum, where,
                   zeros, tan, tanh, dot)

from scipy.misc import factorial
from .go_benchmark import Benchmark


class Easom(Benchmark):

    """
    Easom objective function.

    This class defines the Easom global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Easom}}(\\mathbf{x}) = a - \\frac{a}{e^{b \\sqrt{\\frac{\\sum_{i=1}^{n} x_i^{2}}{n}}}} + e - e^{\\frac{\\sum_{i=1}^{n} \\cos\\left(c x_i\\right)}{n}}

    Where, in this exercise, :math:`a = 20, b = 0.2` and :math:`c = 2\\pi`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        a = 20.0
        b = 0.2
        c = 2 * pi

        return (-a * exp(-b * sqrt(sum(x ** 2) / self.N))
                - exp(sum(cos(c * x)) / self.N) + a + exp(1))


class Eckerle4(Benchmark):
    """
    Eckerle4 objective function.
    Eckerle, K., NIST (1979).
    Circular Interference Transmittance Study.
    """

    # TODO, this is a NIST regression standard dataset
    def __init__(self, dimensions=3):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0., 1., 10.],
                           [20, 20., 600.])
        self.global_optimum = [[1.5543827178, 4.0888321754, 4.5154121844e2]]
        self.fglob = 1.4635887487E-03

        self.a = asarray([1.5750000E-04, 1.6990000E-04, 2.3500000E-04,
                          3.1020000E-04, 4.9170000E-04, 8.7100000E-04,
                          1.7418000E-03, 4.6400000E-03, 6.5895000E-03,
                          9.7302000E-03, 1.4900200E-02, 2.3731000E-02,
                          4.0168300E-02, 7.1255900E-02, 1.2644580E-01,
                          2.0734130E-01, 2.9023660E-01, 3.4456230E-01,
                          3.6980490E-01, 3.6685340E-01, 3.1067270E-01,
                          2.0781540E-01, 1.1643540E-01, 6.1676400E-02,
                          3.3720000E-02, 1.9402300E-02, 1.1783100E-02,
                          7.4357000E-03, 2.2732000E-03, 8.8000000E-04,
                          4.5790000E-04, 2.3450000E-04, 1.5860000E-04,
                          1.1430000E-04, 7.1000000E-05])

        self.b = asarray([4.0000000E+02, 4.0500000E+02, 4.1000000E+02,
                          4.1500000E+02, 4.2000000E+02, 4.2500000E+02,
                          4.3000000E+02, 4.3500000E+02, 4.3650000E+02,
                          4.3800000E+02, 4.3950000E+02, 4.4100000E+02,
                          4.4250000E+02, 4.4400000E+02, 4.4550000E+02,
                          4.4700000E+02, 4.4850000E+02, 4.5000000E+02,
                          4.5150000E+02, 4.5300000E+02, 4.5450000E+02,
                          4.5600000E+02, 4.5750000E+02, 4.5900000E+02,
                          4.6050000E+02, 4.6200000E+02, 4.6350000E+02,
                          4.6500000E+02, 4.7000000E+02, 4.7500000E+02,
                          4.8000000E+02, 4.8500000E+02, 4.9000000E+02,
                          4.9500000E+02, 5.0000000E+02])

    def fun(self, x, *args):
        self.nfev += 1

        vec = x[0] / x[1] * exp(-(self.b - x[2]) ** 2 / (2 * x[1] ** 2))
        return sum((self.a - vec) ** 2)


class EggCrate(Benchmark):

    """
    Egg Crate objective function.

    This class defines the Egg Crate global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{EggCrate}}(\\mathbf{x}) = x_1^2 + x_2^2 + 25 \\left[ \\sin^2(x_1) + \\sin^2(x_2) \\right]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.0] * self.N, [5.0] * self.N)

        self.global_optimum = [[0.0, 0.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return x[0] ** 2 + x[1] ** 2 + 25 * (sin(x[0]) ** 2 + sin(x[1]) ** 2)


class EggHolder(Benchmark):

    """
    Egg Holder objective function.

    This class defines the Egg Holder global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{EggHolder}}(\\mathbf{x}) = - x_{1} \\sin\\left(\\sqrt{\\lvert{x_{1} - x_{2} -47}\\rvert}\\right) - \\left(x_{2} + 47\\right) \\sin\\left(\\sqrt{\\left|{\\frac{1}{2} x_{1} + x_{2} + 47}\\right|}\\right)

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-512, 512]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -959.640662711` for :math:`\\mathbf{x} = [512, 404.2319]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-512.1] * self.N,
                           [512.0] * self.N)

        self.global_optimum = [[512.0, 404.2319]]
        self.fglob = -959.640662711
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        vec = (-(x[1:] + 47) * sin(sqrt(abs(x[1:] + x[:-1] / 2. + 47)))
               - x[:-1] * sin(sqrt(abs(x[:-1] - (x[1:] + 47)))))
        return sum(vec)


class ElAttarVidyasagarDutta(Benchmark):

    """
    El-Attar-Vidyasagar-Dutta objective function.

    This class defines the El-Attar-Vidyasagar-Dutta function global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{ElAttarVidyasagarDutta}}(\\mathbf{x}) = (x_1^2 + x_2 - 10)^2 + (x_1 + x_2^2 - 7)^2 + (x_1^2 + x_2^3 - 1)^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 1.712780354` for :math:`\\mathbf{x} = [3.40918683, -2.17143304]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = [(-4, 4), (-4, 4)]

        self.global_optimum = [[3.40918683, -2.17143304]]
        self.fglob = 1.712780354

    def fun(self, x, *args):
        self.nfev += 1

        return ((x[0] ** 2 + x[1] - 10) ** 2 + (x[0] + x[1] ** 2 - 7) ** 2
                + (x[0] ** 2 + x[1] ** 3 - 1) ** 2)


class Exp2(Benchmark):

    """
    Exp2 objective function.

    This class defines the Exp2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Exp2}}(\\mathbf{x}) = \\sum_{i=0}^9 \\left ( e^{-ix_1/10} - 5e^{-ix_2/10} -e^{-i/10} + 5e^{-i} \\right )^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 20]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = [1, 10.]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N, [20.0] * self.N)
        self.custom_bounds = [(0, 2), (0, 20)]

        self.global_optimum = [[1.0, 10.]]
        self.fglob = 0

    def fun(self, x, *args):
        self.nfev += 1

        i = arange(10.)
        vec = (exp(-i * x[0] / 10.) - 5 * exp(-i * x[1] / 10.) - exp(-i / 10)
               + 5 * exp(-i)) ** 2

        return sum(vec)


class Exponential(Benchmark):

    """
    Exponential objective function.

    This class defines the Exponential global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Exponential}}(\\mathbf{x}) = -e^{-0.5 \\sum_{i=1}^n x_i^2}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = -1` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-1.0] * self.N, [1.0] * self.N)

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = -1.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return -exp(-0.5 * sum(x ** 2.0))
