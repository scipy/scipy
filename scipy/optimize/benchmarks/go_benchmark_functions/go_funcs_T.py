# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import
import numpy as np
from numpy import (abs, arctan2, asarray, cos, exp, floor, log, log10,
                   arange, pi, prod, roll, sign, sin, sqrt, sum, where,
                   zeros, tan, tanh, dot)

from scipy.misc import factorial
from .go_benchmark import Benchmark


class TestTubeHolder(Benchmark):
    """
    TestTubeHolder objective function.

    This class defines the TestTubeHolder global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{TestTubeHolder}}(\\mathbf{x}) = - 4 \\left | {e^{\left|{\\cos 
        \\left(\\frac{1}{200} x_{1}^{2} + \\frac{1}{200} x_{2}^{2}\\right)}
        \\right|}\\sin\\left(x_{1}\\right) \\cos\\left(x_{2}\\right)}\\right|

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in 
    [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -10.872299901558` for :math:`\\mathbf{x}=
    [-\\pi/2, 0]`
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[-pi / 2, 0.0]]
        self.fglob = -10.87229990155800

    def fun(self, x, *args):
        self.nfev += 1

        u = sin(x[0]) * cos(x[1])
        v = (x[0] ** 2 + x[1] ** 2) / 200
        return -4 * abs(u * exp(abs(cos(v))))


class Thurber(Benchmark):

    """
    Thurber objective function.

    """

    # TODO, this is a NIST regression standard dataset
    def __init__(self, dimensions=7):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip(
            [500., 500., 100., 10., 0.1, 0.1, 0.],
            [2000., 2000., 1000., 150., 2., 1., 0.2])
        self.global_optimum = [[1.288139680e3, 1.4910792535e3, 5.8323836877e2,
                                75.416644291, 0.96629502864, 0.39797285797,
                                4.9727297349e-2]]
        self.fglob = 5642.7082397
        self.a = asarray([80.574, 84.248, 87.264, 87.195, 89.076, 89.608,
                          89.868, 90.101, 92.405, 95.854, 100.696, 101.06,
                          401.672, 390.724, 567.534, 635.316, 733.054, 759.087,
                          894.206, 990.785, 1090.109, 1080.914, 1122.643,
                          1178.351, 1260.531, 1273.514, 1288.339, 1327.543,
                          1353.863, 1414.509, 1425.208, 1421.384, 1442.962,
                          1464.350, 1468.705, 1447.894, 1457.628])
        self.b = np.asarray([-3.067, -2.981, -2.921, -2.912, -2.840, -2.797,
                             -2.702, -2.699, -2.633, -2.481, -2.363, -2.322,
                             -1.501, -1.460, -1.274, -1.212, -1.100, -1.046,
                             -0.915, -0.714, -0.566, -0.545, -0.400, -0.309,
                             -0.109, -0.103, 0.010, 0.119, 0.377, 0.790, 0.963,
                             1.006, 1.115, 1.572, 1.841, 2.047, 2.200])

    def fun(self, x, *args):
        self.nfev += 1

        vec = x[0] + x[1] * self.b + x[2] * self.b ** 2 + x[3] * self.b ** 3
        vec /= 1 + x[4] * self.b + x[5] * self.b ** 2 + x[6] * self.b ** 3

        return sum((self.a - vec) ** 2)


class Treccani(Benchmark):

    """
    Treccani objective function.

    This class defines the Treccani global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Treccani}}(\\mathbf{x}) = x_1^4 + 4x_1^3 + 4x_1^2 + x_2^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in
    [-5, 5]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [-2, 0]`
    or :math:`\\mathbf{x} = [0, 0]`.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.0] * self.N, [5.0] * self.N)
        self.custom_bounds = [(-2, 2), (-2, 2)]

        self.global_optimum = [[-2.0, 0.0]]
        self.fglob = 0

    def fun(self, x, *args):
        self.nfev += 1

        return x[0] ** 4 + 4.0 * x[0] ** 3 + 4.0 * x[0] ** 2 + x[1] ** 2


class Trefethen(Benchmark):

    """
    Trefethen objective function.

    This class defines the Trefethen global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Trefethen}}(\\mathbf{x}) = 0.25 x_{1}^{2} + 0.25 x_{2}^{2} + e^{\\sin\\left(50 x_{1}\\right)} - \\sin\\left(10 x_{1} + 10 x_{2}\\right) + \\sin\\left(60 e^{x_{2}}\\right) + \\sin\\left[70 \\sin\\left(x_{1}\\right)\\right] + \\sin\\left[\\sin\\left(80 x_{2}\\right)\\right]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -3.3068686474` for :math:`\\mathbf{x} = [-0.02440307923, 0.2106124261]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.custom_bounds = [(-5, 5), (-5, 5)]

        self.global_optimum = [[-0.02440307923, 0.2106124261]]
        self.fglob = -3.3068686474

    def fun(self, x, *args):
        self.nfev += 1

        val = 0.25 * x[0] ** 2 + 0.25 * x[1] ** 2
        val += exp(sin(50. * x[0])) - sin(10 * x[0] + 10 * x[1])
        val += sin(60 * exp(x[1]))
        val += sin(70 * sin(x[0]))
        val += sin(sin(80 * x[1]))
        return val


class ThreeHumpCamel(Benchmark):

    """
    Three Hump Camel objective function.

    This class defines the Three Hump Camel global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{ThreeHumpCamel}}(\\mathbf{x}) = 2x_1^2 - 1.05x_1^4 + \\frac{x_1^6}{6} + x_1x_2 + x_2^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [0, 0]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.0] * self.N, [5.0] * self.N)
        self.custom_bounds = [(-2, 2), (-1.5, 1.5)]

        self.global_optimum = [[0.0, 0.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return (2.0 * x[0] ** 2.0 - 1.05 * x[0] ** 4.0 + x[0] ** 6 / 6.0
                + x[0] * x[1] + x[1] ** 2.0)


class Trid(Benchmark):

    """
    Trid objective function.

    This class defines the Trid global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Trid}}(\\mathbf{x}) = \\sum_{i=1}^{n}(x_i - 1)^2 - \\sum_{i=2}^{n} x_ix_{i-1}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-20, 20]` for :math:`i=1,...,6`.

    *Global optimum*: :math:`f(x_i) = -50` for :math:`\\mathbf{x} = [6, 10, 12, 12, 10, 6]`

    """

    def __init__(self, dimensions=6):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N, [20.0] * self.N)

        self.global_optimum = [[6, 10, 12, 12, 10, 6]]
        self.fglob = -50.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum((x - 1.0) ** 2.0) - sum(x[1:] * x[:-1])


class Trigonometric01(Benchmark):

    """
    Trigonometric 1 objective function.

    This class defines the Trigonometric 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Trigonometric01}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\left [n - \\sum_{j=1}^{n} \\cos(x_j) + i \\left(1 - cos(x_i) - sin(x_i) \\right ) \\right]^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, \\pi]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N, [pi] * self.N)

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        i = arange(1.0, self.N + 1)
        return sum((self.N - sum(cos(x) + i * (1 - cos(x) - sin(x)))) ** 2.0)


class Trigonometric02(Benchmark):

    """
    Trigonometric 2 objective function.

    This class defines the Trigonometric 2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Trigonometric2}}(\\mathbf{x}) = 1 + \\sum_{i=1}^{n} 8 \\sin^2 \\left[7(x_i - 0.9)^2 \\right] + 6 \\sin^2 \\left[14(x_i - 0.9)^2 \\right] + (x_i - 0.9)^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-500, 500]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 1` for :math:`x_i = 0.9` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-500.0] * self.N,
                           [500.0] * self.N)
        self.custom_bounds = [(0, 2), (0, 2)]

        self.global_optimum = [[0.9 for _ in range(self.N)]]
        self.fglob = 1.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        vec = (8 * sin(7 * (x - 0.9) ** 2) ** 2
               + 6 * sin(14 * (x - 0.9) ** 2) ** 2
               + (x - 0.9) ** 2)
        return 1.0 + sum(vec)


class Tripod(Benchmark):

    """
    Tripod objective function.

    This class defines the Tripod global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Tripod}}(\\mathbf{x}) = p(x_2) \\left[1 + p(x_1) \\right] + \\lvert x_1 + 50p(x_2) \\left[1 - 2p(x_1) \\right] \\rvert + \\lvert x_2 + 50\\left[1 - 2p(x_2)\\right] \\rvert

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [0, -50]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)

        self.global_optimum = [[0.0, -50.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        p1 = float(x[0] >= 0)
        p2 = float(x[1] >= 0)

        return (p2 * (1.0 + p1) + abs(x[0] + 50.0 * p2 * (1.0 - 2.0 * p1))
                + abs(x[1] + 50.0 * (1.0 - 2.0 * p2)))
