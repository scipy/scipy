# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import
import numpy as np
from numpy import (abs, arctan2, asarray, cos, exp, floor, log, log10,
                   arange, pi, prod, roll, sign, sin, sqrt, sum, where,
                   zeros, tan, tanh, dot)

from scipy.misc import factorial
from .go_benchmark import Benchmark


class CarromTable(Benchmark):

    """
    CarromTable objective function.

    The CarromTable [1]_ global optimization problem is a multimodal
    minimization problem defined as follows:

    .. math::

        f_{\text{CarromTable}}(\mathbf{x}) = - \frac{1}{30}\left(\cos(x_1)
        cos(x_2) e^{\left|1 - \frac{\sqrt{x_1^2 + x_2^2}}{\pi}\right|}\right)^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \in [-10,
    10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -24.15681551650653` for :math:`x_i = \pm
    9.646157266348881` for :math:`i=1,...,n`

    ..[1] S. K. Mishra, Global Optimization By Differential Evolution and
     Particle Swarm Methods: Evaluation On Some Benchmark Functions, Munich
     Research Papers in Economics

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.global_optimum = [(9.646157266348881, 9.646134286497169),
                               (-9.646157266348881, 9.646134286497169),
                               (9.646157266348881, -9.646134286497169),
                               (-9.646157266348881, -9.646134286497169)]
        self.fglob = -24.15681551650653

    def fun(self, x, *args):
        self.nfev += 1

        u = cos(x[0]) * cos(x[1])
        v = sqrt(x[0] ** 2 + x[1] ** 2)
        return -((u * exp(abs(1 - v / pi))) ** 2) / 30


class Chichinadze(Benchmark):

    """
    Chichinadze objective function.

    This class defines the Chichinadze global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Chichinadze}}(\\mathbf{x}) = x_{1}^{2} - 12 x_{1} + 8 \\sin\\left(\\frac{5}{2} \\pi x_{1}\\right) + 10 \\cos\\left(\\frac{1}{2} \\pi x_{1}\\right) + 11 - 0.2 \\frac{\\sqrt{5}}{e^{\\frac{1}{2} \\left(x_{2} -0.5\\right)^{2}}}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-30, 30]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -42.94438701899098` for :math:`\\mathbf{x} = [6.189866586965680, 0.5]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-30.0] * self.N, [30.0] * self.N)
        self.custom_bounds = [(-10, 10), (-10, 10)]

        self.global_optimum = [[6.189866586965680, 0.5]]
        self.fglob = -42.94438701899098

    def fun(self, x, *args):
        self.nfev += 1

        return (x[0] ** 2 - 12 * x[0] + 11 + 10 * cos(pi * x[0] / 2)
                + 8 * sin(5 * pi * x[0] / 2)
                - 1.0 / sqrt(5) * exp(-((x[1] - 0.5) ** 2) / 2))


class Cigar(Benchmark):

    """
    Cigar objective function.

    This class defines the Cigar global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Cigar}}(\\mathbf{x}) = x_1^2 + 10^6\\sum_{i=2}^{n} x_i^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = [(-5, 5), (-5, 5)]

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return x[0] ** 2 + 1e6 * sum(x[1:] ** 2)


class Cola(Benchmark):

    """
    Cola objective function.

    This class defines the Cola global optimization problem. The 17-dimensional function computes
    indirectly the formula :math:`f(n, u)` by setting :math:`x_0 = y_0, x_1 = u_0, x_i = u_{2(i2)}, y_i = u_{2(i2)+1}` :

    .. math::

        f_{\\text{Cola}}(\\mathbf{x}) = \\sum_{i<j}^{n} \\left (r_{i,j} - d_{i,j} \\right )^2

    Where :math:`r_{i,j}` is given by:

    .. math::

        r_{i,j} = \\sqrt{(x_i - x_j)^2 + (y_i - y_j)^2}

    And :math:`d` is a symmetric matrix given by:

    .. math::

        \\mathbf{d} = \\left [ d_{ij} \\right ] = \\begin{pmatrix}
        1.27 &  &  &  &  &  &  &  & \\\\
        1.69 & 1.43 &  &  &  &  &  &  & \\\\
        2.04 & 2.35 & 2.43 &  &  &  &  &  & \\\\
        3.09 & 3.18 & 3.26 & 2.85  &  &  &  &  & \\\\
        3.20 & 3.22 & 3.27 & 2.88 & 1.55 &  &  &  & \\\\
        2.86 & 2.56 & 2.58 & 2.59 & 3.12 & 3.06  &  &  & \\\\
        3.17 & 3.18 & 3.18 & 3.12 & 1.31 & 1.64 & 3.00  & \\\\
        3.21 & 3.18 & 3.18 & 3.17 & 1.70 & 1.36 & 2.95 & 1.32  & \\\\
        2.38 & 2.31 & 2.42 & 1.94 & 2.85 & 2.81 & 2.56 & 2.91 & 2.97
        \\end{pmatrix}

    This function has bounds :math:`0 \\leq x_0 \\leq 4` and :math:`-4 \\leq x_i
    \\leq 4` for :math:`i = 1,...,n-1`. It has a global minimum of 11.7464.
    """

    def __init__(self, dimensions=17):
        Benchmark.__init__(self, dimensions)

        self._bounds = [[0.0, 4.0]] + \
            list(zip([-4.0] * (self.N - 1),
                 [4.0] * (self.N - 1)))

        self.global_optimum = [
            [0.651906, 1.30194, 0.099242, -0.883791, -0.8796,
             0.204651, -3.28414, 0.851188, -3.46245, 2.53245,
             -0.895246, 1.40992, -3.07367, 1.96257, -2.97872,
             -0.807849, -1.68978]]
        self.fglob = 11.7464

        self.d = asarray([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                 [1.27, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                 [1.69, 1.43, 0, 0, 0, 0, 0, 0, 0, 0],
                 [2.04, 2.35, 2.43, 0, 0, 0, 0, 0, 0, 0],
                 [3.09, 3.18, 3.26, 2.85, 0, 0, 0, 0, 0, 0],
                 [3.20, 3.22, 3.27, 2.88, 1.55, 0, 0, 0, 0, 0],
                 [2.86, 2.56, 2.58, 2.59, 3.12, 3.06, 0, 0, 0, 0],
                 [3.17, 3.18, 3.18, 3.12, 1.31, 1.64, 3.00, 0, 0, 0],
                 [3.21, 3.18, 3.18, 3.17, 1.70, 1.36, 2.95, 1.32, 0, 0],
                 [2.38, 2.31, 2.42, 1.94, 2.85, 2.81, 2.56, 2.91, 2.97, 0.]])

    def fun(self, x, *args):
        self.nfev += 1

        xi = np.atleast_2d(asarray([0.0, x[0]] + list(x[1::2])))
        xj = np.repeat(xi, np.size(xi, 1), axis=0)
        xi = xi.T

        yi = np.atleast_2d(asarray([0.0, 0.0] + list(x[2::2])))
        yj = np.repeat(yi, np.size(yi, 1), axis=0)
        yi = yi.T

        inner = (sqrt(((xi - xj) ** 2 + (yi - yj) ** 2)) - self.d) ** 2
        inner = np.tril(inner, -1)
        return sum(sum(inner, axis=1))


class Colville(Benchmark):

    """
    Colville objective function.

    This class defines the Colville global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Colville}}(\\mathbf{x}) = \\left(x_{1} -1\\right)^{2} + 100 \\left(x_{1}^{2} - x_{2}\\right)^{2} + 10.1 \\left(x_{2} -1\\right)^{2} + \\left(x_{3} -1\\right)^{2} + 90 \\left(x_{3}^{2} - x_{4}\\right)^{2} + 10.1 \\left(x_{4} -1\\right)^{2} + 19.8 \\frac{x_{4} -1}{x_{2}}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,4`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 1` for :math:`i=1,...,4`

    """

    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[1 for _ in range(self.N)]]
        self.fglob = 0.0

    def fun(self, x, *args):

        self.nfev += 1
        return (100 * (x[0] ** 2 - x[1]) ** 2
                + (x[0] - 1) ** 2 + (x[2] - 1) ** 2
                + 90 * (x[2] ** 2 - x[3]) ** 2
                + 10.1 * ((x[1] - 1) ** 2 + (x[3] - 1) ** 2)
                + 19.8 * (1 / x[1]) * (x[3] - 1))


class Corana(Benchmark):

    """
    Corana objective function.

    This class defines the Corana global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Corana}}(\\mathbf{x}) = \\begin{cases} \\sum_{i=1}^n 0.15 d_i
        [z_i - 0.05\\textrm{sgn}(z_i)]^2 & \\textrm{if}|x_i-z_i| < 0.05 \\\\
        d_ix_i^2 & \\textrm{otherwise}\\end{cases}

    Where, in this exercise:

    .. math::

        z_i = 0.2 \\lfloor |x_i/s_i|+0.49999\\rfloor\\textrm{sgn}(x_i), d_i=(1,1000,10,100, ...)


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,...,4`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,4`

    """

    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.0] * self.N, [5.0] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0

    def fun(self, x, *args):

        self.nfev += 1

        d = [1., 1000., 10., 100.]
        r = 0
        for j in range(4):
            zj = floor(abs(x[j] / 0.2) + 0.49999) * sign(x[j]) * 0.2
            if abs(x[j] - zj) < 0.05:
                r += 0.15 * ((zj - 0.05 * sign(zj)) ** 2) * d[j]
            else:
                r += d[j] * x[j] * x[j]
        return r


class CosineMixture(Benchmark):

    """
    Cosine Mixture objective function.

    This class defines the Cosine Mixture global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{CosineMixture}}(\\mathbf{x}) = -0.1 \\sum_{i=1}^n \\cos(5 \\pi x_i) - \\sum_{i=1}^n x_i^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,...,N`.

    *Global optimum*: :math:`f(x_i) = -0.1N` for :math:`x_i = 0` for :math:`i=1,...,N`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.change_dimensionality = True
        self._bounds = zip([-1.0] * self.N, [1.0] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = -0.1 * self.N

    def fun(self, x, *args):
        self.nfev += 1

        return -0.1 * sum(cos(5.0 * pi * x)) - sum(x ** 2.0)


class CrossInTray(Benchmark):

    """
    Cross-in-Tray objective function.

    This class defines the Cross-in-Tray global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{CrossInTray}}(\\mathbf{x}) = - 0.0001 \\left(\\left|{e^{\\left|{100 - \\frac{\\sqrt{x_{1}^{2} + x_{2}^{2}}}{\\pi}}\\right|} \\sin\\left(x_{1}\\right) \\sin\\left(x_{2}\\right)}\\right| + 1\\right)^{0.1}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-15, 15]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -2.062611870822739` for :math:`x_i = \\pm 1.349406608602084` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [(1.349406685353340, 1.349406608602084),
                               (-1.349406685353340, 1.349406608602084),
                               (1.349406685353340, -1.349406608602084),
                               (-1.349406685353340, -1.349406608602084)]
        self.fglob = -2.062611870822739

    def fun(self, x, *args):

        self.nfev += 1
        return (-0.0001 * (abs(sin(x[0]) * sin(x[1])
                           * exp(abs(100 - sqrt(x[0] ** 2 + x[1] ** 2) / pi)))
                           + 1) ** (0.1))


class CrossLegTable(Benchmark):

    """
    Cross-Leg-Table objective function.

    This class defines the Cross-Leg-Table global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{CrossLegTable}}(\\mathbf{x}) = - \\frac{1}{\\left(\\left|{e^{\\left|{100 - \\frac{\\sqrt{x_{1}^{2} + x_{2}^{2}}}{\\pi}}\\right|} \\sin\\left(x_{1}\\right) \\sin\\left(x_{2}\\right)}\\right| + 1\\right)^{0.1}}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -1`. The global minimum is found on the planes :math:`x_1 = 0` and :math:`x_2 = 0`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[0., 0.]]
        self.fglob = -1.0

    def fun(self, x, *args):

        self.nfev += 1
        u = 100 - sqrt(x[0] ** 2 + x[1] ** 2) / pi
        v = sin(x[0]) * sin(x[1])
        return -(abs(v * exp(abs(u))) + 1) ** (-0.1)


class CrownedCross(Benchmark):

    """
    Crowned Cross objective function.

    This class defines the Crowned Cross global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{CrownedCross}}(\\mathbf{x}) = 0.0001 \\left(\\left|{e^{\\left|{100- \\frac{\\sqrt{x_{1}^{2} + x_{2}^{2}}}{\\pi}}\\right|} \\sin\\left(x_{1}\\right) \\sin\\left(x_{2}\\right)}\\right| + 1\\right)^{0.1}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0.0001`. The global minimum is found on the planes :math:`x_1 = 0` and :math:`x_2 = 0`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[0, 0]]
        self.fglob = 0.0001

    def fun(self, x, *args):

        self.nfev += 1
        u = 100 - sqrt(x[0] ** 2 + x[1] ** 2) / pi
        v = sin(x[0]) * sin(x[1])
        return 0.0001 * (abs(v * exp(abs(u))) + 1) ** (0.1)


class Csendes(Benchmark):

    """
    Csendes objective function.

    This class defines the Csendes global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Csendes}}(\\mathbf{x}) = \\sum_{i=1}^n x_i^6 \\left[ 2 + \\sin \\left( \\frac{1}{x_i} \\right ) \\right]


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,...,N`.

    *Global optimum*: :math:`f(x_i) = 0.0` for :math:`x_i = 0` for :math:`i=1,...,N`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.change_dimensionality = True
        self._bounds = zip([-1.0] * self.N, [1.0] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = np.nan

    def fun(self, x, *args):
        self.nfev += 1

        try:
            return sum((x ** 6.0) * (2.0 + sin(1.0 / x)))
        except ZeroDivisionError, FloatingPointError:
            return np.nan

    def success(self, x):
        """Is a candidate solution at the global minimum"""
        val = self.fun(asarray(x))
        if np.isnan(val):
            return True
        try:
            np.testing.assert_almost_equal(val, 0., 4)
            return True
        except AssertionError:
            return False

        return False


class Cube(Benchmark):

    """
    Cube objective function.

    This class defines the Cube global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Cube}}(\\mathbf{x}) = 100(x_2 - x_1^3)^2 + (1 - x1)^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,N`.

    *Global optimum*: :math:`f(x_i) = 0.0` for :math:`\\mathbf{x} = [1, 1]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.custom_bounds = ([0, 2], [0, 2])

        self.global_optimum = [[1.0, 1.0]]
        self.fglob = 0.0

    def fun(self, x, *args):

        self.nfev += 1
        return 100.0 * (x[1] - x[0] ** 3.0) ** 2.0 + (1.0 - x[0]) ** 2.0
