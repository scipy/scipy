# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import
import numpy as np
from numpy import (abs, arctan2, asarray, cos, exp, floor, log, log10,
                   arange, pi, prod, roll, sign, sin, sqrt, sum, where,
                   zeros, tan, tanh, dot)

from scipy.misc import factorial
from .go_benchmark import Benchmark


class Matyas(Benchmark):

    """
    Matyas objective function.

    This class defines the Matyas global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Matyas}}(\\mathbf{x}) = 0.26(x_1^2 + x_2^2) - 0.48x_1x_2


    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return 0.26 * (x[0] ** 2 + x[1] ** 2) - 0.48 * x[0] * x[1]


class McCormick(Benchmark):

    """
    McCormick objective function.

    This class defines the McCormick global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{McCormick}}(\\mathbf{x}) = - x_{1} + 2 x_{2} + \\left(x_{1} - x_{2}\\right)^{2} + \\sin\\left(x_{1} + x_{2}\\right) + 1


    Here, :math:`n` represents the number of dimensions and :math:`x_1 \\in [-1.5, 4]`, :math:`x_2 \\in [-3, 4]`.

    *Global optimum*: :math:`f(x_i) = -1.913222954981037` for :math:`\\mathbf{x} = [-0.5471975602214493, -1.547197559268372]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = [(-1.5, 4.0), (-3.0, 4.0)]

        self.global_optimum = [[-0.5471975602214493, -1.547197559268372]]
        self.fglob = -1.913222954981037

    def fun(self, x, *args):
        self.nfev += 1

        return (sin(x[0] + x[1]) + (x[0] - x[1]) ** 2 - 1.5 * x[0]
                + 2.5 * x[1] + 1)


class Meyer(Benchmark):

    """
    Meyer objective function.

    """

    # TODO, this is a NIST regression standard dataset
    def __init__(self, dimensions=3):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0., 100., 100.],
                           [1, 1000., 500.])
        self.global_optimum = [[5.6096364710e-3, 6.1813463463e3,
                                3.4522363462e2]]
        self.fglob = 8.7945855171e1
        self.a = asarray([3.478E+04, 2.861E+04, 2.365E+04, 1.963E+04, 1.637E+04,
                          1.372E+04, 1.154E+04, 9.744E+03, 8.261E+03, 7.030E+03,
                          6.005E+03, 5.147E+03, 4.427E+03, 3.820E+03, 3.307E+03,
                          2.872E+03])
        self.b = asarray([5.000E+01, 5.500E+01, 6.000E+01, 6.500E+01, 7.000E+01,
                          7.500E+01, 8.000E+01, 8.500E+01, 9.000E+01, 9.500E+01,
                          1.000E+02, 1.050E+02, 1.100E+02, 1.150E+02, 1.200E+02,
                          1.250E+02])

    def fun(self, x, *args):
        self.nfev += 1

        vec = x[0] * exp(x[1] / (self.b + x[2]))

        return sum((self.a - vec) ** 2)


class Michalewicz(Benchmark):

    """
    Michalewicz objective function.

    This class defines the Michalewicz global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Michalewicz}}(\\mathbf{x}) = - \\sum_{i=1}^{2} \\sin\\left(x_i\\right) \\sin^{2 m}\\left(\\frac{i x_i^{2}}{\\pi}\\right)


    Where, in this exercise, :math:`m = 10`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, \\pi]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -1.8013` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N, [pi] * self.N)

        self.global_optimum = [[2.20290555, 1.570796]]
        self.fglob = -1.8013

    def fun(self, x, *args):
        self.nfev += 1

        m = 10.0
        i = arange(1, self.N + 1)
        return -sum(sin(x) * sin(i * x ** 2 / pi) ** (2 * m))


class MieleCantrell(Benchmark):

    """
    Miele-Cantrell objective function.

    This class defines the Miele-Cantrell global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{MieleCantrell}}(\\mathbf{x}) = (e^{-x_1} - x_2)^4 + 100(x_2 - x_3)^6 + \\tan^4(x_3 - x_4) + x_1^8


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,...,4`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [0, 1, 1, 1]`

    """

    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-1.0] * self.N, [1.0] * self.N)

        self.global_optimum = [[0.0, 1.0, 1.0, 1.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return ((exp(-x[0]) - x[1]) ** 4 + 100 * (x[1] - x[2]) ** 6
                + tan(x[2] - x[3]) ** 4 + x[0] ** 8)


class Mishra01(Benchmark):

    """
    Mishra 1 objective function.

    This class defines the Mishra 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra01}}(\\mathbf{x}) = (1 + x_n)^{x_n} \\hspace{10pt} ; \\hspace{10pt} x_n = n - \\sum_{i=1}^{n-1} x_i


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 1]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 2` for :math:`x_i = 1` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N,
                           [1.0 + 1e-9] * self.N)

        self.global_optimum = [[1.0 for _ in range(self.N)]]
        self.fglob = 2.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        xn = self.N - sum(x[0:-1])
        return (1 + xn) ** xn


class Mishra02(Benchmark):

    """
    Mishra 2 objective function.

    This class defines the Mishra 2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra02}}(\\mathbf{x}) = (1 + x_n)^{x_n} \\hspace{10pt} ; \\hspace{10pt} x_n = n - \\sum_{i=1}^{n-1} \\frac{(x_i + x_{i+1})}{2}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 1]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 2` for :math:`x_i = 1` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N,
                           [1.0 + 1e-9] * self.N)

        self.global_optimum = [[1.0 for _ in range(self.N)]]
        self.fglob = 2.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        xn = self.N - sum((x[:-1] + x[1:]) / 2.0)
        return (1 + xn) ** xn


class Mishra03(Benchmark):

    """
    Mishra 3 objective function.

    This class defines the Mishra 3 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra03}}(\\mathbf{x}) = \\sqrt{\\lvert \\cos{\\sqrt{\\lvert x_1^2 + x_2^2 \\rvert}} \\rvert} + 0.01(x_1 + x_2)


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -0.1999` for :math:`x_i = {-9.99378322, -9.99918927}`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[-9.99378322, -9.99918927]]
        self.fglob = -0.19990562

    def fun(self, x, *args):
        self.nfev += 1

        return ((0.01 * (x[0] + x[1])
                + sqrt(abs(cos(sqrt(abs(x[0] ** 2 + x[1] ** 2)))))))


class Mishra04(Benchmark):

    """
    Mishra 4 objective function.

    This class defines the Mishra 4 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra04}}(\\mathbf{x}) = \\sqrt{\\lvert \\sin{\\sqrt{\\lvert
       x_1^2 + x_2^2 \\rvert}} \\rvert} + 0.01(x_1 + x_2)


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in
    [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -0.17767` for :math:`x_i = {-8.71499636, -9.0533148}`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[-8.92011797, -8.85128046]]
        self.fglob = -0.17768386538

    def fun(self, x, *args):
        self.nfev += 1

        return ((0.01 * (x[0] + x[1])
                + sqrt(abs(sin(sqrt(abs(x[0] ** 2 + x[1] ** 2)))))))


class Mishra05(Benchmark):

    """
    Mishra 5 objective function.

    This class defines the Mishra 5 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra05}}(\\mathbf{x}) = \\left [ \\sin^2 ((\\cos(x_1) + \\cos(x_2))^2) + \\cos^2 ((\\sin(x_1) + \\sin(x_2))^2) + x_1 \\right ]^2 + 0.01(x_1 + x_2)


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -0.119829` for :math:`\\mathbf{x} = [-1.98682, -10]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[-1.98682, -10.0]]
        self.fglob = -0.119829

    def fun(self, x, *args):
        self.nfev += 1

        return (0.01 * (x[0] + x[1])
                + (sin((cos(x[0]) + cos(x[1])) ** 2) ** 2
                   + cos((sin(x[0]) + sin(x[1])) ** 2) ** 2 + x[0]) ** 2)


class Mishra06(Benchmark):

    """
    Mishra 6 objective function.

    This class defines the Mishra 6 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra06}}(\\mathbf{x}) = -\\log{\\left [ \\sin^2 ((\\cos(x_1) + \\cos(x_2))^2) - \\cos^2 ((\\sin(x_1) + \\sin(x_2))^2) + x_1 \\right ]^2} + 0.01 \\left[(x_1 -1)^2 + (x_2 - 1)^2 \\right]


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -2.28395` for :math:`\\mathbf{x} = [2.88631, 1.82326]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[2.88631, 1.82326]]
        self.fglob = -2.28395

    def fun(self, x, *args):
        self.nfev += 1

        a = 0.1 * ((x[0] - 1) ** 2 + (x[1] - 1) ** 2)
        u = (cos(x[0]) + cos(x[1])) ** 2
        v = (sin(x[0]) + sin(x[1])) ** 2
        return a - log((sin(u) ** 2 - cos(v) ** 2 + x[0]) ** 2)


class Mishra07(Benchmark):

    """
    Mishra 7 objective function.

    This class defines the Mishra 7 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra07}}(\\mathbf{x}) = \\left [\\prod_{i=1}^{n} x_i - n! \\right]^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = \\sqrt{n}` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.custom_bounds = [(-2, 2), (-2, 2)]
        self.global_optimum = [[sqrt(self.N)
                               for i in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return (prod(x) - factorial(self.N)) ** 2.0


class Mishra08(Benchmark):

    """
    Mishra 8 objective function.

    This class defines the Mishra 8 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra08}}(\\mathbf{x}) = 0.001 \\left[\\lvert x_1^{10} - 20x_1^9 + 180x_1^8 - 960 x_1^7 + 3360x_1^6 - 8064x_1^5 + 13340x_1^4 - 15360x_1^3 + 11520x_1^2 - 5120x_1 + 2624 \\rvert \\lvert x_2^4 + 12x_2^3 + 54x_2^2 + 108x_2 + 81 \\rvert \\right]^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [2, -3]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.custom_bounds = [(1.0, 2.0), (-4.0, 1.0)]
        self.global_optimum = [[2.0, 3.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        val = abs(x[0] ** 10 - 20 * x[0] ** 9 + 180 * x[0] ** 8
                  - 960 * x[0] ** 7 + 3360 * x[0] ** 6 - 8064 * x[0] ** 5
                  + 13340 * x[0] ** 4 - 15360 * x[0] ** 3 + 11520 * x[0] ** 2
                  - 5120 * x[0] + 2624)
        val *= abs(x[1] ** 2 + 12 * x[1] ** 2 +
                   54 * x[1] ** 2 + 108 * x[1] + 81)
        return 0.001 * val ** 2


class Mishra09(Benchmark):

    """
    Mishra 9 objective function.

    This class defines the Mishra 9 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra09}}(\\mathbf{x}) = \\left[ ab^2c + abc^2 + b^2 + (x_1 + x_2 - x_3)^2 \\right]^2


    Where, in this exercise:

    .. math::

        \\begin{cases} a = 2x_1^3 + 5x_1x_2 + 4x_3 - 2x_1^2x_3 - 18 \\\\
        b = x_1 + x_2^3 + x_1x_2^2 + x_1x_3^2 - 22 \\\\
        c = 8x_1^2 + 2x_2x_3 + 2x_2^2 + 3x_2^3 - 52 \\end{cases}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2,3`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [1, 2, 3]`

    """

    def __init__(self, dimensions=3):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.global_optimum = [[1.0, 2.0, 3.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        a = (2 * x[0] ** 3 + 5 * x[0] * x[1]
             + 4 * x[2] - 2 * x[0] ** 2 * x[2] - 18)
        b = x[0] + x[1] ** 3 + x[0] * x[1] ** 2 + x[0] * x[2] ** 2 - 22.0
        c = 8 * x[0] ** 2 + 2 * x[1] * x[2] + \
            2 * x[1] ** 2 + 3 * x[1] ** 3 - 52

        return (a * c * b ** 2 + a * b * c ** 2 + b ** 2
                + (x[0] + x[1] - x[2]) ** 2) ** 2


class Mishra10(Benchmark):

    """
    Mishra 10 objective function.

    This class defines the Mishra 10 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::
    TODO - int(x) should be used instead of floor(x)!!!!!
       f_{\\text{Mishra10}}(\\mathbf{x}) = \\left[ \\lfloor x_1 \\perp x_2 \\rfloor - \\lfloor x_1 \\rfloor - \\lfloor x_2 \\rfloor \\right]^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [2, 2]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.global_optimum = [[2.0, 2.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        x1, x2 = int(x[0]), int(x[1])
        # TODO rewrite equation above with nint(x)
        f1 = x1 + x2
        f2 = x1 * x2
        return (f1 - f2) ** 2.0


class Mishra11(Benchmark):

    """
    Mishra 11 objective function.

    This class defines the Mishra 11 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra11}}(\\mathbf{x}) = \\left [ \\frac{1}{n} \\sum_{i=1}^{n} \\lvert x_i \\rvert - \\left(\\prod_{i=1}^{n} \\lvert x_i \\rvert \\right )^{\\frac{1}{n}} \\right]^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.custom_bounds = [(-3, 3), (-3, 3)]

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        N = self.N
        return ((1.0 / N) * sum(abs(x)) - (prod(abs(x))) ** 1.0 / N) ** 2.0


class MultiModal(Benchmark):

    """
    MultiModal objective function.

    This class defines the MultiModal global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{MultiModal}}(\\mathbf{x}) = \\left( \\sum_{i=1}^n \\lvert x_i \\rvert \\right) \\left( \\prod_{i=1}^n \\lvert x_i \\rvert \\right)


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.custom_bounds = [(-5, 5), (-5, 5)]

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum(abs(x)) * prod(abs(x))
