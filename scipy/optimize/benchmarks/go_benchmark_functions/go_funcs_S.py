# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import
import numpy as np
from numpy import (abs, arctan2, asarray, cos, exp, floor, log, log10,
                   arange, pi, prod, roll, sign, sin, sqrt, sum, where,
                   zeros, tan, tanh, dot)

from scipy.misc import factorial
from .go_benchmark import Benchmark


class Salomon(Benchmark):

    """
    Salomon objective function.

    This class defines the Salomon global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Salomon}}(\\mathbf{x}) = 1 - \\cos \\left (2 \\pi
       \\sqrt{\\sum_{i=1}^{n} x_i^2} \\right) + 0.1 \\sqrt{\\sum_{i=1}^n x_i^2}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in
    [-100, 100]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for
    :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = [(-50, 50), (-50, 50)]

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        u = sum(x ** 2)
        return 1 - cos(2 * pi * sqrt(u)) + 0.1 * sqrt(u)


class Sargan(Benchmark):

    """
    Sargan objective function.

    This class defines the Sargan global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Sargan}}(\\mathbf{x}) = \\sum_{i=1}^{n} n \\left (x_i^2 + 0.4 \\sum_{i \\neq j}^{n} x_ix_j \\right)

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = [(-5, 5), (-5, 5)]

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        x0 = x[:-1]
        x1 = roll(x, -1)[:-1]

        return sum(self.N * (x ** 2 + 0.4 * sum(x0 * x1)))


class Schaffer01(Benchmark):

    """
    Schaffer 1 objective function.

    This class defines the Schaffer 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schaffer01}}(\\mathbf{x}) = 0.5 + \\frac{\\sin^2 (x_1^2 + x_2^2)^2 - 0.5}{1 + 0.001(x_1^2 + x_2^2)^2}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = [(-10, 10), (-10, 10)]

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        num = sin((x[0] ** 2 + x[1] ** 2)) ** 2 - 0.5
        den = (1 + 0.001 * (x[0] ** 2 + x[1] ** 2)) ** 2
        return 0.5 + num / den


class Schaffer02(Benchmark):

    """
    Schaffer 2 objective function.

    This class defines the Schaffer 2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schaffer02}}(\\mathbf{x}) = 0.5 + \\frac{\\sin^2 (x_1^2 - x_2^2)^2 - 0.5}{1 + 0.001(x_1^2 + x_2^2)^2}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = [(-10, 10), (-10, 10)]

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        num = sin((x[0] ** 2 - x[1] ** 2)) ** 2 - 0.5
        den = (1 + 0.001 * (x[0] ** 2 + x[1] ** 2)) ** 2
        return 0.5 + num / den


class Schaffer03(Benchmark):

    """
    Schaffer 3 objective function.

    This class defines the Schaffer 3 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schaffer03}}(\\mathbf{x}) = 0.5 + \\frac{\\sin^2 \\left( \\cos \\lvert x_1^2 - x_2^2 \\rvert \\right ) - 0.5}{1 + 0.001(x_1^2 + x_2^2)^2}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0.00156685` for :math:`\\mathbf{x} = [0, 1.253115]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = [(-10, 10), (-10, 10)]

        self.global_optimum = [[0.0, 1.253115]]
        self.fglob = 0.00156685

    def fun(self, x, *args):
        self.nfev += 1

        num = sin(cos(abs(x[0] ** 2 - x[1] ** 2))) ** 2 - 0.5
        den = (1 + 0.001 * (x[0] ** 2 + x[1] ** 2)) ** 2
        return 0.5 + num / den


class Schaffer04(Benchmark):

    """
    Schaffer 4 objective function.

    This class defines the Schaffer 4 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schaffer04}}(\\mathbf{x}) = 0.5 + \\frac{\\cos^2 \\left( \\sin(x_1^2 - x_2^2) \\right ) - 0.5}{1 + 0.001(x_1^2 + x_2^2)^2}^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0.292579` for :math:`\\mathbf{x} = [0, 1.253115]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = [(-10, 10), (-10, 10)]

        self.global_optimum = [[0.0, 1.253115]]
        self.fglob = 0.292579

    def fun(self, x, *args):
        self.nfev += 1

        num = cos(sin(abs(x[0] ** 2 - x[1] ** 2))) ** 2 - 0.5
        den = (1 + 0.001 * (x[0] ** 2 + x[1] ** 2)) ** 2
        return 0.5 + num / den


class SchmidtVetters(Benchmark):

    """
    Schmidt-Vetters objective function.

    This class defines the Schmidt-Vetters global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{SchmidtVetters}}(\\mathbf{x}) = \\frac{1}{1 + (x_1 - x_2)^2} + \\sin \\left(\\frac{\\pi x_2 + x_3}{2} \\right) + e^{\\left(\\frac{x_1+x_2}{x_2} - 2\\right)^2}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 10]` for :math:`i=1,2,3`.

    *Global optimum*: :math:`f(x_i) = 2.99643266` for :math:`x_i = [0.79876108,  0.79962581,  0.79848824]`

    """

    def __init__(self, dimensions=3):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([0.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[0.79876108, 0.79962581, 0.79848824]]
        self.fglob = 2.99643266

    def fun(self, x, *args):
        self.nfev += 1

        return (1 / (1 + (x[0] - x[1]) ** 2) + sin((pi * x[1] + x[2]) / 2)
                + exp(((x[0] + x[1]) / x[1] - 2) ** 2))


class Schwefel01(Benchmark):

    """
    Schwefel 1 objective function.

    This class defines the Schwefel 1 global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel01}}(\\mathbf{x}) = \\left(\\sum_{i=1}^n x_i^2 \\right)^{\\alpha}

    Where, in this exercise, :math:`\\alpha = \\sqrt{\\pi}`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = ([-4.0, 4.0], [-4.0, 4.0])

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        alpha = sqrt(pi)
        return (sum(x ** 2.0)) ** alpha


class Schwefel02(Benchmark):

    """
    Schwefel 2 objective function.

    This class defines the Schwefel 2 global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel02}}(\\mathbf{x}) = \\sum_{i=1}^n \\left(\\sum_{j=1}^i x_i \\right)^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = ([-4.0, 4.0], [-4.0, 4.0])

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        mat = np.repeat(np.atleast_2d(x), self.N, axis=0)
        inner = sum(np.tril(mat), axis=1)
        return sum(inner ** 2)


class Schwefel04(Benchmark):

    """
    Schwefel 4 objective function.

    This class defines the Schwefel 4 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel04}}(\\mathbf{x}) = \\sum_{i=1}^n \\left[(x_i - 1)^2 + (x_1 - x_i^2)^2 \\right]


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 1` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([0.0] * self.N, [10.0] * self.N)
        self.custom_bounds = ([0.0, 2.0], [0.0, 2.0])

        self.global_optimum = [[1.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum((x - 1.0) ** 2.0 + (x[0] - x ** 2.0) ** 2.0)


class Schwefel06(Benchmark):

    """
    Schwefel 6 objective function.

    This class defines the Schwefel 6 global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel06}}(\\mathbf{x}) = \\max(\\lvert x_1 + 2x_2 - 7 \\rvert, \\lvert 2x_1 + x_2 - 5 \\rvert)


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [1, 3]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = ([-10.0, 10.0], [-10.0, 10.0])

        self.global_optimum = [[1.0, 3.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return max(abs(x[0] + 2 * x[1] - 7), abs(2 * x[0] + x[1] - 5))


class Schwefel20(Benchmark):

    """
    Schwefel 20 objective function.

    This class defines the Schwefel 20 global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel20}}(\\mathbf{x}) = \\sum_{i=1}^n \\lvert x_i \\rvert


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum(abs(x))


class Schwefel21(Benchmark):

    """
    Schwefel 21 objective function.

    This class defines the Schwefel 21 global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel21}}(\\mathbf{x}) = \\smash{\\displaystyle\\max_{1 \leq i \leq n}} \\lvert x_i \\rvert


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return max(abs(x))


class Schwefel22(Benchmark):

    """
    Schwefel 22 objective function.

    This class defines the Schwefel 22 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel22}}(\\mathbf{x}) = \\sum_{i=1}^n \\lvert x_i \\rvert + \\prod_{i=1}^n \\lvert x_i \\rvert


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = ([-10.0, 10.0], [-10.0, 10.0])

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum(abs(x)) + prod(abs(x))


class Schwefel26(Benchmark):

    """
    Schwefel 26 objective function.

    This class defines the Schwefel 26 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel26}}(\\mathbf{x}) = 418.9829n - \\sum_{i=1}^n x_i \\sin(\\sqrt{|x_i|})

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-500, 500]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 420.968746` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([-500.0] * self.N,
                           [500.0] * self.N)

        self.global_optimum = [[420.968746 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return 418.982887 * self.N - sum(x * sin(sqrt(abs(x))))


class Schwefel36(Benchmark):

    """
    Schwefel 36 objective function.

    This class defines the Schwefel 36 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel36}}(\\mathbf{x}) = -x_1x_2(72 - 2x_1 - 2x_2)


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 500]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -3456` for :math:`\\mathbf{x} = [12, 12]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([0.0] * self.N, [500.0] * self.N)
        self.custom_bounds = ([0.0, 20.0], [0.0, 20.0])

        self.global_optimum = [[12.0, 12.0]]
        self.fglob = -3456.0

    def fun(self, x, *args):
        self.nfev += 1

        return -x[0] * x[1] * (72 - 2 * x[0] - 2 * x[1])


class Shekel05(Benchmark):

    """
    Shekel 5 objective function.

    This class defines the Shekel 5 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Shekel05}}(\\mathbf{x}) = \\sum_{i=1}^{m} \\frac{1}{c_{i} + \\sum_{j=1}^{n} (x_{j} - a_{ij})^2 }`

    Where, in this exercise:

    .. math::

        \\mathbf{a} = \\begin{bmatrix} 4.0 & 4.0 & 4.0 & 4.0 \\\\ 1.0 & 1.0 & 1.0 & 1.0 \\\\ 8.0 & 8.0 & 8.0 & 8.0 \\\\ 6.0 & 6.0 & 6.0 & 6.0 \\\\ 3.0 & 7.0 & 3.0 & 7.0 \\end{bmatrix}

    .. math::

        \\mathbf{c} = \\begin{bmatrix} 0.1 \\\\ 0.2 \\\\ 0.2 \\\\ 0.4 \\\\ 0.4 \\end{bmatrix}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 10]` for :math:`i=1,...,4`.

    *Global optimum*: :math:`f(x_i) = -10.15319585` for :math:`x_i = 4` for :math:`i=1,...,4`

    """

    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[4.0 for _ in range(self.N)]]
        self.fglob = -10.15319585
        self.A = asarray([[4.0, 4.0, 4.0, 4.0],
                          [1.0, 1.0, 1.0, 1.0],
                          [8.0, 8.0, 8.0, 8.0],
                          [6.0, 6.0, 6.0, 6.0],
                          [3.0, 7.0, 3.0, 7.0]])

        self.C = asarray([0.1, 0.2, 0.2, 0.4, 0.4])

    def fun(self, x, *args):
        self.nfev += 1

        return -sum(1.0 / (dot(x - a, x - a) + c) for a, c
                    in zip(self.A, self.C))


class Shekel07(Benchmark):

    """
    Shekel 7 objective function.

    This class defines the Shekel 7 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Shekel07}}(\\mathbf{x}) = \\sum_{i=1}^{m} \\frac{1}{c_{i} + \\sum_{j=1}^{n} (x_{j} - a_{ij})^2 }`

    Where, in this exercise:

    .. math::

        \\mathbf{a} = \\begin{bmatrix} 4.0 & 4.0 & 4.0 & 4.0 \\\\ 1.0 & 1.0 & 1.0 & 1.0 \\\\ 8.0 & 8.0 & 8.0 & 8.0 \\\\
        6.0 & 6.0 & 6.0 & 6.0 \\\\ 3.0 & 7.0 & 3.0 & 7.0 \\\\ 2.0 & 9.0 & 2.0 & 9.0 \\\\ 5.0 & 5.0 & 3.0 & 3.0 \\end{bmatrix}

    .. math::

        \\mathbf{c} = \\begin{bmatrix} 0.1 \\\\ 0.2 \\\\ 0.2 \\\\ 0.4 \\\\ 0.4 \\\\ 0.6 \\\\ 0.3 \\end{bmatrix}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 10]` for :math:`i=1,...,4`.

    *Global optimum*: :math:`f(x_i) = -10.4028188` for :math:`x_i = 4` for :math:`i=1,...,4`

    """

    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[4.0 for _ in range(self.N)]]
        self.fglob = -10.4028188
        self.A = asarray([[4.0, 4.0, 4.0, 4.0],
                          [1.0, 1.0, 1.0, 1.0],
                          [8.0, 8.0, 8.0, 8.0],
                          [6.0, 6.0, 6.0, 6.0],
                          [3.0, 7.0, 3.0, 7.0],
                          [2.0, 9.0, 2.0, 9.0],
                          [5.0, 5.0, 3.0, 3.0]])

        self.C = asarray([0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3])

    def fun(self, x, *args):
        self.nfev += 1

        return -sum(1.0 / (dot(x - a, x - a) + c) for a, c
                    in zip(self.A, self.C))


class Shekel10(Benchmark):

    """
    Shekel 10 objective function.

    This class defines the Shekel 10 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Shekel10}}(\\mathbf{x}) = \\sum_{i=1}^{m} \\frac{1}{c_{i} + \\sum_{j=1}^{n} (x_{j} - a_{ij})^2 }`

    Where, in this exercise:

    .. math::

        \\mathbf{a} = \\begin{bmatrix} 4.0 & 4.0 & 4.0 & 4.0 \\\\ 1.0 & 1.0 & 1.0 & 1.0 \\\\ 8.0 & 8.0 & 8.0 & 8.0 \\\\
        6.0 & 6.0 & 6.0 & 6.0 \\\\ 3.0 & 7.0 & 3.0 & 7.0 \\\\ 2.0 & 9.0 & 2.0 & 9.0 \\\\ 5.0 & 5.0 & 3.0 & 3.0 \\\\
        8.0 & 1.0 & 8.0 & 1.0 \\\\ 6.0 & 2.0 & 6.0 & 2.0 \\\\ 7.0 & 3.6 & 7.0 & 3.6 \\end{bmatrix}

    .. math::

        \\mathbf{c} = \\begin{bmatrix} 0.1 \\\\ 0.2 \\\\ 0.2 \\\\ 0.4 \\\\ 0.4 \\\\ 0.6 \\\\ 0.3 \\\\ 0.7 \\\\ 0.5 \\\\ 0.5 \\end{bmatrix}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 10]` for :math:`i=1,...,4`.

    *Global optimum*: :math:`f(x_i) = -10.5362837` for :math:`x_i = 4` for :math:`i=1,...,4`

    """

    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[4.0 for _ in range(self.N)]]
        self.fglob = -10.5362837262
        self.A = asarray([[4.0, 4.0, 4.0, 4.0],
                          [1.0, 1.0, 1.0, 1.0],
                          [8.0, 8.0, 8.0, 8.0],
                          [6.0, 6.0, 6.0, 6.0],
                          [3.0, 7.0, 3.0, 7.0],
                          [2.0, 9.0, 2.0, 9.0],
                          [5.0, 5.0, 3.0, 3.0],
                          [8.0, 1.0, 8.0, 1.0],
                          [6.0, 2.0, 6.0, 2.0],
                          [7.0, 3.6, 7.0, 3.6]])

        self.C = asarray([0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3, 0.7, 0.5, 0.5])

    def fun(self, x, *args):
        self.nfev += 1

        return -sum(1.0 / (dot(x - a, x - a) + c) for a, c
                    in zip(self.A, self.C))


class Shubert01(Benchmark):

    """
    Shubert 1 objective function.

    This class defines the Shubert 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Shubert01}}(\\mathbf{x}) = \\left( \\sum\\limits_{i=1}^{5} i\\cos[(i+1)x_1 + i] \\right) \\left( \\sum\\limits_{i=1}^{5} i\\cos[(i+1)x_2 + i] \\right)

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = -186.7309` for :math:`\\mathbf{x} = [-7.0835, 4.8580]` (and many others).

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.global_optimum = [[-7.0835, 4.8580]]

        self.fglob = -186.7309

        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        j = np.atleast_2d(arange(1, 6)).T
        y = j * cos((j + 1) * x + j)
        return prod(sum(y, axis=0))
        # TODO change equation to reflect higher dimensions are possible


class Shubert03(Benchmark):

    """
    Shubert 3 objective function.

    This class defines the Shubert 3 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Shubert03}}(\\mathbf{x}) = \\sum_{i=1}^n \\sum_{j=1}^5 j \\sin \\left[(j+1)x_i \\right] + j

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = -24.062499` for :math:`\\mathbf{x} = [5.791794, 5.791794]` (and many others).

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[5.791794, 5.791794]]
        self.fglob = -24.062499

        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        j = np.atleast_2d(arange(1, 6)).T
        y = -j * sin((j + 1) * x + j)
        # TODO change equation to reflect higher dimensions are possible
        return sum(sum(y))


class Shubert04(Benchmark):

    """
    Shubert 4 objective function.

    This class defines the Shubert 4 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Shubert04}}(\\mathbf{x}) = \\sum_{i=1}^n \\sum_{j=1}^5 j \\cos \\left[(j+1)x_i \\right] + j

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = -29.016015` for :math:`\\mathbf{x} = [-0.80032121, -7.08350592]` (and many others).

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[-0.80032121, -7.08350592]]
        self.fglob = -29.016015

        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        j = np.atleast_2d(arange(1, 6)).T
        y = -j * cos((j + 1) * x + j)
        # TODO change equation to reflect higher dimensions are possible
        return sum(sum(y))


class SineEnvelope(Benchmark):

    """
    SineEnvelope objective function.

    This class defines the SineEnvelope global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{SineEnvelope}}(\\mathbf{x}) = -\\sum_{i=1}^{n-1}\\left[\\frac{\\sin^2(\\sqrt{x_{i+1}^2+x_{i}^2}-0.5)}{(0.001(x_{i+1}^2+x_{i}^2)+1)^2}+0.5\\right]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = [(-20, 20), (-20, 20)]

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        X0 = x[:-1]
        X1 = x[1:]
        X02X12 = X0 ** 2 + X1 ** 2
        return sum((sin(sqrt(X02X12)) ** 2 - 0.5) / (1 + 0.001 * X02X12) ** 2
                   + 0.5)


class SixHumpCamel(Benchmark):

    """
    Six Hump Camel objective function.

    This class defines the Six Hump Camel global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{SixHumpCamel}}(\\mathbf{x}) = 4x_1^2+x_1x_2-4x_2^2-2.1x_1^4+4x_2^4+\\frac{1}{3}x_1^6

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -1.031628453489877` for :math:`\\mathbf{x} = [0.08984201368301331 , -0.7126564032704135]`
    or :math:`\\mathbf{x} = [-0.08984201368301331, 0.7126564032704135]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.0] * self.N, [5.0] * self.N)
        self.custom_bounds = [(-2, 2), (-1.5, 1.5)]

        self.global_optimum = [(0.08984201368301331, -0.7126564032704135),
                               (-0.08984201368301331, 0.7126564032704135)]
        self.fglob = -1.031628

    def fun(self, x, *args):
        self.nfev += 1
        return ((4 - 2.1 * x[0] ** 2 + x[0] ** 4 / 3) * x[0] ** 2 + x[0] * x[1]
                + (4 * x[1] ** 2 - 4) * x[1] ** 2)


class Sodp(Benchmark):

    """
    Sodp objective function.

    This class defines the Sum Of Different Powers global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Sodp}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\lvert{x_{i}}\\rvert^{i + 1}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-1.0] * self.N, [1.0] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        i = arange(1, self.N + 1)
        return sum(abs(x) ** (i + 1))


class Sphere(Benchmark):

    """
    Sphere objective function.

    This class defines the Sphere global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Sphere}}(\\mathbf{x}) = \\sum_{i=1}^{n} x_i^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,...,n`.

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

        return sum(x ** 2)


class Step(Benchmark):

    """
    Step objective function.

    This class defines the Step global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Step}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\left ( \\lfloor x_i  + 0.5 \\rfloor \\right )^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0.5` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = ([-5, 5], [-5, 5])

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum((floor(x + 0.5)) ** 2.0)


class Stochastic(Benchmark):

    """
    Stochastic objective function.

    This class defines a Stochastic global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Stochastic}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\epsilon_i \\left | {x_i - \\frac{1}{i}} \\right |

    The variable :math:`\\epsilon_i, (i=1,...,n)` is a random variable uniformly distributed in :math:`[0, 1]`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = [1/n]` for :math:`i=1,...,n`
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.0] * self.N, [5.0] * self.N)

        self.global_optimum = [[1.0 / _ for _ in range(1, self.N + 1)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        rnd = np.random.uniform(0.0, 1.0, size=(self.N, ))
        i = arange(1, self.N + 1)

        return sum(rnd * abs(x - 1.0 / i))


class StretchedV(Benchmark):

    """
    StretchedV objective function.

    This class defines the Stretched V global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{StretchedV}}(\\mathbf{x}) = \sum_{i=1}^{n-1} t^{1/4} [\sin (50t^{0.1}) + 1]^2

    Where, in this exercise:

    .. math::

       t = x_{i+1}^2 + x_i^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [0., 0.]` when :math:`n = 2`.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10] * self.N, [10] * self.N)

        self.global_optimum = [[0, 0]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        t = x[1:] ** 2 + x[: -1] ** 2
#         TODO: fix equation in docs
        return sum(t ** 0.25 * (sin(50.0 * t ** 0.1) + 1) ** 2)


class StyblinskiTang(Benchmark):

    """
    StyblinskiTang objective function.

    This class defines the Styblinski-Tang global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{StyblinskiTang}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\left(x_i^4 - 16x_i^2 + 5x_i \\right)

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = -39.16616570377142n` for :math:`x_i = -2.903534018185960` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.0] * self.N, [5.0] * self.N)

        self.global_optimum = [[-2.903534018185960 for _ in range(self.N)]]
        self.fglob = -39.16616570377142 * self.N
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum(x ** 4 - 16 * x ** 2 + 5 * x) / 2
