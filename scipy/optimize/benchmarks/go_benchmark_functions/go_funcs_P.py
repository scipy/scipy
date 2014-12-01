# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import
import numpy as np
from numpy import (abs, arctan2, asarray, cos, exp, floor, log, log10,
                   arange, pi, prod, roll, sign, sin, sqrt, sum, where,
                   zeros, tan, tanh, dot)

from scipy.misc import factorial
from .go_benchmark import Benchmark


class Parsopoulos(Benchmark):

    """
    Parsopoulos objective function.

    This class defines the Parsopoulos global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Parsopoulos}}(\\mathbf{x}) = \\cos(x_1)^2 + \\sin(x_2)^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,2`.

    *Global optimum*: This function has innite number of global minima in R2, at points :math:`\\left(k\\frac{\\pi}{2}, \\lambda \\pi \\right)`,
    where :math:`k = \\pm1, \\pm3, ...` and :math:`\\lambda = 0, \\pm1, \\pm2, ...`

    In the given domain problem, function has 12 global minima all equal to zero.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.0] * self.N, [5.0] * self.N)

        self.global_optimum = [[pi / 2.0, pi]]
        self.fglob = 0

    def fun(self, x, *args):
        self.nfev += 1

        return cos(x[0]) ** 2.0 + sin(x[1]) ** 2.0


class Pathological(Benchmark):

    """
    Pathological objective function.

    This class defines the Pathological global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Pathological}}(\\mathbf{x}) = \\sum_{i=1}^{n -1} \\frac{\\sin^{2}\\left(\\sqrt{100 x_{i+1}^{2} + x_{i}^{2}}\\right) -0.5}{0.001 \\left(x_{i}^{2} - 2x_{i}x_{i+1} + x_{i+1}^{2}\\right)^{2} + 0.50}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0.` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.

    def fun(self, x, *args):
        self.nfev += 1

        vec = (0.5 + (sin(sqrt(100 * x[: -1] ** 2 + x[1:] ** 2)) ** 2 - 0.5) /
               (1. + 0.001 * (x[: -1] ** 2 - 2 * x[: -1] * x[1:]
                              + x[1:] ** 2) ** 2))
        return sum(vec)


class Paviani(Benchmark):

    """
    Paviani objective function.

    This class defines the Paviani global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Paviani}}(\\mathbf{x}) = \\sum_{i=1}^{10} \\left[\\log^{2}\\left(10 - x_i\\right) + \\log^{2}\\left(x_i -2\\right)\\right] - \\left(\\prod_{i=1}^{10} x_i^{10} \\right)^{0.2}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [2.001, 9.999]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = -45.7784684040686` for :math:`x_i = 9.350266` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=10):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([2.001] * self.N, [9.999] * self.N)

        self.global_optimum = [[9.350266 for _ in range(self.N)]]
        self.fglob = -45.7784684040686

    def fun(self, x, *args):
        self.nfev += 1

        return sum(log(x - 2) ** 2.0 + log(10.0 - x) ** 2.0) - prod(x) ** 0.2


class Penalty01(Benchmark):

    """
    Penalty 1 objective function.

    This class defines the Penalty 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Penalty01}}(\\mathbf{x}) = \\frac{\\pi}{30} \\left\\{10 \\sin^2(\\pi y_1) + \\sum_{i=1}^{n-1} (y_i - 1)^2 \\left[1 + 10 \\sin^2(\\pi y_{i+1}) \\right ] + (y_n - 1)^2 \\right \\} + \\sum_{i=1}^n u(x_i, 10, 100, 4)

    Where, in this exercise:

    .. math::

       y_i = 1 + \\frac{1}{4}(x_i + 1)

    And:

    .. math::

       u(x_i, a, k, m) = \\begin{cases} k(x_i - a)^m & \\textrm{if} \\hspace{5pt} x_i > a \\\\
       0 & \\textrm{if} \\hspace{5pt} -a \\leq x_i \\leq a \\\\
       k(-x_i - a)^m & \\textrm{if} \\hspace{5pt} x_i < -a \\end{cases}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-50, 50]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = -1` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-50.0] * self.N, [50.0] * self.N)
        self.custom_bounds = ([-5.0, 5.0], [-5.0, 5.0])

        self.global_optimum = [[-1.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        a, b, c = 10.0, 100.0, 4.0

        u = where(x > a, b * (x - a) ** c, 0.0)
        v = where(x < -a, b * (-x - a) ** c, 0.0)
        w = sum(u + v)

        y = 1.0 + (x + 1.0) / 4.0

        return (w + (pi / 30.0) * (10.0 * sin(pi * y[0]) ** 2.0
                + sum((y[: -1] - 1.0) ** 2.0
                      * (1.0 + 10.0 * sin(pi * y[1:]) ** 2.0))
                + (y[-1] - 1) ** 2.0))


class Penalty02(Benchmark):

    """
    Penalty 2 objective function.

    This class defines the Penalty 2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Penalty02}}(\\mathbf{x}) = 0.1 \\left\\{\\sin^2(3\\pi x_1) + \\sum_{i=1}^{n-1} (x_i - 1)^2 \\left[1 + \\sin^2(3\\pi x_{i+1}) \\right ] + (x_n - 1)^2 \\left [1 + \\sin^2(2 \\pi x_n) \\right ]\\right \\} + \\sum_{i=1}^n u(x_i, 5, 100, 4)

    Where, in this exercise:

    .. math::

       u(x_i, a, k, m) = \\begin{cases} k(x_i - a)^m & \\textrm{if} \\hspace{5pt} x_i > a \\\\
       0 & \\textrm{if} \\hspace{5pt} -a \\leq x_i \\leq a \\\\
       k(-x_i - a)^m & \\textrm{if} \\hspace{5pt} x_i < -a \\end{cases}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-50, 50]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 1` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-50.0] * self.N, [50.0] * self.N)
        self.custom_bounds = ([-4.0, 4.0], [-4.0, 4.0])

        self.global_optimum = [[1.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        a, b, c = 5.0, 100.0, 4.0

        u = where(x > a, b * (x - a) ** c, 0.0)
        v = where(x < -a, b * (-x - a) ** c, 0.0)
        w = sum(u + v)
        # TODO: is equation correct? If it's not then a factor of 10 needs to appear
        # before first sin**2 term.
        return (w + 0.1 * (sin(3.0 * pi * x[0]) ** 2.0
                + sum((x[:-1] - 1.0) ** 2.0
                      * (1.0 + sin(3 * pi * x[1:]) ** 2.0))
                + (x[-1] - 1) ** 2.0 * (1 + sin(2 * pi * x[-1]) ** 2.0)))


class PenHolder(Benchmark):

    """
    PenHolder objective function.

    This class defines the PenHolder global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{PenHolder}}(\\mathbf{x}) = -e^{\\left|{e^{-\\left|{- \\frac{\\sqrt{x_{1}^{2} + x_{2}^{2}}}{\\pi} + 1}\\right|} \\cos\\left(x_{1}\\right) \\cos\\left(x_{2}\\right)}\\right|^{-1}}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-11, 11]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -0.9635348327265058` for :math:`x_i = \\pm 9.646167671043401` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-11.0] * self.N, [11.0] * self.N)

        self.global_optimum = [[-9.646167708023526, 9.646167671043401]]
        self.fglob = -0.9635348327265058

    def fun(self, x, *args):
        self.nfev += 1

        a = abs(1. - (sqrt(x[0] ** 2 + x[1] ** 2) / pi))
        b = cos(x[0]) * cos(x[1]) * exp(a)
        return -exp(-abs(b) ** -1)


class PermFunction01(Benchmark):

    """
    PermFunction 1 objective function.

    This class defines the Perm Function 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{PermFunction01}}(\\mathbf{x}) = \\sum_{k=1}^n \\left\\{ \\sum_{j=1}^n (j^k + \\beta) \\left[ \\left(\\frac{x_j}{j}\\right)^k - 1 \\right] \\right\\}^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-n, n+1]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = i` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-self.N] * self.N,
                           [self.N + 1] * self.N)

        self.global_optimum = [range(1, self.N + 1)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        b = 0.5
        s_out = 0.0
        for k in range(1, self.N + 1):
            j = arange(1, self.N + 1)
            s_in = (j ** k + b) * ((x[j - 1] / j) ** k - 1)
            s_out += sum(s_in ** 2)

        return s_out


class PermFunction02(Benchmark):

    """
    PermFunction 2 objective function.

    This class defines the Perm Function 2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{PermFunction02}}(\\mathbf{x}) = \\sum_{k=1}^n \\left\\{ \\sum_{j=1}^n (j + \\beta) \\left[ \\left(x_j^k - {\\frac{1}{j}}^{k} \\right ) \\right] \\right\\}^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-n, n+1]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = \\frac{1}{i}` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-self.N] * self.N,
                           [self.N + 1] * self.N)
        self.custom_bounds = ([0, 1.5], [0, 1.0])

        self.global_optimum = [1. / arange(1, self.N + 1)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        b = 10
        outer = 0
        j = arange(1, self.N + 1)
        for k in range(1, self.N + 1):
            inner = (j + b) * (x[j - 1] ** k - (1. / j) ** k)
            outer += sum(inner ** 2)
        return outer


class Pinter(Benchmark):

    """
    Pinter objective function.

    This class defines the Pinter global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Pinter}}(\\mathbf{x}) = \\sum_{i=1}^n ix_i^2 + \\sum_{i=1}^n 20i \\sin^2 A + \\sum_{i=1}^n i \\log_{10} (1 + iB^2)


    Where, in this exercise:

    .. math::

        \\begin{cases} A = x_{i-1} \\sin x_i + \\sin x_{i+1} \\\\
        B = x_{i-1}^2 - 2x_i + 3x_{i+1} - \\cos x_i + 1 \\end{cases}

    Where :math:`x_0 = x_n` and :math:`x_{n+1} = x_1`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        f = 0.0

        for i in range(self.N):
            x_i = x[i]

            if i == 0:
                x_mi = x[-1]
                x_pi = x[i + 1]
            elif i == self.N - 1:
                x_mi = x[i - 1]
                x_pi = x[0]
            else:
                x_mi = x[i - 1]
                x_pi = x[i + 1]

            A = x_mi * sin(x_i) + sin(x_pi)
            B = x_mi ** 2.0 - 2 * x_i + 3 * x_pi - cos(x_i) + 1.0
# TODO: where does this come from, I only find different references fr
# this.

            f += (i + 1.0) * x_i ** 2.0 + 20.0 * (i + 1.0) * sin(A) ** 2.0 + \
                (i + 1.0) * log10(1.0 + (i + 1.0) * B ** 2.0)

        return f


class Plateau(Benchmark):

    """
    Plateau objective function.

    This class defines the Plateau global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Plateau}}(\\mathbf{x}) = 30 + \\sum_{i=1}^n \\lfloor x_i \\rfloor


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5.12, 5.12]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 30` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.12] * self.N, [5.12] * self.N)

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 30.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return 30.0 + sum(floor(abs(x)))


class Powell(Benchmark):

    """
    Powell objective function.

    This class defines the Powell global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Powell}}(\\mathbf{x}) = (x_3+10x_1)^2+5(x_2-x_4)^2+(x_1-2x_2)^4+10(x_3-x_4)^4


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-4, 5]` for :math:`i=1,...,4`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,4`

    """

    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-4.0] * self.N, [5.0] * self.N)
        self.global_optimum = [[0, 0, 0, 0]]
        self.fglob = 0

    def fun(self, x, *args):
        self.nfev += 1

        return ((x[0] + 10 * x[1]) ** 2 + 5 * (x[2] - x[3]) ** 2
                + (x[1] - 2 * x[2]) ** 4 + 10 * (x[0] - x[3]) ** 4)


class PowerSum(Benchmark):

    """
    Power sum objective function.

    This class defines the Power Sum global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{PowerSum}}(\\mathbf{x}) = \\sum_{k=1}^n\\left[\\left(\\sum_{i=1}^n x_i^k \\right) - b_k \\right]^2

    Where, in this exercise, :math:`\\mathbf{b} = [8, 18, 44, 114]`

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 4]` for :math:`i=1,...,4`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [1, 2, 2, 3]`

    """

    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N,
                           [float(self.N)] * self.N)

        self.global_optimum = [[1.0, 2.0, 2.0, 3.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        b = [8.0, 18.0, 44.0, 114.0]
        y = 0.0

        for k in range(1, self.N + 1):
            s_in = sum(x ** k)
            y += (s_in - b[k - 1]) ** 2.0

        return y


class Price01(Benchmark):

    """
    Price 1 objective function.

    This class defines the Price 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Price01}}(\\mathbf{x}) = (\\lvert x_1 \\rvert - 5)^2 + (\\lvert x_2 \\rvert - 5)^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-500, 500]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0.0` for :math:`\\mathbf{x} = [5, 5]` or :math:`\\mathbf{x} = [5, -5]`
    or :math:`\\mathbf{x} = [-5, 5]` or :math:`\\mathbf{x} = [-5, -5]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-500.0] * self.N,
                           [500.0] * self.N)
        self.custom_bounds = ([-10.0, 10.0], [-10.0, 10.0])

        self.global_optimum = [[5.0, 5.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return (abs(x[0]) - 5.0) ** 2.0 + (abs(x[1]) - 5.0) ** 2.0


class Price02(Benchmark):

    """
    Price 2 objective function.

    This class defines the Price 2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Price02}}(\\mathbf{x}) = 1 + \\sin^2(x_1) + \\sin^2(x_2) - 0.1e^{(-x_1^2 - x_2^2)}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0.9` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[0.0, 0.0]]
        self.fglob = 0.9

    def fun(self, x, *args):
        self.nfev += 1

        return 1.0 + sum(sin(x) ** 2) - 0.1 * exp(-x[0] ** 2.0 - x[1] ** 2.0)


class Price03(Benchmark):

    """
    Price 3 objective function.

    This class defines the Price 3 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Price03}}(\\mathbf{x}) = 100(x_2 - x_1^2)^2 + \\left[6.4(x_2 - 0.5)^2 - x_1 - 0.6 \\right]^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-50, 50]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [-5, -5]`, :math:`\\mathbf{x} = [-5, 5]`,
    :math:`\\mathbf{x} = [5, -5]`, :math:`\\mathbf{x} = [5, 5]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-50.0] * self.N, [50.0] * self.N)
        self.custom_bounds = ([0, 2], [0, 2])

        self.global_optimum = [[1.0, 1.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return (100 * (x[1] - x[0] ** 2) ** 2
                + (6.4 * (x[1] - 0.5) ** 2 - x[0] - 0.6) ** 2)


class Price04(Benchmark):

    """
    Price 4 objective function.

    This class defines the Price 4 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Price04}}(\\mathbf{x}) = (2x_1^3x_2 - x_2^3)^2 + (6x_1 - x_2^2 + x_2)^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-50, 50]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [0, 0]`, :math:`\\mathbf{x} = [2, 4]` and
    :math:`\\mathbf{x} = [1.464, -2.506]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-50.0] * self.N, [50.0] * self.N)
        self.custom_bounds = ([0, 2], [0, 2])

        self.global_optimum = [[2.0, 4.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return ((2.0 * x[1] * x[0] ** 3.0 - x[1] ** 3.0) ** 2.0
                + (6.0 * x[0] - x[1] ** 2.0 + x[1]) ** 2.0)
