# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import
import numpy as np
from numpy import (abs, arctan2, asarray, cos, exp, floor, log, log10,
                   arange, pi, prod, roll, sign, sin, sqrt, sum, where,
                   zeros, tan, tanh, dot)

from scipy.misc import factorial
from .go_benchmark import Benchmark


class NeedleEye(Benchmark):

    """
    NeedleEye objective function.

    This class defines the Needle-Eye global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{NeedleEye}}(\\mathbf{x}) = \\begin{cases} 1 & \\textrm{if} \\hspace{5pt} \\lvert x_i \\rvert  <  eye \\hspace{5pt} \\forall i \\\\
               \\sum_{i=1}^n (100 + \\lvert x_i \\rvert) & \\textrm{if} \\hspace{5pt} \\lvert x_i \\rvert > eye \\\\
               0 & \\textrm{otherwise} \\end{cases}

    Where, in this exercise, :math:`eye = 0.0001`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 1` for :math:`x_i = 0.` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 1.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        f = fp = 0.0
        eye = 0.0001

        for i in range(self.N):
            if abs(x[i]) >= eye:
                fp = 1.0
                f += 100.0 + abs(x[i])
            else:
                f += 1.0

        if fp < 1e-6:
            f = f / self.N

        return f


class NewFunction01(Benchmark):

    """
    NewFunction01 objective function.

    This class defines the NewFunction01 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{NewFunction01}}(\\mathbf{x}) = \\left | {\\cos\\left(\\sqrt{\\left|{x_{1}^{2} + x_{2}}\\right|}\\right)} \\right |^{0.5} + (x_{1} + x_{2})/100


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -0.18459899925` for :math:`\\mathbf{x} = [-8.46669057, -9.99982177]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[-8.46669057, -9.99982177]]
        self.fglob = -0.18459899925

    def fun(self, x, *args):
        self.nfev += 1

        return ((abs(cos(sqrt(abs(x[0] ** 2 + x[1]))))) ** 0.5
                + 0.01 * (x[0] + x[1]))


class NewFunction02(Benchmark):

    """
    NewFunction02 objective function.

    This class defines the NewFunction02 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{NewFunction02}}(\\mathbf{x}) = \\left | {\\sin\\left(\\sqrt{\\lvert{x_{1}^{2} + x_{2}}\\rvert}\\right)} \\right |^{0.5} + (x_{1} + x_{2})/100


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -0.19933159253` for :math:`\\mathbf{x} = [-9.94103375, -9.99771235]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[-9.94103375, -9.99771235]]
        self.fglob = -0.19933159253

    def fun(self, x, *args):
        self.nfev += 1

        return ((abs(sin(sqrt(abs(x[0] ** 2 + x[1]))))) ** 0.5
                + 0.01 * (x[0] + x[1]))


class NewFunction03(Benchmark):

    """
    NewFunction03 objective function.

    This class defines the NewFunction03 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{NewFunction03}}(\\mathbf{x}) = 0.01 x_{1} + 0.1 x_{2} + \\left\{x_{1} + \\sin^{2}\\left[\\left(\\cos\\left(x_{1}\\right) + \\cos\\left(x_{2}\\right)\\right)^{2}\\right] + \\cos^{2}\\left[\\left(\\sin\\left(x_{1}\\right) + \\sin\\left(x_{2}\\right)\\right)^{2}\\right]\\right\}^{2}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -1.019829` for :math:`\\mathbf{x} = [-1.98682, -10]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[-1.98682, -10.0]]
        self.fglob = -1.019829

    def fun(self, x, *args):
        self.nfev += 1

        f1 = sin((cos(x[0]) + cos(x[1])) ** 2) ** 2
        f2 = cos((sin(x[0]) + sin(x[1])) ** 2) ** 2
        f = (f1 + f2 + x[0]) ** 2
        f = f + 0.01 * x[0] + 0.1 * x[1]

        return f
