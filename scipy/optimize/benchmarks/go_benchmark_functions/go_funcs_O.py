# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import
import numpy as np
from numpy import (abs, arctan2, asarray, cos, exp, floor, log, log10,
                   arange, pi, prod, roll, sign, sin, sqrt, sum, where,
                   zeros, tan, tanh, dot)

from scipy.misc import factorial
from .go_benchmark import Benchmark


class OddSquare(Benchmark):

    """
    Odd Square objective function.

    This class defines the Odd Square global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{OddSquare}}(\\mathbf{x}) = -e^{-\\frac{d}{2\\pi}} \\cos(\\pi d) \\left( 1 + \\frac{0.02h}{d + 0.01} \\right )

    Where, in this exercise:

    .. math::

        \\begin{cases} d = n \\cdot \\smash{\\displaystyle\\max_{1 \leq i \leq n}} \\left[ (x_i - b_i)^2 \\right ] \\\\
        \\\\
        h = \\sum_{i=1}^{n} (x_i - b_i)^2 \\end{cases}

    And :math:`\\mathbf{b} = [1, 1.3, 0.8, -0.4, -1.3, 1.6, -0.2, -0.6, 0.5, 1.4, 1, 1.3, 0.8, -0.4, -1.3, 1.6, -0.2, -0.6, 0.5, 1.4]`

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5 \\pi, 5 \\pi]` for :math:`i=1,...,n` and :math:`n \\leq 20`.

    *Global optimum*: :math:`f(x_i) = -1.0084` for :math:`\\mathbf{x} \\approx b`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.0 * pi] * self.N,
                           [5.0 * pi] * self.N)
        self.custom_bounds = ([-2.0, 4.0], [-2.0, 4.0])
        self.a = asarray([1, 1.3, 0.8, -0.4, -1.3, 1.6, -0.2, -0.6, 0.5, 1.4]
                         * 2)
        self.global_optimum = [[1.09263477, 1.39263477]]

        self.fglob = -1.0084

    def fun(self, x, *args):
        self.nfev += 1
        b = self.a[0: self.N]
        d = self.N * max((x - b) ** 2.0)
        h = sum((x - b) ** 2.0)
        # TODO for n = 2, optimum = [1.09263477,  1.39263477].  The best solution
        # changes on dimensionality, but not the minimum energy
        return (-exp(-d / (2.0 * pi)) * cos(pi * d)
                * (1.0 + 0.02 * h / (d + 0.01)))
