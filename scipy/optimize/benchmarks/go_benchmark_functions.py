#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
=======================================================
:mod:`go_benchmark` -- Benchmark optimization functions
=======================================================

This module provides a set of benchmark problems for global optimization.

.. Copyright 2013 Andrea Gavana

.. module:: go_benchmark
.. moduleauthor:: Andrea Gavana <andrea.gavana@gmail.com>

"""
from __future__ import division
import os

import numpy as np

from numpy import abs, arange, arctan2, asarray, atleast_1d, cos, exp, floor, inf, log, ones, log10, arange
from numpy import pi, prod, roll, seterr, sign, sin, sqrt, sum, where, zeros, zeros_like, tan, tanh
from numpy import linspace, meshgrid, dot

from numpy.linalg import norm
from numpy.random import uniform

from scipy.misc import factorial

# Plotting stuff
import matplotlib
# matplotlib.use('AGG')

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import colors

import matplotlib.pyplot as plt

seterr(all='ignore')

# -------------------------------------------------------------------------------- #


class Benchmark(object):

    """
    Defines a global optimization benchmark problem.

    This abstract class defines the basic structure of a global
    optimization problem. Subclasses should implement the ``evaluator`` method
    for a particular optimization problem.

    Public Attributes:

    - *dimensions* -- the number of inputs to the problem
    - *fun_evals* -- stores the number of function evaluations, as some crappy
      optimization frameworks (i.e., `nlopt`) do not return this value
    - *change_dimensionality* -- whether we can change the benchmark function `x`
      variable length (i.e., the dimensionality of the problem)
    - *custom_bounds* -- a set of lower/upper bounds for plot purposes (if needed).
    - *spacing* -- the spacing to use to generate evenly spaced samples across the
      lower/upper bounds on the variables, for plotting purposes
    """

    def __init__(self, dimensions):
        self.dimensions = dimensions
        self.fun_evals = 0
        self.change_dimensionality = False
        self.custom_bounds = None

        if dimensions == 1:
            self.spacing = 1001
        else:
            self.spacing = 201

    def __str__(self):
        return '{0} ({1} dimensions)'.format(self.__class__.__name__, self.dimensions)

    def __repr__(self):
        return self.__class__.__name__

    def generator(self):
        """The generator function for the benchmark problem."""
        return [uniform(l, u) for l, u in self.bounds]

    def evaluator(self, candidates):
        """The evaluator function for the benchmark problem."""
        raise NotImplementedError

    def set_dimensions(self, ndim):
        self.dimensions = ndim

    def lower_bounds_constraints(self, x):

        lower = asarray([b[0] for b in self.bounds])
        return asarray(x) - lower

    def upper_bounds_constraints(self, x):

        upper = asarray([b[1] for b in self.bounds])
        return upper - asarray(x)

    def plot(self, set_title=False):

        if self.dimensions > 2:
            return

##        plt.rcParams['text.usetex'] = True
        fig = plt.figure()

        if self.custom_bounds:
            bounds = self.custom_bounds
        else:
            bounds = self.bounds

        xmin, xmax = bounds[0]

        if xmin < -1e4:
            xmin, xmax = -100, 100

        name = self.__class__.__name__
        spacing = self.spacing

        X = linspace(xmin, xmax, spacing)

        if self.dimensions == 2:
            # 3D functions
            ax = Axes3D(fig)
            ymin, ymax = bounds[1]

            if ymin < -1e4:
                ymin, ymax = -100, 100

            Y = linspace(ymin, ymax, spacing)

            X, Y = meshgrid(X, Y)
            Z = zeros(X.shape)

            for i in range(X.shape[0]):
                for j in range(X.shape[1]):
                    Z[i, j] = self.evaluator(asarray([X[i, j], Y[i, j]]))

            min_index = Z.argmin()
            cmap = matplotlib.cm.jet

            if np.any(numpy.isnan(Z)):
                Z = numpy.ma.array(Z, mask=numpy.isnan(Z))
                lev = linspace(Z.min(), Z.max(), self.spacing)
                norml = colors.BoundaryNorm(lev, 256)
                ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
                                cmap=cmap, linewidth=0.0, shade=True, norm=norml)
            else:
                ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
                                cmap=cmap, linewidth=0.0, shade=True)

            cset = plt.contour(X, Y, Z, zdir='z',
                               offset=ax.get_zlim()[0], alpha=0.3)

            ax.set_xlabel(r'$x_1$', fontsize=16)
            ax.set_ylabel(r'$x_2$', fontsize=16)
            ax.set_zlabel(r'$f(x_1, x_2)$', fontsize=16)

        else:
            # 2D functions
            ax = fig.add_subplot(111)
            Y = zeros(X.shape)

            for i in range(X.shape[0]):
                Y[i] = self.evaluator(atleast_1d(X[i]))

            ax.plot(X, Y, 'b-', lw=1.5, zorder=30)
            xf, yf = self.global_optimum, self.fglob

            ax.plot(xf, yf, 'r.', ms=11, zorder=40)

            ax.grid()
            ax.set_xlabel(r'$x$', fontsize=16)
            ax.set_ylabel(r'$f(x)$', fontsize=16)
            zlabels = []

        if set_title:
            ax.set_title(name + ' Test Function', fontweight='bold')

        out_folder = os.path.join(os.getcwd(), 'docs', 'figures')
        if not os.path.isdir(out_folder):
            os.makedirs(out_folder)

        filename = os.path.join(out_folder, '%s.png' % name)
        fig.savefig(filename)
        plt.close(fig)
        del fig


#-----------------------------------------------------------------------
#                     SINGLE-OBJECTIVE PROBLEMS
#-----------------------------------------------------------------------

class Ackley(Benchmark):

    """
    Ackley test objective function.

    This class defines the Ackley global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Ackley}}(\\mathbf{x}) = -20e^{-0.2 \\sqrt{\\frac{1}{n} \\sum_{i=1}^n x_i^2}} - e^{ \\frac{1}{n} \\sum_{i=1}^n \\cos(2 \\pi x_i)} + 20 + e

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-32, 32]` for :math:`i=1,...,n`.

    .. figure:: figures/Ackley.png
        :alt: Ackley function
        :align: center

        **Two-dimensional Ackley function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self.bounds = zip([-30.0] * self.dimensions, [30.0] * self.dimensions)

        self.global_optimum = [0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        a = 20.0
        b = 0.2
        c = 2.0 * pi
        return (-a * exp(-b * sqrt(1. / self.dimensions * sum(x ** 2)))
                - exp(1. / self.dimensions * sum(cos(c * x)))
                + a + exp(1.))


#-----------------------------------------------------------------------

class Adjiman(Benchmark):

    """
    Adjiman test objective function.

    This class defines the Adjiman global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Adjiman}}(\\mathbf{x}) = \\cos(x_1)\\sin(x_2) - \\frac{x_1}{(x_2^2 + 1)}

    Here, :math:`n` represents the number of dimensions and :math:`x_1 \\in [-1, 2]` and :math:`x_2 \\in [-1, 1]`.

    .. figure:: figures/Adjiman.png
        :alt: Adjiman function
        :align: center

        **Two-dimensional Adjiman function**


    *Global optimum*: :math:`f(x_i) = -2.02181` for :math:`\\mathbf{x} = [2, 0.10578]`
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self.bounds = ([-1.0, 2.0], [-1.0, 1.0])

        self.global_optimum = [2.0, 0.10578]
        self.fglob = -2.02180678

    def evaluator(self, x, *args):
        self.fun_evals += 1
        return cos(x[0]) * sin(x[1]) - x[0] / (x[1] ** 2 + 1)

# -------------------------------------------------------------------------------- #


class Alpine01(Benchmark):

    """
    Alpine 1 test objective function.

    This class defines the Alpine 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Alpine01}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\lvert {x_i \\sin \\left( x_i \\right) + 0.1 x_i} \\rvert

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    .. figure:: figures/Alpine01.png
        :alt: Alpine 1 function
        :align: center

        **Two-dimensional Alpine 1 function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return sum(abs(x * sin(x) + 0.1 * x))

# -------------------------------------------------------------------------------- #


class Alpine02(Benchmark):

    """
    Alpine 2 test objective function.

    This class defines the Alpine 2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Alpine02}}(\\mathbf{x}) = \\prod_{i=1}^{n} \\sqrt{x_i} \\sin(x_i)

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 10]` for :math:`i=1,...,n`.

    .. figure:: figures/Alpine02.png
        :alt: Alpine 2 function
        :align: center

        **Two-dimensional Alpine 2 function**


    *Global optimum*: :math:`f(x_i) = -6.1295` for :math:`x_i = 7.917` for :math:`i=1,...,n`
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([0.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [7.91705268, 4.81584232]
        self.fglob = -6.12950
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return prod(sqrt(x) * sin(x))

# -------------------------------------------------------------------------------- #


class AMGM(Benchmark):

    """
    AMGM test objective function.

    This class defines the Arithmetic Mean - Geometric Mean Equality global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{AMGM}}(\\mathbf{x}) = \\left ( \\frac{1}{n} \\sum_{i=1}^{n} x_i - \\sqrt[n]{ \\prod_{i=1}^{n} x_i} \\right )^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 10]` for :math:`i=1,...,n`.

    .. figure:: figures/AMGM.png
        :alt: AMGM function
        :align: center

        **Two-dimensional Arithmetic Mean - Geometric Mean Equality function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_1 = x_2 = ... = x_n` for :math:`i=1,...,n`
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([0.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [1, 1]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        n = self.dimensions

        f1 = sum(x)
        f2 = prod(x)
        f1 = f1 / n
        f2 = f2 ** (1.0 / n)
        f = (f1 - f2) ** 2

        return f

# -------------------------------------------------------------------------------- #


class BartelsConn(Benchmark):

    """
    Bartels-Conn test objective function.

    This class defines the Bartels-Conn global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{BartelsConn}}(\\mathbf{x}) = \\lvert {x_1^2 + x_2^2 + x_1x_2} \\rvert + \\lvert {\\sin(x_1)} \\rvert + \\lvert {\\cos(x_2)} \\rvert

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-50, 50]` for :math:`i=1,...,n`.

    .. figure:: figures/BartelsConn.png
        :alt: Bartels-Conn function
        :align: center

        **Two-dimensional Bartels-Conn function**


    *Global optimum*: :math:`f(x_i) = 1` for :math:`x_i = 0` for :math:`i=1,...,n`
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-5.0] * self.dimensions, [5.0] * self.dimensions)

        self.global_optimum = [0 for _ in range(self.dimensions)]
        self.fglob = 1.0

    def evaluator(self, x, *args):
        self.fun_evals += 1

        return (abs(x[0] ** 2.0 + x[1] ** 2.0 + x[0] * x[1]) + abs(sin(x[1]))
                + abs(cos(x[1])))

# -------------------------------------------------------------------------------- #


class Beale(Benchmark):

    """
    Beale test objective function.

    This class defines the Beale global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Beale}}(\\mathbf{x}) = \\left(x_1 x_2 - x_1 + 1.5\\right)^{2} + \\left(x_1 x_2^{2} - x_1 + 2.25\\right)^{2} + \\left(x_1 x_2^{3} - x_1 + 2.625\\right)^{2}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    .. figure:: figures/Beale.png
        :alt: Beale function
        :align: center

        **Two-dimensional Beale function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [3, 0.5]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-4.5] * self.dimensions, [4.5] * self.dimensions)

        self.global_optimum = [3.0, 0.5]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return ((1.5 - x[0] + x[0] * x[1]) ** 2
                + (2.25 - x[0] + x[0] * x[1] ** 2) ** 2
                + (2.625 - x[0] + x[0] * x[1] ** 3) ** 2)

# -------------------------------------------------------------------------------- #


class Bird(Benchmark):

    """
    Bird test objective function.

    This class defines the Bird global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Bird}}(\\mathbf{x}) = \\left(x_1 - x_2\\right)^{2} + e^{\left[1 - \\sin\\left(x_1\\right) \\right]^{2}} \\cos\\left(x_2\\right) + e^{\left[1 - \\cos\\left(x_2\\right)\\right]^{2}} \\sin\\left(x_1\\right)

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-2\\pi, 2\\pi]` for :math:`i=1,2`.

    .. figure:: figures/Bird.png
        :alt: Bird function
        :align: center

        **Two-dimensional Bird function**


    *Global optimum*: :math:`f(x_i) = -106.7645367198034` for :math:`\\mathbf{x} = [4.701055751981055 , 3.152946019601391]` or
    :math:`\\mathbf{x} = [-1.582142172055011, -3.130246799635430]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-2.0 * pi] * self.dimensions,
                          [2.0 * pi] * self.dimensions)

        self.global_optimum = ([4.701055751981055, 3.152946019601391],
                               [-1.582142172055011, -3.130246799635430])
        self.fglob = -106.7645367198034

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return (sin(x[0]) * exp((1 - cos(x[1])) ** 2)
                + cos(x[1]) * exp((1 - sin(x[0])) ** 2) + (x[0] - x[1]) ** 2)

# -------------------------------------------------------------------------------- #


class Bohachevsky(Benchmark):

    """
    Bohachevsky test objective function.

    This class defines the Bohachevsky global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Bohachevsky}}(\\mathbf{x}) = \\sum_{i=1}^{n-1}\\left[x_i^2 + 2x_{i+1}^2 - 0.3\\cos(3\\pi x_i) - 0.4\\cos(4\\pi x_{i+1}) + 0.7\\right]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-15, 15]` for :math:`i=1,...,n`.

    .. figure:: figures/Bohachevsky.png
        :alt: Bohachevsky function
        :align: center

        **Two-dimensional Bohachevsky function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-15.0] * self.dimensions, [15.0] * self.dimensions)
        self.custom_bounds = [(-2, 2), (-2, 2)]

        self.global_optimum = [0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x0 = x[:-1]
        x1 = roll(x, -1)[:-1]

        return sum(x0 ** 2 + 2 * x1 ** 2 - 0.3 * cos(3 * pi * x0)
                   - 0.4 * cos(4 * pi * x1) + 0.7)

# -------------------------------------------------------------------------------- #


class BoxBetts(Benchmark):

    """
    BoxBetts test objective function.

    This class defines the Box-Betts global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{BoxBetts}}(\\mathbf{x}) = \\sum_{i=1}^k g(x_i)^2

    Where, in this exercise:

    .. math:: g(x) = e^{-0.1(i+1)x_1} - e^{-0.1(i+1)x_2} - \\left[(e^{-0.1(i+1)}) - e^{-(i+1)}x_3\\right]


    And :math:`k = 10`.

    Here, :math:`n` represents the number of dimensions and :math:`x_1 \\in [0.9, 1.2], x_2 \\in [9, 11.2], x_3 \\in [0.9, 1.2]`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [1, 10, 1]`

    """

    def __init__(self, dimensions=3):
        Benchmark.__init__(self, dimensions)

        self.bounds = ([0.9, 1.2], [9.0, 11.2], [0.9, 1.2])

        self.global_optimum = [1.0, 10.0, 1.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        y = 0.0

        for i in range(1, 11):
            y += (exp(-0.1 * i * x[0]) - exp(-0.1 * i * x[1])
                  - (exp(-0.1 * i) - exp(-1.0 * i)) * x[2]) ** 2.0

        return y

# -------------------------------------------------------------------------------- #


class Branin01(Benchmark):

    """
    Branin 1 test objective function.

    This class defines the Branin 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Branin01}}(\\mathbf{x}) = \\left(- 1.275 \\frac{x_1^{2}}{\pi^{2}} + 5 \\frac{x_1}{\pi} + x_2 -6\\right)^{2} + \\left(10 - \\frac{5}{4 \\pi} \\right) \\cos\\left(x_1\\right) + 10

    Here, :math:`n` represents the number of dimensions and :math:`x_1 \\in [-5, 10], x_2 \\in [0, 15]`

    .. figure:: figures/Branin01.png
        :alt: Branin 1 function
        :align: center

        **Two-dimensional Branin 1 function**


    *Global optimum*: :math:`f(x_i) = 0.39788735772973816` for :math:`\\mathbf{x} = [-\\pi, 12.275]` or
    :math:`\\mathbf{x} = [\\pi, 2.275]` or :math:`\\mathbf{x} = [9.42478, 2.475]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(-5., 10.), (0., 15.)]

        self.global_optimum = [(-pi, 12.275), (pi, 2.275), (9.42478, 2.475)]
        self.fglob = 0.39788735772973816

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return (x[1] - (5.1 / (4 * pi ** 2)) * x[0] ** 2 + 5 * x[0] / pi - 6) ** 2 + 10 * (1 - 1 / (8 * pi)) * cos(x[0]) + 10

# -------------------------------------------------------------------------------- #


class Branin02(Benchmark):

    """
    Branin 2 test objective function.

    This class defines the Branin 2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Branin02}}(\\mathbf{x}) = \\left(- 1.275 \\frac{x_1^{2}}{\pi^{2}} + 5 \\frac{x_1}{\pi} + x_2 -6\\right)^{2} + \\left(10 - \\frac{5}{4 \\pi} \\right) \\cos\\left(x_1\\right) \\cos\\left(x_2\\right) + \\log(x_1^2+x_2^2 +1) + 10

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 15]` for :math:`i=1,2`.

    .. figure:: figures/Branin02.png
        :alt: Branin 2 function
        :align: center

        **Two-dimensional Branin 2 function**


    *Global optimum*: :math:`f(x_i) = 5.559037` for :math:`\\mathbf{x} = [-3.2, 12.53]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(-5.0, 15.0), (-5.0, 15.0)]

        self.global_optimum = [-3.2, 12.53]
        self.fglob = 5.559037

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return (x[1] - (5.1 / (4 * pi ** 2)) * x[0] ** 2 + 5 * x[0] / pi - 6) ** 2 + 10 * (1 - 1 / (8 * pi)) * cos(x[0]) * cos(x[1]) + log(x[0] ** 2.0 + x[1] ** 2.0 + 1.0) + 10

# -------------------------------------------------------------------------------- #


class Brent(Benchmark):

    """
    Brent test objective function.

    This class defines the Brent global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Brent}}(\\mathbf{x}) = (x_1 + 10)^2 + (x_2 + 10)^2 + e^{(-x_1^2-x_2^2)}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    .. figure:: figures/Brent.png
        :alt: Brent function
        :align: center

        **Two-dimensional Brent function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [-10, -10]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)
        self.custom_bounds = ([-10, 2], [-10, 2])

        self.global_optimum = [-10.0, -10.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return (x[0] + 10.0) ** 2.0 + (x[1] + 10.0) ** 2.0 + exp(-x[0] ** 2.0 - x[1] ** 2.0)

# -------------------------------------------------------------------------------- #


class Brown(Benchmark):

    """
    Brown test objective function.

    This class defines the Brown global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Brown}}(\\mathbf{x}) = \\sum_{i=1}^{n-1}\\left[ \\left(x_i^2\\right)^{x_{i+1}^2+1} + \\left(x_{i+1}^2\\right)^{x_i^2+1} \\right]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 4]` for :math:`i=1,...,n`.

    .. figure:: figures/Brown.png
        :alt: Brown function
        :align: center

        **Two-dimensional Brown function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-1.0] * self.dimensions, [4.0] * self.dimensions)
        self.custom_bounds = ([-1.0, 1.0], [-1.0, 1.0])

        self.global_optimum = [0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x0 = x[:-1]
        x1 = x[1:]
        return sum((x0 ** 2.0) ** (x1 ** 2.0 + 1.0) + (x1 ** 2.0) ** (x0 ** 2.0 + 1.0))

# -------------------------------------------------------------------------------- #


class Bukin02(Benchmark):

    """
    Bukin 2 test objective function.

    This class defines the Bukin 2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Bukin02}}(\\mathbf{x}) = 100 (x_2 - 0.01x_1^2 + 1) + 0.01(x_1 + 10)^2

    Here, :math:`n` represents the number of dimensions and :math:`x_1 \\in [-15, -5], x_2 \\in [-3, 3]`

    .. figure:: figures/Bukin02.png
        :alt: Bukin 2 function
        :align: center

        **Two-dimensional Bukin 2 function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [-10, 0]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(-15.0, -5.0), (-3.0, 3.0)]

        self.global_optimum = [-10.0, 0.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return 100 * (x[1] ** 2 - 0.01 * x[0] ** 2 + 1.0) + 0.01 * (x[0] + 10.0) ** 2.0

# -------------------------------------------------------------------------------- #


class Bukin04(Benchmark):

    """
    Bukin 4 test objective function.

    This class defines the Bukin 4 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Bukin04}}(\\mathbf{x}) = 100 x_2^{2} + 0.01 \\lvert{x_1 + 10} \\rvert

    Here, :math:`n` represents the number of dimensions and :math:`x_1 \\in [-15, -5], x_2 \\in [-3, 3]`

    .. figure:: figures/Bukin04.png
        :alt: Bukin 4 function
        :align: center

        **Two-dimensional Bukin 4 function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [-10, 0]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(-15.0, -5.0), (-3.0, 3.0)]

        self.global_optimum = [-10.0, 0.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return 100 * x[1] ** 2 + 0.01 * abs(x[0] + 10)

# -------------------------------------------------------------------------------- #


class Bukin06(Benchmark):

    """
    Bukin 6 test objective function.

    This class defines the Bukin 6 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Bukin06}}(\\mathbf{x}) = 100 \\sqrt{ \\lvert{x_2 - 0.01 x_1^{2}} \\rvert} + 0.01 \\lvert{x_1 + 10} \\rvert

    Here, :math:`n` represents the number of dimensions and :math:`x_1 \\in [-15, -5], x_2 \\in [-3, 3]`

    .. figure:: figures/Bukin06.png
        :alt: Bukin 6 function
        :align: center

        **Two-dimensional Bukin 6 function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [-10, 1]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(-15.0, -5.0), (-3.0, 3.0)]

        self.global_optimum = [-10.0, 1.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return 100 * sqrt(abs(x[1] - 0.01 * x[0] ** 2)) + 0.01 * abs(x[0] + 10)

# -------------------------------------------------------------------------------- #


class CarromTable(Benchmark):

    """
    CarromTable test objective function.

    This class defines the CarromTable global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{CarromTable}}(\\mathbf{x}) = - \\frac{1}{30} e^{2 \\left|{1 - \\frac{\\sqrt{x_{1}^{2} + x_{2}^{2}}}{\pi}}\\right|} \\cos^{2}\\left(x_{1}\\right) \\cos^{2}\\left(x_{2}\\right)

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    .. figure:: figures/CarromTable.png
        :alt: CarromTable function
        :align: center

        **Two-dimensional CarromTable function**


    *Global optimum*: :math:`f(x_i) = -24.15681551650653` for :math:`x_i = \\pm 9.646157266348881` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [(9.646157266348881, 9.646134286497169),
                               (-9.646157266348881, 9.646134286497169),
                               (9.646157266348881, -9.646134286497169),
                               (-9.646157266348881, -9.646134286497169)]
        self.fglob = -24.15681551650653

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return -((cos(x[0]) * cos(x[1]) * exp(abs(1 - sqrt(x[0] ** 2 + x[1] ** 2) / pi))) ** 2) / 30

# -------------------------------------------------------------------------------- #


class Chichinadze(Benchmark):

    """
    Chichinadze test objective function.

    This class defines the Chichinadze global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Chichinadze}}(\\mathbf{x}) = x_{1}^{2} - 12 x_{1} + 8 \\sin\\left(\\frac{5}{2} \\pi x_{1}\\right) + 10 \\cos\\left(\\frac{1}{2} \\pi x_{1}\\right) + 11 - 0.2 \\frac{\\sqrt{5}}{e^{\\frac{1}{2} \\left(x_{2} -0.5\\right)^{2}}}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-30, 30]` for :math:`i=1,2`.

    .. figure:: figures/Chichinadze.png
        :alt: Chichinadze function
        :align: center

        **Two-dimensional Chichinadze function**


    *Global optimum*: :math:`f(x_i) = -42.94438701899098` for :math:`\\mathbf{x} = [6.189866586965680, 0.5]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-30.0] * self.dimensions, [30.0] * self.dimensions)
        self.custom_bounds = [(-10, 10), (-10, 10)]

        self.global_optimum = [6.189866586965680, 0.5]
        self.fglob = -42.94438701899098

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return x[0] ** 2 - 12 * x[0] + 11 + 10 * cos(pi * x[0] / 2) + 8 * sin(5 * pi * x[0] / 2) - 1.0 / sqrt(5) * exp(-((x[1] - 0.5) ** 2) / 2)

# -------------------------------------------------------------------------------- #


class Cigar(Benchmark):

    """
    Cigar test objective function.

    This class defines the Cigar global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Cigar}}(\\mathbf{x}) = x_1^2 + 10^6\\sum_{i=2}^{n} x_i^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    .. figure:: figures/Cigar.png
        :alt: Cigar function
        :align: center

        **Two-dimensional Cigar function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-100.0] * self.dimensions,
                          [100.0] * self.dimensions)
        self.custom_bounds = [(-5, 5), (-5, 5)]

        self.global_optimum = [0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return x[0] ** 2 + 1e6 * sum(x[1:] ** 2)

# -------------------------------------------------------------------------------- #


class Cola(Benchmark):

    """
    Cola test objective function.

    This class defines the Cola global optimization problem. The 17-dimensional function computes
    indirectly the formula :math:`f(n, u)` by setting :math:`x_0 = y_0, x_1 = u_0, x_i = u_{2(iâˆ’2)}, y_i = u_{2(iâˆ’2)+1}` :

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

    This function has bounds :math:`0 \\leq x_0 \\leq 4` and :math:`-4 \\leq x_i \\leq 4` for :math:`i = 1,...,n-1`. It
    has a global minimum of 11.7464.
    """

    def __init__(self, dimensions=17):
        Benchmark.__init__(self, dimensions)

        self.bounds = [[0.0, 4.0]] + \
            list(zip([-4.0] * (self.dimensions - 1),
                 [4.0] * (self.dimensions - 1)))

        self.global_optimum = [0.651906, 1.30194, 0.099242, -0.883791, -0.8796,
                               0.204651, -3.28414, 0.851188, -
                               3.46245, 2.53245, -0.895246,
                               1.40992, -3.07367, 1.96257, -2.97872, -0.807849, -1.68978]
        self.fglob = 11.7464

    def evaluator(self, x, *args):

        self.fun_evals += 1

        # C implementation - doesn't work
# dis = [1.27,
##               1.69, 1.43,
##               2.04, 2.35, 2.43,
##               3.09, 3.18, 3.26, 2.85,
##               3.20, 3.22, 3.27, 2.88, 1.55,
##               2.86, 2.56, 2.58, 2.59, 3.12, 3.06,
##               3.17, 3.18, 3.18, 3.12, 1.31, 1.64, 3.00,
##               3.21, 3.18, 3.18, 3.17, 1.70, 1.36, 2.95, 1.32,
# 2.38, 2.31, 2.42, 1.94, 2.85, 2.81, 2.56, 2.91, 2.97]
##
##        s = 0.0
##        k = 1
##        mt = zeros((20, ))
##        mt[4:] = x[1:]
##
# for i in range(1, 10):
# for j in range(i):
##                temp = 0.0
# for t in range(2):
##                    temp += (mt[i*2+t] - mt[j*2+t])**2.0
##                s += (dis[k-1] - sqrt(temp))**2.0
##                k += 1
##
# return s

        # Scilab implementation

        d = asarray([[0, 0,  0,  0,  0,  0,  0,  0,  0],
                     [1.27, 0,  0,  0,  0,  0,  0,  0,  0],
                     [1.69, 1.43, 0,  0,  0,  0,  0,  0,  0],
                     [2.04, 2.35, 2.43, 0,    0,    0,    0,    0,    0],
                     [3.09, 3.18, 3.26, 2.85, 0,    0,    0,    0,    0],
                     [3.20, 3.22, 3.27, 2.88, 1.55, 0,    0,    0,    0],
                     [2.86, 2.56, 2.58, 2.59, 3.12, 3.06, 0,    0,    0],
                     [3.17, 3.18, 3.18, 3.12, 1.31, 1.64, 3.00, 0,    0],
                     [3.21, 3.18, 3.18, 3.17, 1.70, 1.36, 2.95, 1.32, 0],
                     [2.38, 2.31, 2.42, 1.94, 2.85, 2.81, 2.56, 2.91, 2.97]])

        x1 = asarray([0.0, x[0]] + list(x[1::2]))
        x2 = asarray([0.0, 0.0] + list(x[2::2]))
        y = 0.0

        for i in range(1, len(x1)):
            y += sum((sqrt((x1[i] - x1[0:i]) ** 2.0 +
                     (x2[i] - x2[0:i]) ** 2.0) - d[i, 0:i]) ** 2.0)

# for j in range(i):
##                y += (sqrt((x1[i] - x1[j])**2.0 + (x2[i] - x2[j])**2.0) - d[i, j])**2.0

        return y

# -------------------------------------------------------------------------------- #


class Colville(Benchmark):

    """
    Colville test objective function.

    This class defines the Colville global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Colville}}(\\mathbf{x}) = \\left(x_{1} -1\\right)^{2} + 100 \\left(x_{1}^{2} - x_{2}\\right)^{2} + 10.1 \\left(x_{2} -1\\right)^{2} + \\left(x_{3} -1\\right)^{2} + 90 \\left(x_{3}^{2} - x_{4}\\right)^{2} + 10.1 \\left(x_{4} -1\\right)^{2} + 19.8 \\frac{x_{4} -1}{x_{2}}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,4`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 1` for :math:`i=1,...,4`

    """

    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [1 for _ in range(self.dimensions)]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return 100 * (x[0] ** 2 - x[1]) ** 2 + (x[0] - 1) ** 2 + (x[2] - 1) ** 2 + 90 * (x[2] ** 2 - x[3]) ** 2 + 10.1 * ((x[1] - 1) ** 2 + (x[3] - 1) ** 2) + 19.8 * (1 / x[1]) * (x[3] - 1)

# -------------------------------------------------------------------------------- #


class Corana(Benchmark):

    """
    Corana test objective function.

    This class defines the Corana global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Corana}}(\\mathbf{x}) = \\begin{cases} \\sum_{i=1}^n 0.15 d_i [z_i - 0.05\\textrm{sgn}(z_i)]^2 & \\textrm{if}|x_i-z_i| < 0.05 \\\\
               d_ix_i^2 & \\textrm{otherwise}\\end{cases}

    Where, in this exercise:

    .. math::

        z_i = 0.2 \\lfloor |x_i/s_i|+0.49999\\rfloor\\textrm{sgn}(x_i), d_i=(1,1000,10,100, ...)


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,...,4`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,4`

    """

    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-5.0] * self.dimensions, [5.0] * self.dimensions)

        self.global_optimum = [0 for _ in range(self.dimensions)]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        d = [1., 1000., 10., 100.]
        r = 0
        for j in range(4):
            zj = floor(abs(x[j] / 0.2) + 0.49999) * sign(x[j]) * 0.2
            if abs(x[j] - zj) < 0.05:
                r += 0.15 * ((zj - 0.05 * sign(zj)) ** 2) * d[j]
            else:
                r += d[j] * x[j] * x[j]
        return r

# -------------------------------------------------------------------------------- #


class CosineMixture(Benchmark):

    """
    Cosine Mixture test objective function.

    This class defines the Cosine Mixture global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{CosineMixture}}(\\mathbf{x}) = -0.1 \\sum_{i=1}^n \\cos(5 \\pi x_i) - \\sum_{i=1}^n x_i^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,...,N`.

    .. figure:: figures/CosineMixture.png
        :alt: Cosine Mixture function
        :align: center

        **Two-dimensional Cosine Mixture function**

    *Global optimum*: :math:`f(x_i) = -0.1N` for :math:`x_i = 0` for :math:`i=1,...,N`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.change_dimensionality = True
        self.bounds = zip([-1.0] * self.dimensions, [1.0] * self.dimensions)

        self.global_optimum = [0 for _ in range(self.dimensions)]
        self.fglob = -0.1 * self.dimensions

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return -0.1 * sum(cos(5.0 * pi * x)) - sum(x ** 2.0)

# -------------------------------------------------------------------------------- #


class CrossInTray(Benchmark):

    """
    Cross-in-Tray test objective function.

    This class defines the Cross-in-Tray global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{CrossInTray}}(\\mathbf{x}) = - 0.0001 \\left(\\left|{e^{\\left|{100 - \\frac{\\sqrt{x_{1}^{2} + x_{2}^{2}}}{\\pi}}\\right|} \\sin\\left(x_{1}\\right) \\sin\\left(x_{2}\\right)}\\right| + 1\\right)^{0.1}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-15, 15]` for :math:`i=1,2`.

    .. figure:: figures/CrossInTray.png
        :alt: Cross-in-Tray function
        :align: center

        **Two-dimensional Cross-in-Tray function**


    *Global optimum*: :math:`f(x_i) = -2.062611870822739` for :math:`x_i = \\pm 1.349406608602084` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [(1.349406685353340, 1.349406608602084),
                               (-1.349406685353340, 1.349406608602084),
                               (1.349406685353340, -1.349406608602084),
                               (-1.349406685353340, -1.349406608602084)]
        self.fglob = -2.062611870822739

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return -0.0001 * (abs(sin(x[0]) * sin(x[1]) * exp(abs(100 - sqrt(x[0] ** 2 + x[1] ** 2) / pi))) + 1) ** (0.1)

# -------------------------------------------------------------------------------- #


class CrossLegTable(Benchmark):

    """
    Cross-Leg-Table test objective function.

    This class defines the Cross-Leg-Table global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{CrossLegTable}}(\\mathbf{x}) = - \\frac{1}{\\left(\\left|{e^{\\left|{100 - \\frac{\\sqrt{x_{1}^{2} + x_{2}^{2}}}{\\pi}}\\right|} \\sin\\left(x_{1}\\right) \\sin\\left(x_{2}\\right)}\\right| + 1\\right)^{0.1}}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    .. figure:: figures/CrossLegTable.png
        :alt: Cross-Leg-Table function
        :align: center

        **Two-dimensional Cross-Leg-Table function**


    *Global optimum*: :math:`f(x_i) = -1`. The global minimum is found on the planes :math:`x_1 = 0` and :math:`x_2 = 0`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [0., 0.]
        self.fglob = -1.0

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return -(abs(sin(x[0]) * sin(x[1]) * exp(abs(100 - sqrt(x[0] ** 2 + x[1] ** 2) / pi))) + 1) ** (-0.1)

# -------------------------------------------------------------------------------- #


class CrownedCross(Benchmark):

    """
    Crowned Cross test objective function.

    This class defines the Crowned Cross global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{CrownedCross}}(\\mathbf{x}) = 0.0001 \\left(\\left|{e^{\\left|{100- \\frac{\\sqrt{x_{1}^{2} + x_{2}^{2}}}{\\pi}}\\right|} \\sin\\left(x_{1}\\right) \\sin\\left(x_{2}\\right)}\\right| + 1\\right)^{0.1}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    .. figure:: figures/CrownedCross.png
        :alt: Crowned Cross function
        :align: center

        **Two-dimensional Crowned Cross function**


    *Global optimum*: :math:`f(x_i) = 0.0001`. The global minimum is found on the planes :math:`x_1 = 0` and :math:`x_2 = 0`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [0, 0]
        self.fglob = 0.0001

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return 0.0001 * (abs(sin(x[0]) * sin(x[1]) * exp(abs(100 - sqrt(x[0] ** 2 + x[1] ** 2) / pi))) + 1) ** (0.1)

# -------------------------------------------------------------------------------- #


class Csendes(Benchmark):

    """
    Csendes test objective function.

    This class defines the Csendes global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Csendes}}(\\mathbf{x}) = \\sum_{i=1}^n x_i^6 \\left[ 2 + \\sin \\left( \\frac{1}{x_i} \\right ) \\right]


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,...,N`.

    .. figure:: figures/Csendes.png
        :alt: Csendes function
        :align: center

        **Two-dimensional Csendes function**

    *Global optimum*: :math:`f(x_i) = 0.0` for :math:`x_i = 0` for :math:`i=1,...,N`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.change_dimensionality = True
        self.bounds = zip([-1.0] * self.dimensions, [1.0] * self.dimensions)

        self.global_optimum = [0 for _ in range(self.dimensions)]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return sum((x ** 6.0) * (2.0 + sin(1.0 / x)))

# -------------------------------------------------------------------------------- #


class Cube(Benchmark):

    """
    Cube test objective function.

    This class defines the Cube global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Cube}}(\\mathbf{x}) = 100(x_2 - x_1^3)^2 + (1 - x1)^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,N`.

    .. figure:: figures/Cube.png
        :alt: Cube function
        :align: center

        **Two-dimensional Cube function**

    *Global optimum*: :math:`f(x_i) = 0.0` for :math:`\\mathbf{x} = [1, 1]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)
        self.custom_bounds = ([0, 2], [0, 2])

        self.global_optimum = [1.0, 1.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return 100.0 * (x[1] - x[0] ** 3.0) ** 2.0 + (1.0 - x[0]) ** 2.0

# -------------------------------------------------------------------------------- #


class Damavandi(Benchmark):

    """
    Damavandi test objective function.

    This class defines the Damavandi global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Damavandi}}(\\mathbf{x}) = \\left[ 1 - \\lvert{\\frac{\\sin[\\pi(x_1-2)]\\sin[\\pi(x2-2)]}{\\pi^2(x_1-2)(x_2-2)}} \\rvert^5 \\right] \\left[2 + (x_1-7)^2 + 2(x_2-7)^2 \\right]


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 14]` for :math:`i=1,...,n`.

    .. figure:: figures/Damavandi.png
        :alt: Damavandi function
        :align: center

        **Two-dimensional Damavandi function**

    *Global optimum*: :math:`f(x_i) = 0.0` for :math:`x_i = 2` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.change_dimensionality = True
        self.bounds = zip([0.0] * self.dimensions, [14.0] * self.dimensions)

        self.global_optimum = [2 for _ in range(self.dimensions)]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        num = sin(pi * (x[0] - 2.0)) * sin(pi * (x[1] - 2.0))
        den = (pi ** 2) * (x[0] - 2.0) * (x[1] - 2.0)
        factor1 = 1.0 - (abs(num / den)) ** 5.0
        factor2 = 2 + (x[0] - 7.0) ** 2.0 + 2 * (x[1] - 7.0) ** 2.0

        return factor1 * factor2

# -------------------------------------------------------------------------------- #


class Deb01(Benchmark):

    """
    Deb 1 test objective function.

    This class defines the Deb 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Deb01}}(\\mathbf{x}) = - \\frac{1}{N} \\sum_{i=1}^n \\sin^6(5 \\pi x_i)


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,...,n`.

    .. figure:: figures/Deb01.png
        :alt: Deb 1 function
        :align: center

        **Two-dimensional Deb 1 function**

    *Global optimum*: :math:`f(x_i) = 0.0`. The number of global minima is :math:`5^n` that are evenly spaced
    in the function landscape, where :math:`n` represents the dimension of the problem.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.change_dimensionality = True

        self.bounds = zip([-1.0] * self.dimensions, [1.0] * self.dimensions)

        self.global_optimum = [0.3, -0.3]
        self.fglob = -1.0

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return -(1.0 / self.dimensions) * sum(sin(5 * pi * x) ** 6.0)

# -------------------------------------------------------------------------------- #


class Deb02(Benchmark):

    """
    Deb 2 test objective function.

    This class defines the Deb 2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Deb02}}(\\mathbf{x}) = - \\frac{1}{N} \\sum_{i=1}^n \\sin^6 \\left[ 5 \\pi \\left ( x_i^{3/4} - 0.05 \\right) \\right ]


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 1]` for :math:`i=1,...,n`.

    .. figure:: figures/Deb02.png
        :alt: Deb 2 function
        :align: center

        **Two-dimensional Deb 2 function**

    *Global optimum*: :math:`f(x_i) = 0.0`. The number of global minima is :math:`5^n` that are evenly spaced
    in the function landscape, where :math:`n` represents the dimension of the problem.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.change_dimensionality = True

        self.bounds = zip([0.0] * self.dimensions, [1.0] * self.dimensions)

        self.global_optimum = [0.93388314, 0.68141781]
        self.fglob = -1.0

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return -(1.0 / self.dimensions) * sum(sin(5 * pi * (x ** 0.75 - 0.05)) ** 6.0)

# -------------------------------------------------------------------------------- #


class Decanomial(Benchmark):

    """
    Decanomial test objective function.

    This class defines the Decanomial function global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Decanomial}}(\\mathbf{x}) = 0.001 \\left(\\lvert{x_{2}^{4} + 12 x_{2}^{3} + 54 x_{2}^{2} + 108 x_{2} + 81.0}\\rvert + \\lvert{x_{1}^{10} - 20 x_{1}^{9} + 180 x_{1}^{8} - 960 x_{1}^{7} + 3360 x_{1}^{6} - 8064 x_{1}^{5} + 13340 x_{1}^{4} - 15360 x_{1}^{3} + 11520 x_{1}^{2} - 5120 x_{1} + 2624.0}\\rvert\\right)^{2}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    .. figure:: figures/Decanomial.png
        :alt: Decanomial function
        :align: center

        **Two-dimensional Decanomial function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [2, -3]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)
        self.custom_bounds = [(0, 2.5), (-2, -4)]

        self.global_optimum = [2.0, -3.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        F1 = abs(x[0] ** 10 - 20 * x[0] ** 9 + 180 * x[0] ** 8 - 960 * x[0] ** 7 + 3360 * x[0] ** 6 - 8064 * x[0] ** 5 +
                 13340 * x[0] ** 4 - 15360 * x[0] ** 3 + 11520 * x[0] ** 2 - 5120 * x[0] + 2624.0)
        F2 = abs(x[1] ** 4 + 12 * x[1] ** 3 + 54 *
                 x[1] ** 2 + 108 * x[1] + 81.0)
        F = 0.001 * (F1 + F2) ** 2

        return F

# -------------------------------------------------------------------------------- #


class Deceptive(Benchmark):

    """
    Deceptive test objective function.

    This class defines the Deceptive global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Deceptive}}(\\mathbf{x}) = - \\left [\\frac{1}{n} \\sum_{i=1}^{n} g_i(x_i) \\right ]^{\\beta}


    Where :math:`\\beta` is a fixed non-linearity factor; in this exercise, :math:`\\beta = 2`. The function :math:`g_i(x_i)`
    is given by:

    .. math::

        g_i(x_i) = \\begin{cases} - \\frac{x}{\\alpha_i} + \\frac{4}{5} & \\textrm{if} \\hspace{5pt} 0 \\leq x_i \\leq \\frac{4}{5} \\alpha_i \\\\
           \\frac{5x}{\\alpha_i} -4 & \\textrm{if} \\hspace{5pt} \\frac{4}{5} \\alpha_i \\le x_i \\leq \\alpha_i \\\\
           \\frac{5(x - \\alpha_i)}{\\alpha_i-1} & \\textrm{if} \\hspace{5pt} \\alpha_i \\le x_i \\leq \\frac{1 + 4\\alpha_i}{5} \\\\
           \\frac{x - 1}{1 - \\alpha_i} & \\textrm{if} \\hspace{5pt} \\frac{1 + 4\\alpha_i}{5} \\le x_i \\leq 1 \\end{cases}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 1]` for :math:`i=1,...,n`.

    .. figure:: figures/Deceptive.png
        :alt: Deceptive function
        :align: center

        **Two-dimensional Deceptive function**

    *Global optimum*: :math:`f(x_i) = -1` for :math:`x_i = \\alpha_i` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([0.0] * self.dimensions, [1.0] * self.dimensions)

        n = self.dimensions
        alpha = arange(1.0, n + 1.0) / (n + 1.0)

        self.global_optimum = alpha
        self.fglob = -1.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        n = self.dimensions
        alpha = arange(1.0, n + 1.0) / (n + 1.0)
        beta = 2.0

        g = zeros((n, ))

        for i in range(n):
            if x[i] <= 0.0:
                g[i] = x[i]
            elif x[i] < 0.8 * alpha[i]:
                g[i] = -x[i] / alpha[i] + 0.8
            elif x[i] < alpha[i]:
                g[i] = 5.0 * x[i] / alpha[i] - 4.0
            elif x[i] < (1.0 + 4 * alpha[i]) / 5.0:
                g[i] = 5.0 * (x[i] - alpha[i]) / (alpha[i] - 1.0) + 1.0
            elif x[i] <= 1.0:
                g[i] = (x[i] - 1.0) / (1.0 - alpha[i]) + 4.0 / 5.0
            else:
                g[i] = x[i] - 1.0

        return -((1.0 / n) * sum(g)) ** beta

# -------------------------------------------------------------------------------- #


class DeckkersAarts(Benchmark):

    """
    Deckkers-Aarts test objective function.

    This class defines the Deckkers-Aarts global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{DeckkersAarts}}(\\mathbf{x}) = 10^5x_1^2 + x_2^2 - (x_1^2 + x_2^2)^2 + 10^{-5}(x_1^2 + x_2^2)^4


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-20, 20]` for :math:`i=1,2`.

    .. figure:: figures/DeckkersAarts.png
        :alt: DeckkersAarts function
        :align: center

        **Two-dimensional Deckkers-Aarts function**

    *Global optimum*: :math:`f(x_i) = -24776.518242168` for :math:`\\mathbf{x} = [0, \\pm 14.9451209]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-20.0] * self.dimensions, [20.0] * self.dimensions)
        self.custom_bounds = ([-1, 1], [14, 16])

        self.global_optimum = [0.0, 14.9451209]
        self.fglob = -24776.518342168

    def evaluator(self, x, *args):
        self.fun_evals += 1
        return (1.e5 * x[0] ** 2 + x[1] ** 2 - (x[0] ** 2 + x[1] ** 2) ** 2
                + 1.e-5 * (x[0] ** 2 + x[1] ** 2) ** 4)

# -------------------------------------------------------------------------------- #


class DeflectedCorrugatedSpring(Benchmark):

    """
    DeflectedCorrugatedSpring test objective function.

    This class defines the Deflected Corrugated Spring function global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{DeflectedCorrugatedSpring}}(\\mathbf{x}) = 0.1\\sum_{i=1}^n \\left[ (x_i - \\alpha)^2 - \\cos \\left( K \\sqrt {\\sum_{i=1}^n (x_i - \\alpha)^2} \\right ) \\right ]


    Where, in this exercise, :math:`K = 5` and :math:`\\alpha = 5`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 2\\alpha]` for :math:`i=1,...,n`.

    .. figure:: figures/DeflectedCorrugatedSpring.png
        :alt: Deflected Corrugated Spring function
        :align: center

        **Two-dimensional Deflected Corrugated Spring function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = \\alpha` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        alpha = 5.0
        self.bounds = zip([0] * self.dimensions, [2 * alpha] * self.dimensions)

        self.global_optimum = [alpha for _ in range(self.dimensions)]
        self.fglob = -1.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        K, alpha = 5.0, 5.0

        return -cos(K * sqrt(sum((x - alpha) ** 2))) + 0.1 * sum((x - alpha) ** 2)

# -------------------------------------------------------------------------------- #


class DeVilliersGlasser01(Benchmark):

    """
    DeVilliers-Glasser 1 test objective function.

    This class defines the DeVilliers-Glasser 1 function global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{DeVilliersGlasser01}}(\\mathbf{x}) = \\sum_{i=1}^{24} \\left[ x_1x_2^{t_i} \\sin(x_3t_i + x_4) - y_i \\right ]^2


    Where, in this exercise, :math:`t_i = 0.1(i-1)` and :math:`y_i = 60.137(1.371^{t_i}) \\sin(3.112t_i + 1.761)`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [1, 100]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`\\mathbf{x} = [60.137, 1.371, 3.112, 1.761]`.

    """

    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([1.0] * self.dimensions, [100.0] * self.dimensions)

        self.global_optimum = [60.137, 1.371, 3.112, 1.761]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        t_i = 0.1 * arange(24)
        y_i = 60.137 * (1.371 ** t_i) * sin(3.112 * t_i + 1.761)

        return sum((x[0] * (x[1] ** t_i) * sin(x[2] * t_i + x[3]) - y_i) ** 2.0)

# -------------------------------------------------------------------------------- #


class DeVilliersGlasser02(Benchmark):

    """
    DeVilliers-Glasser 2 test objective function.

    This class defines the DeVilliers-Glasser 2 function global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{DeVilliersGlasser01}}(\\mathbf{x}) = \\sum_{i=1}^{24} \\left[ x_1x_2^{t_i} \\tanh \\left [x_3t_i + \\sin(x_4t_i) \\right] \\cos(t_ie^{x_5}) - y_i \\right ]^2


    Where, in this exercise, :math:`t_i = 0.1(i-1)` and :math:`y_i = 53.81(1.27^{t_i}) \\tanh (3.012t_i + \\sin(2.13t_i)) \\cos(e^{0.507}t_i)`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [1, 60]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`\\mathbf{x} = [53.81, 1.27, 3.012, 2.13, 0.507]`.

    """

    def __init__(self, dimensions=5):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([1.0] * self.dimensions, [60.0] * self.dimensions)

        self.global_optimum = [53.81, 1.27, 3.012, 2.13, 0.507]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        t_i = 0.1 * arange(16)
        y_i = (53.81 * 1.27 ** t_i * tanh(3.012 * t_i + sin(2.13 * t_i))
               * cos(exp(0.507) * t_i))

        return sum((x[0] * (x[1] ** t_i) * tanh(x[2] * t_i + sin(x[3] * t_i))
                   * cos(t_i * exp(x[4])) - y_i) ** 2.0)

# -------------------------------------------------------------------------------- #


class DixonPrice(Benchmark):

    """
    Dixon and Price test objective function.

    This class defines the Dixon and Price global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{DixonPrice}}(\\mathbf{x}) = (x_i - 1)^2 + \\sum_{i=2}^n i(2x_i^2 - x_{i-1})^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    .. figure:: figures/DixonPrice.png
        :alt: Dixon and Price function
        :align: center

        **Two-dimensional Dixon and Price function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 2^{- \\frac{(2^i-2)}{2^i}}` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)
        self.custom_bounds = [(-2, 3), (-2, 3)]

        self.global_optimum = [2.0 ** (-(2.0 ** i - 2.0) / 2.0 ** i)
                               for i in range(1, self.dimensions + 1)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        s = 0.0
        for i in range(1, self.dimensions):
            s += i * (2.0 * x[i] ** 2.0 - x[i - 1]) ** 2.0

        y = s + (x[0] - 1.0) ** 2.0
        return y


# -------------------------------------------------------------------------------- #

class Dolan(Benchmark):

    """
    Dolan test objective function.

    This class defines the Dolan global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Dolan}}(\\mathbf{x}) = \\lvert (x_1 + 1.7x_2)\\sin(x_1) - 1.5x_3 - 0.1x_4\\cos(x_5 + x_5 - x_1) + 0.2x_5^2 - x_2 - 1 \\rvert

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 10^{-5}` for :math:`\\mathbf{x} = [8.39045925, 4.81424707, 7.34574133, 68.88246895, 3.85470806]`

    """

    def __init__(self, dimensions=5):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-100.0] * self.dimensions,
                          [100.0] * self.dimensions)

        self.global_optimum = [8.39045925, 4.81424707, 7.34574133, 68.88246895,
                               3.85470806]
        self.fglob = 1e-5

    def evaluator(self, x, *args):

        self.fun_evals += 1

        return abs((x[0] + 1.7 * x[1]) * sin(x[0]) - 1.5 * x[2] - 0.1 * x[3] * cos(x[3] + x[4] - x[0]) + 0.2 * x[4] ** 2.0 - x[1] - 1.0)

# -------------------------------------------------------------------------------- #


class DropWave(Benchmark):

    """
    DropWave test objective function.

    This class defines the DropWave global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{DropWave}}(\\mathbf{x}) = - \\frac{1 + \\cos\\left(12 \\sqrt{\\sum_{i=1}^{n} x_i^{2}}\\right)}{2 + 0.5 \\sum_{i=1}^{n} x_i^{2}}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5.12, 5.12]` for :math:`i=1,2`.

    .. figure:: figures/DropWave.png
        :alt: DropWave function
        :align: center

        **Two-dimensional DropWave function**


    *Global optimum*: :math:`f(x_i) = -1` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-5.12] * self.dimensions, [5.12] * self.dimensions)

        self.global_optimum = [0 for _ in range(self.dimensions)]
        self.fglob = -1.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        norm_x = sum(x ** 2)
        return -(1 + cos(12 * sqrt(norm_x))) / (0.5 * norm_x + 2)

# -------------------------------------------------------------------------------- #


class Easom(Benchmark):

    """
    Easom test objective function.

    This class defines the Easom global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Easom}}(\\mathbf{x}) = a - \\frac{a}{e^{b \\sqrt{\\frac{\\sum_{i=1}^{n} x_i^{2}}{n}}}} + e - e^{\\frac{\\sum_{i=1}^{n} \\cos\\left(c x_i\\right)}{n}}

    Where, in this exercise, :math:`a = 20, b = 0.2` and :math:`c = 2\\pi`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    .. figure:: figures/Easom.png
        :alt: Easom function
        :align: center

        **Two-dimensional Easom function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-100.0] * self.dimensions,
                          [100.0] * self.dimensions)

        self.global_optimum = [0 for _ in range(self.dimensions)]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        a = 20.0
        b = 0.2
        c = 2 * pi
        n = self.dimensions

        return -a * exp(-b * sqrt(sum(x ** 2) / n)) - exp(sum(cos(c * x)) / n) + a + exp(1)

# -------------------------------------------------------------------------------- #


class EggCrate(Benchmark):

    """
    Egg Crate test objective function.

    This class defines the Egg Crate global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{EggCrate}}(\\mathbf{x}) = x_1^2 + x_2^2 + 25 \\left[ \\sin^2(x_1) + \\sin^2(x_2) \\right]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,2`.

    .. figure:: figures/EggCrate.png
        :alt: Egg Crate function
        :align: center

        **Two-dimensional Egg Crate function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-5.0] * self.dimensions, [5.0] * self.dimensions)

        self.global_optimum = [0.0, 0.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2 = x
        return x1 ** 2.0 + x2 ** 2.0 + 25.0 * (sin(x1) ** 2.0 + sin(x2) ** 2.0)

# -------------------------------------------------------------------------------- #


class EggHolder(Benchmark):

    """
    Egg Holder test objective function.

    This class defines the Egg Holder global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{EggHolder}}(\\mathbf{x}) = - x_{1} \\sin\\left(\\sqrt{\\lvert{x_{1} - x_{2} -47}\\rvert}\\right) - \\left(x_{2} + 47\\right) \\sin\\left(\\sqrt{\\left|{\\frac{1}{2} x_{1} + x_{2} + 47}\\right|}\\right)

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-512, 512]` for :math:`i=1,2`.

    .. figure:: figures/EggHolder.png
        :alt: Egg Holder function
        :align: center

        **Two-dimensional Egg Holder function**


    *Global optimum*: :math:`f(x_i) = -959.640662711` for :math:`\\mathbf{x} = [512, 404.2319]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-512.1] * self.dimensions,
                          [512.0] * self.dimensions)

        self.global_optimum = [512.0, 404.2319]
        self.fglob = -959.640662711

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return -(x[1] + 47) * sin(sqrt(abs(x[1] + x[0] / 2 + 47))) - x[0] * sin(sqrt(abs(x[0] - (x[1] + 47))))

# -------------------------------------------------------------------------------- #


class ElAttarVidyasagarDutta(Benchmark):

    """
    El-Attar-Vidyasagar-Dutta test objective function.

    This class defines the El-Attar-Vidyasagar-Dutta function global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{ElAttarVidyasagarDutta}}(\\mathbf{x}) = (x_1^2 + x_2 - 10)^2 + (x_1 + x_2^2 - 7)^2 + (x_1^2 + x_2^3 - 1)^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    .. figure:: figures/ElAttarVidyasagarDutta.png
        :alt: El-Attar-Vidyasagar-Dutta function
        :align: center

        **Two-dimensional El-Attar-Vidyasagar-Dutta function**


    *Global optimum*: :math:`f(x_i) = 1.712780354` for :math:`\\mathbf{x} = [3.40918683, -2.17143304]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-100.0] * self.dimensions,
                          [100.0] * self.dimensions)
        self.custom_bounds = [(-4, 4), (-4, 4)]

        self.global_optimum = [3.40918683, -2.17143304]
        self.fglob = 1.712780354

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2 = x
        return (x1 ** 2.0 + x2 - 10) ** 2.0 + (x1 + x2 ** 2.0 - 7) ** 2.0 + (x1 ** 2.0 + x2 ** 3.0 - 1) ** 2.0

# -------------------------------------------------------------------------------- #


class Exp2(Benchmark):

    """
    Exp2 test objective function.

    This class defines the Exp2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Exp2}}(\\mathbf{x}) = \\sum_{i=0}^9 \\left ( e^{-ix_1/10} - 5e^{-ix_2/10} -e^{-i/10} + 5e^{-i} \\right )^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 20]` for :math:`i=1,2`.

    .. figure:: figures/Exp2.png
        :alt: Exp2 function
        :align: center

        **Two-dimensional Exp2 function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = [1, 0.1]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([0.0] * self.dimensions, [20.0] * self.dimensions)
        self.custom_bounds = [(0, 2), (0, 2)]

        self.global_optimum = [1.0, 0.1]
        self.fglob = 0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        y = 0.0
        for i in range(10):
            y += (exp(-i * x[0] / 10.0) - 5 * exp(-i * x[1] * 10)
                  - exp(-i / 10.0) + 5 * exp(-i)) ** 2.0

        return y

# -------------------------------------------------------------------------------- #


class Exponential(Benchmark):

    """
    Exponential test objective function.

    This class defines the Exponential global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Exponential}}(\\mathbf{x}) = -e^{-0.5 \\sum_{i=1}^n x_i^2}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,...,n`.

    .. figure:: figures/Exponential.png
        :alt: Exponential function
        :align: center

        **Two-dimensional Exponential function**


    *Global optimum*: :math:`f(x_i) = -1` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-1.0] * self.dimensions, [1.0] * self.dimensions)

        self.global_optimum = [0.0 for _ in range(self.dimensions)]
        self.fglob = -1.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return -exp(-0.5 * sum(x ** 2.0))

# -------------------------------------------------------------------------------- #


class FreudensteinRoth(Benchmark):

    """
    FreudensteinRoth test objective function.

    This class defines the Freudenstein & Roth global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{FreudensteinRoth}}(\\mathbf{x}) =  \\left\{x_1 - 13 + \\left[(5 - x_2)x_2 - 2 \\right] x_2 \\right\}^2 + \\left \{x_1 - 29 + \\left[(x_2 + 1)x_2 - 14 \\right] x_2 \\right\}^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    .. figure:: figures/FreudensteinRoth.png
        :alt: FreudensteinRoth function
        :align: center

        **Two-dimensional FreudensteinRoth function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [5, 4]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)
        self.custom_bounds = [(-3, 3), (-5, 5)]

        self.global_optimum = [5.0, 4.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1
        f1 = (-13.0 + x[0] + ((5.0 - x[1]) * x[1] - 2.0) * x[1]) ** 2
        f2 = (-29.0 + x[0] + ((x[1] + 1.0) * x[1] - 14.0) * x[1]) ** 2

        return f1 + f2

# -------------------------------------------------------------------------------- #


class Gear(Benchmark):

    """
    Gear test objective function.

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

        self.bounds = zip([12.0] * self.dimensions, [60.0] * self.dimensions)
        self.global_optimum = [16, 19, 43, 49]
        self.fglob = 2.7e-12

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return (1.0 / 6.931 - floor(x[0]) * floor(x[1]) / (floor(x[2]) * floor(x[3]))) ** 2

# -------------------------------------------------------------------------------- #


class Giunta(Benchmark):

    """
    Giunta test objective function.

    This class defines the Giunta global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Giunta}}(\\mathbf{x}) = 0.6 + \\sum_{i=1}^{n} \\left[\\sin^{2}\\left(1 - \\frac{16}{15} x_i\\right) - \\frac{1}{50} \\sin\\left(4 - \\frac{64}{15} x_i\\right) - \\sin\\left(1 - \\frac{16}{15} x_i\\right)\\right]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,2`.

    .. figure:: figures/Giunta.png
        :alt: Giunta function
        :align: center

        **Two-dimensional Giunta function**


    *Global optimum*: :math:`f(x_i) = 0.06447042053690566` for :math:`\\mathbf{x} = [0.4673200277395354, 0.4673200169591304]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-1.0] * self.dimensions, [1.0] * self.dimensions)

        self.global_optimum = [0.4673200277395354, 0.4673200169591304]
        self.fglob = 0.06447042053690566

    def evaluator(self, x, *args):

        self.fun_evals += 1

        arg = 16 * x / 15.0 - 1
        return 0.6 + sum(sin(arg) + sin(arg) ** 2 + sin(4 * arg) / 50)

# -------------------------------------------------------------------------------- #


class GoldsteinPrice(Benchmark):

    """
    Goldstein-Price test objective function.

    This class defines the Goldstein-Price global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{GoldsteinPrice}}(\\mathbf{x}) = \\left[ 1+(x_1+x_2+1)^2(19-14x_1+3x_1^2-14x_2+6x_1x_2+3x_2^2) \\right] \\left[ 30+(2x_1-3x_2)^2(18-32x_1+12x_1^2+48x_2-36x_1x_2+27x_2^2) \\right]


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-2, 2]` for :math:`i=1,2`.

    .. figure:: figures/GoldsteinPrice.png
        :alt: Goldstein-Price function
        :align: center

        **Two-dimensional Goldstein-Price function**


    *Global optimum*: :math:`f(x_i) = 3` for :math:`\\mathbf{x} = [0, -1]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-2.0] * self.dimensions, [2.0] * self.dimensions)

        self.global_optimum = [0., -1.]
        self.fglob = 3.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        a = 1 + (x[0] + x[1] + 1) ** 2 * \
            (19 - 14 * x[0] + 3 * x[0] ** 2 -
             14 * x[1] + 6 * x[0] * x[1] + 3 * x[1] ** 2)
        b = 30 + (2 * x[0] - 3 * x[1]) ** 2 * \
            (18 - 32 * x[0] + 12 * x[0] ** 2 +
             48 * x[1] - 36 * x[0] * x[1] + 27 * x[1] ** 2)
        return a * b

# -------------------------------------------------------------------------------- #


class Griewank(Benchmark):

    """
    Griewank test objective function.

    This class defines the Griewank global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Griewank}}(\\mathbf{x}) = \\frac{1}{4000}\\sum_{i=1}^n x_i^2 - \\prod_{i=1}^n\\cos\\left(\\frac{x_i}{\\sqrt{i}}\\right) + 1

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-600, 600]` for :math:`i=1,...,n`.

    .. figure:: figures/Griewank.png
        :alt: Griewank function
        :align: center

        **Two-dimensional Griewank function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-600.0] * self.dimensions,
                          [600.0] * self.dimensions)
        self.custom_bounds = [(-50, 50), (-50, 50)]

        self.global_optimum = [0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return sum(x ** 2) / 4000.0 - prod(cos(x / sqrt(1.0 + arange(len(x))))) + 1.0

# -------------------------------------------------------------------------------- #


class Gulf(Benchmark):

    """
    Gulf test objective function.

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

        self.bounds = zip([0.0] * self.dimensions, [50.0] * self.dimensions)

        self.global_optimum = [50.0, 25.0, 1.5]
        self.fglob = 0.0

    def evaluator(self, x, *args):
        self.fun_evals += 1

        m = 99.
        i = arange(1., m + 1)
        y = 25 + (-50 * log(i / 100.)) ** (2 / 3.)
        vec = (exp(-((abs(y - x[1])) ** x[2] / x[0])) - i / 100.)
        return sum(vec ** 2)

# -------------------------------------------------------------------------------- #


class Hansen(Benchmark):

    """
    Hansen test objective function.

    This class defines the Hansen global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Hansen}}(\\mathbf{x}) = \\left[ \\sum_{i=0}^4(i+1)\\cos(ix_1+i+1)\\right ] \\left[\\sum_{j=0}^4(j+1)\\cos[(j+2)x_2+j+1])\\right ]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    .. figure:: figures/Hansen.png
        :alt: Hansen function
        :align: center

        **Two-dimensional Hansen function**


    *Global optimum*: :math:`f(x_i) = -176.54179` for :math:`\\mathbf{x} = [-7.58989583, -7.70831466]`.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [-7.58989583, -7.70831466]
        self.fglob = -176.54179

    def evaluator(self, x, *args):

        self.fun_evals += 1

        f1 = f2 = 0.0
        for i in range(5):
            f1 += (i + 1) * cos(i * x[0] + i + 1)
            f2 += (i + 1) * cos((i + 2) * x[1] + i + 1)

        return f1 * f2

# -------------------------------------------------------------------------------- #


class Hartmann3(Benchmark):

    """
    Hartmann3 test objective function.

    This class defines the Hartmann3 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Hartmann3}}(\\mathbf{x}) = -\\sum\\limits_{i=1}^{4} c_i e^{-\\sum\\limits_{j=1}^{n}a_{ij}(x_j - p_{ij})^2}

    Where, in this exercise:

    .. math::

        \\begin{array}{l|ccc|c|ccr}
        \\hline
        i & & a_{ij}&  & c_i & & p_{ij} &  \\\\
        \\hline
        1 & 3.0 & 10.0 & 30.0 & 1.0 & 0.3689  & 0.1170 & 0.2673 \\\\
        2 & 0.1 & 10.0 & 35.0 & 1.2 & 0.4699 & 0.4387 & 0.7470 \\\\
        3 & 3.0 & 10.0 & 30.0 & 3.0 & 0.1091 & 0.8732 & 0.5547 \\\\
        4 & 0.1 & 10.0 & 35.0 & 3.2 & 0.03815 & 0.5743 & 0.8828 \\\\
        \\hline
        \\end{array}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 1]` for :math:`i=1,2,3`.

    *Global optimum*: :math:`f(x_i) = -3.8627821478` for :math:`\\mathbf{x} = [0.11461292,  0.55564907,  0.85254697]`

    """

    def __init__(self, dimensions=3):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([0.0] * self.dimensions, [1.0] * self.dimensions)

        self.global_optimum = [0.11461292,  0.55564907,  0.85254697]
        self.fglob = -3.8627821478

    def evaluator(self, x, *args):

        self.fun_evals += 1

        a = asarray([[3.0,  0.1,  3.0,  0.1],
                     [10.0, 10.0, 10.0, 10.0],
                     [30.0, 35.0, 30.0, 35.0]])
        p = asarray([[0.36890, 0.46990, 0.10910, 0.03815],
                     [0.11700, 0.43870, 0.87320, 0.57430],
                     [0.26730, 0.74700, 0.55470, 0.88280]])
        c = asarray([1.0, 1.2, 3.0, 3.2])
        d = zeros_like(c)

        for i in range(4):
            d[i] = sum(a[:, i] * (x - p[:, i]) ** 2)

        return -sum(c * exp(-d))

# -------------------------------------------------------------------------------- #


class Hartmann6(Benchmark):

    """
    Hartmann6 test objective function.

    This class defines the Hartmann6 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Hartmann6}}(\\mathbf{x}) = -\\sum\\limits_{i=1}^{4} c_i e^{-\\sum\\limits_{j=1}^{n}a_{ij}(x_j - p_{ij})^2}

    Where, in this exercise:

    .. math::

        \\begin{array}{l|cccccc|r}
        \\hline
        i & &   &   a_{ij} &  &  & & c_i  \\\\
        \\hline
        1 & 10.0  & 3.0  & 17.0 & 3.50  & 1.70  & 8.00  & 1.0 \\\\
        2 & 0.05  & 10.0 & 17.0 & 0.10  & 8.00  & 14.00 & 1.2 \\\\
        3 & 3.00  & 3.50 & 1.70 & 10.0  & 17.00 & 8.00  & 3.0 \\\\
        4 & 17.00 & 8.00 & 0.05 & 10.00 & 0.10  & 14.00 & 3.2 \\\\
        \\hline
        \\end{array}

        \\newline
        \\\\
        \\newline

        \\begin{array}{l|cccccr}
        \\hline
        i &  &   & p_{ij} &  & & \\\\
        \\hline
        1 & 0.1312 & 0.1696 & 0.5569 & 0.0124 & 0.8283 & 0.5886 \\\\
        2 & 0.2329 & 0.4135 & 0.8307 & 0.3736 & 0.1004 & 0.9991 \\\\
        3 & 0.2348 & 0.1451 & 0.3522 & 0.2883 & 0.3047 & 0.6650 \\\\
        4 & 0.4047 & 0.8828 & 0.8732 & 0.5743 & 0.1091 & 0.0381 \\\\
        \\hline
        \\end{array}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 1]` for :math:`i=1,...,6`.

    *Global optimum*: :math:`f(x_i) = -3.32236801141551` for :math:`\\mathbf{x} = [0.20168952, 0.15001069, 0.47687398, 0.27533243, 0.31165162, 0.65730054]`

    """

    def __init__(self, dimensions=6):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([0.0] * self.dimensions, [1.0] * self.dimensions)

        self.global_optimum = [0.20168952, 0.15001069,
                               0.47687398, 0.27533243, 0.31165162, 0.65730054]
        self.fglob = -3.32236801141551

    def evaluator(self, x, *args):

        self.fun_evals += 1

        a = asarray([[10.00,  0.05,  3.00, 17.00],
                     [3.00, 10.00,  3.50,  8.00],
                     [17.00, 17.00,  1.70,  0.05],
                     [3.50,  0.10, 10.00, 10.00],
                     [1.70,  8.00, 17.00,  0.10],
                     [8.00, 14.00,  8.00, 14.00]])

        p = asarray([[0.1312, 0.2329, 0.2348, 0.4047],
                     [0.1696, 0.4135, 0.1451, 0.8828],
                     [0.5569, 0.8307, 0.3522, 0.8732],
                     [0.0124, 0.3736, 0.2883, 0.5743],
                     [0.8283, 0.1004, 0.3047, 0.1091],
                     [0.5886, 0.9991, 0.6650, 0.0381]])

        c = asarray([1.0, 1.2, 3.0, 3.2])
        d = zeros_like(c)

        for i in range(4):
            d[i] = sum(a[:, i] * (x - p[:, i]) ** 2)

        return -sum(c * exp(-d))

# -------------------------------------------------------------------------------- #


class HelicalValley(Benchmark):

    """
    HelicalValley test objective function.

    This class defines the HelicalValley global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{HelicalValley}}(\\mathbf{x}) = 100{[z-10\\Psi(x_1,x_2)]^2+(\\sqrt{x_1^2+x_2^2}-1)^2}+x_3^2

    Where, in this exercise:

    .. math::

        2\\pi\\Psi(x,y) =  \\begin{cases} \\arctan(y/x) & \\textrm{for} x > 0 \\\\
        \\pi + \\arctan(y/x) & \\textrm{for} x < 0 \\end{cases}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-\infty, \\infty]` for :math:`i=1,2,3`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [1, 0, 0]`

    """

    def __init__(self, dimensions=3):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-100.0] * self.dimensions, [100] * self.dimensions)

        self.global_optimum = [1.0, 0.0, 0.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return 100 * ((x[2] - 10 * arctan2(x[1], x[0]) / 2 / pi) ** 2 + (sqrt(x[0] ** 2 + x[1] ** 2) - 1) ** 2) + x[2] ** 2

# -------------------------------------------------------------------------------- #


class HimmelBlau(Benchmark):

    """
    HimmelBlau test objective function.

    This class defines the HimmelBlau global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{HimmelBlau}}(\\mathbf{x}) = (x_1^2 + x_2 - 11)^2 + (x_1 + x_2^2 -7)^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-6, 6]` for :math:`i=1,2`.

    .. figure:: figures/HimmelBlau.png
        :alt: HimmelBlau function
        :align: center

        **Two-dimensional HimmelBlau function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [3, 2]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-6] * self.dimensions, [6] * self.dimensions)

        self.global_optimum = [3.0, 2.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return (x[0] ** 2 + x[1] - 11) ** 2 + (x[0] + x[1] ** 2 - 7) ** 2

# -------------------------------------------------------------------------------- #


class HolderTable(Benchmark):

    """
    HolderTable test objective function.

    This class defines the HolderTable global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{HolderTable}}(\\mathbf{x}) = - \\left|{e^{\\left|{1 - \\frac{\\sqrt{x_{1}^{2} + x_{2}^{2}}}{\\pi} }\\right|} \\sin\\left(x_{1}\\right) \\cos\\left(x_{2}\\right)}\\right|

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    .. figure:: figures/HolderTable.png
        :alt: HolderTable function
        :align: center

        **Two-dimensional HolderTable function**


    *Global optimum*: :math:`f(x_i) = -19.20850256788675` for :math:`x_i = \\pm 9.664590028909654` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [(8.055023472141116, 9.664590028909654),
                               (-8.055023472141116, 9.664590028909654),
                               (8.055023472141116, -9.664590028909654),
                               (-8.055023472141116, -9.664590028909654)]
        self.fglob = -19.20850256788675

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return -abs(sin(x[0]) * cos(x[1]) * exp(abs(1 - sqrt(x[0] ** 2 + x[1] ** 2) / pi)))

# -------------------------------------------------------------------------------- #


class Holzman(Benchmark):

    """
    Holzman test objective function.

    This class defines the Holzman global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Holzman}}(\\mathbf{x}) = \\sum_{i=1}^{100} \\left [ e^{\\frac{1}{x_1} (u_i-x_2)^{x_3}} -0.01(i) \\right ]

    Where, in this exercise:

    .. math::

        u_i = 25 + (-50 \\log{[0.01i]})^{2/3}


    Here, :math:`n` represents the number of dimensions and :math:`x_1 \\in [0, 100], x_2 \\in [0, 25.6], x_3 \\in [0, 5]`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [50, 25, 1.5]`

    """

    def __init__(self, dimensions=3):
        Benchmark.__init__(self, dimensions)

        self.bounds = ([0.0, 100.0], [0.0, 25.6], [0.0, 5.0])

        self.global_optimum = [50.0, 25.0, 1.5]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        val = 0
        i = arange(1, 101)
        t = 2 / 3.
        u = 25 + (-50 * log(0.01 * i)) ** t
        v = (u - x[1]) ** x[2]
        w = exp(-v / x[0])
        return sum(-0.01 * i + w)

# -------------------------------------------------------------------------------- #


class Hosaki(Benchmark):

    """
    Hosaki test objective function.

    This class defines the Hosaki global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Hosaki}}(\\mathbf{x}) = \\left ( 1 - 8x_1 + 7x_1^2 - \\frac{7}{3}x_1^3 + \\frac{1}{4}x_1^4 \\right )x_2^2e^{-x_1}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 10]` for :math:`i=1,2`.

    .. figure:: figures/Hosaki.png
        :alt: Hosaki function
        :align: center

        **Two-dimensional Hosaki function**


    *Global optimum*: :math:`f(x_i) = -2.3458` for :math:`\\mathbf{x} = [4, 2]`.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([0.0] * self.dimensions, [10.0] * self.dimensions)
        self.custom_bounds = [(0, 5), (0, 5)]

        self.global_optimum = [4, 2]
        self.fglob = -2.3458

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return (1 + x[0] * (-8 + x[0] * (7 + x[0] * (-7.0 / 3.0 + x[0] * 1.0 / 4.0)))) * x[1] * x[1] * exp(-x[1])

# -------------------------------------------------------------------------------- #


class Infinity(Benchmark):

    """
    Infinity test objective function.

    This class defines the Infinity global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Infinity}}(\\mathbf{x}) = \\sum_{i=1}^{n} x_i^{6} \\left [ \\sin\\left ( \\frac{1}{x_i} \\right )+2 \\right ]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,...,n`.

    .. figure:: figures/Infinity.png
        :alt: Infinity function
        :align: center

        **Two-dimensional Infinity function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-1.0] * self.dimensions, [1.0] * self.dimensions)

        self.global_optimum = [1e-16 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return sum(x ** 6.0 * (sin(1.0 / x) + 2.0))

# -------------------------------------------------------------------------------- #


class JennrichSampson(Benchmark):

    """
    Jennrich-Sampson test objective function.

    This class defines the Jennrich-Sampson global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{JennrichSampson}}(\\mathbf{x}) = \\sum_{i=1}^{10} \\left [2 + 2i - (e^{ix_1} + e^{ix_2}) \\right ]^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,2`.

    .. figure:: figures/JennrichSampson.png
        :alt: Jennrich-Sampson function
        :align: center

        **Two-dimensional Jennrich-Sampson function**


    *Global optimum*: :math:`f(x_i) = 124.3621824` for :math:`\\mathbf{x} = [0.257825, 0.257825]`.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-1.0] * self.dimensions, [1.0] * self.dimensions)

        self.global_optimum = [0.257825, 0.257825]
        self.custom_bounds = [(-1, 0.34), (-1, 0.34)]
        self.fglob = 124.3621824

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2 = x
        rng = arange(1.0, 11.0)
        return sum((2.0 + 2.0 * rng - (exp(rng * x1) + exp(rng * x2))) ** 2.0)

# -------------------------------------------------------------------------------- #


class Judge(Benchmark):

    """
    Judge test objective function.

    This class defines the Judge global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Judge}}(\\mathbf{x}) = \\sum_{i=1}^{20} \\left [ \\left (x_1 + A_i x_2 + B x_2^2 \\right ) - C_i \\right ]^2

    Where, in this exercise:

    .. math::

        \\begin{cases} A = [4.284, 4.149, 3.877, 0.533, 2.211, 2.389, 2.145,  3.231, 1.998, 1.379, 2.106, 1.428, 1.011, 2.179, 2.858, 1.388, 1.651, 1.593, 1.046, 2.152] \\\\
        B = [0.286, 0.973, 0.384, 0.276, 0.973, 0.543, 0.957, 0.948, 0.543, 0.797, 0.936, 0.889, 0.006, 0.828, 0.399, 0.617, 0.939, 0.784, 0.072, 0.889] \\\\
        C = [0.645, 0.585, 0.310, 0.058, 0.455, 0.779, 0.259, 0.202, 0.028, 0.099, 0.142, 0.296, 0.175, 0.180, 0.842, 0.039, 0.103, 0.620, 0.158, 0.704] \\end{cases}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    .. figure:: figures/Judge.png
        :alt: Judge function
        :align: center

        **Two-dimensional Judge function**


    *Global optimum*: :math:`f(x_i) = 16.0817307` for :math:`\\mathbf{x} = [0.86479, 1.2357]`.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [0.86479, 1.2357]
        self.custom_bounds = [(-2.0, 2.0), (-2.0, 2.0)]
        self.fglob = 16.0817307

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2 = x

        Y = asarray(
            [4.284, 4.149, 3.877, 0.533, 2.211, 2.389, 2.145,  3.231, 1.998, 1.379,
             2.106, 1.428, 1.011, 2.179, 2.858, 1.388, 1.651, 1.593, 1.046, 2.152])

        X2 = asarray(
            [0.286, 0.973, 0.384, 0.276, 0.973, 0.543, 0.957, 0.948, 0.543, 0.797,
             0.936, 0.889, 0.006, 0.828, 0.399, 0.617, 0.939, 0.784, 0.072, 0.889])

        X3 = asarray(
            [0.645, 0.585, 0.310, 0.058, 0.455, 0.779, 0.259, 0.202, 0.028, 0.099,
             0.142, 0.296, 0.175, 0.180, 0.842, 0.039, 0.103, 0.620, 0.158, 0.704])

        return sum(((x1 + x2 * X2 + (x2 ** 2.0) * X3) - Y) ** 2.0)

# -------------------------------------------------------------------------------- #


class Katsuura(Benchmark):

    """
    Katsuura test objective function.

    This class defines the Katsuura global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Katsuura}}(\\mathbf{x}) = \\prod_{i=0}^{n-1} \\left [ 1 + (i+1) \\sum_{k=1}^{d} \\lfloor (2^k x_i) \\rfloor 2^{-k} \\right ]

    Where, in this exercise, :math:`d = 32`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 100]` for :math:`i=1,...,n`.

    .. figure:: figures/Katsuura.png
        :alt: Katsuura function
        :align: center

        **Two-dimensional Katsuura function**


    *Global optimum*: :math:`f(x_i) = 1` for :math:`x_i = 0` for :math:`i=1,...,n`.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([0.0] * self.dimensions, [100.0] * self.dimensions)

        self.global_optimum = [0.0 for _ in range(self.dimensions)]
        self.custom_bounds = [(0, 1), (0, 1)]
        self.fglob = 1.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        d = 32

        prod = 1.0
        for i in range(self.dimensions):
            s = 0.0
            for k in range(1, d + 1):
                pow2 = 2.0 ** k
                s += round(pow2 * x[i]) / pow2
            prod = prod * (1.0 + (i + 1.0) * s)

        return prod

# -------------------------------------------------------------------------------- #


class Keane(Benchmark):

    """
    Keane test objective function.

    This class defines the Keane global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Keane}}(\\mathbf{x}) = \\frac{\\sin^2(x_1 - x_2)\\sin^2(x_1 + x_2)}{\\sqrt{x_1^2 + x_2^2}}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 10]` for :math:`i=1,2`.

    .. figure:: figures/Keane.png
        :alt: Keane function
        :align: center

        **Two-dimensional Keane function**


    *Global optimum*: :math:`f(x_i) = 0.673668` for :math:`\\mathbf{x} = [0.0, 1.39325]`.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([0.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [0.0, 1.39325]
        self.custom_bounds = [(-1, 0.34), (-1, 0.34)]
        self.fglob = 0.673668

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2 = x

        return (sin(x1 - x2) ** 2.0 * sin(x1 + x2) ** 2.0) / sqrt(x1 ** 2.0 + x2 ** 2.0)

# -------------------------------------------------------------------------------- #


class Kowalik(Benchmark):

    """
    Kowalik test objective function.

    This class defines the Kowalik global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Kowalik}}(\\mathbf{x}) = \\sum_{i=0}^{10} \\left [ a_i - \\frac{x_1(b_i^2+b_ix_2)}{b_i^2 + b_ix_3 + x_4} \\right ]^2

    Where:

    .. math::

       \\mathbf{a} = [4, 2, 1, 1/2, 1/4 1/8, 1/10, 1/12, 1/14, 1/16] \\\\
       \\mathbf{b} = [0.1957, 0.1947, 0.1735, 0.1600, 0.0844, 0.0627, 0.0456, 0.0342, 0.0323, 0.0235, 0.0246]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,...,4`.

    *Global optimum*: :math:`f(x_i) = 0.00030748610` for :math:`\\mathbf{x} = [0.192833, 0.190836, 0.123117, 0.135766]`.

    """

    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-5.0] * self.dimensions, [5.0] * self.dimensions)
        self.global_optimum = [0.192833, 0.190836, 0.123117, 0.135766]
        self.fglob = 0.00030748610

    def evaluator(self, x, *args):

        self.fun_evals += 1

        b = asarray([4.0,   2.0,   1.0,    1 / 2.0,  1 / 4.0,
                     1 / 6.0, 1 / 8.0, 1 / 10.0, 1 / 12.0, 1 / 14.0,
                     1 / 16.0])
        a = asarray([0.1957, 0.1947, 0.1735, 0.1600, 0.0844,
                     0.0627, 0.0456, 0.0342, 0.0323, 0.0235,
                     0.0246])

        y = 0.0
        for i in range(11):
            bb = b[i] * b[i]
            t = a[i] - (x[0] * (bb + b[i] * x[1]) / (bb + b[i] * x[2] + x[3]))
            y += t * t

        return y


# -------------------------------------------------------------------------------- #

class Langermann(Benchmark):

    """
    Langermann test objective function.

    This class defines the Langermann global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Langermann}}(\\mathbf{x}) = - \\sum_{i=1}^{5} \\frac{c_i \\cos\\left\{\\pi \\left[\\left(x_{1}- a_i\\right)^{2} + \\left(x_{2} - b_i \\right)^{2}\\right]\\right\}}{e^{\\frac{\\left( x_{1} - a_i\\right)^{2} + \\left( x_{2} - b_i\\right)^{2}}{\\pi}}}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 10]` for :math:`i=1,2`.

    .. figure:: figures/Langermann.png
        :alt: Langermann function
        :align: center

        **Two-dimensional Langermann function**


    *Global optimum*: :math:`f(x_i) = -5.1621259` for :math:`\\mathbf{x} = [2.00299219, 1.006096]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([0.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [2.00299219, 1.006096]
        self.fglob = -5.1621259

    def evaluator(self, x, *args):

        self.fun_evals += 1
        a = [3, 5, 2, 1, 7]
        b = [5, 2, 1, 4, 9]
        c = [1, 2, 5, 2, 3]

        return -sum(c * exp(-(1 / pi) * ((x[0] - a) ** 2 +
                    (x[1] - b) ** 2)) * cos(pi * ((x[0] - a) ** 2 + (x[1] - b) ** 2)))

# -------------------------------------------------------------------------------- #


class LennardJones(Benchmark):

    """
    LennardJones test objective function.

    This class defines the Lennard-Jones global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{LennardJones}}(\\mathbf{x}) = \\sum_{i=0}^{n-2}\\sum_{j>1}^{n-1}\\frac{1}{r_{ij}^{12}} - \\frac{1}{r_{ij}^{6}}


    Where, in this exercise:

    .. math::

        r_{ij} = \\sqrt{(x_{3i}-x_{3j})^2 + (x_{3i+1}-x_{3j+1})^2) + (x_{3i+2}-x_{3j+2})^2}


    Valid for any dimension, :math:`n = 3*k, k=2,3,4,...,20`. :math:`k` is the number of atoms in 3-D space
    constraints: unconstrained type: multi-modal with one global minimum; non-separable

    Value-to-reach: :math:`minima[k-2] + 0.0001`. See array of minima below; additional minima available at
    the Cambridge cluster database:

    http://www-wales.ch.cam.ac.uk/~jon/structures/LJ/tables.150.html

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-4, 4]` for :math:`i=1,...,n`.

    *Global optimum*:

    .. math::

       minima = [-1.,-3.,-6.,-9.103852,-12.712062,-16.505384,-19.821489,-24.113360, \\\\
       -28.422532,-32.765970,-37.967600,-44.326801,-47.845157,-52.322627, \\\\
       -56.815742,-61.317995, -66.530949,-72.659782,-77.1777043]

    """

    def __init__(self, dimensions=6):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-4.0] * self.dimensions, [4.0] * self.dimensions)

        self.global_optimum = []

        minima = [
            -1.0, -3.0, -6.0, -9.103852, -12.712062, -
            16.505384, -19.821489, -24.113360, -28.422532,
            -32.765970, -37.967600, -44.326801, -
            47.845157, -52.322627, -56.815742, -61.317995,
            -66.530949, -72.659782, -77.1777043]

        k = dimensions / 3
        self.fglob = minima[k - 2]
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        k = self.dimensions / 3
        s = 0.0

        for i in range(k - 1):
            for j in range(i + 1, k):
                a = 3 * i
                b = 3 * j
                xd = x[a] - x[b]
                yd = x[a + 1] - x[b + 1]
                zd = x[a + 2] - x[b + 2]
                ed = xd * xd + yd * yd + zd * zd
                ud = ed * ed * ed
                if ed > 0.0:
                    s += (1.0 / ud - 2.0) / ud

        return s

# -------------------------------------------------------------------------------- #


class Leon(Benchmark):

    """
    Leon test objective function.

    This class defines the Leon global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Leon}}(\\mathbf{x}) = \\left(1 - x_{1}\\right)^{2} + 100 \\left(x_{2} - x_{1}^{2} \\right)^{2}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1.2, 1.2]` for :math:`i=1,2`.

    .. figure:: figures/Leon.png
        :alt: Leon function
        :align: center

        **Two-dimensional Leon function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 1` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-1.2] * self.dimensions, [1.2] * self.dimensions)

        self.global_optimum = [1 for _ in range(self.dimensions)]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return 100 * (x[1] - x[0] ** 2.0) ** 2.0 + (1 - x[0]) ** 2.0

# -------------------------------------------------------------------------------- #


class Levy03(Benchmark):

    """
    Levy 3 test objective function.

    This class defines the Levy 3 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Levy03}}(\\mathbf{x}) = \\sin^2(\\pi y_1)+\\sum_{i=1}^{n-1}(y_i-1)^2[1+10\\sin^2(\\pi y_{i+1})]+(y_n-1)^2

    Where, in this exercise:

    .. math::

        y_i=1+\\frac{x_i-1}{4}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    .. figure:: figures/Levy03.png
        :alt: Levy 3 function
        :align: center

        **Two-dimensional Levy 3 function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 1` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)
        self.custom_bounds = [(-5, 5), (-5, 5)]

        self.global_optimum = [1 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        n = len(x)
        z = zeros_like(x)

        for i in range(n):
            z[i] = 1 + (x[i] - 1) / 4

        s = sin(pi * z[0]) ** 2
        for i in range(n - 1):
            s = s + (z[i] - 1) ** 2 * (1 + 10 * (sin(pi * z[i] + 1)) ** 2)

        y = s + (z[n - 1] - 1) ** 2 * (1 + (sin(2 * pi * z[n - 1])) ** 2)
        return y

# -------------------------------------------------------------------------------- #


class Levy05(Benchmark):

    """
    Levy 5 test objective function.

    This class defines the Levy 5 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Levy05}}(\\mathbf{x}) = \\sum_{i=1}^{5} i \\cos \\left[(i-1)x_1 + i \\right] \\times \\sum_{j=1}^{5} j \\cos \\left[(j+1)x_2 + j \\right] + (x_1 + 1.42513)^2 + (x_2 + 0.80032)^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    .. figure:: figures/Levy05.png
        :alt: Levy 5 function
        :align: center

        **Two-dimensional Levy 5 function**


    *Global optimum*: :math:`f(x_i) = -176.1375779` for :math:`\\mathbf{x} = [-1.30685, -1.42485]`.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)
        self.custom_bounds = ([-2.0, 2.0], [-2.0, 2.0])

        self.global_optimum = [-1.30685, -1.42485]
        self.fglob = -176.1375779

    def evaluator(self, x, *args):

        self.fun_evals += 1
        i = arange(1, 6)
        a = i * cos((i - 1) * x[0] + i)
        b = i * cos((i + 1) * x[1] + i)

        return sum(a) * sum(b) + (x[0] + 1.42513) ** 2 + (x[1] + 0.80032) ** 2

# -------------------------------------------------------------------------------- #


class Levy13(Benchmark):

    """
    Levy13 test objective function.

    This class defines the Levy13 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Levy13}}(\\mathbf{x}) = \\left(x_{1} -1\\right)^{2} \\left[\sin^{2}\\left(3 \\pi x_{2}\\right) + 1\\right] + \\left(x_{2} -1\\right)^{2} \\left[\\sin^{2}\\left(2 \\pi x_{2}\\right) + 1\\right] + \\sin^{2}\\left(3 \\pi x_{1}\\right)


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    .. figure:: figures/Levy13.png
        :alt: Levy13 function
        :align: center

        **Two-dimensional Levy13 function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 1` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)
        self.custom_bounds = [(-5, 5), (-5, 5)]

        self.global_optimum = [1 for _ in range(self.dimensions)]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return (sin(3 * pi * x[0])) ** 2 + ((x[0] - 1) ** 2) * (1 + (sin(3 * pi * x[1])) ** 2) + ((x[1] - 1) ** 2) * (1 + (sin(2 * pi * x[1])) ** 2)

# -------------------------------------------------------------------------------- #


class Matyas(Benchmark):

    """
    Matyas test objective function.

    This class defines the Matyas global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Matyas}}(\\mathbf{x}) = 0.26(x_1^2 + x_2^2) - 0.48x_1x_2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    .. figure:: figures/Matyas.png
        :alt: Matyas function
        :align: center

        **Two-dimensional Matyas function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [0 for _ in range(self.dimensions)]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return 0.26 * (x[0] ** 2 + x[1] ** 2) - 0.48 * x[0] * x[1]

# -------------------------------------------------------------------------------- #


class McCormick(Benchmark):

    """
    McCormick test objective function.

    This class defines the McCormick global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{McCormick}}(\\mathbf{x}) = - x_{1} + 2 x_{2} + \\left(x_{1} - x_{2}\\right)^{2} + \\sin\\left(x_{1} + x_{2}\\right) + 1


    Here, :math:`n` represents the number of dimensions and :math:`x_1 \\in [-1.5, 4]`, :math:`x_2 \\in [-3, 4]`.

    .. figure:: figures/McCormick.png
        :alt: McCormick function
        :align: center

        **Two-dimensional McCormick function**


    *Global optimum*: :math:`f(x_i) = -1.913222954981037` for :math:`\\mathbf{x} = [-0.5471975602214493, -1.547197559268372]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(-1.5, 4.0), (-3.0, 4.0)]

        self.global_optimum = [-0.5471975602214493, -1.547197559268372]
        self.fglob = -1.913222954981037

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return sin(x[0] + x[1]) + (x[0] - x[1]) ** 2 - 1.5 * x[0] + 2.5 * x[1] + 1

# -------------------------------------------------------------------------------- #


class Michalewicz(Benchmark):

    """
    Michalewicz test objective function.

    This class defines the Michalewicz global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Michalewicz}}(\\mathbf{x}) = - \\sum_{i=1}^{2} \\sin\\left(x_i\\right) \\sin^{2 m}\\left(\\frac{i x_i^{2}}{\\pi}\\right)


    Where, in this exercise, :math:`m = 10`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, \\pi]` for :math:`i=1,2`.

    .. figure:: figures/Michalewicz.png
        :alt: Michalewicz function
        :align: center

        **Two-dimensional Michalewicz function**


    *Global optimum*: :math:`f(x_i) = -1.8013` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([0.0] * self.dimensions, [pi] * self.dimensions)

        self.global_optimum = [2.20290555, 1.570796]
        self.fglob = -1.8013

    def evaluator(self, x, *args):

        self.fun_evals += 1

        m = 10.0
        i = arange(1, self.dimensions + 1)
        return -sum(sin(x) * sin(i * x ** 2 / pi) ** (2 * m))

# -------------------------------------------------------------------------------- #


class MieleCantrell(Benchmark):

    """
    Miele-Cantrell test objective function.

    This class defines the Miele-Cantrell global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{MieleCantrell}}(\\mathbf{x}) = (e^{-x_1} - x_2)^4 + 100(x_2 - x_3)^6 + \\tan^4(x_3 - x_4) + x_1^8


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,...,4`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [0, 1, 1, 1]`

    """

    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-1.0] * self.dimensions, [1.0] * self.dimensions)

        self.global_optimum = [0.0, 1.0, 1.0, 1.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2, x3, x4 = x
        return (exp(-x1) - x2) ** 4.0 + 100.0 * (x2 - x3) ** 6.0 + (tan(x3 - x4)) ** 4.0 + x1 ** 8.0

# -------------------------------------------------------------------------------- #


class Mishra01(Benchmark):

    """
    Mishra 1 test objective function.

    This class defines the Mishra 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra01}}(\\mathbf{x}) = (1 + x_n)^{x_n} \\hspace{10pt} ; \\hspace{10pt} x_n = n - \\sum_{i=1}^{n-1} x_i


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 1]` for :math:`i=1,...,n`.

    .. figure:: figures/Mishra01.png
        :alt: Mishra 1 function
        :align: center

        **Two-dimensional Mishra 1 function**


    *Global optimum*: :math:`f(x_i) = 2` for :math:`x_i = 1` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([0.0] * self.dimensions,
                          [1.0 + 1e-9] * self.dimensions)

        self.global_optimum = [1.0 for _ in range(self.dimensions)]
        self.fglob = 2.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        n = self.dimensions

        xn = n - sum(x[0:-1])
        return (1 + xn) ** xn

# -------------------------------------------------------------------------------- #


class Mishra02(Benchmark):

    """
    Mishra 2 test objective function.

    This class defines the Mishra 2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra02}}(\\mathbf{x}) = (1 + x_n)^{x_n} \\hspace{10pt} ; \\hspace{10pt} x_n = n - \\sum_{i=1}^{n-1} \\frac{(x_i + x_{i+1})}{2}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 1]` for :math:`i=1,...,n`.

    .. figure:: figures/Mishra02.png
        :alt: Mishra 2 function
        :align: center

        **Two-dimensional Mishra 2 function**


    *Global optimum*: :math:`f(x_i) = 2` for :math:`x_i = 1` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([0.0] * self.dimensions,
                          [1.0 + 1e-9] * self.dimensions)

        self.global_optimum = [1.0 for _ in range(self.dimensions)]
        self.fglob = 2.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        n = self.dimensions

        xn = n - sum((x[0:-1] + x[1:]) / 2.0)
        return (1 + xn) ** xn

# -------------------------------------------------------------------------------- #


class Mishra03(Benchmark):

    """
    Mishra 3 test objective function.

    This class defines the Mishra 3 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra03}}(\\mathbf{x}) = \\sqrt{\\lvert \\cos{\\sqrt{\\lvert x_1^2 + x_2^2 \\rvert}} \\rvert} + 0.01(x_1 + x_2)


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    .. figure:: figures/Mishra03.png
        :alt: Mishra 3 function
        :align: center

        **Two-dimensional Mishra 3 function**


    *Global optimum*: :math:`f(x_i) = -0.1999` for :math:`x_i = {-9.99378322, -9.99918927}`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [-9.99378322, -9.99918927]
        self.fglob = -0.19990562

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return ((0.01 * (x[0] + x[1])
                + sqrt(abs(cos(sqrt(abs(x[0] ** 2 + x[1] ** 2)))))))

# -------------------------------------------------------------------------------- #


class Mishra04(Benchmark):

    """
    Mishra 4 test objective function.

    This class defines the Mishra 4 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra04}}(\\mathbf{x}) = \\sqrt{\\lvert \\sin{\\sqrt{\\lvert x_1^2 + x_2^2 \\rvert}} \\rvert} + 0.01(x_1 + x_2)


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    .. figure:: figures/Mishra04.png
        :alt: Mishra 4 function
        :align: center

        **Two-dimensional Mishra 4 function**


    *Global optimum*: :math:`f(x_i) = -0.17767` for :math:`x_i = {-8.71499636, -9.0533148}`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [-8.71499636, -9.0533148]
        self.fglob = -0.17767

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return ((0.01 * (x[0] + x[1])
                + sqrt(abs(sin(sqrt(abs(x[0] ** 2 + x[1] ** 2)))))))

# -------------------------------------------------------------------------------- #


class Mishra05(Benchmark):

    """
    Mishra 5 test objective function.

    This class defines the Mishra 5 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra05}}(\\mathbf{x}) = \\left [ \\sin^2 ((\\cos(x_1) + \\cos(x_2))^2) + \\cos^2 ((\\sin(x_1) + \\sin(x_2))^2) + x_1 \\right ]^2 + 0.01(x_1 + x_2)


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    .. figure:: figures/Mishra05.png
        :alt: Mishra 5 function
        :align: center

        **Two-dimensional Mishra 5 function**


    *Global optimum*: :math:`f(x_i) = -0.119829` for :math:`\\mathbf{x} = [-1.98682, -10]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [-1.98682, -10.0]
        self.fglob = -0.119829

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2 = x
        return ((sin((cos(x1) + cos(x2)) ** 2.0) ** 2.0) + (cos((sin(x1) + sin(x2)) ** 2.0) ** 2.0) + x1) ** 2.0 + 0.01 * (x1 + x2)

# -------------------------------------------------------------------------------- #


class Mishra06(Benchmark):

    """
    Mishra 6 test objective function.

    This class defines the Mishra 6 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra06}}(\\mathbf{x}) = -\\log{\\left [ \\sin^2 ((\\cos(x_1) + \\cos(x_2))^2) - \\cos^2 ((\\sin(x_1) + \\sin(x_2))^2) + x_1 \\right ]^2} + 0.01 \\left[(x_1 -1)^2 + (x_2 - 1)^2 \\right]


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    .. figure:: figures/Mishra06.png
        :alt: Mishra 6 function
        :align: center

        **Two-dimensional Mishra 6 function**


    *Global optimum*: :math:`f(x_i) = -2.28395` for :math:`\\mathbf{x} = [2.88631, 1.82326]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [2.88631, 1.82326]
        self.fglob = -2.28395

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2 = x
        return -log(((sin((cos(x1) + cos(x2)) ** 2.0) ** 2.0) - (cos((sin(x1) + sin(x2)) ** 2.0) ** 2.0) + x1) ** 2.0) + 0.1 * ((x1 - 1.0) ** 2.0 + (x2 - 1.0) ** 2.0)

# -------------------------------------------------------------------------------- #


class Mishra07(Benchmark):

    """
    Mishra 7 test objective function.

    This class defines the Mishra 7 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra07}}(\\mathbf{x}) = \\left [\\prod_{i=1}^{n} x_i - n! \\right]^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    .. figure:: figures/Mishra07.png
        :alt: Mishra 7 function
        :align: center

        **Two-dimensional Mishra 7 function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = \\sqrt{n}` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)
        self.custom_bounds = [(-2, 2), (-2, 2)]
        self.global_optimum = [sqrt(self.dimensions)
                               for i in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return (prod(x) - factorial(self.dimensions)) ** 2.0

# -------------------------------------------------------------------------------- #


class Mishra08(Benchmark):

    """
    Mishra 8 test objective function.

    This class defines the Mishra 8 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra08}}(\\mathbf{x}) = 0.001 \\left[\\lvert x_1^{10} - 20x_1^9 + 180x_1^8 - 960 x_1^7 + 3360x_1^6 - 8064x_1^5 + 13340x_1^4 - 15360x_1^3 + 11520x_1^2 - 5120x_1 + 2624 \\rvert \\lvert x_2^4 + 12x_2^3 + 54x_2^2 + 108x_2 + 81 \\rvert \\right]^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    .. figure:: figures/Mishra08.png
        :alt: Mishra 8 function
        :align: center

        **Two-dimensional Mishra 8 function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [2, -3]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)
        self.custom_bounds = [(1.0, 2.0), (-4.0, 1.0)]
        self.global_optimum = [2.0, -3.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        F1 = abs(x[0] ** 10 - 20 * x[0] ** 9 + 180 * x[0] ** 8 - 960 * x[0] ** 7 + 3360 * x[0] ** 6 -
                 8064 * x[0] ** 5 + 13340 * x[0] ** 4 - 15360 * x[0] ** 3 + 11520 * x[0] ** 2 - 5120 * x[0] + 2624)
        F2 = abs(x[1] ** 4 + 12 * x[1] ** 3 + 54 *
                 x[1] ** 2 + 108 * x[1] + 81.0)
        return 0.001 * (F1 + F2) ** 2

# -------------------------------------------------------------------------------- #


class Mishra09(Benchmark):

    """
    Mishra 9 test objective function.

    This class defines the Mishra 9 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra09}}(\\mathbf{x}) = \\left[ ab^2c + abc^2 + b^2 + (x_1 + x_2 - x_3)^2 \\right]^2


    Where, in this exercise:

    .. math::

        \\begin{cases} a = 2x_1^3 + 5x_1x_2^2 + 4x_3 - 2x_1^2x_3 - 18 \\\\
        b = x_1 + x_2^3 + x_1x_3^2 - 22 \\\\
        c = 8x_1^2 + 2x_2x_3 + 2x_2^2 + 3x_2^3 - 52 \\end{cases}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2,3`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [1, 2, 3]`

    """

    def __init__(self, dimensions=3):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)
        self.global_optimum = [1.0, 2.0, 3.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2, x3 = x
        F1 = 2 * x1 ** 3 + 5 * x1 * x2 + 4 * x3 - 2 * x1 ** 2 * x3 - 18.0
        F2 = x1 + x2 ** 3 + x1 * x2 ** 2 + x1 * x3 ** 2 - 22.0
        F3 = 8 * x1 ** 2 + 2 * x2 * x3 + 2 * x2 ** 2 + 3 * x2 ** 3 - 52.0
        return (F1 * F3 * F2 ** 2 + F1 * F2 * F3 ** 2 + F2 ** 2 + (x1 + x2 - x3) ** 2) ** 2

# -------------------------------------------------------------------------------- #


class Mishra10(Benchmark):

    """
    Mishra 10 test objective function.

    This class defines the Mishra 10 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra10}}(\\mathbf{x}) = \\left[ \\lfloor x_1 \\perp x_2 \\rfloor - \\lfloor x_1 \\rfloor - \\lfloor x_2 \\rfloor \\right]^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    .. figure:: figures/Mishra10.png
        :alt: Mishra 10 function
        :align: center

        **Two-dimensional Mishra 10 function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [2, 2]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)
        self.global_optimum = [2.0, 2.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1
        x1, x2 = int(x[0]), int(x[1])

        f1 = x1 + x2
        f2 = x1 * x2
        return (f1 - f2) ** 2.0

# -------------------------------------------------------------------------------- #


class Mishra11(Benchmark):

    """
    Mishra 11 test objective function.

    This class defines the Mishra 11 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra11}}(\\mathbf{x}) = \\left [ \\frac{1}{n} \\sum_{i=1}^{n} \\lvert x_i \\rvert - \\left(\\prod_{i=1}^{n} \\lvert x_i \\rvert \\right )^{\\frac{1}{n}} \\right]^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    .. figure:: figures/Mishra11.png
        :alt: Mishra 11 function
        :align: center

        **Two-dimensional Mishra 11 function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)
        self.custom_bounds = [(-3, 3), (-3, 3)]

        self.global_optimum = [0.0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        n = self.dimensions

        return ((1.0 / n) * sum(abs(x)) - (prod(abs(x))) ** 1.0 / n) ** 2.0

# -------------------------------------------------------------------------------- #


class MultiModal(Benchmark):

    """
    MultiModal test objective function.

    This class defines the MultiModal global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{MultiModal}}(\\mathbf{x}) = \\left( \\sum_{i=1}^n \\lvert x_i \\rvert \\right) \\left( \\prod_{i=1}^n \\lvert x_i \\rvert \\right)


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    .. figure:: figures/MultiModal.png
        :alt: MultiModal function
        :align: center

        **Two-dimensional MultiModal function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)
        self.custom_bounds = [(-5, 5), (-5, 5)]

        self.global_optimum = [0.0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return sum(abs(x)) * prod(abs(x))

# -------------------------------------------------------------------------------- #


class NeedleEye(Benchmark):

    """
    NeedleEye test objective function.

    This class defines the Needle-Eye global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{NeedleEye}}(\\mathbf{x}) = \\begin{cases} 1 & \\textrm{if} \\hspace{5pt} \\lvert x_i \\rvert  <  eye \\hspace{5pt} \\forall i \\\\
               \\sum_{i=1}^n (100 + \\lvert x_i \\rvert) & \\textrm{if} \\hspace{5pt} \\lvert x_i \\rvert > eye \\\\
               0 & \\textrm{otherwise} \\end{cases}

    Where, in this exercise, :math:`eye = 0.0001`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    .. figure:: figures/NeedleEye.png
        :alt: NeedleEye function
        :align: center

        **Two-dimensional NeedleEye function**


    *Global optimum*: :math:`f(x_i) = 1` for :math:`x_i = -1` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [-1.0 for _ in range(self.dimensions)]
        self.fglob = 1.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        f = fp = 0.0
        eye = 0.0001

        for i in range(self.dimensions):
            if abs(x[i]) >= eye:
                fp = 1.0
                f += 100.0 + abs(x[i])
            else:
                f += 1.0

        if fp < 1e-6:
            f = f / self.dimensions

        return f

# -------------------------------------------------------------------------------- #


class NewFunction01(Benchmark):

    """
    NewFunction01 test objective function.

    This class defines the NewFunction01 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{NewFunction01}}(\\mathbf{x}) = \\left | {\\cos\\left(\\sqrt{\\left|{x_{1}^{2} + x_{2}}\\right|}\\right)} \\right |^{0.5} + (x_{1} + x_{2})/100


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    .. figure:: figures/NewFunction01.png
        :alt: NewFunction01 function
        :align: center

        **Two-dimensional NewFunction01 function**


    *Global optimum*: :math:`f(x_i) = -0.184642678` for :math:`\\mathbf{x} = [-8.46669057, -9.99982177]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [-8.46669057, -9.99982177]
        self.fglob = -0.184642678

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return (abs(cos(sqrt(abs(x[0] ** 2 + x[1]))))) ** 0.5 + 0.01 * (x[0] + x[1])

# -------------------------------------------------------------------------------- #


class NewFunction02(Benchmark):

    """
    NewFunction02 test objective function.

    This class defines the NewFunction02 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{NewFunction02}}(\\mathbf{x}) = \\left | {\\sin\\left(\\sqrt{\\lvert{x_{1}^{2} + x_{2}}\\rvert}\\right)} \\right |^{0.5} + (x_{1} + x_{2})/100


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    .. figure:: figures/NewFunction02.png
        :alt: NewFunction02 function
        :align: center

        **Two-dimensional NewFunction02 function**


    *Global optimum*: :math:`f(x_i) = -0.19937167` for :math:`\\mathbf{x} = [-9.94103375, -9.99771235]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [-9.94103375, -9.99771235]
        self.fglob = -0.19937167547710213

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return (abs(sin(sqrt(abs(x[0] ** 2 + x[1]))))) ** 0.5 + 0.01 * (x[0] + x[1])

# -------------------------------------------------------------------------------- #


class NewFunction03(Benchmark):

    """
    NewFunction03 test objective function.

    This class defines the NewFunction03 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{NewFunction03}}(\\mathbf{x}) = 0.01 x_{1} + 0.1 x_{2} + \\left\{x_{1} + \\sin^{2}\\left[\\left(\\cos\\left(x_{1}\\right) + \\cos\\left(x_{2}\\right)\\right)^{2}\\right] + \\cos^{2}\\left[\\left(\\sin\\left(x_{1}\\right) + \\sin\\left(x_{2}\\right)\\right)^{2}\\right]\\right\}^{2}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    .. figure:: figures/NewFunction03.png
        :alt: NewFunction03 function
        :align: center

        **Two-dimensional NewFunction03 function**


    *Global optimum*: :math:`f(x_i) = -1.019829` for :math:`\\mathbf{x} = [-1.98682, -10]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [-1.98682, -10.0]
        self.fglob = -1.019829

    def evaluator(self, x, *args):

        self.fun_evals += 1
        f1 = sin((cos(x[0]) + cos(x[1])) ** 2) ** 2
        f2 = cos((sin(x[0]) + sin(x[1])) ** 2) ** 2
        f = (f1 + f2 + x[0]) ** 2
        f = f + 0.01 * x[0] + 0.1 * x[1]

        return f

# -------------------------------------------------------------------------------- #


class OddSquare(Benchmark):

    """
    Odd Square test objective function.

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

    .. figure:: figures/OddSquare.png
        :alt: Odd Square function
        :align: center

        **Two-dimensional Odd Square function**


    *Global optimum*: :math:`f(x_i) = -1.0084` for :math:`\\mathbf{x} \\approx b`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-5.0 * pi] * self.dimensions,
                          [5.0 * pi] * self.dimensions)
        self.custom_bounds = ([-2.0, 4.0], [-2.0, 4.0])

        a = asarray([1, 1.3, 0.8, -0.4, -1.3, 1.6, -0.2, -0.6, 0.5, 1.4] * 2)

        self.global_optimum = a[0:self.dimensions]
        self.fglob = -1.0084

    def evaluator(self, x, *args):

        self.fun_evals += 1

        c = 0.02
        a = asarray([1, 1.3, 0.8, -0.4, -1.3, 1.6, -0.2, -0.6, 0.5, 1.4] * 2)

        b = a[0:self.dimensions]
        d = self.dimensions * max((x - b) ** 2.0)
        h = sum((x - b) ** 2.0)

        return -exp(-d / (2.0 * pi)) * cos(pi * d) * (1.0 + 0.02 * h / (d + 0.01))

# -------------------------------------------------------------------------------- #


class Parsopoulos(Benchmark):

    """
    Parsopoulos test objective function.

    This class defines the Parsopoulos global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Parsopoulos}}(\\mathbf{x}) = \\cos(x_1)^2 + \\sin(x_2)^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,2`.

    .. figure:: figures/Parsopoulos.png
        :alt: Parsopoulos function
        :align: center

        **Two-dimensional Parsopoulos function**


    *Global optimum*: This function has inï¬nite number of global minima in R2, at points :math:`\\left(k\\frac{\\pi}{2}, \\lambda \\pi \\right)`,
    where :math:`k = \\pm1, \\pm3, ...` and :math:`\\lambda = 0, \\pm1, \\pm2, ...`

    In the given domain problem, function has 12 global minima all equal to zero.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-5.0] * self.dimensions, [5.0] * self.dimensions)

        self.global_optimum = [pi / 2.0, pi]
        self.fglob = 0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2 = x
        return cos(x1) ** 2.0 + sin(x2) ** 2.0

# -------------------------------------------------------------------------------- #


class Pathological(Benchmark):

    """
    Pathological test objective function.

    This class defines the Pathological global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Pathological}}(\\mathbf{x}) = \\sum_{i=1}^{n -1} \\frac{\\sin^{2}\\left(\\sqrt{100 x_{i+1}^{2} + x_{i}^{2}}\\right) -0.5}{0.001 \\left(x_{i}^{2} - 2x_{i}x_{i+1} + x_{i+1}^{2}\\right)^{2} + 0.50}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    .. figure:: figures/Pathological.png
        :alt: Pathological function
        :align: center

        **Two-dimensional Pathological function**


    *Global optimum*: :math:`f(x_i) = 0.` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-100.0] * self.dimensions,
                          [100.0] * self.dimensions)

        self.global_optimum = [0 for _ in range(self.dimensions)]
        self.fglob = 0.

    def evaluator(self, x, *args):

        self.fun_evals += 1

        vec = (0.5 + (sin(sqrt(100 * x[: -1] ** 2 + x[1:] ** 2)) ** 2 - 0.5) /
               (1. + 0.001 * (x[: -1] ** 2 - 2 * x[: -1] * x[1:]
                              + x[1:] ** 2) ** 2))
        return sum(vec)

# -------------------------------------------------------------------------------- #


class Paviani(Benchmark):

    """
    Paviani test objective function.

    This class defines the Paviani global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Paviani}}(\\mathbf{x}) = \\sum_{i=1}^{10} \\left[\\log^{2}\\left(10 - x_i\\right) + \\log^{2}\\left(x_i -2\\right)\\right] - \\left(\\prod_{i=1}^{10} x_i^{10} \\right)^{0.2}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [2.001, 9.999]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = -45.7784684040686` for :math:`x_i = 9.350266` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=10):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([2.001] * self.dimensions, [9.999] * self.dimensions)

        self.global_optimum = [9.350266 for _ in range(self.dimensions)]
        self.fglob = -45.7784684040686

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return sum(log(x - 2) ** 2.0 + log(10.0 - x) ** 2.0) - prod(x) ** 0.2

# -------------------------------------------------------------------------------- #


class Penalty01(Benchmark):

    """
    Penalty 1 test objective function.

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

    .. figure:: figures/Penalty01.png
        :alt: Penalty 1 function
        :align: center

        **Two-dimensional Penalty 1 function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = -1` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-50.0] * self.dimensions, [50.0] * self.dimensions)
        self.custom_bounds = ([-5.0, 5.0], [-5.0, 5.0])

        self.global_optimum = [-1.0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        a, b, c = 10.0, 100.0, 4.0
        xx = abs(x)
        u = where(xx > a, b * (xx - a) ** c, 0.0)

        y = 1.0 + (x + 1.0) / 4.0
        return sum(u) + (pi / 30.0) * (10.0 * sin(pi * y[0]) ** 2.0 + sum((y[0:-1] - 1.0) ** 2.0 * (1.0 + 10.0 * sin(pi * y[1:]) ** 2.0)) + (y[-1] - 1) ** 2.0)


# -------------------------------------------------------------------------------- #

class Penalty02(Benchmark):

    """
    Penalty 2 test objective function.

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

    .. figure:: figures/Penalty02.png
        :alt: Penalty 2 function
        :align: center

        **Two-dimensional Penalty 2 function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 1` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-50.0] * self.dimensions, [50.0] * self.dimensions)
        self.custom_bounds = ([-4.0, 4.0], [-4.0, 4.0])

        self.global_optimum = [1.0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        a, b, c = 5.0, 100.0, 4.0
        xx = abs(x)
        u = where(xx > a, b * (xx - a) ** c, 0.0)

        return sum(u) + 0.1 * (10.0 * sin(3.0 * pi * x[0]) ** 2.0 + sum((x[0:-1] - 1.0) ** 2.0 * (1.0 + sin(pi * x[1:]) ** 2.0)) + (x[-1] - 1) ** 2.0 * (1 + sin(2 * pi * x[-1]) ** 2.0))

# -------------------------------------------------------------------------------- #


class PenHolder(Benchmark):

    """
    PenHolder test objective function.

    This class defines the PenHolder global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{PenHolder}}(\\mathbf{x}) = -e^{\\left|{e^{\\left|{- \\frac{\\sqrt{x_{1}^{2} + x_{2}^{2}}}{\\pi} + 1}\\right|} \\cos\\left(x_{1}\\right) \\cos\\left(x_{2}\\right)}\\right|^{-1}}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-11, 11]` for :math:`i=1,2`.

    .. figure:: figures/PenHolder.png
        :alt: PenHolder function
        :align: center

        **Two-dimensional PenHolder function**


    *Global optimum*: :math:`f(x_i) = -0.9635348327265058` for :math:`x_i = \\pm 9.646167671043401` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-11.0] * self.dimensions, [11.0] * self.dimensions)

        self.global_optimum = [-9.646167708023526, 9.646167671043401]
        self.fglob = -0.9635348327265058

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return -exp(-(abs(cos(x[0]) * cos(x[1]) * exp(abs(1 - sqrt(x[0] ** 2 + x[1] ** 2) / pi)))) ** (-1))

# -------------------------------------------------------------------------------- #


class PermFunction01(Benchmark):

    """
    PermFunction 1 test objective function.

    This class defines the Perm Function 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{PermFunction01}}(\\mathbf{x}) = \\sum_{k=1}^n \\left\\{ \\sum_{j=1}^n (j^k + \\beta) \\left[ \\left(\\frac{x_j}{j}\\right)^k - 1 \\right] \\right\\}^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-n, n+1]` for :math:`i=1,...,n`.

    .. figure:: figures/PermFunction01.png
        :alt: PermFunction 1 function
        :align: center

        **Two-dimensional PermFunction 1 function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = i` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-self.dimensions] * self.dimensions,
                          [self.dimensions + 1] * self.dimensions)

        self.global_optimum = range(1, self.dimensions + 1)
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        b = 0.5
        s_out = 0.0
        for k in range(1, self.dimensions + 1):
            j = arange(1, self.dimensions + 1)
            s_in = (j ** k + b) * ((x[j - 1] / j) ** k - 1)
            s_out += sum(s_in ** 2)

        return s_out

# -------------------------------------------------------------------------------- #


class PermFunction02(Benchmark):

    """
    PermFunction 2 test objective function.

    This class defines the Perm Function 2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{PermFunction02}}(\\mathbf{x}) = \\sum_{k=1}^n \\left\\{ \\sum_{j=1}^n (j + \\beta) \\left[ \\left(x_j^k - {\\frac{1}{j}}^{k} \\right ) \\right] \\right\\}^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-n, n+1]` for :math:`i=1,...,n`.

    .. figure:: figures/PermFunction02.png
        :alt: PermFunction 2 function
        :align: center

        **Two-dimensional PermFunction 2 function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = \\frac{1}{i}` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-self.dimensions] * self.dimensions,
                          [self.dimensions + 1] * self.dimensions)
        self.custom_bounds = ([0, 1.5], [0, 1.0])

        self.global_optimum = 1. / arange(1, self.dimensions + 1)
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):
        self.fun_evals += 1
        b = 10
        outer = 0
        j = arange(1, self.dimensions + 1)
        for k in range(1, self.dimensions + 1):
            inner = (j + b) * (x[j - 1] ** k - (1. / j) ** k)
            outer += sum(inner ** 2)
        return outer
# -------------------------------------------------------------------------------- #


class Pinter(Benchmark):

    """
    Pinter test objective function.

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

    .. figure:: figures/Pinter.png
        :alt: Pinter function
        :align: center

        **Two-dimensional Pinter function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [0.0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        n = self.dimensions

        f = 0.0

        for i in range(n):
            x_i = x[i]

            if i == 0:
                x_mi = x[-1]
                x_pi = x[i + 1]
            elif i == n - 1:
                x_mi = x[i - 1]
                x_pi = x[0]
            else:
                x_mi = x[i - 1]
                x_pi = x[i + 1]

            A = x_mi * sin(x_i) + sin(x_pi)
            B = x_mi ** 2.0 - 2 * x_i + 3 * x_pi - cos(x_i) + 1.0

            f += (i + 1.0) * x_i ** 2.0 + 20.0 * (i + 1.0) * sin(A) ** 2.0 + \
                (i + 1.0) * log10(1.0 + (i + 1.0) * B ** 2.0)

        return f

# -------------------------------------------------------------------------------- #


class Plateau(Benchmark):

    """
    Plateau test objective function.

    This class defines the Plateau global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Plateau}}(\\mathbf{x}) = 30 + \\sum_{i=1}^n \\lfloor x_i \\rfloor


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5.12, 5.12]` for :math:`i=1,...,n`.

    .. figure:: figures/Plateau.png
        :alt: Plateau function
        :align: center

        **Two-dimensional Plateau function**


    *Global optimum*: :math:`f(x_i) = 30` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-5.12] * self.dimensions, [5.12] * self.dimensions)

        self.global_optimum = [0.0 for _ in range(self.dimensions)]
        self.fglob = 30.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return 30.0 + sum(floor(abs(x)))

# -------------------------------------------------------------------------------- #


class Powell(Benchmark):

    """
    Powell test objective function.

    This class defines the Powell global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Powell}}(\\mathbf{x}) = (x_3+10x_1)^2+5(x_2-x_4)^2+(x_1-2x_2)^4+10(x_3-x_4)^4


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-4, 5]` for :math:`i=1,...,4`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,4`

    """

    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-4.0] * self.dimensions, [5.0] * self.dimensions)
        self.global_optimum = [0, 0, 0, 0]
        self.fglob = 0

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return (x[0] + 10 * x[1]) ** 2 + 5 * (x[2] - x[3]) ** 2 + (x[1] - 2 * x[2]) ** 4 + 10 * (x[0] - x[3]) ** 4

# -------------------------------------------------------------------------------- #


class PowerSum(Benchmark):

    """
    Power sum test objective function.

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

        self.bounds = zip([0.0] * self.dimensions,
                          [float(self.dimensions)] * self.dimensions)

        self.global_optimum = [1.0, 2.0, 2.0, 3.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        b = [8.0, 18.0, 44.0, 114.0]
        y = 0.0

        for k in range(1, self.dimensions + 1):
            s_in = 0.0
            for i in range(self.dimensions):
                s_in = s_in + x[i] ** k

            y = y + (s_in - b[k - 1]) ** 2.0

        return y

# -------------------------------------------------------------------------------- #


class Price01(Benchmark):

    """
    Price 1 test objective function.

    This class defines the Price 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Price01}}(\\mathbf{x}) = (\\lvert x_1 \\rvert - 5)^2 + (\\lvert x_2 \\rvert - 5)^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-500, 500]` for :math:`i=1,2`.

    .. figure:: figures/Price01.png
        :alt: Price 1 function
        :align: center

        **Two-dimensional Price 1 function**


    *Global optimum*: :math:`f(x_i) = 0.0` for :math:`\\mathbf{x} = [5, 5]` or :math:`\\mathbf{x} = [5, -5]`
    or :math:`\\mathbf{x} = [-5, 5]` or :math:`\\mathbf{x} = [-5, -5]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-500.0] * self.dimensions,
                          [500.0] * self.dimensions)
        self.custom_bounds = ([-10.0, 10.0], [-10.0, 10.0])

        self.global_optimum = [5.0, 5.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2 = x
        return (abs(x1) - 5.0) ** 2.0 + (abs(x2) - 5.0) ** 2.0

# -------------------------------------------------------------------------------- #


class Price02(Benchmark):

    """
    Price 2 test objective function.

    This class defines the Price 2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Price02}}(\\mathbf{x}) = 1 + \\sin^2(x_1) + \\sin^2(x_2) - 0.1e^{(-x_1^2 - x_2^2)}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    .. figure:: figures/Price02.png
        :alt: Price 2 function
        :align: center

        **Two-dimensional Price 2 function**


    *Global optimum*: :math:`f(x_i) = 0.9` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [0.0, 0.0]
        self.fglob = 0.9

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2 = x
        return 1.0 + sin(x1) ** 2.0 + sin(x2) ** 2.0 - 0.1 * exp(-x1 ** 2.0 - x2 ** 2.0)

# -------------------------------------------------------------------------------- #


class Price03(Benchmark):

    """
    Price 3 test objective function.

    This class defines the Price 3 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Price03}}(\\mathbf{x}) = 100(x_2 - x_1^2)^2 + \\left[6.4(x_2 - 0.5)^2 - x_1 - 0.6 \\right]^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-50, 50]` for :math:`i=1,2`.

    .. figure:: figures/Price03.png
        :alt: Price 3 function
        :align: center

        **Two-dimensional Price 3 function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [-5, -5]`, :math:`\\mathbf{x} = [-5, 5]`,
    :math:`\\mathbf{x} = [5, -5]`, :math:`\\mathbf{x} = [5, 5]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-50.0] * self.dimensions, [50.0] * self.dimensions)
        self.custom_bounds = ([0, 2], [0, 2])

        self.global_optimum = [1.0, 1.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):
        self.fun_evals += 1
        return (100 * (x[1] - x[0] ** 2) ** 2
                + (6.4 * (x[1] - 0.5) ** 2 - x[0] - 0.6) ** 2)

# -------------------------------------------------------------------------------- #


class Price04(Benchmark):

    """
    Price 4 test objective function.

    This class defines the Price 4 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Price04}}(\\mathbf{x}) = (2x_1^3x_2 - x_2^3)^2 + (6x_1 - x_2^2 + x_2)^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-50, 50]` for :math:`i=1,2`.

    .. figure:: figures/Price04.png
        :alt: Price 4 function
        :align: center

        **Two-dimensional Price 4 function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [0, 0]`, :math:`\\mathbf{x} = [2, 4]` and
    :math:`\\mathbf{x} = [1.464, -2.506]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-50.0] * self.dimensions, [50.0] * self.dimensions)
        self.custom_bounds = ([0, 2], [0, 2])

        self.global_optimum = [2.0, 4.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2 = x
        return (2.0 * x2 * x1 ** 3.0 - x2 ** 3.0) ** 2.0 + (6.0 * x1 - x2 ** 2.0 + x2) ** 2.0

# -------------------------------------------------------------------------------- #


class Qing(Benchmark):

    """
    Qing test objective function.

    This class defines the Qing global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Qing}}(\\mathbf{x}) = \\sum_{i=1}^{n} (x_i^2 - i)^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-500, 500]` for :math:`i=1,...,n`.

    .. figure:: figures/Qing.png
        :alt: Qing function
        :align: center

        **Two-dimensional Qing function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = \\pm \\sqrt(i)` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-500.0] * self.dimensions,
                          [500.0] * self.dimensions)
        self.custom_bounds = [(-2, 2), (-2, 2)]

        self.global_optimum = [sqrt(_) for _ in range(1, self.dimensions + 1)]
        self.fglob = 0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        rng = arange(1, self.dimensions + 1)
        return sum((x ** 2.0 - rng) ** 2.0)

# -------------------------------------------------------------------------------- #


class Quadratic(Benchmark):

    """
    Quadratic test objective function.

    This class defines the Quadratic global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Quadratic}}(\\mathbf{x}) = -3803.84 - 138.08x_1 - 232.92x_2 + 128.08x_1^2 + 203.64x_2^2 + 182.25x_1x_2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    .. figure:: figures/Quadratic.png
        :alt: Quadratic function
        :align: center

        **Two-dimensional Quadratic function**


    *Global optimum*: :math:`f(x_i) = -3873.72418` for :math:`\\mathbf{x} = [0.19388, 0.48513]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)
        self.custom_bounds = [(0, 1), (0, 1)]

        self.global_optimum = [0.19388, 0.48513]
        self.fglob = -3873.72418
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2 = x
        return -3803.84 - 138.08 * x1 - 232.92 * x2 + 128.08 * x1 ** 2.0 + 203.64 * x2 ** 2.0 + 182.25 * x1 * x2

# -------------------------------------------------------------------------------- #


class Quintic(Benchmark):

    """
    Quintic test objective function.

    This class defines the Quintic global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Quintic}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\left|{x_{i}^{5} - 3 x_{i}^{4} + 4 x_{i}^{3} + 2 x_{i}^{2} - 10 x_{i} -4}\\right|


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    .. figure:: figures/Quintic.png
        :alt: Quintic function
        :align: center

        **Two-dimensional Quintic function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = -1` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)
        self.custom_bounds = [(-2, 2), (-2, 2)]

        self.global_optimum = [-1.0 for _ in range(self.dimensions)]
        self.fglob = 0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return sum(abs(x ** 5 - 3 * x ** 4 + 4 * x ** 3 + 2 * x ** 2 - 10 * x - 4))

# -------------------------------------------------------------------------------- #


class Rana(Benchmark):

    """
    Rana test objective function.

    This class defines the Rana global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Rana}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\left[x_{i} \\sin\\left(\\sqrt{\\lvert{x_{1} - x_{i} + 1}\\rvert}\\right) \\cos\\left(\\sqrt{\\lvert{x_{1} + x_{i} + 1}\\rvert}\\right) + \\left(x_{1} + 1\\right) \\sin\\left(\\sqrt{\\lvert{x_{1} + x_{i} + 1}\\rvert}\\right) \\cos\\left(\\sqrt{\\lvert{x_{1} - x_{i} + 1}\\rvert}\\right)\\right]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-500.000001, 500.000001]` for :math:`i=1,...,n`.

    .. figure:: figures/Rana.png
        :alt: Rana function
        :align: center

        **Two-dimensional Rana function**

    *Global optimum*: :math:`f(x_i) = -928.5478` for :math:`x_i = -500` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-500.000001] * self.dimensions,
                          [500.000001] * self.dimensions)

        self.global_optimum = [-300.3376, 500.]
        self.fglob = -500.8021602966615
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        t1 = sqrt(abs(x[1:] + x[: -1] + 1))
        t2 = sqrt(abs(x[1:] - x[: -1] + 1))
        return sum((x[1:] + 1) * cos(t2) * sin(t1) + x[:-1] * cos(t1) * sin(t2))


# -------------------------------------------------------------------------------- #

class Rastrigin(Benchmark):

    """
    Rastrigin test objective function.

    This class defines the Rastrigin global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Rastrigin}}(\\mathbf{x}) = 10n \\sum_{i=1}^n \\left[ x_i^2 - 10 \\cos(2\\pi x_i) \\right]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5.12, 5.12]` for :math:`i=1,...,n`.

    .. figure:: figures/Rastrigin.png
        :alt: Rastrigin function
        :align: center

        **Two-dimensional Rastrigin function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self.bounds = zip([-5.12] * self.dimensions, [5.12] * self.dimensions)

        self.global_optimum = [0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return 10.0 * self.dimensions + sum(x ** 2.0 - 10.0 * cos(2.0 * pi * x))

# -------------------------------------------------------------------------------- #


class Ripple01(Benchmark):

    """
    Ripple 1 test objective function.

    This class defines the Ripple 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Ripple01}}(\\mathbf{x}) = \\sum_{i=1}^2 -e^{-2 \\log 2 (\\frac{x_i-0.1}{0.8})^2} \\left[\\sin^6(5 \\pi x_i) + 0.1\\cos^2(500 \\pi x_i) \\right]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 1]` for :math:`i=1,...,n`.

    .. figure:: figures/Ripple01.png
        :alt: Ripple 1 function
        :align: center

        **Two-dimensional Ripple 1 function**

    *Global optimum*: :math:`f(x_i) = -2.2` for :math:`x_i = 0.1` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self.bounds = zip([0.0] * self.dimensions, [1.0] * self.dimensions)

        self.global_optimum = [0.1 for _ in range(self.dimensions)]
        self.fglob = -2.2

    def evaluator(self, x, *args):

        self.fun_evals += 1

        return sum(-exp(-2.0 * log(2.0) * ((x - 0.1) / 0.8) ** 2.0) * (sin(5.0 * pi * x) ** 6.0 + 0.1 * cos(500.0 * pi * x) ** 2.0))

# -------------------------------------------------------------------------------- #


class Ripple25(Benchmark):

    """
    Ripple 25 test objective function.

    This class defines the Ripple 25 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Ripple25}}(\\mathbf{x}) = \\sum_{i=1}^2 -e^{-2 \\log 2 (\\frac{x_i-0.1}{0.8})^2} \\left[\\sin^6(5 \\pi x_i) \\right]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 1]` for :math:`i=1,...,n`.

    .. figure:: figures/Ripple25.png
        :alt: Ripple 25 function
        :align: center

        **Two-dimensional Ripple 25 function**

    *Global optimum*: :math:`f(x_i) = -2` for :math:`x_i = 0.1` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self.bounds = zip([0.0] * self.dimensions, [1.0] * self.dimensions)

        self.global_optimum = [0.1 for _ in range(self.dimensions)]
        self.fglob = -2.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        return sum(-exp(-2.0 * log(2.0) * ((x - 0.1) / 0.8) ** 2.0) * (sin(5.0 * pi * x) ** 6.0))

# -------------------------------------------------------------------------------- #


class Rosenbrock(Benchmark):

    """
    Rosenbrock test objective function.

    This class defines the Rosenbrock global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Rosenbrock}}(\\mathbf{x}) = \\sum_{i=1}^{n-1} [100(x_i^2 - x_{i+1})^2 + (x_i - 1)^2]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 10]` for :math:`i=1,...,n`.

    .. figure:: figures/Rosenbrock.png
        :alt: Rosenbrock function
        :align: center

        **Two-dimensional Rosenbrock function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 1` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-5.0] * self.dimensions, [10.0] * self.dimensions)
        self.custom_bounds = [(-2, 2), (-2, 2)]

        self.global_optimum = [1 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return sum(100.0 * (x[1:] - x[:-1] ** 2.0) ** 2.0 + (1 - x[:-1]) ** 2.0)

# -------------------------------------------------------------------------------- #


class RosenbrockModified(Benchmark):

    """
    Modified Rosenbrock test objective function.

    This class defines the Modified Rosenbrock global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{RosenbrockModified}}(\\mathbf{x}) = 74 + 100(x_2 - x_1^2)^2 + (1 - x_1)^2 - 400 e^{-\\frac{(x_1+1)^2 + (x_2 + 1)^2}{0.1}}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-2, 2]` for :math:`i=1,2`.

    .. figure:: figures/RosenbrockModified.png
        :alt: Modified Rosenbrock function
        :align: center

        **Two-dimensional Modified Rosenbrock function**

    *Global optimum*: :math:`f(x_i) = 34.04024310` for :math:`\\mathbf{x} = [-0.90955374, -0.95057172]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-2.0] * self.dimensions, [2.0] * self.dimensions)
        self.custom_bounds = ([-1.0, 0.5], [-1.0, 1.0])

        self.global_optimum = [-0.90955374, -0.95057172]
        self.fglob = 34.040243106640844

    def evaluator(self, x, *args):

        self.fun_evals += 1
        a = 74 + 100. * (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2
        a -= 400 * exp(-((x[0] + 1.) ** 2 + (x[1] + 1.) ** 2) / 0.1)
        return a

# -------------------------------------------------------------------------------- #


class RotatedEllipse01(Benchmark):

    """
    Rotated Ellipse 1 test objective function.

    This class defines the Rotated Ellipse 1 global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{RotatedEllipse01}}(\\mathbf{x}) = 7x_1^2 - 6 \\sqrt{3} x_1x_2 + 13x_2^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-500, 500]` for :math:`i=1,2`.

    .. figure:: figures/RotatedEllipse01.png
        :alt: Rotated Ellipse 1 function
        :align: center

        **Two-dimensional Rotated Ellipse 1 function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [0, 0]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-500.0] * self.dimensions,
                          [500.0] * self.dimensions)
        self.custom_bounds = ([-2.0, 2.0], [-2.0, 2.0])

        self.global_optimum = [0.0, 0.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2 = x
        return 7.0 * x1 ** 2.0 - 6.0 * sqrt(3) * x1 * x2 + 13 * x2 ** 2.0

# -------------------------------------------------------------------------------- #


class RotatedEllipse02(Benchmark):

    """
    Rotated Ellipse 2 test objective function.

    This class defines the Rotated Ellipse 2 global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{RotatedEllipse02}}(\\mathbf{x}) = x_1^2 - x_1x_2 + x_2^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-500, 500]` for :math:`i=1,2`.

    .. figure:: figures/RotatedEllipse02.png
        :alt: Rotated Ellipse 2 function
        :align: center

        **Two-dimensional Rotated Ellipse 2 function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [0, 0]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-500.0] * self.dimensions,
                          [500.0] * self.dimensions)
        self.custom_bounds = ([-2.0, 2.0], [-2.0, 2.0])

        self.global_optimum = [0.0, 0.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2 = x
        return x1 ** 2.0 - x1 * x2 + x2 ** 2.0

# -------------------------------------------------------------------------------- #


class Salomon(Benchmark):

    """
    Salomon test objective function.

    This class defines the Salomon global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Salomon}}(\\mathbf{x}) = 1 - \\cos \\left (2 \\pi \\sqrt{\\sum_{i=1}^{n} x_i^2} \\right) + 0.1 \\sqrt{\\sum_{i=1}^n x_i^2}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    .. figure:: figures/Salomon.png
        :alt: Salomon function
        :align: center

        **Two-dimensional Salomon function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self.bounds = zip([-100.0] * self.dimensions,
                          [100.0] * self.dimensions)
        self.custom_bounds = [(-50, 50), (-50, 50)]

        self.global_optimum = [0.0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return 1.0 - cos(2.0 * pi * sqrt(sum(x ** 2.0))) + 0.1 * sqrt(sum(x ** 2.0))

# -------------------------------------------------------------------------------- #


class Sargan(Benchmark):

    """
    Sargan test objective function.

    This class defines the Sargan global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Sargan}}(\\mathbf{x}) = \\sum_{i=1}^{n} n \\left (x_i^2 + 0.4 \\sum_{i \\neq j}^{n} x_ix_j \\right)

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    .. figure:: figures/Sargan.png
        :alt: Sargan function
        :align: center

        **Two-dimensional Sargan function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self.bounds = zip([-100.0] * self.dimensions,
                          [100.0] * self.dimensions)
        self.custom_bounds = [(-5, 5), (-5, 5)]

        self.global_optimum = [0.0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x0 = x[:-1]
        x1 = roll(x, -1)[:-1]

        return sum(self.dimensions * (x ** 2 + 0.4 * sum(x0 * x1)))

# -------------------------------------------------------------------------------- #


class Schaffer01(Benchmark):

    """
    Schaffer 1 test objective function.

    This class defines the Schaffer 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schaffer01}}(\\mathbf{x}) = 0.5 + \\frac{\\sin^2 (x_1^2 + x_2^2)^2 - 0.5}{1 + 0.001(x_1^2 + x_2^2)^2}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    .. figure:: figures/Schaffer01.png
        :alt: Schaffer 1 function
        :align: center

        **Two-dimensional Schaffer 1 function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-100.0] * self.dimensions,
                          [100.0] * self.dimensions)
        self.custom_bounds = [(-10, 10), (-10, 10)]

        self.global_optimum = [0.0 for _ in range(self.dimensions)]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2 = x
        return 0.5 + (sin(x1 ** 2.0 + x2 ** 2.0) ** 2.0 - 0.5) / (1 + 0.001 * (x1 ** 2.0 + x2 ** 2.0) ** 2.0)

# -------------------------------------------------------------------------------- #


class Schaffer02(Benchmark):

    """
    Schaffer 2 test objective function.

    This class defines the Schaffer 2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schaffer02}}(\\mathbf{x}) = 0.5 + \\frac{\\sin^2 (x_1^2 - x_2^2)^2 - 0.5}{1 + 0.001(x_1^2 + x_2^2)^2}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    .. figure:: figures/Schaffer02.png
        :alt: Schaffer 2 function
        :align: center

        **Two-dimensional Schaffer 2 function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-100.0] * self.dimensions,
                          [100.0] * self.dimensions)
        self.custom_bounds = [(-10, 10), (-10, 10)]

        self.global_optimum = [0.0 for _ in range(self.dimensions)]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2 = x
        return 0.5 + (sin(x1 ** 2.0 - x2 ** 2.0) ** 2.0 - 0.5) / (1 + 0.001 * (x1 ** 2.0 + x2 ** 2.0) ** 2.0)

# -------------------------------------------------------------------------------- #


class Schaffer03(Benchmark):

    """
    Schaffer 3 test objective function.

    This class defines the Schaffer 3 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schaffer03}}(\\mathbf{x}) = 0.5 + \\frac{\\sin^2 \\left( \\cos \\lvert x_1^2 - x_2^2 \\rvert \\right ) - 0.5}{1 + 0.001(x_1^2 + x_2^2)^2}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    .. figure:: figures/Schaffer03.png
        :alt: Schaffer 3 function
        :align: center

        **Two-dimensional Schaffer 3 function**

    *Global optimum*: :math:`f(x_i) = 0.00156685` for :math:`\\mathbf{x} = [0, 1.253115]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-100.0] * self.dimensions,
                          [100.0] * self.dimensions)
        self.custom_bounds = [(-10, 10), (-10, 10)]

        self.global_optimum = [0.0, 1.253115]
        self.fglob = 0.00156685

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2 = x

        num = sin(cos(abs(x[0] ** 2 - x[1] ** 2))) ** 2 - 0.5
        den = (1 + 0.001 * (x[0] ** 2 + x[1] ** 2)) ** 2
        return 0.5 + num / den

# -------------------------------------------------------------------------------- #


class Schaffer04(Benchmark):

    """
    Schaffer 4 test objective function.

    This class defines the Schaffer 4 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schaffer04}}(\\mathbf{x}) = 0.5 + \\frac{\\cos^2 \\left( \\sin(x_1^2 - x_2^2) \\right ) - 0.5}{1 + 0.001(x_1^2 + x_2^2)^2}^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    .. figure:: figures/Schaffer04.png
        :alt: Schaffer 4 function
        :align: center

        **Two-dimensional Schaffer 4 function**

    *Global optimum*: :math:`f(x_i) = 0.292579` for :math:`\\mathbf{x} = [0, 1.253115]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-100.0] * self.dimensions,
                          [100.0] * self.dimensions)
        self.custom_bounds = [(-10, 10), (-10, 10)]

        self.global_optimum = [0.0, 1.253115]
        self.fglob = 0.292579

    def evaluator(self, x, *args):

        self.fun_evals += 1

        num = cos(sin(abs(x[0] ** 2 - x[1] ** 2))) ** 2 - 0.5
        den = (1 + 0.001 * (x[0] ** 2 + x[1] ** 2)) ** 2
        return 0.5 + num / den

# -------------------------------------------------------------------------------- #


class SchmidtVetters(Benchmark):

    """
    Schmidt-Vetters test objective function.

    This class defines the Schmidt-Vetters global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{SchmidtVetters}}(\\mathbf{x}) = \\frac{1}{1 + (x_1 - x_2)^2} + \\sin \\left(\\frac{\\pi x_2 + x_3}{2} \\right) + e^{\\left(\\frac{x_1+x_2}{x_2} - 2\\right)^2}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 10]` for :math:`i=1,2,3`.

    *Global optimum*: :math:`f(x_i) = 2.99643266` for :math:`x_i = [0.79876108,  0.79962581,  0.79848824]`

    """

    def __init__(self, dimensions=3):
        Benchmark.__init__(self, dimensions)
        self.bounds = zip([0.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [0.79876108,  0.79962581,  0.79848824]
        self.fglob = 2.99643266

    def evaluator(self, x, *args):

        self.fun_evals += 1

        return (1 / (1 + (x[0] - x[1]) ** 2) + sin((pi * x[1] + x[2]) / 2)
                + exp(((x[0] + x[1]) / x[1] - 2) ** 2))


# -------------------------------------------------------------------------------- #

class Schwefel01(Benchmark):

    """
    Schwefel 1 test objective function.

    This class defines the Schwefel 1 global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel01}}(\\mathbf{x}) = \\left(\\sum_{i=1}^n x_i^2 \\right)^{\\alpha}

    Where, in this exercise, :math:`\\alpha = \\sqrt{\\pi}`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    .. figure:: figures/Schwefel01.png
        :alt: Schwefel 1 function
        :align: center

        **Two-dimensional Schwefel 1 function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self.bounds = zip([-100.0] * self.dimensions,
                          [100.0] * self.dimensions)
        self.custom_bounds = ([-4.0, 4.0], [-4.0, 4.0])

        self.global_optimum = [0.0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        alpha = sqrt(pi)
        return (sum(x ** 2.0)) ** alpha

# -------------------------------------------------------------------------------- #


class Schwefel02(Benchmark):

    """
    Schwefel 2 test objective function.

    This class defines the Schwefel 2 global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel02}}(\\mathbf{x}) = \\sum_{i=1}^n \\left(\\sum_{j=1}^i x_i \\right)^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    .. figure:: figures/Schwefel02.png
        :alt: Schwefel 2 function
        :align: center

        **Two-dimensional Schwefel 2 function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self.bounds = zip([-100.0] * self.dimensions,
                          [100.0] * self.dimensions)
        self.custom_bounds = ([-4.0, 4.0], [-4.0, 4.0])

        self.global_optimum = [0.0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        s = 0.0

        for i in range(self.dimensions):
            temp = 0.0
            for j in range(i):
                temp += x[j]
            s += temp ** 2.0

        return s

# -------------------------------------------------------------------------------- #


class Schwefel04(Benchmark):

    """
    Schwefel 4 test objective function.

    This class defines the Schwefel 4 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel04}}(\\mathbf{x}) = \\sum_{i=1}^n \\left[(x_i - 1)^2 + (x_1 - x_i^2)^2 \\right]


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 10]` for :math:`i=1,...,n`.

    .. figure:: figures/Schwefel04.png
        :alt: Schwefel 4 function
        :align: center

        **Two-dimensional Schwefel 4 function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 1` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self.bounds = zip([0.0] * self.dimensions, [10.0] * self.dimensions)
        self.custom_bounds = ([0.0, 2.0], [0.0, 2.0])

        self.global_optimum = [1.0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        return sum((x - 1.0) ** 2.0 + (x[0] - x ** 2.0) ** 2.0)

# -------------------------------------------------------------------------------- #


class Schwefel06(Benchmark):

    """
    Schwefel 6 test objective function.

    This class defines the Schwefel 6 global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel06}}(\\mathbf{x}) = \\max(\\lvert x_1 + 2x_2 - 7 \\rvert, \\lvert 2x_1 + x_2 - 5 \\rvert)


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    .. figure:: figures/Schwefel06.png
        :alt: Schwefel 6 function
        :align: center

        **Two-dimensional Schwefel 6 function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [1, 3]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self.bounds = zip([-100.0] * self.dimensions,
                          [100.0] * self.dimensions)
        self.custom_bounds = ([-10.0, 10.0], [-10.0, 10.0])

        self.global_optimum = [1.0, 3.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2 = x
        vector = [abs(x1 + 2 * x2 - 7), abs(2 * x1 + x2 - 5)]
        return max(vector)

# -------------------------------------------------------------------------------- #


class Schwefel20(Benchmark):

    """
    Schwefel 20 test objective function.

    This class defines the Schwefel 20 global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel20}}(\\mathbf{x}) = \\sum_{i=1}^n \\lvert x_i \\rvert


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    .. figure:: figures/Schwefel20.png
        :alt: Schwefel 20 function
        :align: center

        **Two-dimensional Schwefel 20 function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self.bounds = zip([-100.0] * self.dimensions,
                          [100.0] * self.dimensions)

        self.global_optimum = [0.0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        return sum(abs(x))

# -------------------------------------------------------------------------------- #


class Schwefel21(Benchmark):

    """
    Schwefel 21 test objective function.

    This class defines the Schwefel 21 global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel21}}(\\mathbf{x}) = \\smash{\\displaystyle\\max_{1 \leq i \leq n}} \\lvert x_i \\rvert


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    .. figure:: figures/Schwefel21.png
        :alt: Schwefel 21 function
        :align: center

        **Two-dimensional Schwefel 21 function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self.bounds = zip([-100.0] * self.dimensions,
                          [100.0] * self.dimensions)

        self.global_optimum = [0.0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        return max(abs(x))

# -------------------------------------------------------------------------------- #


class Schwefel22(Benchmark):

    """
    Schwefel 22 test objective function.

    This class defines the Schwefel 22 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel22}}(\\mathbf{x}) = \\sum_{i=1}^n \\lvert x_i \\rvert + \\prod_{i=1}^n \\lvert x_i \\rvert


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    .. figure:: figures/Schwefel22.png
        :alt: Schwefel 22 function
        :align: center

        **Two-dimensional Schwefel 22 function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self.bounds = zip([-100.0] * self.dimensions,
                          [100.0] * self.dimensions)
        self.custom_bounds = ([-10.0, 10.0], [-10.0, 10.0])

        self.global_optimum = [0.0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        return sum(abs(x)) + prod(abs(x))

# -------------------------------------------------------------------------------- #


class Schwefel26(Benchmark):

    """
    Schwefel 26 test objective function.

    This class defines the Schwefel 26 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel26}}(\\mathbf{x}) = 418.9829n - \\sum_{i=1}^n x_i \\sin(\\sqrt{|x_i|})

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-500, 500]` for :math:`i=1,...,n`.

    .. figure:: figures/Schwefel26.png
        :alt: Schwefel 26 function
        :align: center

        **Two-dimensional Schwefel 26 function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 420.968746` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self.bounds = zip([-500.0] * self.dimensions,
                          [500.0] * self.dimensions)

        self.global_optimum = [420.968746 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return 418.982887 * self.dimensions - sum([x * sin(sqrt(abs(x)))])

# -------------------------------------------------------------------------------- #


class Schwefel36(Benchmark):

    """
    Schwefel 36 test objective function.

    This class defines the Schwefel 36 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel36}}(\\mathbf{x}) = -x_1x_2(72 - 2x_1 - 2x_2)


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 500]` for :math:`i=1,2`.

    .. figure:: figures/Schwefel36.png
        :alt: Schwefel 36 function
        :align: center

        **Two-dimensional Schwefel 36 function**

    *Global optimum*: :math:`f(x_i) = -3456` for :math:`\\mathbf{x} = [12, 12]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self.bounds = zip([0.0] * self.dimensions, [500.0] * self.dimensions)
        self.custom_bounds = ([0.0, 20.0], [0.0, 20.0])

        self.global_optimum = [12.0, 12.0]
        self.fglob = -3456.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2 = x
        return -x1 * x2 * (72.0 - 2.0 * x1 - 2.0 * x2)

# -------------------------------------------------------------------------------- #


class Shekel05(Benchmark):

    """
    Shekel 5 test objective function.

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

        self.bounds = zip([0.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [4.0 for _ in range(self.dimensions)]
        self.fglob = -10.15319585

    def evaluator(self, x, *args):

        self.fun_evals += 1
        m = 5

        A = asarray([[4.0, 4.0, 4.0, 4.0],
                     [1.0, 1.0, 1.0, 1.0],
                     [8.0, 8.0, 8.0, 8.0],
                     [6.0, 6.0, 6.0, 6.0],
                     [3.0, 7.0, 3.0, 7.0]])

        C = asarray([0.1, 0.2, 0.2, 0.4, 0.4])

        return -sum(1.0 / (dot(x - a, x - a) + c) for a, c in zip(A, C))

# -------------------------------------------------------------------------------- #


class Shekel07(Benchmark):

    """
    Shekel 7 test objective function.

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

        self.bounds = zip([0.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [4.0 for _ in range(self.dimensions)]
        self.fglob = -10.4028188

    def evaluator(self, x, *args):

        self.fun_evals += 1
        m = 7

        A = asarray([[4.0, 4.0, 4.0, 4.0],
                     [1.0, 1.0, 1.0, 1.0],
                     [8.0, 8.0, 8.0, 8.0],
                     [6.0, 6.0, 6.0, 6.0],
                     [3.0, 7.0, 3.0, 7.0],
                     [2.0, 9.0, 2.0, 9.0],
                     [5.0, 5.0, 3.0, 3.0]])

        C = asarray([0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3])

        return -sum(1.0 / (dot(x - a, x - a) + c) for a, c in zip(A, C))

# -------------------------------------------------------------------------------- #


class Shekel10(Benchmark):

    """
    Shekel 10 test objective function.

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

        self.bounds = zip([0.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [4.0 for _ in range(self.dimensions)]
        self.fglob = -10.5362837262

    def evaluator(self, x, *args):

        self.fun_evals += 1
        m = 10

        A = asarray([[4.0, 4.0, 4.0, 4.0],
                     [1.0, 1.0, 1.0, 1.0],
                     [8.0, 8.0, 8.0, 8.0],
                     [6.0, 6.0, 6.0, 6.0],
                     [3.0, 7.0, 3.0, 7.0],
                     [2.0, 9.0, 2.0, 9.0],
                     [5.0, 5.0, 3.0, 3.0],
                     [8.0, 1.0, 8.0, 1.0],
                     [6.0, 2.0, 6.0, 2.0],
                     [7.0, 3.6, 7.0, 3.6]])

        C = asarray([0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3, 0.7, 0.5, 0.5])

        return -sum(1.0 / (dot(x - a, x - a) + c) for a, c in zip(A, C))

# -------------------------------------------------------------------------------- #


class Shubert01(Benchmark):

    """
    Shubert 1 test objective function.

    This class defines the Shubert 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Shubert01}}(\\mathbf{x}) = \\left( \\sum\\limits_{i=1}^{5} i\\cos[(i+1)x_1 + i] \\right) \\left( \\sum\\limits_{i=1}^{5} i\\cos[(i+1)x_2 + i] \\right)

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    .. figure:: figures/Shubert01.png
        :alt: Shubert 1 function
        :align: center

        **Two-dimensional Shubert 1 function**

    *Global optimum*: :math:`f(x_i) = -186.7309` for :math:`\\mathbf{x} = [-7.0835, 4.8580]` (and many others).

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [-7.0835, 4.8580]
        self.fglob = -186.7309

    def evaluator(self, x, *args):

        self.fun_evals += 1

        s1 = s2 = 0.0
        for i in range(1, 6):
            s1 = s1 + i * cos((i + 1) * x[0] + i)
            s2 = s2 + i * cos((i + 1) * x[1] + i)

        y = s1 * s2
        return y

# -------------------------------------------------------------------------------- #


class Shubert03(Benchmark):

    """
    Shubert 3 test objective function.

    This class defines the Shubert 3 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Shubert03}}(\\mathbf{x}) = \\sum_{i=1}^n \\sum_{j=1}^5 j \\sin \\left[(j+1)x_i \\right] + j

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    .. figure:: figures/Shubert03.png
        :alt: Shubert 3 function
        :align: center

        **Two-dimensional Shubert 3 function**

    *Global optimum*: :math:`f(x_i) = -24.062499` for :math:`\\mathbf{x} = [5.791794, 5.791794]` (and many others).

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [5.791794, 5.791794]
        self.fglob = -24.062499

    def evaluator(self, x, *args):

        self.fun_evals += 1

        return -sin(2.0 * x[0] + 1.0) - 2.0 * sin(3.0 * x[0] + 2.0) - 3.0 * sin(4.0 * x[0] + 3.0) - 4.0 * sin(5.0 * x[0] + 4.0) \
               - 5.0 * sin(6.0 * x[0] + 5.0) - sin(2.0 * x[1] + 1.0) - 2.0 * sin(3.0 * x[1] + 2.0) - 3.0 * sin(4.0 * x[1] + 3.0) \
               - 4.0 * sin(5.0 * x[1] + 4.0) - 5.0 * sin(6.0 * x[1] + 5.0)

# -------------------------------------------------------------------------------- #


class Shubert04(Benchmark):

    """
    Shubert 4 test objective function.

    This class defines the Shubert 4 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Shubert04}}(\\mathbf{x}) = \\sum_{i=1}^n \\sum_{j=1}^5 j \\cos \\left[(j+1)x_i \\right] + j

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    .. figure:: figures/Shubert04.png
        :alt: Shubert 4 function
        :align: center

        **Two-dimensional Shubert 4 function**

    *Global optimum*: :math:`f(x_i) = -29.016015` for :math:`\\mathbf{x} = [-0.80032121, -7.08350592]` (and many others).

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [-0.80032121, -7.08350592]
        self.fglob = -29.016015

    def evaluator(self, x, *args):

        self.fun_evals += 1

        return -cos(2.0 * x[0] + 1.0) - 2.0 * cos(3.0 * x[0] + 2.0) - 3.0 * cos(4.0 * x[0] + 3.0) - 4.0 * cos(5.0 * x[0] + 4.0) \
               - 5.0 * cos(6.0 * x[0] + 5.0) - cos(2.0 * x[1] + 1.0) - 2.0 * cos(3.0 * x[1] + 2.0) - 3.0 * cos(4.0 * x[1] + 3.0) \
               - 4.0 * cos(5.0 * x[1] + 4.0) - 5.0 * cos(6.0 * x[1] + 5.0)

# -------------------------------------------------------------------------------- #


class SineEnvelope(Benchmark):

    """
    SineEnvelope test objective function.

    This class defines the SineEnvelope global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{SineEnvelope}}(\\mathbf{x}) = -\\sum_{i=1}^{n-1}\\left[\\frac{\\sin^2(\\sqrt{x_{i+1}^2+x_{i}^2}-0.5)}{(0.001(x_{i+1}^2+x_{i}^2)+1)^2}+0.5\\right]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    .. figure:: figures/SineEnvelope.png
        :alt: SineEnvelope function
        :align: center

        **Two-dimensional SineEnvelope function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-100.0] * self.dimensions,
                          [100.0] * self.dimensions)
        self.custom_bounds = [(-20, 20), (-20, 20)]

        self.global_optimum = [0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        X1 = x[0:-1]
        X2 = x[1:]
        X12X22 = X1 ** 2 + X2 ** 2
        return sum((sin(sqrt(X12X22)) ** 2 - 0.5) / (1 + 0.001 * X12X22) ** 2 + 0.5)

# -------------------------------------------------------------------------------- #


class SixHumpCamel(Benchmark):

    """
    Six Hump Camel test objective function.

    This class defines the Six Hump Camel global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{SixHumpCamel}}(\\mathbf{x}) = 4x_1^2+x_1x_2-4x_2^2-2.1x_1^4+4x_2^4+\\frac{1}{3}x_1^6

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,2`.

    .. figure:: figures/SixHumpCamel.png
        :alt: Six Hump Camel function
        :align: center

        **Two-dimensional Six Hump Camel function**

    *Global optimum*: :math:`f(x_i) = -1.031628453489877` for :math:`\\mathbf{x} = [0.08984201368301331 , -0.7126564032704135]`
    or :math:`\\mathbf{x} = [-0.08984201368301331, 0.7126564032704135]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-5.0] * self.dimensions, [5.0] * self.dimensions)
        self.custom_bounds = [(-2, 2), (-1.5, 1.5)]

        self.global_optimum = [(0.08984201368301331, -0.7126564032704135),
                               (-0.08984201368301331, 0.7126564032704135)]
        self.fglob = -1.031628

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return (4 - 2.1 * x[0] ** 2 + x[0] ** 4 / 3) * x[0] ** 2 + x[0] * x[1] + (4 * x[1] ** 2 - 4) * x[1] ** 2

# -------------------------------------------------------------------------------- #


class Sodp(Benchmark):

    """
    Sodp test objective function.

    This class defines the Sum Of Different Powers global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Sodp}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\lvert{x_{i}}\\rvert^{i + 1}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,...,n`.

    .. figure:: figures/Sodp.png
        :alt: Sodp function
        :align: center

        **Two-dimensional Sum Of Different Powers function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-1.0] * self.dimensions, [1.0] * self.dimensions)

        self.global_optimum = [0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        i = arange(1, self.dimensions + 1)
        return sum(abs(x) ** (i + 1))

# -------------------------------------------------------------------------------- #


class Sphere(Benchmark):

    """
    Sphere test objective function.

    This class defines the Sphere global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Sphere}}(\\mathbf{x}) = \\sum_{i=1}^{n} x_i^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,...,n`.

    .. figure:: figures/Sphere.png
        :alt: Sphere function
        :align: center

        **Two-dimensional Sphere function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self.bounds = zip([-5.12] * self.dimensions, [5.12] * self.dimensions)

        self.global_optimum = [0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return sum(x ** 2)

# -------------------------------------------------------------------------------- #


class Step(Benchmark):

    """
    Step test objective function.

    This class defines the Step global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Step}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\left ( \\lfloor x_i  + 0.5 \\rfloor \\right )^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    .. figure:: figures/Step.png
        :alt: Step function
        :align: center

        **Two-dimensional Step function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0.5` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self.bounds = zip([-100.0] * self.dimensions,
                          [100.0] * self.dimensions)
        self.custom_bounds = ([-5, 5], [-5, 5])

        self.global_optimum = [0.0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return sum((floor(x + 0.5)) ** 2.0)

# -------------------------------------------------------------------------------- #


class Stochastic(Benchmark):

    """
    Stochastic test objective function.

    This class defines a Stochastic global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Stochastic}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\epsilon_i \\left | {x_i - \\frac{1}{i}} \\right |

    The variable :math:`\\epsilon_i, (i=1,...,n)` is a random variable uniformly distributed in :math:`[0, 1]`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,...,n`.

    .. figure:: figures/Stochastic.png
        :alt: Stochastic function
        :align: center

        **Two-dimensional Stochastic function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = [1/n]` for :math:`i=1,...,n`
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-5.0] * self.dimensions, [5.0] * self.dimensions)

        self.global_optimum = [1.0 / _ for _ in range(1, self.dimensions + 1)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        rnd = uniform(0.0, 1.0, size=(self.dimensions, ))
        rng = arange(1, self.dimensions + 1)

        return sum(rnd * abs(x - 1.0 / rng))

# -------------------------------------------------------------------------------- #


class StretchedV(Benchmark):

    """
    StretchedV test objective function.

    This class defines the Stretched V global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{StretchedV}}(\\mathbf{x}) = \sum_{i=1}^{n-1} t^{1/4} [\sin (50t^{0.1}) + 1]^2

    Where, in this exercise:

    .. math::

       t = x_{i+1}^2 + x_i^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    .. figure:: figures/StretchedV.png
        :alt: StretchedV function
        :align: center

        **Two-dimensional StretchedV function**


    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [-9.38723188, 9.34026753]` when :math:`n = 2`.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10] * self.dimensions, [10] * self.dimensions)

        self.global_optimum = [-9.38723188, 9.34026753]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        s = 0.0
        for i in range(self.dimensions - 1):
            t = x[i + 1] * x[i + 1] + x[i] * x[i]
            s += t ** 0.25 * (sin(50.0 * t ** 0.1 + 1.0)) ** 2.0

        return s

# -------------------------------------------------------------------------------- #


class StyblinskiTang(Benchmark):

    """
    StyblinskiTang test objective function.

    This class defines the Styblinski-Tang global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{StyblinskiTang}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\left(x_i^4 - 16x_i^2 + 5x_i \\right)

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,...,n`.

    .. figure:: figures/StyblinskiTang.png
        :alt: StyblinskiTang function
        :align: center

        **Two-dimensional Styblinski-Tang function**

    *Global optimum*: :math:`f(x_i) = -39.16616570377142n` for :math:`x_i = -2.903534018185960` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-5.0] * self.dimensions, [5.0] * self.dimensions)

        self.global_optimum = [
            -2.903534018185960 for _ in range(self.dimensions)]
        self.fglob = -39.16616570377142 * self.dimensions
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return sum(x ** 4 - 16 * x ** 2 + 5 * x) / 2

# -------------------------------------------------------------------------------- #


class TestTubeHolder(Benchmark):

    """
    TestTubeHolder test objective function.

    This class defines the TestTubeHolder global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{TestTubeHolder}}(\\mathbf{x}) = - 4 \\left | {e^{\\left|{\\cos\\left(\\frac{1}{200} x_{1}^{2} + \\frac{1}{200} x_{2}^{2}\\right)}\\right|} \\sin\\left(x_{1}\\right) \\cos\\left(x_{2}\\right)}\\right |

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    .. figure:: figures/TestTubeHolder.png
        :alt: TestTubeHolder function
        :align: center

        **Two-dimensional TestTubeHolder function**

    *Global optimum*: :math:`f(x_i) = -10.872299901558` for :math:`\\mathbf{x} = [-\\pi/2, 0]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [-pi / 2, 0.0]
        self.fglob = -10.87229990155800

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return -4 * abs(sin(x[0]) * cos(x[1]) * exp(abs(cos((x[0] ** 2 + x[1] ** 2) / 200))))

# -------------------------------------------------------------------------------- #


class Treccani(Benchmark):

    """
    Treccani test objective function.

    This class defines the Treccani global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Treccani}}(\\mathbf{x}) = x_1^4 + 4x_1^3 + 4x_1^2 + x_2^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,2`.

    .. figure:: figures/Treccani.png
        :alt: Treccani function
        :align: center

        **Two-dimensional Treccani function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [-2, 0]` or :math:`\\mathbf{x} = [0, 0]`.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-5.0] * self.dimensions, [5.0] * self.dimensions)
        self.custom_bounds = [(-2, 2), (-2, 2)]

        self.global_optimum = [-2.0, 0.0]
        self.fglob = 0

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return x[0] ** 4 + 4.0 * x[0] ** 3 + 4.0 * x[0] ** 2 + x[1] ** 2

# -------------------------------------------------------------------------------- #


class Trefethen(Benchmark):

    """
    Trefethen test objective function.

    This class defines the Trefethen global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Trefethen}}(\\mathbf{x}) = 0.25 x_{1}^{2} + 0.25 x_{2}^{2} + e^{\\sin\\left(50 x_{1}\\right)} - \\sin\\left(10 x_{1} + 10 x_{2}\\right) + \\sin\\left(60 e^{x_{2}}\\right) + \\sin\\left[70 \\sin\\left(x_{1}\\right)\\right] + \\sin\\left[\\sin\\left(80 x_{2}\\right)\\right]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    .. figure:: figures/Trefethen.png
        :alt: Trefethen function
        :align: center

        **Two-dimensional Trefethen function**

    *Global optimum*: :math:`f(x_i) = -3.3068686474` for :math:`\\mathbf{x} = [-0.02440307923, 0.2106124261]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)
        self.custom_bounds = [(-5, 5), (-5, 5)]

        self.global_optimum = [-0.02440307923, 0.2106124261]
        self.fglob = -3.3068686474

    def evaluator(self, x, *args):

        self.fun_evals += 1
        F = exp(sin(50 * x[0])) + sin(60 * exp(x[1])) + sin(70 * sin(x[0])) + \
            sin(sin(80 * x[1])) - sin(10 * (x[0] + x[1])) + \
            1.0 / 4 * (x[0] ** 2 + x[1] ** 2)

        return F

# -------------------------------------------------------------------------------- #


class ThreeHumpCamel(Benchmark):

    """
    Three Hump Camel test objective function.

    This class defines the Three Hump Camel global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{ThreeHumpCamel}}(\\mathbf{x}) = 2x_1^2 - 1.05x_1^4 + \\frac{x_1^6}{6} + x_1x_2 + x_2^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,2`.

    .. figure:: figures/ThreeHumpCamel.png
        :alt: Three Hump Camel function
        :align: center

        **Two-dimensional Three Hump Camel function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [0, 0]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-5.0] * self.dimensions, [5.0] * self.dimensions)
        self.custom_bounds = [(-2, 2), (-1.5, 1.5)]

        self.global_optimum = [0.0, 0.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return (2.0 * x[0] ** 2.0 - 1.05 * x[0] ** 4.0 + x[0] ** 6 / 6.0
                + x[0] * x[1] + x[1] ** 2.0)

# -------------------------------------------------------------------------------- #


class Trid(Benchmark):

    """
    Trid test objective function.

    This class defines the Trid global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Trid}}(\\mathbf{x}) = \\sum_{i=1}^{n}(x_i - 1)^2 - \\sum_{i=2}^{n} x_ix_{i-1}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-20, 20]` for :math:`i=1,...,6`.

    *Global optimum*: :math:`f(x_i) = -50` for :math:`\\mathbf{x} = [6, 10, 12, 12, 10, 6]`

    """

    def __init__(self, dimensions=6):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([0.0] * self.dimensions, [20.0] * self.dimensions)

        self.global_optimum = [6, 10, 12, 12, 10, 6]
        self.fglob = -50.0

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return sum((x - 1.0) ** 2.0) - sum(x[1:] * x[0:-1])

# -------------------------------------------------------------------------------- #


class Trigonometric01(Benchmark):

    """
    Trigonometric 1 test objective function.

    This class defines the Trigonometric 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Trigonometric01}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\left [n - \\sum_{j=1}^{n} \\cos(x_j) + i \\left(1 - cos(x_i) - sin(x_i) \\right ) \\right]^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, \\pi]` for :math:`i=1,...,n`.

    .. figure:: figures/Trigonometric01.png
        :alt: Trigonometric 1 function
        :align: center

        **Two-dimensional Trigonometric 1 function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([0.0] * self.dimensions, [pi] * self.dimensions)

        self.global_optimum = [0.0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        rng = arange(1.0, self.dimensions + 1)
        return sum((self.dimensions - sum(cos(x) + rng * (1 - cos(x) - sin(x)))) ** 2.0)

# -------------------------------------------------------------------------------- #


class Trigonometric02(Benchmark):

    """
    Trigonometric 2 test objective function.

    This class defines the Trigonometric 2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Trigonometric2}}(\\mathbf{x}) = 1 + \\sum_{i=1}^{n} 8 \\sin^2 \\left[7(x_i - 0.9)^2 \\right] + 6 \\sin^2 \\left[14(x_i - 0.9)^2 \\right] + (x_i - 0.9)^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-500, 500]` for :math:`i=1,...,n`.

    .. figure:: figures/Trigonometric02.png
        :alt: Trigonometric 2 function
        :align: center

        **Two-dimensional Trigonometric 2 function**

    *Global optimum*: :math:`f(x_i) = 1` for :math:`x_i = 0.9` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-500.0] * self.dimensions,
                          [500.0] * self.dimensions)
        self.custom_bounds = [(0, 2), (0, 2)]

        self.global_optimum = [0.9 for _ in range(self.dimensions)]
        self.fglob = 1.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        return 1.0 + sum(8.0 * (sin(7.0 * (x - 0.9) ** 2.0) ** 2.0) + 6.0 * (sin(14.0 * (x - 0.9) ** 2.0) ** 2.0) + (x - 0.9) ** 2.0)

# -------------------------------------------------------------------------------- #


class Tripod(Benchmark):

    """
    Tripod test objective function.

    This class defines the Tripod global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Tripod}}(\\mathbf{x}) = p(x_2) \\left[1 + p(x_1) \\right] + \\lvert x_1 + 50p(x_2) \\left[1 - 2p(x_1) \\right] \\rvert + \\lvert x_2 + 50\\left[1 - 2p(x_2)\\right] \\rvert

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    .. figure:: figures/Tripod.png
        :alt: Tripod function
        :align: center

        **Two-dimensional Tripod function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [0, -50]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-100.0] * self.dimensions,
                          [100.0] * self.dimensions)

        self.global_optimum = [0.0, -50.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        p1 = float(x[0] >= 0)
        p2 = float(x[1] >= 0)

        return (p2 * (1.0 + p1) + abs(x[0] + 50.0 * p2 * (1.0 - 2.0 * p1))
                + abs(x[1] + 50.0 * (1.0 - 2.0 * p2)))

# -------------------------------------------------------------------------------- #


class Ursem01(Benchmark):

    """
    Ursem 1 test objective function.

    This class defines the Ursem 1 global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Ursem01}}(\\mathbf{x}) = - \\sin(2x_1 - 0.5 \\pi) - 3 \\cos(x_2) - 0.5x_1

    Here, :math:`n` represents the number of dimensions and :math:`x_1 \\in [-2.5, 3]`, :math:`x_2 \\in [-2, 2]`.

    .. figure:: figures/Ursem01.png
        :alt: Ursem 1 function
        :align: center

        **Two-dimensional Ursem 1 function**

    *Global optimum*: :math:`f(x_i) = -4.8168` for :math:`\\mathbf{x} = [1.69714, 0.0]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(-2.5, 3.0), (-2.0, 2.0)]

        self.global_optimum = [1.69714, 0.0]
        self.fglob = -4.8168

    def evaluator(self, x, *args):
        self.fun_evals += 1
        return (-sin(2 * x[0] - 0.5 * pi) - 3.0 * cos(x[1]) - 0.5 * x[0])

# -------------------------------------------------------------------------------- #


class Ursem03(Benchmark):

    """
    Ursem 3 test objective function.

    This class defines the Ursem 3 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Ursem03}}(\\mathbf{x}) = - \\sin(2.2 \\pi x_1 + 0.5 \\pi) \\frac{2 - \\lvert x_1 \\rvert}{2} \\frac{3 - \\lvert x_1 \\rvert}{2} - \\sin(2.2 \\pi x_2 + 0.5 \\pi) \\frac{2 - \\lvert x_2 \\rvert}{2} \\frac{3 - \\lvert x_2 \\rvert}{2}

    Here, :math:`n` represents the number of dimensions and :math:`x_1 \\in [-2, 2]`, :math:`x_2 \\in [-1.5, 1.5]`.

    .. figure:: figures/Ursem03.png
        :alt: Ursem 3 function
        :align: center

        **Two-dimensional Ursem 3 function**

    *Global optimum*: :math:`f(x_i) = -3` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(-2, 2), (-1.5, 1.5)]

        self.global_optimum = [0.0 for _ in range(self.dimensions)]
        self.fglob = -3.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        return (-sin(2.2 * pi * x[0] + 0.5 * pi)
                * ((2.0 - abs(x[0])) / 2.0) * ((3.0 - abs(x[0])) / 2)
                - sin(2.2 * pi * x[1] + 0.5 * pi)
                * ((2.0 - abs(x[1])) / 2) * ((3.0 - abs(x[1])) / 2))

# -------------------------------------------------------------------------------- #


class Ursem04(Benchmark):

    """
    Ursem 4 test objective function.

    This class defines the Ursem 4 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Ursem04}}(\\mathbf{x}) = -3 \\sin(0.5 \\pi x_1 + 0.5 \\pi) \\frac{2 - \\sqrt{x_1^2 + x_2 ^ 2}}{4}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-2, 2]` for :math:`i=1,2`.

    .. figure:: figures/Ursem04.png
        :alt: Ursem 4 function
        :align: center

        **Two-dimensional Ursem 4 function**

    *Global optimum*: :math:`f(x_i) = -1.5` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-2.0] * self.dimensions, [2.0] * self.dimensions)

        self.global_optimum = [0.0 for _ in range(self.dimensions)]
        self.fglob = -1.5

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2 = x
        return -3 * sin(0.5 * pi * x1 + 0.5 * pi) * (2.0 - sqrt(x1 ** 2.0 + x2 ** 2.0)) / 4.0

# -------------------------------------------------------------------------------- #


class UrsemWaves(Benchmark):

    """
    Ursem Waves test objective function.

    This class defines the Ursem Waves global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{UrsemWaves}}(\\mathbf{x}) = -0.9x_1^2 + (x_2^2 - 4.5x_2^2)x_1x_2 + 4.7 \\cos \\left[ 2x_1 - x_2^2(2 + x_1) \\right ] \\sin(2.5 \\pi x_1)

    Here, :math:`n` represents the number of dimensions and :math:`x_1 \\in [-0.9, 1.2]`, :math:`x_2 \\in [-1.2, 1.2]`.

    .. figure:: figures/UrsemWaves.png
        :alt: Ursem Waves function
        :align: center

        **Two-dimensional Ursem Waves function**

    *Global optimum*: :math:`f(x_i) = -8.5536` for :math:`x_i = 1.2` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(-0.9, 1.2), (-1.2, 1.2)]

        self.global_optimum = [1.2 for _ in range(self.dimensions)]
        self.fglob = -8.5536

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2 = x
        return -0.9 * x1 ** 2.0 + (x2 ** 2.0 - 4.5 * x2 ** 2.0) * x1 * x2 + 4.7 * cos(3 * x1 - x2 ** 2.0 * (2.0 + x1)) * sin(2.5 * pi * x1)

# -------------------------------------------------------------------------------- #


class VenterSobiezcczanskiSobieski(Benchmark):

    """
    Venter Sobiezcczanski-Sobieski test objective function.

    This class defines the Venter Sobiezcczanski-Sobieski global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{VenterSobiezcczanskiSobieski}}(\\mathbf{x}) = x_1^2 - 100 \\cos^2(x_1) - 100 \\cos(x_1^2/30) + x_2^2 - 100 \\cos^2(x_2) - 100 \\cos(x_2^2/30)

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-50, 50]` for :math:`i=1,2`.

    .. figure:: figures/VenterSobiezcczanskiSobieski.png
        :alt: Venter Sobiezcczanski-Sobieski function
        :align: center

        **Two-dimensional Venter Sobiezcczanski-Sobieski function**

    *Global optimum*: :math:`f(x_i) = -400` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-50.0] * self.dimensions, [50.0] * self.dimensions)
        self.custom_bounds = ([-10, 10], [-10, 10])

        self.global_optimum = [0.0 for _ in range(self.dimensions)]
        self.fglob = -400

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2 = x
        return x1 ** 2.0 - 100.0 * cos(x1) ** 2.0 - 100.0 * cos(x1 ** 2.0 / 30.0) + x2 ** 2.0 - 100.0 * cos(x2) ** 2.0 - 100.0 * cos(x2 ** 2.0 / 30.0)

# -------------------------------------------------------------------------------- #


class Vincent(Benchmark):

    """
    Vincent test objective function.

    This class defines the Vincent global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Vincent}}(\\mathbf{x}) = - \\sum_{i=1}^{n} \\sin(10 \\log(x))

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0.25, 10]` for :math:`i=1,...,n`.

    .. figure:: figures/Vincent.png
        :alt: Vincent function
        :align: center

        **Two-dimensional Vincent function**

    *Global optimum*: :math:`f(x_i) = -n` for :math:`x_i = 7.70628098` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([0.25] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [7.70628098 for _ in range(self.dimensions)]
        self.fglob = -float(self.dimensions)
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        return -sum(sin(10.0 * log(x)))

# -------------------------------------------------------------------------------- #


class Watson(Benchmark):

    """
    Watson test objective function.

    This class defines the Watson global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Watson}}(\\mathbf{x}) = \\sum_{i=0}^{29} \\left\\{ \\sum_{j=0}^4 ((j - 1)a_i^j x_{j+1}) - \\left[ \\sum_{j=0}^5 a_i^j x_{j+1} \\right ]^2 - 1 \\right\\}^2 + x_1^2


    Where, in this exercise, :math:`a_i = i/29`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,...,6`.

    *Global optimum*: :math:`f(x_i) = 0.002288` for :math:`\\mathbf{x} = [-0.0158, 1.012, -0.2329, 1.260, -1.513, 0.9928]`

    """

    def __init__(self, dimensions=6):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-5.0] * self.dimensions, [5.0] * self.dimensions)

        self.global_optimum = [-0.0158, 1.012, -0.2329, 1.260, -1.513, 0.9928]
        self.fglob = 0.002288

    def evaluator(self, x, *args):

        self.fun_evals += 1

        vec = zeros((31, ))
        div = (arange(29) + 1.0) / 29.0
        s1 = 0.0
        dx = 1.0

        for j in range(1, self.dimensions):
            s1 += j * dx * x[j]
            dx *= div

        s2 = 0.0
        dx = 1.0

        for j in range(self.dimensions):
            s2 += dx * x[j]
            dx *= div

        vec[:29] = s1 - s2 ** 2.0 - 1.0
        vec[29] = x[0]
        vec[30] = x[1] - x[0] ** 2 - 1

        return sum(vec ** 2.0)

# -------------------------------------------------------------------------------- #


class Wavy(Benchmark):

    """
    W / Wavy test objective function.

    This class defines the W / Wavy global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Wavy}}(\\mathbf{x}) = 1 - \\frac{1}{n} \\sum_{i=1}^{n} \\cos(kx_i)e^{-\\frac{x_i^2}{2}}


    Where, in this exercise, :math:`k = 10`. The number of local minima is :math:`kn` and :math:`(k + 1)n` for odd and even :math:`k` respectively.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-\\pi, \\pi]` for :math:`i=1,2`.

    .. figure:: figures/Wavy.png
        :alt: W / Wavy function
        :align: center

        **Two-dimensional W / Wavy function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-pi] * self.dimensions, [pi] * self.dimensions)

        self.global_optimum = [0.0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return 1.0 - (1.0 / self.dimensions) * sum(cos(10 * x) * exp(-x ** 2.0 / 2.0))

# -------------------------------------------------------------------------------- #


class WayburnSeader01(Benchmark):

    """
    Wayburn and Seader 1 test objective function.

    This class defines the Wayburn and Seader 1 global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{WayburnSeader01}}(\\mathbf{x}) = (x_1^6 + x_2^4 - 17)^2 + (2x_1 + x_2 - 4)^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,2`.

    .. figure:: figures/WayburnSeader01.png
        :alt: Wayburn and Seader 1 function
        :align: center

        **Two-dimensional Wayburn and Seader 1 function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [1, 2]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-5.0] * self.dimensions, [5.0] * self.dimensions)
        self.custom_bounds = ([-2, 2], [-2, 2])

        self.global_optimum = [1.0, 2.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2 = x
        return (x1 ** 6.0 + x2 ** 4.0 - 17.0) ** 2.0 + (2 * x1 + x2 - 4.0) ** 2.0

# -------------------------------------------------------------------------------- #


class WayburnSeader02(Benchmark):

    """
    Wayburn and Seader 2 test objective function.

    This class defines the Wayburn and Seader 2 global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{WayburnSeader02}}(\\mathbf{x}) = \\left[ 1.613 - 4(x_1 - 0.3125)^2 - 4(x_2 - 1.625)^2 \\right]^2 + (x_2 - 1)^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-500, 500]` for :math:`i=1,2`.

    .. figure:: figures/WayburnSeader02.png
        :alt: Wayburn and Seader 2 function
        :align: center

        **Two-dimensional Wayburn and Seader 2 function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [0.2, 1]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-500.0] * self.dimensions,
                          [500.0] * self.dimensions)
        self.custom_bounds = ([-1, 2], [-1, 2])

        self.global_optimum = [0.2, 1.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2 = x
        return (1.613 - 4.0 * (x1 - 0.3125) ** 2.0 - 4.0 * (x2 - 1.625) ** 2.0) ** 2.0 + (x2 - 1.0) ** 2.0

# -------------------------------------------------------------------------------- #


class Weierstrass(Benchmark):

    """
    Weierstrass test objective function.

    This class defines the Weierstrass global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Weierstrass}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\left [ \\sum_{k=0}^{kmax} a^k \\cos \\left( 2 \\pi b^k (x_i + 0.5) \\right) - n \\sum_{k=0}^{kmax} a^k \\cos(\\pi b^k) \\right ]


    Where, in this exercise, :math:`kmax = 20`, :math:`a = 0.5` and :math:`b = 3`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-0.5, 0.5]` for :math:`i=1,...,n`.

    .. figure:: figures/Weierstrass.png
        :alt: Weierstrass function
        :align: center

        **Two-dimensional Weierstrass function**

    *Global optimum*: :math:`f(x_i) = 4` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-0.5] * self.dimensions, [0.5] * self.dimensions)

        self.global_optimum = [0.0 for _ in range(self.dimensions)]
        self.fglob = 4.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        kmax = 20
        a, b = 0.5, 3.0

        y = zeros((self.dimensions, ))

        ak = a ** (arange(0, kmax + 1))
        bk = b ** (arange(0, kmax + 1))
        for i in range(self.dimensions):
            y[i] = sum(ak * cos(2 * pi * bk * (x[i] + 0.5))) - \
                self.dimensions * sum(ak * cos(pi * bk))

        return sum(y)

# -------------------------------------------------------------------------------- #


class Whitley(Benchmark):

    """
    Whitley test objective function.

    This class defines the Whitley global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Whitley}}(\\mathbf{x}) = \\sum_{i=1}^n \\sum_{j=1}^n \\left[\\frac{(100(x_i^2-x_j)^2 + (1-x_j)^2)^2}{4000} - \\cos(100(x_i^2-x_j)^2 + (1-x_j)^2)+1 \\right]


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10.24, 10.24]` for :math:`i=1,...,n`.

    .. figure:: figures/Whitley.png
        :alt: Whitley function
        :align: center

        **Two-dimensional Whitley function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 1` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.24] * self.dimensions,
                          [10.24] * self.dimensions)
        self.custom_bounds = ([-1, 2], [-1, 2])

        self.global_optimum = [1.0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        n = self.dimensions

        s_out = 0.0
        for i in xrange(n):
            for j in range(n):
                temp = 100.0 * ((x[i] ** 2.0) - x[j]) + (1.0 - x[j]) ** 2.0
                s_out += (float(temp ** 2.0) / 4000.0) - cos(temp) + 1.0

        return s_out

# -------------------------------------------------------------------------------- #


class Wolfe(Benchmark):

    """
    Wolfe test objective function.

    This class defines the Wolfe global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Wolfe}}(\\mathbf{x}) = \\frac{4}{3}(x_1^2 + x_2^2 - x_1x_2)^{0.75} + x_3


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 2]` for :math:`i=1,2,3`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,2,3`

    """

    def __init__(self, dimensions=3):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([0.0] * self.dimensions, [2.0] * self.dimensions)

        self.global_optimum = [0.0 for _ in range(self.dimensions)]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2, x3 = x

        return (4.0 / 3.0) * (x1 ** 2.0 + x2 ** 2.0 - x1 * x2) ** 0.75 + x3

# -------------------------------------------------------------------------------- #


class XinSheYang01(Benchmark):

    """
    Xin-She Yang 1 test objective function.

    This class defines the Xin-She Yang 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{XinSheYang01}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\epsilon_i \\lvert x_i \\rvert^i


    The variable :math:`\\epsilon_i, (i=1,...,n)` is a random variable uniformly distributed in :math:`[0, 1]`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,...,n`.

    .. figure:: figures/XinSheYang01.png
        :alt: Xin-She Yang 1 function
        :align: center

        **Two-dimensional Xin-She Yang 1 function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-5.0] * self.dimensions, [5.0] * self.dimensions)
        self.custom_bounds = ([-2, 2], [-2, 2])

        self.global_optimum = [0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        epsilon = uniform(0.0, 1.0, size=self.dimensions)

        rng = arange(1.0, self.dimensions + 1.0)
        return sum(epsilon * (abs(x) ** rng))

# -------------------------------------------------------------------------------- #


class XinSheYang02(Benchmark):

    """
    Xin-She Yang 2 test objective function.

    This class defines the Xin-She Yang 2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{XinSheYang02}}(\\mathbf{x}) = \\frac{\\sum_{i=1}^{n} \\lvert{x_{i}}\\rvert}{e^{\\sum_{i=1}^{n} \\sin\\left(x_{i}^{2.0}\\right)}}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-2\\pi, 2\\pi]` for :math:`i=1,...,n`.

    .. figure:: figures/XinSheYang02.png
        :alt: Xin-She Yang 2 function
        :align: center

        **Two-dimensional Xin-She Yang 2 function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-2 * pi] * self.dimensions,
                          [2 * pi] * self.dimensions)

        self.global_optimum = [0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return sum(abs(x)) * exp(-sum(sin(x ** 2.0)))

# -------------------------------------------------------------------------------- #


class XinSheYang03(Benchmark):

    """
    Xin-She Yang 3 test objective function.

    This class defines the Xin-She Yang 3 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{XinSheYang03}}(\\mathbf{x}) = e^{-\\sum_{i=1}^{n} (x_i/\\beta)^{2m}} - 2e^{-\\sum_{i=1}^{n} x_i^2} \\prod_{i=1}^{n} \\cos^2(x_i)


    Where, in this exercise, :math:`\\beta = 15` and :math:`m = 3`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-20, 20]` for :math:`i=1,...,n`.

    .. figure:: figures/XinSheYang03.png
        :alt: Xin-She Yang 3 function
        :align: center

        **Two-dimensional Xin-She Yang 3 function**

    *Global optimum*: :math:`f(x_i) = -1` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-20.0] * self.dimensions, [20.0] * self.dimensions)

        self.global_optimum = [0 for _ in range(self.dimensions)]
        self.fglob = -1.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        beta, m = 15.0, 5.0
        return exp(-sum((x / beta) ** (2 * m))) - 2.0 * exp(-sum(x ** 2.0)) * prod(cos(x) ** 2.0)

# -------------------------------------------------------------------------------- #


class XinSheYang04(Benchmark):

    """
    Xin-She Yang 4 test objective function.

    This class defines the Xin-She Yang 4 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{XinSheYang04}}(\\mathbf{x}) = \\left[ \\sum_{i=1}^{n} \\sin^2(x_i) - e^{-\\sum_{i=1}^{n} x_i^2} \\right ] e^{-\\sum_{i=1}^{n} \\sin^2 \\sqrt{ \\lvert x_i \\rvert }}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    .. figure:: figures/XinSheYang04.png
        :alt: Xin-She Yang 4 function
        :align: center

        **Two-dimensional Xin-She Yang 4 function**

    *Global optimum*: :math:`f(x_i) = -1` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [0 for _ in range(self.dimensions)]
        self.fglob = -1.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        return (sum(sin(x) ** 2.0) - exp(-sum(x ** 2.0))) * exp(-sum(sin(sqrt(abs(x))) ** 2.0))

# -------------------------------------------------------------------------------- #


class Xor(Benchmark):

    """
    Xor test objective function.

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

        self.bounds = zip([-1.0] * self.dimensions, [1.0] * self.dimensions)

        self.global_optimum = [1.0, -1.0, 1.0,
                               -1.0, -1.0, 1.0, 1.0, -1.0, 0.421134]
        self.fglob = 0.9597588

    def evaluator(self, x, *args):

        self.fun_evals += 1

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

# -------------------------------------------------------------------------------- #


class YaoLiu04(Benchmark):

    """
    Yao-Liu 4 test objective function.

    This class defines the Yao-Liu function 4 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{YaoLiu04}}(\\mathbf{x}) = {max}_i \\left\{ \\left | x_i \\right | , 1 \\leq i \\leq n \\right\}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    .. figure:: figures/YaoLiu04.png
        :alt: YaoLiu04 function
        :align: center

        **Two-dimensional Yao-Liu 4 function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = [0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return abs(x.max())

# -------------------------------------------------------------------------------- #


class YaoLiu09(Benchmark):

    """
    Yao-Liu 9 test objective function.

    This class defines the Yao-Liu function 9 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{YaoLiu09}}(\\mathbf{x}) = \\sum_{i=1}^n \\left [ x_i^2 - 10 \\cos(2 \\pi x_i ) + 10 \\right ]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5.12, 5.12]` for :math:`i=1,...,n`.

    .. figure:: figures/YaoLiu09.png
        :alt: Yao-Liu 9 function
        :align: center

        **Two-dimensional Yao-Liu 9 function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-5.12] * self.dimensions, [5.12] * self.dimensions)

        self.global_optimum = [0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return sum(x ** 2.0 - 10.0 * cos(2 * pi * x) + 10)

# -------------------------------------------------------------------------------- #


class Zacharov(Benchmark):

    """
    Zacharov test objective function.

    This class defines the Zacharov global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

         f_{\\text{Zacharov}}(\\mathbf{x}) = \\sum_{i=1}^{n} x_i^2 + \\left ( \\frac{1}{2} \\sum_{i=1}^{n} i x_i \\right )^2 + \\left ( \\frac{1}{2} \\sum_{i=1}^{n} i x_i \\right )^4

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 10]` for :math:`i=1,...,n`.

    .. figure:: figures/Zacharov.png
        :alt: Zacharov function
        :align: center

        **Two-dimensional Zacharov function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-5.0] * self.dimensions, [10.0] * self.dimensions)
        self.custom_bounds = ([-1, 1], [-1, 1])

        self.global_optimum = [0 for _ in range(self.dimensions)]
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1

        rng = arange(1.0, self.dimensions + 1.0)
        return sum(x ** 2.0) + (0.5 * sum(rng * x)) ** 2.0 + (0.5 * sum(rng * x)) ** 4.0

# -------------------------------------------------------------------------------- #


class ZeroSum(Benchmark):

    """
    ZeroSum test objective function.

    This class defines the ZeroSum global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

         f_{\\text{ZeroSum}}(\\mathbf{x}) = \\begin{cases}0 & \\textrm{if} \\sum_{i=1}^n x_i = 0 \\\\
                1 + \\left(10000 \\left |\\sum_{i=1}^n x_i\\right| \\right)^{0.5} & \\textrm{otherwise}\\end{cases}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    .. figure:: figures/ZeroSum.png
        :alt: ZeroSum function
        :align: center

        **Two-dimensional ZeroSum function**

    *Global optimum*: :math:`f(x_i) = 0` where :math:`\\sum_{i=1}^n x_i = 0`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)

        self.global_optimum = []
        self.fglob = 0.0
        self.change_dimensionality = True

    def evaluator(self, x, *args):

        self.fun_evals += 1
        if abs(sum(x)) < 3e-16:
            return 0.0

        return 1.0 + (10000.0 * abs(sum(x))) ** 0.5

# -------------------------------------------------------------------------------- #


class Zettl(Benchmark):

    """
    Zettl test objective function.

    This class defines the Zettl global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Zettl}}(\\mathbf{x}) = \\frac{1}{4} x_{1} + \\left(x_{1}^{2} - 2 x_{1} + x_{2}^{2}\\right)^{2}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 5]` for :math:`i=1,2`.

    .. figure:: figures/Zettl.png
        :alt: Zettl function
        :align: center

        **Two-dimensional Zettl function**

    *Global optimum*: :math:`f(x_i) = -0.0037912` for :math:`\\mathbf{x} = [-0.029896, 0.0]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-1.0] * self.dimensions, [5.0] * self.dimensions)

        self.global_optimum = [-0.02989597760285287, 0.0]
        self.fglob = -0.003791237220468656

    def evaluator(self, x, *args):

        self.fun_evals += 1
        return (x[0] ** 2 + x[1] ** 2 - 2 * x[0]) ** 2 + x[0] / 4

# -------------------------------------------------------------------------------- #


class Zimmerman(Benchmark):

    """
    Zimmerman test objective function.

    This class defines the Zimmerman global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Zimmerman}}(\\mathbf{x}) = \\max \\left[Zh1(x), Zp(Zh2(x))\\textrm{sgn}(Zh2(x)), Zp(Zh3(x)) \\textrm{sgn}(Zh3(x)), Zp(-x_1)\\textrm{sgn}(x_1), Zp(-x_2)\\textrm{sgn}(x_2) \\right]

    Where, in this exercise:

    .. math::

        \\begin{cases} Zh1(x) = 9 - x_1 - x_2 \\\\
        Zh2(x) = (x_1 - 3)^2 + (x_2 - 2)^2 \\\\
        Zh3(x) = x_1x_2 - 14 \\\\
        Zp(t) = 100(1 + t) \\end{cases}

    Where :math:`x` is a vector and :math:`t` is a scalar.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 100]` for :math:`i=1,2`.

    .. figure:: figures/Zimmerman.png
        :alt: Zimmerman function
        :align: center

        **Two-dimensional Zimmerman function**

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [7, 2]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([0.0] * self.dimensions, [100.0] * self.dimensions)
        self.custom_bounds = ([0.0, 8.0], [0.0, 8.0])

        self.global_optimum = [7.0, 2.0]
        self.fglob = 0.0

    def evaluator(self, x, *args):

        self.fun_evals += 1

        Zh1 = (lambda x: 9.0 - x[0] - x[1])
        Zh2 = (lambda x: (x[0] - 3.0) ** 2.0 + (x[1] - 2.0) ** 2.0 - 16.0)
        Zh3 = (lambda x: x[0] * x[1] - 14.0)
        Zp = (lambda x: 100.0 * (1.0 + x))

        px = [Zh1(x), Zp(Zh2(x)) * sign(Zh2(x)), Zp(Zh3(x))
              * sign(Zh3(x)), Zp(-x[0]) * sign(x[0]), Zp(-x[1]) * sign(x[1])]
        return max(px)

# -------------------------------------------------------------------------------- #


class Zirilli(Benchmark):

    """
    Zettl test objective function.

    This class defines the Zirilli global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Zirilli}}(\\mathbf{x}) = 0.25x_1^4 - 0.5x_1^2 + 0.1x_1 + 0.5x_2^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    .. figure:: figures/Zirilli.png
        :alt: Zirilli function
        :align: center

        **Two-dimensional Zirilli function**

    *Global optimum*: :math:`f(x_i) = -0.3523` for :math:`\\mathbf{x} = [-1.0465, 0]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.bounds = zip([-10.0] * self.dimensions, [10.0] * self.dimensions)
        self.custom_bounds = ([-2.0, 2.0], [-2.0, 2.0])

        self.global_optimum = [-1.0465, 0.0]
        self.fglob = -0.35238603

    def evaluator(self, x, *args):

        self.fun_evals += 1

        x1, x2 = x
        return 0.25 * x1 ** 4.0 - 0.5 * x1 ** 2.0 + 0.1 * x1 + 0.5 * x2 ** 2

# -------------------------------------------------------------------------------- #


#-----------------------------------------------------------------------
#                 UNIVARIATE SINGLE-OBJECTIVE PROBLEMS
#-----------------------------------------------------------------------

class Problem02(Benchmark):

    """
    Univariate Problem02 test objective function.

    This class defines the Univariate Problem02 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Problem02}}(x) = \\sin(x) + \\sin \\left(\\frac{10}{3}x \\right)

    Bound constraints: :math:`x \\in [2.7, 7.5]`

    .. figure:: figures/Problem02.png
        :alt: Univariate Problem02 function
        :align: center

        **Univariate Problem02 function**

    *Global optimum*: :math:`f(x)=-1.899599` for :math:`x = 5.145735`

    """

    def __init__(self, dimensions=1):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(2.7, 7.5)]

        self.global_optimum = 5.145735
        self.fglob = -1.899599

    def evaluator(self, x, *args):

        self.fun_evals += 1
        x = x[0]
        return sin(x) + sin(10.0 / 3.0 * x)

# -------------------------------------------------------------------------------- #


class Problem03(Benchmark):

    """
    Univariate Problem03 test objective function.

    This class defines the Univariate Problem03 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Problem03}}(x) = - \\sum_{k=1}^6 k \\sin[(k+1)x+k]

    Bound constraints: :math:`x \\in [-10, 10]`

    .. figure:: figures/Problem03.png
        :alt: Univariate Problem03 function
        :align: center

        **Univariate Problem03 function**

    *Global optimum*: :math:`f(x)=-12.03124` for :math:`x = -6.7745761`

    """

    def __init__(self, dimensions=1):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(-10, 10)]

        self.global_optimum = -6.7745761
        self.fglob = -12.03124

    def evaluator(self, x, *args):

        self.fun_evals += 1
        x = x[0]

        y = 0.0
        for k in range(1, 6):
            y += k * sin((k + 1) * x + k)

        return -y

# -------------------------------------------------------------------------------- #


class Problem04(Benchmark):

    """
    Univariate Problem04 test objective function.

    This class defines the Univariate Problem04 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Problem04}}(x) = - \\left(16x^2 - 24x + 5 \\right) e^{-x}

    Bound constraints: :math:`x \\in [1.9, 3.9]`

    .. figure:: figures/Problem04.png
        :alt: Univariate Problem04 function
        :align: center

        **Univariate Problem04 function**

    *Global optimum*: :math:`f(x)=-3.85045` for :math:`x = 2.868034`

    """

    def __init__(self, dimensions=1):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(1.9, 3.9)]

        self.global_optimum = 2.868034
        self.fglob = -3.85045

    def evaluator(self, x, *args):

        self.fun_evals += 1
        x = x[0]
        return -(16 * x ** 2 - 24 * x + 5) * exp(-x)

# -------------------------------------------------------------------------------- #


class Problem05(Benchmark):

    """
    Univariate Problem05 test objective function.

    This class defines the Univariate Problem05 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Problem05}}(x) = - \\left(1.4 - 3x \\right) \\sin(18x)

    Bound constraints: :math:`x \\in [0, 1.2]`

    .. figure:: figures/Problem05.png
        :alt: Univariate Problem05 function
        :align: center

        **Univariate Problem05 function**

    *Global optimum*: :math:`f(x)=-1.48907` for :math:`x = 0.96609`

    """

    def __init__(self, dimensions=1):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(0.0, 1.2)]

        self.global_optimum = 0.96609
        self.fglob = -1.48907

    def evaluator(self, x, *args):

        self.fun_evals += 1
        x = x[0]
        return -(1.4 - 3 * x) * sin(18.0 * x)

# -------------------------------------------------------------------------------- #


class Problem06(Benchmark):

    """
    Univariate Problem06 test objective function.

    This class defines the Univariate Problem06 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Problem06}}(x) = - \\left[x + \\sin(x) \\right] e^{-x^2}

    Bound constraints: :math:`x \\in [-10, 10]`

    .. figure:: figures/Problem06.png
        :alt: Univariate Problem06 function
        :align: center

        **Univariate Problem06 function**

    *Global optimum*: :math:`f(x)=-0.824239` for :math:`x = 0.67956`

    """

    def __init__(self, dimensions=1):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(-10.0, 10.0)]

        self.global_optimum = 0.67956
        self.fglob = -0.824239

    def evaluator(self, x, *args):

        self.fun_evals += 1
        x = x[0]
        return -(x + sin(x)) * exp(-x ** 2.0)

# -------------------------------------------------------------------------------- #


class Problem07(Benchmark):

    """
    Univariate Problem07 test objective function.

    This class defines the Univariate Problem07 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Problem07}}(x) = \\sin(x) + \\sin \\left(\\frac{10}{3}x \\right) + \\log(x) - 0.84x + 3

    Bound constraints: :math:`x \\in [2.7, 7.5]`

    .. figure:: figures/Problem07.png
        :alt: Univariate Problem07 function
        :align: center

        **Univariate Problem07 function**

    *Global optimum*: :math:`f(x)=-1.6013` for :math:`x = 5.19978`

    """

    def __init__(self, dimensions=1):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(2.7, 7.5)]

        self.global_optimum = 5.19978
        self.fglob = -1.6013

    def evaluator(self, x, *args):

        self.fun_evals += 1
        x = x[0]
        return sin(x) + sin(10.0 / 3.0 * x) + log(x) - 0.84 * x + 3

# -------------------------------------------------------------------------------- #


class Problem08(Benchmark):

    """
    Univariate Problem08 test objective function.

    This class defines the Univariate Problem08 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Problem08}}(x) = - \\sum_{k=1}^6 k \\cos[(k+1)x+k]

    Bound constraints: :math:`x \\in [-10, 10]`

    .. figure:: figures/Problem08.png
        :alt: Univariate Problem08 function
        :align: center

        **Univariate Problem08 function**

    *Global optimum*: :math:`f(x)=-14.508` for :math:`x = -7.083506`

    """

    def __init__(self, dimensions=1):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(-10, 10)]

        self.global_optimum = -7.083506
        self.fglob = -14.508

    def evaluator(self, x, *args):

        self.fun_evals += 1
        x = x[0]

        y = 0.0
        for k in range(1, 6):
            y += k * cos((k + 1) * x + k)

        return -y

# -------------------------------------------------------------------------------- #


class Problem09(Benchmark):

    """
    Univariate Problem09 test objective function.

    This class defines the Univariate Problem09 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Problem09}}(x) = \\sin(x) + \\sin \\left(\\frac{2}{3} x \\right)

    Bound constraints: :math:`x \\in [3.1, 20.4]`

    .. figure:: figures/Problem09.png
        :alt: Univariate Problem09 function
        :align: center

        **Univariate Problem09 function**

    *Global optimum*: :math:`f(x)=-1.90596` for :math:`x = 17.039`

    """

    def __init__(self, dimensions=1):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(3.1, 20.4)]

        self.global_optimum = 17.039
        self.fglob = -1.90596

    def evaluator(self, x, *args):

        self.fun_evals += 1
        x = x[0]
        return sin(x) + sin(2.0 / 3.0 * x)

# -------------------------------------------------------------------------------- #


class Problem10(Benchmark):

    """
    Univariate Problem10 test objective function.

    This class defines the Univariate Problem10 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Problem10}}(x) = -x\\sin(x)

    Bound constraints: :math:`x \\in [0, 10]`

    .. figure:: figures/Problem10.png
        :alt: Univariate Problem10 function
        :align: center

        **Univariate Problem10 function**

    *Global optimum*: :math:`f(x)=-7.916727` for :math:`x = 7.9787`

    """

    def __init__(self, dimensions=1):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(0, 10)]

        self.global_optimum = 7.9787
        self.fglob = -7.916727

    def evaluator(self, x, *args):

        self.fun_evals += 1
        x = x[0]
        return -x * sin(x)

# -------------------------------------------------------------------------------- #


class Problem11(Benchmark):

    """
    Univariate Problem11 test objective function.

    This class defines the Univariate Problem11 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Problem11}}(x) = 2\\cos(x) + \\cos(2x)

    Bound constraints: :math:`x \\in [-\\pi/2, 2\\pi]`

    .. figure:: figures/Problem11.png
        :alt: Univariate Problem11 function
        :align: center

        **Univariate Problem11 function**

    *Global optimum*: :math:`f(x)=-1.5` for :math:`x = 2.09439`

    """

    def __init__(self, dimensions=1):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(-pi / 2, 2 * pi)]

        self.global_optimum = 2.09439
        self.fglob = -1.5

    def evaluator(self, x, *args):

        self.fun_evals += 1
        x = x[0]
        return 2 * cos(x) + cos(2 * x)

# -------------------------------------------------------------------------------- #


class Problem12(Benchmark):

    """
    Univariate Problem12 test objective function.

    This class defines the Univariate Problem12 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Problem12}}(x) = \\sin^3(x) + \\cos^3(x)

    Bound constraints: :math:`x \\in [0, 2\\pi]`

    .. figure:: figures/Problem12.png
        :alt: Univariate Problem12 function
        :align: center

        **Univariate Problem12 function**

    *Global optimum*: :math:`f(x)=-1` for :math:`x = \\pi`

    """

    def __init__(self, dimensions=1):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(0, 2 * pi)]

        self.global_optimum = pi
        self.fglob = -1

    def evaluator(self, x, *args):

        self.fun_evals += 1
        x = x[0]
        return (sin(x)) ** 3.0 + (cos(x)) ** 3.0

# -------------------------------------------------------------------------------- #


class Problem13(Benchmark):

    """
    Univariate Problem13 test objective function.

    This class defines the Univariate Problem13 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Problem13}}(x) = -x^{2/3} - (1 - x^2)^{1/3}

    Bound constraints: :math:`x \\in [0.001, 0.99]`

    .. figure:: figures/Problem13.png
        :alt: Univariate Problem13 function
        :align: center

        **Univariate Problem13 function**

    *Global optimum*: :math:`f(x)=-1.5874` for :math:`x = 1/\\sqrt(2)`

    """

    def __init__(self, dimensions=1):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(0.001, 0.99)]

        self.global_optimum = 1.0 / sqrt(2)
        self.fglob = -1.5874

    def evaluator(self, x, *args):

        self.fun_evals += 1
        x = x[0]
        return -x ** (2.0 / 3.0) - (1.0 - x ** 2) ** (1.0 / 3.0)

# -------------------------------------------------------------------------------- #


class Problem14(Benchmark):

    """
    Univariate Problem14 test objective function.

    This class defines the Univariate Problem14 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Problem14}}(x) = -e^{-x} \\sin(2\\pi x)

    Bound constraints: :math:`x \\in [0, 4]`

    .. figure:: figures/Problem14.png
        :alt: Univariate Problem14 function
        :align: center

        **Univariate Problem14 function**

    *Global optimum*: :math:`f(x)=-0.788685` for :math:`x = 0.224885`

    """

    def __init__(self, dimensions=1):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(0.0, 4.0)]

        self.global_optimum = 0.224885
        self.fglob = -0.788685

    def evaluator(self, x, *args):

        self.fun_evals += 1
        x = x[0]
        return -exp(-x) * sin(2.0 * pi * x)

# -------------------------------------------------------------------------------- #


class Problem15(Benchmark):

    """
    Univariate Problem15 test objective function.

    This class defines the Univariate Problem15 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Problem15}}(x) = \\frac{x^{2} - 5 x + 6}{x^{2} + 1}

    Bound constraints: :math:`x \\in [-5, 5]`

    .. figure:: figures/Problem15.png
        :alt: Univariate Problem15 function
        :align: center

        **Univariate Problem15 function**

    *Global optimum*: :math:`f(x)=-0.03553` for :math:`x = 2.41422`

    """

    def __init__(self, dimensions=1):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(-5.0, 5.0)]

        self.global_optimum = 2.41422
        self.fglob = -0.03553

    def evaluator(self, x, *args):

        self.fun_evals += 1
        x = x[0]
        return -(-x ** 2.0 + 5 * x - 6) / (x ** 2 + 1)

# -------------------------------------------------------------------------------- #


class Problem18(Benchmark):

    """
    Univariate Problem18 test objective function.

    This class defines the Univariate Problem18 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

         f_{\\text{Problem18}}(x) = \\begin{cases}(x-2)^2 & \\textrm{if} \\hspace{5pt} x \\leq 3 \\\\
                                                  2\\log(x-2)+1&\\textrm{otherwise}\\end{cases}

    Bound constraints: :math:`x \\in [0, 6]`

    .. figure:: figures/Problem18.png
        :alt: Univariate Problem18 function
        :align: center

        **Univariate Problem18 function**

    *Global optimum*: :math:`f(x)=0` for :math:`x = 2`

    """

    def __init__(self, dimensions=1):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(0.0, 6.0)]

        self.global_optimum = 2
        self.fglob = 0

    def evaluator(self, x, *args):

        self.fun_evals += 1
        x = x[0]

        if x <= 3:
            return (x - 2.0) ** 2.0

        return 2 * log(x - 2.0) + 1

# -------------------------------------------------------------------------------- #


class Problem20(Benchmark):

    """
    Univariate Problem20 test objective function.

    This class defines the Univariate Problem20 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Problem20}}(x) = -[x-\\sin(x)]e^{-x^2}

    Bound constraints: :math:`x \\in [-10, 10]`

    .. figure:: figures/Problem20.png
        :alt: Univariate Problem20 function
        :align: center

        **Univariate Problem20 function**

    *Global optimum*: :math:`f(x)=-0.0634905` for :math:`x = 1.195137`

    """

    def __init__(self, dimensions=1):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(-10, 10)]

        self.global_optimum = 1.195137
        self.fglob = -0.0634905

    def evaluator(self, x, *args):

        self.fun_evals += 1
        x = x[0]
        return -(x - sin(x)) * exp(-x ** 2.0)

# -------------------------------------------------------------------------------- #


class Problem21(Benchmark):

    """
    Univariate Problem21 test objective function.

    This class defines the Univariate Problem21 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Problem21}}(x) = x \\sin(x) + x \\cos(2x)

    Bound constraints: :math:`x \\in [0, 10]`

    .. figure:: figures/Problem21.png
        :alt: Univariate Problem21 function
        :align: center

        **Univariate Problem21 function**

    *Global optimum*: :math:`f(x)=-9.50835` for :math:`x = 4.79507`

    """

    def __init__(self, dimensions=1):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(0, 10)]

        self.global_optimum = 4.79507
        self.fglob = -9.50835

    def evaluator(self, x, *args):

        self.fun_evals += 1
        x = x[0]
        return x * sin(x) + x * cos(2.0 * x)

# -------------------------------------------------------------------------------- #


class Problem22(Benchmark):

    """
    Univariate Problem22 test objective function.

    This class defines the Univariate Problem22 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Problem22}}(x) = e^{-3x} - \\sin^3(x)

    Bound constraints: :math:`x \\in [0, 20]`

    .. figure:: figures/Problem22.png
        :alt: Univariate Problem22 function
        :align: center

        **Univariate Problem22 function**

    *Global optimum*: :math:`f(x)=e^{-27\\pi/2} - 1` for :math:`x = 9\\pi/2`

    """

    def __init__(self, dimensions=1):
        Benchmark.__init__(self, dimensions)

        self.bounds = [(0, 20)]

        self.global_optimum = 9.0 * pi / 2.0
        self.fglob = exp(-27.0 * pi / 2.0) - 1.0

    def evaluator(self, x, *args):

        self.fun_evals += 1
        x = x[0]
        return exp(-3.0 * x) - (sin(x)) ** 3.0

# -------------------------------------------------------------------------------- #
