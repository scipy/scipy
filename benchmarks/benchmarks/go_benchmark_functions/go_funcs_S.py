# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import

from numpy import (abs, arctan2, asarray, cos, exp, floor, log, log10,
                   arange, pi, prod, roll, sign, sin, sqrt, sum,
                   tan, tanh, dot, repeat, atleast_2d, tril)
from numpy.random import uniform
from .go_benchmark import Benchmark


class Salomon(Benchmark):

    r"""
    Salomon objective function.

    This class defines the Salomon [1]_ global optimization problem. This is a
    multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Salomon}}(x) = 1 - \cos \left (2 \pi
        \sqrt{\sum_{i=1}^{n} x_i^2} \right) + 0.1 \sqrt{\sum_{i=1}^n x_i^2}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \in
    [-100, 100]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x_i = 0` for
    :math:`i = 1, ..., n`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = list(zip([-100.0] * self.N,
                           [100.0] * self.N))
        self.custom_bounds = [(-50, 50), (-50, 50)]

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        u = sqrt(sum(x ** 2))
        return 1 - cos(2 * pi * u) + 0.1 * u


class Sargan(Benchmark):

    r"""
    Sargan objective function.

    This class defines the Sargan [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Sargan}}(x) = \sum_{i=1}^{n} n \left (x_i^2
        + 0.4 \sum_{i \neq j}^{n} x_ix_j \right)

    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-100, 100]` for
    :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x_i = 0` for
    :math:`i = 1, ..., n`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = list(zip([-100.0] * self.N,
                           [100.0] * self.N))
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

    r"""
    Schaffer 1 objective function.

    This class defines the Schaffer 1 [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Schaffer01}}(x) = 0.5 + \frac{\sin^2 (x_1^2 + x_2^2)^2 - 0.5}
        {1 + 0.001(x_1^2 + x_2^2)^2}

    with :math:`x_i \in [-100, 100]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x = [0, 0]` for
    :math:`i = 1, 2`

    .. [1] Mishra, S. Some new test functions for global optimization and
    performance of repulsive particle swarm method.
    Munich Personal RePEc Archive, 2006, 2718
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-100.0] * self.N,
                           [100.0] * self.N))
        self.custom_bounds = [(-10, 10), (-10, 10)]

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        u = (x[0] ** 2 + x[1] ** 2)
        num = sin(u) ** 2 - 0.5
        den = (1 + 0.001 * u) ** 2
        return 0.5 + num / den


class Schaffer02(Benchmark):

    r"""
    Schaffer 2 objective function.

    This class defines the Schaffer 2 [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Schaffer02}}(x) = 0.5 + \frac{\sin^2 (x_1^2 - x_2^2)^2 - 0.5}
        {1 + 0.001(x_1^2 + x_2^2)^2}

    with :math:`x_i \in [-100, 100]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x = [0, 0]` for
    :math:`i = 1, 2`

    .. [1] Mishra, S. Some new test functions for global optimization and
    performance of repulsive particle swarm method.
    Munich Personal RePEc Archive, 2006, 2718
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-100.0] * self.N,
                           [100.0] * self.N))
        self.custom_bounds = [(-10, 10), (-10, 10)]

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        num = sin((x[0] ** 2 - x[1] ** 2)) ** 2 - 0.5
        den = (1 + 0.001 * (x[0] ** 2 + x[1] ** 2)) ** 2
        return 0.5 + num / den


class Schaffer03(Benchmark):

    r"""
    Schaffer 3 objective function.

    This class defines the Schaffer 3 [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\text{Schaffer03}}(x) = 0.5 + \frac{\sin^2 \left( \cos \lvert x_1^2
       - x_2^2 \rvert \right ) - 0.5}{1 + 0.001(x_1^2 + x_2^2)^2}

    with :math:`x_i \in [-100, 100]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = 0.00156685` for :math:`x = [0, 1.253115]`

    .. [1] Mishra, S. Some new test functions for global optimization and
    performance of repulsive particle swarm method.
    Munich Personal RePEc Archive, 2006, 2718
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-100.0] * self.N,
                           [100.0] * self.N))
        self.custom_bounds = [(-10, 10), (-10, 10)]

        self.global_optimum = [[0.0, 1.253115]]
        self.fglob = 0.00156685

    def fun(self, x, *args):
        self.nfev += 1

        num = sin(cos(abs(x[0] ** 2 - x[1] ** 2))) ** 2 - 0.5
        den = (1 + 0.001 * (x[0] ** 2 + x[1] ** 2)) ** 2
        return 0.5 + num / den


class Schaffer04(Benchmark):

    r"""
    Schaffer 4 objective function.

    This class defines the Schaffer 4 [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Schaffer04}}(x) = 0.5 + \frac{\cos^2 \left( \sin(x_1^2 - x_2^2)
        \right ) - 0.5}{1 + 0.001(x_1^2 + x_2^2)^2}^2

    with :math:`x_i \in [-100, 100]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = 0.292579` for :math:`x = [0, 1.253115]`

    .. [1] Mishra, S. Some new test functions for global optimization and
    performance of repulsive particle swarm method.
    Munich Personal RePEc Archive, 2006, 2718
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-100.0] * self.N,
                           [100.0] * self.N))
        self.custom_bounds = [(-10, 10), (-10, 10)]

        self.global_optimum = [[0.0, 1.253115]]
        self.fglob = 0.292579

    def fun(self, x, *args):
        self.nfev += 1

        num = cos(sin(abs(x[0] ** 2 - x[1] ** 2))) ** 2 - 0.5
        den = (1 + 0.001 * (x[0] ** 2 + x[1] ** 2)) ** 2
        return 0.5 + num / den


# class SchmidtVetters(Benchmark):
#
#     r"""
#     Schmidt-Vetters objective function.
#
#     This class defines the Schmidt-Vetters global optimization problem. This
#     is a multimodal minimization problem defined as follows:
#
#     .. math::
#
#         f_{\text{SchmidtVetters}}(x) = \frac{1}{1 + (x_1 - x_2)^2}
#         + \sin \left(\frac{\pi x_2 + x_3}{2} \right)
#         + e^{\left(\frac{x_1+x_2}{x_2} - 2\right)^2}
#
#     with :math:`x_i \in [0, 10]` for :math:`i = 1, 2, 3`.
#
#     *Global optimum*: :math:`f(x) = 2.99643266` for
#     :math:`x = [0.79876108,  0.79962581,  0.79848824]`
#
#     TODO equation seems right, but [7.07083412 , 10., 3.14159293] produces a
#     lower minimum, 0.193973
#     """
#
#     def __init__(self, dimensions=3):
#         Benchmark.__init__(self, dimensions)
#         self._bounds = zip([0.0] * self.N, [10.0] * self.N)
#
#         self.global_optimum = [[0.79876108, 0.79962581, 0.79848824]]
#         self.fglob = 2.99643266
#
#     def fun(self, x, *args):
#         self.nfev += 1
#
#         return (1 / (1 + (x[0] - x[1]) ** 2) + sin((pi * x[1] + x[2]) / 2)
#                 + exp(((x[0] + x[1]) / x[1] - 2) ** 2))


class Schwefel01(Benchmark):

    r"""
    Schwefel 1 objective function.

    This class defines the Schwefel 1 [1]_ global optimization problem. This is a
    unimodal minimization problem defined as follows:

    .. math::

       f_{\text{Schwefel01}}(x) = \left(\sum_{i=1}^n x_i^2 \right)^{\alpha}


    Where, in this exercise, :math:`\alpha = \sqrt{\pi}`.

    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-100, 100]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x_i = 0`
    for :math:`i = 1, ..., n`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = list(zip([-100.0] * self.N,
                           [100.0] * self.N))
        self.custom_bounds = ([-4.0, 4.0], [-4.0, 4.0])

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        alpha = sqrt(pi)
        return (sum(x ** 2.0)) ** alpha


class Schwefel02(Benchmark):

    r"""
    Schwefel 2 objective function.

    This class defines the Schwefel 2 [1]_ global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

        f_{\text{Schwefel02}}(x) = \sum_{i=1}^n \left(\sum_{j=1}^i 
        x_i \right)^2


    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-100, 100]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x_i = 0` for
    :math:`i = 1, ..., n`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = list(zip([-100.0] * self.N,
                           [100.0] * self.N))
        self.custom_bounds = ([-4.0, 4.0], [-4.0, 4.0])

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        mat = repeat(atleast_2d(x), self.N, axis=0)
        inner = sum(tril(mat), axis=1)
        return sum(inner ** 2)


class Schwefel04(Benchmark):

    r"""
    Schwefel 4 objective function.

    This class defines the Schwefel 4 [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Schwefel04}}(x) = \sum_{i=1}^n \left[(x_i - 1)^2
        + (x_1 - x_i^2)^2 \right]


    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [0, 10]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0` for:math:`x_i = 1` for
    :math:`i = 1, ..., n`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = list(zip([0.0] * self.N, [10.0] * self.N))
        self.custom_bounds = ([0.0, 2.0], [0.0, 2.0])

        self.global_optimum = [[1.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum((x - 1.0) ** 2.0 + (x[0] - x ** 2.0) ** 2.0)


class Schwefel06(Benchmark):

    r"""
    Schwefel 6 objective function.

    This class defines the Schwefel 6 [1]_ global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\text{Schwefel06}}(x) = \max(\lvert x_1 + 2x_2 - 7 \rvert,
                                   \lvert 2x_1 + x_2 - 5 \rvert)


    with :math:`x_i \in [-100, 100]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x = [1, 3]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = list(zip([-100.0] * self.N,
                           [100.0] * self.N))
        self.custom_bounds = ([-10.0, 10.0], [-10.0, 10.0])

        self.global_optimum = [[1.0, 3.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return max(abs(x[0] + 2 * x[1] - 7), abs(2 * x[0] + x[1] - 5))


class Schwefel20(Benchmark):

    r"""
    Schwefel 20 objective function.

    This class defines the Schwefel 20 [1]_ global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\text{Schwefel20}}(x) = \sum_{i=1}^n \lvert x_i \rvert


    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-100, 100]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x_i = 0` for
    :math:`i = 1, ..., n`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    TODO: Jamil #122 is incorrect.  There shouldn't be a leading minus sign.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = list(zip([-100.0] * self.N,
                           [100.0] * self.N))

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum(abs(x))


class Schwefel21(Benchmark):

    r"""
    Schwefel 21 objective function.

    This class defines the Schwefel 21 [1]_ global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

        f_{\text{Schwefel21}}(x) = \smash{\displaystyle\max_{1 \leq i \leq n}}
                                   \lvert x_i \rvert


    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-100, 100]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x_i = 0` for
    :math:`i = 1, ..., n`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = list(zip([-100.0] * self.N,
                           [100.0] * self.N))

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return max(abs(x))


class Schwefel22(Benchmark):

    r"""
    Schwefel 22 objective function.

    This class defines the Schwefel 22 [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Schwefel22}}(x) = \sum_{i=1}^n \lvert x_i \rvert
                                  + \prod_{i=1}^n \lvert x_i \rvert


    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-100, 100]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x_i = 0` for
    :math:`i = 1, ..., n`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = list(zip([-100.0] * self.N,
                           [100.0] * self.N))
        self.custom_bounds = ([-10.0, 10.0], [-10.0, 10.0])

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum(abs(x)) + prod(abs(x))


class Schwefel26(Benchmark):

    r"""
    Schwefel 26 objective function.

    This class defines the Schwefel 26 [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Schwefel26}}(x) = 418.9829n - \sum_{i=1}^n x_i
                                  \sin(\sqrt{|x_i|})

    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-500, 500]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x_i = 420.968746` for
    :math:`i = 1, ..., n`

    .. [1] Gavana, A. Global Optimization Benchmarks and AMPGO retrieved 2015
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = list(zip([-500.0] * self.N,
                           [500.0] * self.N))

        self.global_optimum = [[420.968746 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return 418.982887 * self.N - sum(x * sin(sqrt(abs(x))))


class Schwefel36(Benchmark):

    r"""
    Schwefel 36 objective function.

    This class defines the Schwefel 36 [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Schwefel36}}(x) = -x_1x_2(72 - 2x_1 - 2x_2)


    with :math:`x_i \in [0, 500]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = -3456` for :math:`x = [12, 12]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = list(zip([0.0] * self.N, [500.0] * self.N))
        self.custom_bounds = ([0.0, 20.0], [0.0, 20.0])

        self.global_optimum = [[12.0, 12.0]]
        self.fglob = -3456.0

    def fun(self, x, *args):
        self.nfev += 1

        return -x[0] * x[1] * (72 - 2 * x[0] - 2 * x[1])


class Shekel05(Benchmark):

    r"""
    Shekel 5 objective function.

    This class defines the Shekel 5 [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Shekel05}}(x) = \sum_{i=1}^{m} \frac{1}{c_{i}
        + \sum_{j=1}^{n} (x_{j} - a_{ij})^2 }`

    Where, in this exercise:

    .. math::

        a = 
        \begin{bmatrix}
        4.0 & 4.0 & 4.0 & 4.0 \\ 1.0 & 1.0 & 1.0 & 1.0 \\
        8.0 & 8.0 & 8.0 & 8.0 \\ 6.0 & 6.0 & 6.0 & 6.0 \\
        3.0 & 7.0 & 3.0 & 7.0 
        \end{bmatrix}
    .. math::

        c = \begin{bmatrix} 0.1 \\ 0.2 \\ 0.2 \\ 0.4 \\ 0.4 \end{bmatrix}

    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [0, 10]` for :math:`i = 1, ..., 4`.

    *Global optimum*: :math:`f(x) = -10.15319585` for :math:`x_i = 4` for
    :math:`i = 1, ..., 4`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    TODO: this is a different global minimum compared to Jamil#130.  The
    minimum is found by doing lots of optimisations. The solution is supposed
    to be at [4] * N, is there any numerical overflow?
    """

    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([0.0] * self.N, [10.0] * self.N))

        self.global_optimum = [[4.00003715092,
                                4.00013327435,
                                4.00003714871,
                                4.0001332742]]
        self.fglob = -10.1531996791
        self.A = asarray([[4.0, 4.0, 4.0, 4.0],
                          [1.0, 1.0, 1.0, 1.0],
                          [8.0, 8.0, 8.0, 8.0],
                          [6.0, 6.0, 6.0, 6.0],
                          [3.0, 7.0, 3.0, 7.0]])

        self.C = asarray([0.1, 0.2, 0.2, 0.4, 0.4])

    def fun(self, x, *args):
        self.nfev += 1

        return -sum(1 / (sum((x - self.A) ** 2, axis=1) + self.C))


class Shekel07(Benchmark):

    r"""
    Shekel 7 objective function.

    This class defines the Shekel 7 [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Shekel07}}(x) = \sum_{i=1}^{m} \frac{1}{c_{i}
                                 + \sum_{j=1}^{n} (x_{j} - a_{ij})^2 }`

    Where, in this exercise:

    .. math::

        a =
        \begin{bmatrix}
        4.0 & 4.0 & 4.0 & 4.0 \\ 1.0 & 1.0 & 1.0 & 1.0 \\
        8.0 & 8.0 & 8.0 & 8.0 \\ 6.0 & 6.0 & 6.0 & 6.0 \\
        3.0 & 7.0 & 3.0 & 7.0 \\ 2.0 & 9.0 & 2.0 & 9.0 \\
        5.0 & 5.0 & 3.0 & 3.0
        \end{bmatrix}


    .. math::

        c =
        \begin{bmatrix}
        0.1 \\ 0.2 \\ 0.2 \\ 0.4 \\ 0.4 \\ 0.6 \\ 0.3 
        \end{bmatrix}


    with :math:`x_i \in [0, 10]` for :math:`i = 1, ..., 4`.

    *Global optimum*: :math:`f(x) = -10.4028188` for :math:`x_i = 4` for
    :math:`i = 1, ..., 4`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    TODO: this is a different global minimum compared to Jamil#131. This
    minimum is obtained after running lots of minimisations!  Is there any
    numerical overflow that causes the minimum solution to not be [4] * N?
    """

    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([0.0] * self.N, [10.0] * self.N))

        self.global_optimum = [[4.00057291078,
                                4.0006893679,
                                3.99948971076,
                                3.99960615785]]
        self.fglob = -10.4029405668
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

        return -sum(1 / (sum((x - self.A) ** 2, axis=1) + self.C))


class Shekel10(Benchmark):

    r"""
    Shekel 10 objective function.

    This class defines the Shekel 10 [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\text{Shekel10}}(x) = \sum_{i=1}^{m} \frac{1}{c_{i} 
                                + \sum_{j=1}^{n} (x_{j} - a_{ij})^2 }`

    Where, in this exercise:

    .. math::

        a =
        \begin{bmatrix}
        4.0 & 4.0 & 4.0 & 4.0 \\ 1.0 & 1.0 & 1.0 & 1.0 \\
        8.0 & 8.0 & 8.0 & 8.0 \\ 6.0 & 6.0 & 6.0 & 6.0 \\
        3.0 & 7.0 & 3.0 & 7.0 \\ 2.0 & 9.0 & 2.0 & 9.0 \\
        5.0 & 5.0 & 3.0 & 3.0 \\ 8.0 & 1.0 & 8.0 & 1.0 \\
        6.0 & 2.0 & 6.0 & 2.0 \\ 7.0 & 3.6 & 7.0 & 3.6
        \end{bmatrix}


    .. math::

        c =
        \begin{bmatrix}
        0.1 \\ 0.2 \\ 0.2 \\ 0.4 \\ 0.4 \\ 0.6 \\ 0.3 \\ 0.7 \\ 0.5 \\ 0.5
        \end{bmatrix}


    with :math:`x_i \in [0, 10]` for :math:`i = 1, ..., 4`.

    *Global optimum*: :math:`f(x) = -10.5362837` for :math:`x_i = 4` for
    :math:`i = 1, ..., 4`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    TODO Found a lower global minimum than Jamil#132... Is this numerical overflow?
    """

    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([0.0] * self.N, [10.0] * self.N))

        self.global_optimum = [[4.0007465377266271,
                                4.0005929234621407,
                                3.9996633941680968,
                                3.9995098017834123]]
        self.fglob = -10.536409816692023
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

        return -sum(1 / (sum((x - self.A) ** 2, axis=1) + self.C))


class Shubert01(Benchmark):

    r"""
    Shubert 1 objective function.

    This class defines the Shubert 1 [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Shubert01}}(x) = \prod_{i=1}^{n}\left(\sum_{j=1}^{5}
                                  cos(j+1)x_i+j \right )

    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-10, 10]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = -186.7309` for
    :math:`x = [-7.0835, 4.8580]` (and many others).

    .. [1] Gavana, A. Global Optimization Benchmarks and AMPGO retrieved 2015

    TODO: Jamil#133 is missing a prefactor of j before the cos function.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-10.0] * self.N, [10.0] * self.N))
        self.global_optimum = [[-7.0835, 4.8580]]

        self.fglob = -186.7309

        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        j = atleast_2d(arange(1, 6)).T
        y = j * cos((j + 1) * x + j)
        return prod(sum(y, axis=0))


class Shubert03(Benchmark):

    r"""
    Shubert 3 objective function.

    This class defines the Shubert 3 [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Shubert03}}(x) = \sum_{i=1}^n \sum_{j=1}^5 -j 
                                  \sin((j+1)x_i + j)

    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-10, 10]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = -24.062499` for
    :math:`x = [5.791794, 5.791794]` (and many others).

     .. [1] Gavana, A. Global Optimization Benchmarks and AMPGO retrieved 2015

    TODO: Jamil#134 has wrong global minimum value, and is missing a minus sign
    before the whole thing.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-10.0] * self.N, [10.0] * self.N))

        self.global_optimum = [[5.791794, 5.791794]]
        self.fglob = -24.062499

        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        j = atleast_2d(arange(1, 6)).T
        y = -j * sin((j + 1) * x + j)
        return sum(sum(y))


class Shubert04(Benchmark):

    r"""
    Shubert 4 objective function.

    This class defines the Shubert 4 [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Shubert04}}(x) = \left(\sum_{i=1}^n \sum_{j=1}^5 -j
                                  \cos ((j+1)x_i + j)\right)

    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-10, 10]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = -29.016015` for
    :math:`x = [-0.80032121, -7.08350592]` (and many others).

    .. [1] Gavana, A. Global Optimization Benchmarks and AMPGO retrieved 2015

    TODO: Jamil#135 has wrong global minimum value, and is missing a minus sign
    before the whole thing.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-10.0] * self.N, [10.0] * self.N))

        self.global_optimum = [[-0.80032121, -7.08350592]]
        self.fglob = -29.016015

        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        j = atleast_2d(arange(1, 6)).T
        y = -j * cos((j + 1) * x + j)
        return sum(sum(y))


class SineEnvelope(Benchmark):

    r"""
    SineEnvelope objective function.

    This class defines the SineEnvelope [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{SineEnvelope}}(x) = -\sum_{i=1}^{n-1}\left[\frac{\sin^2(
                                       \sqrt{x_{i+1}^2+x_{i}^2}-0.5)}
                                       {(0.001(x_{i+1}^2+x_{i}^2)+1)^2}
                                       + 0.5\right]

    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-100, 100]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x_i = 0` for
    :math:`i = 1, ..., n`

    .. [1] Gavana, A. Global Optimization Benchmarks and AMPGO retrieved 2015

    TODO: Jamil #136
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-100.0] * self.N,
                           [100.0] * self.N))
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

    r"""
    Six Hump Camel objective function.

    This class defines the Six Hump Camel [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{SixHumpCamel}}(x) = 4x_1^2+x_1x_2-4x_2^2-2.1x_1^4+
                                    4x_2^4+\frac{1}{3}x_1^6

    with :math:`x_i \in [-5, 5]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = -1.031628453489877` for
    :math:`x = [0.08984201368301331 , -0.7126564032704135]` or 
    :math:`x = [-0.08984201368301331, 0.7126564032704135]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-5.0] * self.N, [5.0] * self.N))
        self.custom_bounds = [(-2, 2), (-1.5, 1.5)]

        self.global_optimum = [(0.08984201368301331, -0.7126564032704135),
                               (-0.08984201368301331, 0.7126564032704135)]
        self.fglob = -1.031628

    def fun(self, x, *args):
        self.nfev += 1
        return ((4 - 2.1 * x[0] ** 2 + x[0] ** 4 / 3) * x[0] ** 2 + x[0] * x[1]
                + (4 * x[1] ** 2 - 4) * x[1] ** 2)


class Sodp(Benchmark):

    r"""
    Sodp objective function.

    This class defines the Sum Of Different Powers [1]_ global optimization
    problem. This is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Sodp}}(x) = \sum_{i=1}^{n} \lvert{x_{i}}\rvert^{i + 1}

    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-1, 1]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x_i = 0` for
    :math:`i = 1, ..., n`

    .. [1] Gavana, A. Global Optimization Benchmarks and AMPGO retrieved 2015
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-1.0] * self.N, [1.0] * self.N))

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        i = arange(1, self.N + 1)
        return sum(abs(x) ** (i + 1))


class Sphere(Benchmark):

    r"""
    Sphere objective function.

    This class defines the Sphere [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Sphere}}(x) = \sum_{i=1}^{n} x_i^2

    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-1, 1]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x_i = 0` for
    :math:`i = 1, ..., n`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    TODO Jamil has stupid limits
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = list(zip([-5.12] * self.N, [5.12] * self.N))

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum(x ** 2)


class Step(Benchmark):

    r"""
    Step objective function.

    This class defines the Step [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Step}}(x) = \sum_{i=1}^{n} \left ( \lfloor x_i
                             + 0.5 \rfloor \right )^2

    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-100, 100]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x_i = 0.5` for
    :math:`i = 1, ..., n`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = list(zip([-100.0] * self.N,
                           [100.0] * self.N))
        self.custom_bounds = ([-5, 5], [-5, 5])

        self.global_optimum = [[0. for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum(floor(abs(x)))


class Step2(Benchmark):

    r"""
    Step objective function.

    This class defines the Step 2 [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Step}}(x) = \sum_{i=1}^{n} \left ( \lfloor x_i
                             + 0.5 \rfloor \right )^2

    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-100, 100]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x_i = 0.5` for
    :math:`i = 1, ..., n`

    .. [1] Gavana, A. Global Optimization Benchmarks and AMPGO retrieved 2015
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = list(zip([-100.0] * self.N,
                           [100.0] * self.N))
        self.custom_bounds = ([-5, 5], [-5, 5])

        self.global_optimum = [[0.5 for _ in range(self.N)]]
        self.fglob = 0.5
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum((floor(x) + 0.5) ** 2.0)


class Stochastic(Benchmark):

    r"""
    Stochastic objective function.

    This class defines the Stochastic [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Stochastic}}(x) = \sum_{i=1}^{n} \epsilon_i 
                                    \left | {x_i - \frac{1}{i}} \right |

    The variable :math:`\epsilon_i, (i=1,...,n)` is a random variable uniformly
    distributed in :math:`[0, 1]`.

    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-5, 5]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x_i = [1/n]` for
    :math:`i = 1, ..., n`

    .. [1] Gavana, A. Global Optimization Benchmarks and AMPGO retrieved 2015
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-5.0] * self.N, [5.0] * self.N))

        self.global_optimum = [[1.0 / _ for _ in range(1, self.N + 1)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        rnd = uniform(0.0, 1.0, size=(self.N, ))
        i = arange(1, self.N + 1)

        return sum(rnd * abs(x - 1.0 / i))


class StretchedV(Benchmark):

    r"""
    StretchedV objective function.

    This class defines the Stretched V [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{StretchedV}}(x) = \sum_{i=1}^{n-1} t^{1/4}
                                   [\sin (50t^{0.1}) + 1]^2

    Where, in this exercise:

    .. math::

       t = x_{i+1}^2 + x_i^2


    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-10, 10]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x = [0., 0.]` when
    :math:`n = 2`.

    .. [1] Adorio, E. MVF - "Multivariate Test Functions Library in C for
    Unconstrained Global Optimization", 2005

    TODO All the sources disagree on the equation, in some the 1 is in the
    brackets, in others it is outside. In Jamil#142 it's not even 1. Here
    we go with the Adorio option.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-10] * self.N, [10] * self.N))

        self.global_optimum = [[0, 0]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        t = x[1:] ** 2 + x[: -1] ** 2
        return sum(t ** 0.25 * (sin(50.0 * t ** 0.1 + 1) ** 2))


class StyblinskiTang(Benchmark):

    r"""
    StyblinskiTang objective function.

    This class defines the Styblinski-Tang [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\text{StyblinskiTang}}(x) = \sum_{i=1}^{n} \left(x_i^4
                                       - 16x_i^2 + 5x_i \right)

    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-5, 5]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = -39.16616570377142n` for
    :math:`x_i = -2.903534018185960` for :math:`i = 1, ..., n`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-5.0] * self.N, [5.0] * self.N))

        self.global_optimum = [[-2.903534018185960 for _ in range(self.N)]]
        self.fglob = -39.16616570377142 * self.N
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum(x ** 4 - 16 * x ** 2 + 5 * x) / 2
