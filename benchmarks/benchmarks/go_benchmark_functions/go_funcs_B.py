# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import

from numpy import abs, cos, exp, log, arange, pi, roll, sin, sqrt, sum
from .go_benchmark import Benchmark


class BartelsConn(Benchmark):

    r"""
    Bartels-Conn objective function.

    The BartelsConn [1]_ global optimization problem is a multimodal
    minimization problem defined as follows:

    .. math::

        f_{\text{BartelsConn}}(x) = \lvert {x_1^2 + x_2^2 + x_1x_2} \rvert +
         \lvert {\sin(x_1)} \rvert + \lvert {\cos(x_2)} \rvert


    with :math:`x_i \in [-500, 500]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = 1` for :math:`x = [0, 0]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-500.] * self.N, [500.] * self.N)
        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 1.0

    def fun(self, x, *args):
        self.nfev += 1

        return (abs(x[0] ** 2.0 + x[1] ** 2.0 + x[0] * x[1]) + abs(sin(x[0]))
                + abs(cos(x[1])))


class Beale(Benchmark):

    r"""
    Beale objective function.

    The Beale [1]_ global optimization problem is a multimodal
    minimization problem defined as follows:

    .. math::

        f_{\text{Beale}}(x) = \left(x_1 x_2 - x_1 + 1.5\right)^{2} +
        \left(x_1 x_2^{2} - x_1 + 2.25\right)^{2} + \left(x_1 x_2^{3} - x_1 +
        2.625\right)^{2}


    with :math:`x_i \in [-4.5, 4.5]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x=[3, 0.5]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-4.5] * self.N, [4.5] * self.N)
        self.global_optimum = [[3.0, 0.5]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return ((1.5 - x[0] + x[0] * x[1]) ** 2
                + (2.25 - x[0] + x[0] * x[1] ** 2) ** 2
                + (2.625 - x[0] + x[0] * x[1] ** 3) ** 2)


class BiggsExp02(Benchmark):

    r"""
    BiggsExp02 objective function.

    The BiggsExp02 [1]_ global optimization problem is a multimodal minimization
    problem defined as follows

    .. math::

        \begin{matrix}
        f_{\text{BiggsExp02}}(x) = \sum_{i=1}^{10} (e^{-t_i x_1}
           - 5 e^{-t_i x_2} - y_i)^2 \\
        t_i = 0.1 i\\
        y_i = e^{-t_i} - 5 e^{-10t_i}\\
        \end{matrix}


    with :math:`x_i \in [0, 20]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x = [1, 10]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0] * 2,
                           [20] * 2)
        self.global_optimum = [[1., 10.]]
        self.fglob = 0

    def fun(self, x, *args):
        self.nfev += 1

        t = arange(1, 11.) * 0.1
        y = exp(-t) - 5 * exp(-10 * t)
        vec = (exp(-t * x[0]) - 5 * exp(-t * x[1]) - y) ** 2

        return sum(vec)


class BiggsExp03(Benchmark):

    r"""
    BiggsExp03 objective function.

    The BiggsExp03 [1]_ global optimization problem is a multimodal minimization
    problem defined as follows

    .. math::

        \begin{matrix}\ f_{\text{BiggsExp03}}(x) = \sum_{i=1}^{10}
        (e^{-t_i x_1} - x_3e^{-t_i x_2} - y_i)^2\\
        t_i = 0.1i\\
        y_i = e^{-t_i} - 5e^{-10 t_i}\\
        \end{matrix}


    with :math:`x_i \in [0, 20]` for :math:`i = 1, 2, 3`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x = [1, 10, 5]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    """

    def __init__(self, dimensions=3):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0] * 3,
                           [20] * 3)
        self.global_optimum = [[1., 10., 5.]]
        self.fglob = 0

    def fun(self, x, *args):
        self.nfev += 1

        t = arange(1., 11.) * 0.1
        y = exp(-t) - 5 * exp(-10 * t)
        vec = (exp(-t * x[0]) - x[2] * exp(-t * x[1]) - y) ** 2

        return sum(vec)


class BiggsExp04(Benchmark):

    r"""
    BiggsExp04 objective function.

    The BiggsExp04 [1]_ global optimization problem is a multimodal
    minimization problem defined as follows

    .. math::

        \begin{matrix}\ f_{\text{BiggsExp04}}(x) = \sum_{i=1}^{10}
        (x_3 e^{-t_i x_1} - x_4 e^{-t_i x_2} - y_i)^2\\
        t_i = 0.1i\\
        y_i = e^{-t_i} - 5 e^{-10 t_i}\\
        \end{matrix}


    with :math:`x_i \in [0, 20]` for :math:`i = 1, ..., 4`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x = [1, 10, 1, 5]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    """

    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.] * 4,
                           [20.] * 4)
        self.global_optimum = [[1., 10., 1., 5.]]
        self.fglob = 0

    def fun(self, x, *args):
        self.nfev += 1

        t = arange(1, 11.) * 0.1
        y = exp(-t) - 5 * exp(-10 * t)
        vec = (x[2] * exp(-t * x[0]) - x[3] * exp(-t * x[1]) - y) ** 2

        return sum(vec)


class BiggsExp05(Benchmark):

    r"""
    BiggsExp05 objective function.

    The BiggsExp05 [1]_ global optimization problem is a multimodal minimization
    problem defined as follows

    .. math::

        \begin{matrix}\ f_{\text{BiggsExp05}}(x) = \sum_{i=1}^{11}
        (x_3 e^{-t_i x_1} - x_4 e^{-t_i x_2} + 3 e^{-t_i x_5} - y_i)^2\\
        t_i = 0.1i\\
        y_i = e^{-t_i} - 5e^{-10 t_i} + 3e^{-4 t_i}\\
        \end{matrix}


    with :math:`x_i \in [0, 20]` for :math:`i=1, ..., 5`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x = [1, 10, 1, 5, 4]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    """

    def __init__(self, dimensions=5):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.] * 5,
                           [20.] * 5)
        self.global_optimum = [[1., 10., 1., 5., 4.]]
        self.fglob = 0

    def fun(self, x, *args):
        self.nfev += 1
        t = arange(1, 12.) * 0.1
        y = exp(-t) - 5 * exp(-10 * t) + 3 * exp(-4 * t)
        vec = (x[2] * exp(-t * x[0]) - x[3] * exp(-t * x[1])
               + 3 * exp(-t * x[4]) - y) ** 2

        return sum(vec)


class Bird(Benchmark):

    r"""
    Bird objective function.

    The Bird global optimization problem is a multimodal minimization
    problem defined as follows

    .. math::

        f_{\text{Bird}}(x) = \left(x_1 - x_2\right)^{2} + e^{\left[1 -
         \sin\left(x_1\right) \right]^{2}} \cos\left(x_2\right) + e^{\left[1 -
          \cos\left(x_2\right)\right]^{2}} \sin\left(x_1\right)


    with :math:`x_i \in [-2\pi, 2\pi]`

    *Global optimum*: :math:`f(x) = -106.7645367198034` for :math:`x
    = [4.701055751981055, 3.152946019601391]` or :math:`x =
    [-1.582142172055011, -3.130246799635430]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-2.0 * pi] * self.N,
                           [2.0 * pi] * self.N)
        self.global_optimum = [[4.701055751981055, 3.152946019601391],
                               [-1.582142172055011, -3.130246799635430]]
        self.fglob = -106.7645367198034

    def fun(self, x, *args):
        self.nfev += 1

        return (sin(x[0]) * exp((1 - cos(x[1])) ** 2)
                + cos(x[1]) * exp((1 - sin(x[0])) ** 2) + (x[0] - x[1]) ** 2)


class Bohachevsky1(Benchmark):

    r"""
    Bohachevsky 1 objective function.

    The Bohachevsky 1 [1]_ global optimization problem is a multimodal
    minimization problem defined as follows

        .. math::

        f_{\text{Bohachevsky}}(x) = \sum_{i=1}^{n-1}\left[x_i^2 + 2 x_{i+1}^2 -
        0.3 \cos(3 \pi x_i) - 0.4 \cos(4 \pi x_{i + 1}) + 0.7 \right]


    Here, :math:`n` represents the number of dimensions and :math:`x_i \in
    [-15, 15]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x_i = 0` for :math:`i = 1,
    ..., n`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    TODO: equation needs to be fixed up in the docstring. see Jamil#17
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-100.0] * self.N, [100.0] * self.N)
        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return (x[0] ** 2 + 2 * x[1] ** 2 - 0.3 * cos(3 * pi * x[0])
                - 0.4 * cos(4 * pi * x[1]) + 0.7)


class Bohachevsky2(Benchmark):

    r"""
    Bohachevsky 2 objective function.

    The Bohachevsky 2 [1]_ global optimization problem is a multimodal
    minimization problem defined as follows

        .. math::

        f_{\text{Bohachevsky}}(x) = \sum_{i=1}^{n-1}\left[x_i^2 + 2 x_{i+1}^2 -
        0.3 \cos(3 \pi x_i) - 0.4 \cos(4 \pi x_{i + 1}) + 0.7 \right]


    Here, :math:`n` represents the number of dimensions and :math:`x_i \in
    [-15, 15]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x_i = 0` for :math:`i = 1,
    ..., n`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    TODO: equation needs to be fixed up in the docstring. Jamil is also wrong.
    There should be no 0.4 factor in front of the cos term
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-100.0] * self.N, [100.0] * self.N)
        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return (x[0] ** 2 + 2 * x[1] ** 2 - 0.3 * cos(3 * pi * x[0])
                 * cos(4 * pi * x[1]) + 0.3)


class Bohachevsky3(Benchmark):

    r"""
    Bohachevsky 3 objective function.

    The Bohachevsky 3 [1]_ global optimization problem is a multimodal
    minimization problem defined as follows

        .. math::

        f_{\text{Bohachevsky}}(x) = \sum_{i=1}^{n-1}\left[x_i^2 + 2 x_{i+1}^2 -
        0.3 \cos(3 \pi x_i) - 0.4 \cos(4 \pi x_{i + 1}) + 0.7 \right]


    Here, :math:`n` represents the number of dimensions and :math:`x_i \in
    [-15, 15]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x_i = 0` for :math:`i = 1,
    ..., n`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    TODO: equation needs to be fixed up in the docstring. Jamil#19
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-100.0] * self.N, [100.0] * self.N)
        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return (x[0] ** 2 + 2 * x[1] ** 2
                - 0.3 * cos(3 * pi * x[0] + 4 * pi * x[1]) + 0.3)


class BoxBetts(Benchmark):

    r"""
    BoxBetts objective function.

    The BoxBetts global optimization problem is a multimodal
    minimization problem defined as follows

    .. math::

        f_{\text{BoxBetts}}(x) = \sum_{i=1}^k g(x_i)^2


    Where, in this exercise:

    .. math::

        g(x) = e^{-0.1i x_1} - e^{-0.1i x_2} - x_3\left[e^{-0.1i}
        - e^{-i}\right]


    And :math:`k = 10`.

    Here, :math:`x_1 \in [0.9, 1.2], x_2 \in [9, 11.2], x_3 \in [0.9, 1.2]`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x = [1, 10, 1]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=3):
        Benchmark.__init__(self, dimensions)

        self._bounds = ([0.9, 1.2], [9.0, 11.2], [0.9, 1.2])
        self.global_optimum = [[1.0, 10.0, 1.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        i = arange(1, 11)
        g = (exp(-0.1 * i * x[0]) - exp(-0.1 * i * x[1])
             - (exp(-0.1 * i) - exp(-i)) * x[2])
        return sum(g**2)


class Branin01(Benchmark):

    r"""
    Branin01  objective function.

    The Branin01 global optimization problem is a multimodal minimization
    problem defined as follows

    .. math::

        f_{\text{Branin01}}(x) = \left(- 1.275 \frac{x_1^{2}}{\pi^{2}} + 5
        \frac{x_1}{\pi} + x_2 -6\right)^{2} + \left(10 -\frac{5}{4 \pi} \right)
        \cos\left(x_1\right) + 10


    with :math:`x_1 \in [-5, 10], x_2 \in [0, 15]`

    *Global optimum*: :math:`f(x) = 0.39788735772973816` for :math:`x =
    [-\pi, 12.275]` or :math:`x = [\pi, 2.275]` or :math:`x = [3\pi, 2.475]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    TODO: Jamil#22, one of the solutions is different
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = [(-5., 10.), (0., 15.)]

        self.global_optimum = [[-pi, 12.275], [pi, 2.275], [3 * pi, 2.475]]
        self.fglob = 0.39788735772973816

    def fun(self, x, *args):
        self.nfev += 1

        return ((x[1] - (5.1 / (4 * pi ** 2)) * x[0] ** 2
                + 5 * x[0] / pi - 6) ** 2
                + 10 * (1 - 1 / (8 * pi)) * cos(x[0]) + 10)


class Branin02(Benchmark):

    r"""
    Branin02 objective function.

    The Branin02 global optimization problem is a multimodal minimization
    problem defined as follows

    .. math::

        f_{\text{Branin02}}(x) = \left(- 1.275 \frac{x_1^{2}}{\pi^{2}}
        + 5 \frac{x_1}{\pi} + x_2 - 6 \right)^{2} + \left(10 - \frac{5}{4 \pi}
        \right) \cos\left(x_1\right) \cos\left(x_2\right)
        + \log(x_1^2+x_2^2 + 1) + 10


    with :math:`x_i \in [-5, 15]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = 5.559037` for :math:`x = [-3.2, 12.53]`

    .. [1] Gavana, A. Global Optimization Benchmarks and AMPGO retrieved 2015
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = [(-5.0, 15.0), (-5.0, 15.0)]

        self.global_optimum = [[-3.1969884, 12.52625787]]
        self.fglob = 5.5589144038938247

    def fun(self, x, *args):
        self.nfev += 1

        return ((x[1] - (5.1 / (4 * pi ** 2)) * x[0] ** 2
                + 5 * x[0] / pi - 6) ** 2
                + 10 * (1 - 1 / (8 * pi)) * cos(x[0]) * cos(x[1])
                + log(x[0] ** 2.0 + x[1] ** 2.0 + 1.0) + 10)


class Brent(Benchmark):

    r"""
    Brent objective function.

    The Brent [1]_ global optimization problem is a multimodal minimization
    problem defined as follows:

    .. math::

        f_{\text{Brent}}(x) = (x_1 + 10)^2 + (x_2 + 10)^2 + e^{(-x_1^2 -x_2^2)}


    with :math:`x_i \in [-10, 10]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x = [-10, -10]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    TODO solution is different to Jamil#24
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.custom_bounds = ([-10, 2], [-10, 2])

        self.global_optimum = [[-10.0, -10.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1
        return ((x[0] + 10.0) ** 2.0 + (x[1] + 10.0) ** 2.0
                + exp(-x[0] ** 2.0 - x[1] ** 2.0))


class Brown(Benchmark):

    r"""
    Brown objective function.

    The Brown [1]_ global optimization problem is a multimodal minimization
    problem defined as follows:

    .. math::

        f_{\text{Brown}}(x) = \sum_{i=1}^{n-1}\left[
        \left(x_i^2\right)^{x_{i + 1}^2 + 1}
        + \left(x_{i + 1}^2\right)^{x_i^2 + 1}\right]


    with :math:`x_i \in [-1, 4]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for
    :math:`i=1,...,n`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-1.0] * self.N, [4.0] * self.N)
        self.custom_bounds = ([-1.0, 1.0], [-1.0, 1.0])

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        x0 = x[:-1]
        x1 = x[1:]
        return sum((x0 ** 2.0) ** (x1 ** 2.0 + 1.0)
                   + (x1 ** 2.0) ** (x0 ** 2.0 + 1.0))


class Bukin02(Benchmark):

    r"""
    Bukin02 objective function.

    The Bukin02 [1]_ global optimization problem is a multimodal minimization
    problem defined as follows:

    .. math::

        f_{\text{Bukin02}}(x) = 100 (x_2^2 - 0.01x_1^2 + 1)
        + 0.01(x_1 + 10)^2


    with :math:`x_1 \in [-15, -5], x_2 \in [-3, 3]`

    *Global optimum*: :math:`f(x) = -124.75` for :math:`x = [-15, 0]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    TODO: I think that Gavana and Jamil are wrong on this function. In both
    sources the x[1] term is not squared. As such there will be a minimum at
    the smallest value of x[1].
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = [(-15.0, -5.0), (-3.0, 3.0)]

        self.global_optimum = [[-15.0, 0.0]]
        self.fglob = -124.75

    def fun(self, x, *args):

        self.nfev += 1
        return (100 * (x[1] ** 2 - 0.01 * x[0] ** 2 + 1.0)
                + 0.01 * (x[0] + 10.0) ** 2.0)


class Bukin04(Benchmark):

    r"""
    Bukin04 objective function.

    The Bukin04 [1]_ global optimization problem is a multimodal minimization
    problem defined as follows:

    .. math::

        f_{\text{Bukin04}}(x) = 100 x_2^{2} + 0.01 \lvert{x_1 + 10}
        \rvert


    with :math:`x_1 \in [-15, -5], x_2 \in [-3, 3]`

    *Global optimum*: :math:`f(x) = 0` for :math:`x = [-10, 0]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = [(-15.0, -5.0), (-3.0, 3.0)]

        self.global_optimum = [[-10.0, 0.0]]
        self.fglob = 0.0

    def fun(self, x, *args):

        self.nfev += 1
        return 100 * x[1] ** 2 + 0.01 * abs(x[0] + 10)


class Bukin06(Benchmark):

    r"""
    Bukin06 objective function.

    The Bukin06 [1]_ global optimization problem is a multimodal minimization
    problem defined as follows:

    .. math::

        f_{\text{Bukin06}}(x) = 100 \sqrt{ \lvert{x_2 - 0.01 x_1^{2}}
        \rvert} + 0.01 \lvert{x_1 + 10} \rvert


    with :math:`x_1 \in [-15, -5], x_2 \in [-3, 3]`

    *Global optimum*: :math:`f(x) = 0` for :math:`x = [-10, 1]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = [(-15.0, -5.0), (-3.0, 3.0)]
        self.global_optimum = [[-10.0, 1.0]]
        self.fglob = 0.0

    def fun(self, x, *args):

        self.nfev += 1
        return 100 * sqrt(abs(x[1] - 0.01 * x[0] ** 2)) + 0.01 * abs(x[0] + 10)
