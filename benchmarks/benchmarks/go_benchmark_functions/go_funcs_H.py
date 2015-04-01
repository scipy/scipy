# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import

import numpy as np
from numpy import (abs, arctan2, asarray, cos, exp, floor, log, log10,
                   arange, pi, sign, sin, sqrt, sum, tan, tanh)
from .go_benchmark import Benchmark


class Hansen(Benchmark):

    r"""
    Hansen objective function.

    This class defines the Hansen [1]_ global optimization problem. This is a
    multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Hansen}}(x) = \left[ \sum_{i=0}^4(i+1)\cos(ix_1+i+1)\right ]
        \left[\sum_{j=0}^4(j+1)\cos[(j+2)x_2+j+1])\right ]


    with :math:`x_i \in [-10, 10]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = -176.54179` for
    :math:`x = [-7.58989583, -7.70831466]`.

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    TODO Jamil #61 is missing the starting value of i.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[-7.58989583, -7.70831466]]
        self.fglob = -176.54179

    def fun(self, x, *args):
        self.nfev += 1

        i = arange(5.)
        a = (i + 1) * cos(i * x[0] + i + 1)
        b = (i + 1) * cos((i + 2) * x[1] + i + 1)

        return sum(a) * sum(b)


class Hartmann3(Benchmark):

    r"""
    Hartmann3 objective function.

    This class defines the Hartmann3 [1]_ global optimization problem. This is a
    multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Hartmann3}}(x) = -\sum\limits_{i=1}^{4} c_i 
        e^{-\sum\limits_{j=1}^{n}a_{ij}(x_j - p_{ij})^2}


    Where, in this exercise:

    .. math::

        \begin{array}{l|ccc|c|ccr}
        \hline
        i & & a_{ij}&  & c_i & & p_{ij} &  \\
        \hline
        1 & 3.0 & 10.0 & 30.0 & 1.0 & 0.3689  & 0.1170 & 0.2673 \\
        2 & 0.1 & 10.0 & 35.0 & 1.2 & 0.4699 & 0.4387 & 0.7470 \\
        3 & 3.0 & 10.0 & 30.0 & 3.0 & 0.1091 & 0.8732 & 0.5547 \\
        4 & 0.1 & 10.0 & 35.0 & 3.2 & 0.03815 & 0.5743 & 0.8828 \\
        \hline
        \end{array}


    with :math:`x_i \in [0, 1]` for :math:`i = 1, 2, 3`.

    *Global optimum*: :math:`f(x) = -3.8627821478` 
    for :math:`x = [0.11461292,  0.55564907,  0.85254697]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    TODO Jamil #62 has an incorrect coefficient. p[1, 1] should be 0.4387
    """

    def __init__(self, dimensions=3):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N, [1.0] * self.N)

        self.global_optimum = [[0.11461292, 0.55564907, 0.85254697]]
        self.fglob = -3.8627821478

        self.a = asarray([[  3. ,  10. ,  30. ],
                          [  0.1,  10. ,  35. ],
                          [  3. ,  10. ,  30. ],
                          [  0.1,  10. ,  35. ]])

        self.p = asarray([[ 0.3689 ,  0.117  ,  0.2673 ],
                          [ 0.4699 ,  0.4387 ,  0.747  ],
                          [ 0.1091 ,  0.8732 ,  0.5547 ],
                          [ 0.03815,  0.5743 ,  0.8828 ]])

        self.c = asarray([1., 1.2, 3., 3.2])

    def fun(self, x, *args):
        self.nfev += 1

        XX = np.atleast_2d(x)
        d = sum(self.a * (XX - self.p) ** 2, axis=1)
        return -sum(self.c * exp(-d))


class Hartmann6(Benchmark):

    r"""
    Hartmann6 objective function.

    This class defines the Hartmann6 [1]_ global optimization problem. This is a
    multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Hartmann6}}(x) = -\sum\limits_{i=1}^{4} c_i
        e^{-\sum\limits_{j=1}^{n}a_{ij}(x_j - p_{ij})^2}


    Where, in this exercise:

    .. math::

        \begin{array}{l|cccccc|r}
        \hline
        i & &   &   a_{ij} &  &  & & c_i  \\
        \hline
        1 & 10.0  & 3.0  & 17.0 & 3.50  & 1.70  & 8.00  & 1.0 \\
        2 & 0.05  & 10.0 & 17.0 & 0.10  & 8.00  & 14.00 & 1.2 \\
        3 & 3.00  & 3.50 & 1.70 & 10.0  & 17.00 & 8.00  & 3.0 \\
        4 & 17.00 & 8.00 & 0.05 & 10.00 & 0.10  & 14.00 & 3.2 \\
        \hline
        \end{array}

        \newline
        \
        \newline

        \begin{array}{l|cccccr}
        \hline
        i &  &   & p_{ij} &  & & \\
        \hline
        1 & 0.1312 & 0.1696 & 0.5569 & 0.0124 & 0.8283 & 0.5886 \\
        2 & 0.2329 & 0.4135 & 0.8307 & 0.3736 & 0.1004 & 0.9991 \\
        3 & 0.2348 & 0.1451 & 0.3522 & 0.2883 & 0.3047 & 0.6650 \\
        4 & 0.4047 & 0.8828 & 0.8732 & 0.5743 & 0.1091 & 0.0381 \\
        \hline
        \end{array}


    with :math:`x_i \in [0, 1]` for :math:`i = 1, ..., 6`.

    *Global optimum*: :math:`f(x_i) = -3.32236801141551` for
    :math:`{x} = [0.20168952, 0.15001069, 0.47687398, 0.27533243, 0.31165162,
    0.65730054]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=6):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N, [1.0] * self.N)

        self.global_optimum = [[0.20168952, 0.15001069, 0.47687398, 0.27533243,
                                0.31165162, 0.65730054]]

        self.fglob = -3.32236801141551

        self.a = asarray([[ 10.  ,   3.  ,  17.  ,   3.5 ,   1.7 ,   8.  ],
                          [  0.05,  10.  ,  17.  ,   0.1 ,   8.  ,  14.  ],
                          [  3.  ,   3.5 ,   1.7 ,  10.  ,  17.  ,   8.  ],
                          [ 17.  ,   8.  ,   0.05,  10.  ,   0.1 ,  14.  ]])

        self.p = asarray([[0.1312, 0.1696, 0.5569, 0.0124, 0.8283, 0.5886],
                          [0.2329, 0.4135, 0.8307, 0.3736, 0.1004, 0.9991],
                          [0.2348, 0.1451, 0.3522, 0.2883, 0.3047, 0.665 ],
                          [0.4047, 0.8828, 0.8732, 0.5743, 0.1091, 0.0381]])

        self.c = asarray([1.0, 1.2, 3.0, 3.2])

    def fun(self, x, *args):
        self.nfev += 1

        XX = np.atleast_2d(x)
        d = sum(self.a * (XX - self.p) ** 2, axis=1)
        return -sum(self.c * exp(-d))


class HelicalValley(Benchmark):

    r"""
    HelicalValley objective function.

    This class defines the HelicalValley [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{HelicalValley}}({x}) = 100{[z-10\Psi(x_1,x_2)]^2
        +(\sqrt{x_1^2+x_2^2}-1)^2}+x_3^2


    Where, in this exercise:

    .. math::

        2\pi\Psi(x,y) =  \begin{cases} \arctan(y/x) & \textrm{for} x > 0 \\
        \pi + \arctan(y/x) & \textrm{for } x < 0 \end{cases}

    with :math:`x_i \in [-100, 100]` for :math:`i = 1, 2, 3`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x = [1, 0, 0]`

    .. [1] Fletcher, R. & Powell, M. A Rapidly Convergent Descent Method for
    Minimzation, Computer Journal, 1963, 62, 163-168

    TODO: Jamil equation is different to original reference. The above paper
    can be obtained from
    http://galton.uchicago.edu/~lekheng/courses/302/classics/
    fletcher-powell.pdf
    """

    def __init__(self, dimensions=3):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.] * self.N, [10.] * self.N)

        self.global_optimum = [[1.0, 0.0, 0.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        r = sqrt(x[0] ** 2 + x[1] ** 2)
        theta = 1 / (2. * pi) * arctan2(x[1], x[0])

        return x[2] ** 2 + 100 * ((x[2] - 10 * theta) ** 2 + (r - 1) ** 2)


class HimmelBlau(Benchmark):

    r"""
    HimmelBlau objective function.

    This class defines the HimmelBlau [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{HimmelBlau}}({x}) = (x_1^2 + x_2 - 11)^2 + (x_1 + x_2^2 - 7)^2


    with :math:`x_i \in [-6, 6]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x = [3, 2]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.] * self.N, [5.] * self.N)

        self.global_optimum = [[3.0, 2.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return (x[0] ** 2 + x[1] - 11) ** 2 + (x[0] + x[1] ** 2 - 7) ** 2


class HolderTable(Benchmark):

    r"""
    HolderTable objective function.

    This class defines the HolderTable [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{HolderTable}}({x}) = - \left|{e^{\left|{1
        - \frac{\sqrt{x_{1}^{2} + x_{2}^{2}}}{\pi} }\right|}
        \sin\left(x_{1}\right) \cos\left(x_{2}\right)}\right|


    with :math:`x_i \in [-10, 10]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = -19.20850256788675` for
    :math:`x_i = \pm 9.664590028909654` for :math:`i = 1, 2`

    .. [1] Gavana, A. Global Optimization Benchmarks and AMPGO retrieved 2015

    TODO: Jamil #146 equation is wrong - should be squaring the x1 and x2
    terms, but isn't. Gavana does.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [(8.055023472141116, 9.664590028909654),
                               (-8.055023472141116, 9.664590028909654),
                               (8.055023472141116, -9.664590028909654),
                               (-8.055023472141116, -9.664590028909654)]
        self.fglob = -19.20850256788675

    def fun(self, x, *args):
        self.nfev += 1

        return -abs(sin(x[0]) * cos(x[1])
                    * exp(abs(1 - sqrt(x[0] ** 2 + x[1] ** 2) / pi)))


class Hosaki(Benchmark):

    r"""
    Hosaki objective function.

    This class defines the Hosaki [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Hosaki}}(x) = \left ( 1 - 8 x_1 + 7 x_1^2 - \frac{7}{3} x_1^3
        + \frac{1}{4} x_1^4 \right ) x_2^2 e^{-x_1}


    with :math:`x_i \in [0, 10]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = -2.3458115` for :math:`x = [4, 2]`.

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = ([0., 5.], [0., 6.])
        self.custom_bounds = [(0, 5), (0, 5)]

        self.global_optimum = [[4, 2]]
        self.fglob = -2.3458115

    def fun(self, x, *args):
        self.nfev += 1

        val = (1 - 8 * x[0] + 7 * x[0] ** 2 - 7 / 3. * x[0] ** 3
               + 0.25 * x[0] ** 4)
        return val * x[1] ** 2 * exp(-x[1])
