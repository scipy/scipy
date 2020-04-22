# -*- coding: utf-8 -*-
from numpy import abs, asarray, cos, exp, arange, pi, sin, sqrt, sum
from .go_benchmark import Benchmark


class Easom(Benchmark):

    r"""
    Easom objective function.

    This class defines the Easom [1]_ global optimization problem. This is a
    a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Easom}}({x}) = a - \frac{a}{e^{b \sqrt{\frac{\sum_{i=1}^{n}
        x_i^{2}}{n}}}} + e - e^{\frac{\sum_{i=1}^{n} \cos\left(c x_i\right)}
        {n}}


    Where, in this exercise, :math:`a = 20, b = 0.2` and :math:`c = 2 \pi`.

    Here, :math:`x_i \in [-100, 100]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x = [0, 0]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    TODO Gavana website disagrees with Jamil, etc. Gavana equation in docstring is totally wrong.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-100.0] * self.N,
                           [100.0] * self.N))

        self.global_optimum = [[pi for _ in range(self.N)]]
        self.fglob = -1.0

    def fun(self, x, *args):
        self.nfev += 1
        a = (x[0] - pi)**2 + (x[1] - pi)**2
        return -cos(x[0]) * cos(x[1]) * exp(-a)


class Eckerle4(Benchmark):
    r"""
    Eckerle4 objective function.
    Eckerle, K., NIST (1979).
    Circular Interference Transmittance Study.

    ..[1] https://www.itl.nist.gov/div898/strd/nls/data/eckerle4.shtml

    #TODO, this is a NIST regression standard dataset, docstring needs
    improving
    """

    def __init__(self, dimensions=3):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([0., 1., 10.],
                           [20, 20., 600.]))
        self.global_optimum = [[1.5543827178, 4.0888321754, 4.5154121844e2]]
        self.fglob = 1.4635887487E-03

        self.a = asarray([1.5750000E-04, 1.6990000E-04, 2.3500000E-04,
                          3.1020000E-04, 4.9170000E-04, 8.7100000E-04,
                          1.7418000E-03, 4.6400000E-03, 6.5895000E-03,
                          9.7302000E-03, 1.4900200E-02, 2.3731000E-02,
                          4.0168300E-02, 7.1255900E-02, 1.2644580E-01,
                          2.0734130E-01, 2.9023660E-01, 3.4456230E-01,
                          3.6980490E-01, 3.6685340E-01, 3.1067270E-01,
                          2.0781540E-01, 1.1643540E-01, 6.1676400E-02,
                          3.3720000E-02, 1.9402300E-02, 1.1783100E-02,
                          7.4357000E-03, 2.2732000E-03, 8.8000000E-04,
                          4.5790000E-04, 2.3450000E-04, 1.5860000E-04,
                          1.1430000E-04, 7.1000000E-05])

        self.b = asarray([4.0000000E+02, 4.0500000E+02, 4.1000000E+02,
                          4.1500000E+02, 4.2000000E+02, 4.2500000E+02,
                          4.3000000E+02, 4.3500000E+02, 4.3650000E+02,
                          4.3800000E+02, 4.3950000E+02, 4.4100000E+02,
                          4.4250000E+02, 4.4400000E+02, 4.4550000E+02,
                          4.4700000E+02, 4.4850000E+02, 4.5000000E+02,
                          4.5150000E+02, 4.5300000E+02, 4.5450000E+02,
                          4.5600000E+02, 4.5750000E+02, 4.5900000E+02,
                          4.6050000E+02, 4.6200000E+02, 4.6350000E+02,
                          4.6500000E+02, 4.7000000E+02, 4.7500000E+02,
                          4.8000000E+02, 4.8500000E+02, 4.9000000E+02,
                          4.9500000E+02, 5.0000000E+02])

    def fun(self, x, *args):
        self.nfev += 1

        vec = x[0] / x[1] * exp(-(self.b - x[2]) ** 2 / (2 * x[1] ** 2))
        return sum((self.a - vec) ** 2)


class EggCrate(Benchmark):

    r"""
    Egg Crate objective function.

    This class defines the Egg Crate [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{EggCrate}}(x) = x_1^2 + x_2^2 + 25 \left[ \sin^2(x_1)
        + \sin^2(x_2) \right]


    with :math:`x_i \in [-5, 5]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x = [0, 0]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-5.0] * self.N, [5.0] * self.N))

        self.global_optimum = [[0.0, 0.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1
        return x[0] ** 2 + x[1] ** 2 + 25 * (sin(x[0]) ** 2 + sin(x[1]) ** 2)


class EggHolder(Benchmark):

    r"""
    Egg Holder [1]_ objective function.

    This class defines the Egg Holder global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{EggHolder}}=\sum_{1}^{n - 1}\left[-\left(x_{i + 1}
        + 47 \right ) \sin\sqrt{\lvert x_{i+1} + x_i/2 + 47 \rvert}
        - x_i \sin\sqrt{\lvert x_i - (x_{i + 1} + 47)\rvert}\right ]


    Here, :math:`n` represents the number of dimensions and :math:`x_i \in
    [-512, 512]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = -959.640662711` for
    :math:`{x} = [512, 404.2319]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    TODO: Jamil is missing a minus sign on the fglob value
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-512.1] * self.N,
                           [512.0] * self.N))

        self.global_optimum = [[512.0, 404.2319]]
        self.fglob = -959.640662711
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        vec = (-(x[1:] + 47) * sin(sqrt(abs(x[1:] + x[:-1] / 2. + 47)))
               - x[:-1] * sin(sqrt(abs(x[:-1] - (x[1:] + 47)))))
        return sum(vec)


class ElAttarVidyasagarDutta(Benchmark):

    r"""
    El-Attar-Vidyasagar-Dutta [1]_ objective function.

    This class defines the El-Attar-Vidyasagar-Dutta function global
    optimization problem. This is a multimodal minimization problem defined as
    follows:

    .. math::

       f_{\text{ElAttarVidyasagarDutta}}(x) = (x_1^2 + x_2 - 10)^2
       + (x_1 + x_2^2 - 7)^2 + (x_1^2 + x_2^3 - 1)^2


    with :math:`x_i \in [-100, 100]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = 1.712780354` for
    :math:`x= [3.40918683, -2.17143304]`

    .. [1] Gavana, A. Global Optimization Benchmarks and AMPGO retrieved 2015
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-100.0] * self.N,
                           [100.0] * self.N))
        self.custom_bounds = [(-4, 4), (-4, 4)]

        self.global_optimum = [[3.40918683, -2.17143304]]
        self.fglob = 1.712780354

    def fun(self, x, *args):
        self.nfev += 1

        return ((x[0] ** 2 + x[1] - 10) ** 2 + (x[0] + x[1] ** 2 - 7) ** 2
                + (x[0] ** 2 + x[1] ** 3 - 1) ** 2)


class Exp2(Benchmark):

    r"""
    Exp2 objective function.

    This class defines the Exp2 global optimization problem. This is a
    multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Exp2}}(x) = \sum_{i=0}^9 \left ( e^{-ix_1/10} - 5e^{-ix_2/10}
        - e^{-i/10} + 5e^{-i} \right )^2


    with :math:`x_i \in [0, 20]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x = [1, 10.]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([0.0] * self.N, [20.0] * self.N))
        self.custom_bounds = [(0, 2), (0, 20)]

        self.global_optimum = [[1.0, 10.]]
        self.fglob = 0.

    def fun(self, x, *args):
        self.nfev += 1

        i = arange(10.)
        vec = (exp(-i * x[0] / 10.) - 5 * exp(-i * x[1] / 10.) - exp(-i / 10.)
               + 5 * exp(-i)) ** 2

        return sum(vec)


class Exponential(Benchmark):

    r"""
    Exponential [1] objective function.

    This class defines the Exponential global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Exponential}}(x) = -e^{-0.5 \sum_{i=1}^n x_i^2}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \in
    [-1, 1]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x_i) = -1` for :math:`x_i = 0` for
    :math:`i = 1, ..., n`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    TODO Jamil are missing a minus sign on fglob
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-1.0] * self.N, [1.0] * self.N))

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = -1.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return -exp(-0.5 * sum(x ** 2.0))
