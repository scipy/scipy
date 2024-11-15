import numpy as np
from numpy import abs, cos, exp, arange, pi, sin, sqrt, sum, zeros, tanh
from numpy.testing import assert_almost_equal
from .go_benchmark import Benchmark


class Damavandi(Benchmark):
    r"""
    Damavandi objective function.

    This class defines the Damavandi [1]_ global optimization problem. This is a
    multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Damavandi}}(x) = \left[ 1 - \lvert{\frac{
        \sin[\pi (x_1 - 2)]\sin[\pi (x2 - 2)]}{\pi^2 (x_1 - 2)(x_2 - 2)}}
        \rvert^5 \right] \left[2 + (x_1 - 7)^2 + 2(x_2 - 7)^2 \right]


    Here, :math:`n` represents the number of dimensions and :math:`x_i \in
    [0, 14]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0.0` for :math:`x_i = 2` for
    :math:`i = 1, ..., n`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, 2)

        self._bounds = list(zip([0.0] * self.N, [14.0] * self.N))

        self.global_optimum = [[2 for _ in range(self.N)]]
        self.fglob = np.nan

    def fun(self, x, *args):
        self.nfev += 1

        try:
            num = sin(pi * (x[0] - 2.0)) * sin(pi * (x[1] - 2.0))
            den = (pi ** 2) * (x[0] - 2.0) * (x[1] - 2.0)
            factor1 = 1.0 - (abs(num / den)) ** 5.0
            factor2 = 2 + (x[0] - 7.0) ** 2.0 + 2 * (x[1] - 7.0) ** 2.0
            return factor1 * factor2
        except ZeroDivisionError:
            return np.nan

    def success(self, x):
        """Is a candidate solution at the global minimum"""
        val = self.fun(x)
        if np.isnan(val):
            return True
        try:
            assert_almost_equal(val, 0., 4)
            return True
        except AssertionError:
            return False

        return False


class Deb01(Benchmark):
    r"""
    Deb 1 objective function.

    This class defines the Deb 1 [1]_ global optimization problem. This is a
    multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Deb01}}(x) = - \frac{1}{N} \sum_{i=1}^n \sin^6(5 \pi x_i)


    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-1, 1]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x_i) = 0.0`. The number of global minima is
    :math:`5^n` that are evenly spaced in the function landscape, where
    :math:`n` represents the dimension of the problem.

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.change_dimensionality = True

        self._bounds = list(zip([-1.0] * self.N, [1.0] * self.N))

        self.global_optimum = [[0.3, -0.3]]
        self.fglob = -1.0

    def fun(self, x, *args):
        self.nfev += 1
        return -(1.0 / self.N) * sum(sin(5 * pi * x) ** 6.0)


class Deb03(Benchmark):
    r"""
    Deb 3 objective function.

    This class defines the Deb 3 [1]_ global optimization problem. This is a
    multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Deb03}}(x) = - \frac{1}{N} \sum_{i=1}^n \sin^6 \left[ 5 \pi
        \left ( x_i^{3/4} - 0.05 \right) \right ]


    Here, :math:`n` represents the number of dimensions and :math:`x_i \in
    [0, 1]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0.0`. The number of global minima is
    :math:`5^n` that are evenly spaced in the function landscape, where
    :math:`n` represents the dimension of the problem.

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.change_dimensionality = True

        # lower limit changed to zero because of fractional power
        self._bounds = list(zip([0.0] * self.N, [1.0] * self.N))

        self.global_optimum = [[0.93388314, 0.68141781]]
        self.fglob = -1.0

    def fun(self, x, *args):
        self.nfev += 1

        return -(1.0 / self.N) * sum(sin(5 * pi * (x ** 0.75 - 0.05)) ** 6.0)


class Decanomial(Benchmark):
    r"""
    Decanomial objective function.

    This class defines the Decanomial function global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\text{Decanomial}}(x) = 0.001 \left(\lvert{x_{2}^{4} + 12 x_{2}^{3}
       + 54 x_{2}^{2} + 108 x_{2} + 81.0}\rvert + \lvert{x_{1}^{10}
       - 20 x_{1}^{9} + 180 x_{1}^{8} - 960 x_{1}^{7} + 3360 x_{1}^{6}
       - 8064 x_{1}^{5} + 13340 x_{1}^{4} - 15360 x_{1}^{3} + 11520 x_{1}^{2}
       - 5120 x_{1} + 2624.0}\rvert\right)^{2}


    with :math:`x_i \in [-10, 10]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x = [2, -3]`

    .. [1] Gavana, A. Global Optimization Benchmarks and AMPGO retrieved 2015
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-10.0] * self.N, [10.0] * self.N))
        self.custom_bounds = [(0, 2.5), (-2, -4)]

        self.global_optimum = [[2.0, -3.0]]
        self.fglob = 0.0

    def fun(self, x, *args):

        self.nfev += 1

        val = x[1] ** 4 + 12 * x[1] ** 3 + 54 * x[1] ** 2 + 108 * x[1] + 81.0
        val2 = x[0] ** 10. - 20 * x[0] ** 9 + 180 * x[0] ** 8 - 960 * x[0] ** 7
        val2 += 3360 * x[0] ** 6 - 8064 * x[0] ** 5 + 13340 * x[0] ** 4
        val2 += - 15360 * x[0] ** 3 + 11520 * x[0] ** 2 - 5120 * x[0] + 2624
        return 0.001 * (abs(val) + abs(val2)) ** 2.


class Deceptive(Benchmark):
    r"""
    Deceptive objective function.

    This class defines the Deceptive [1]_ global optimization problem. This is a
    multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Deceptive}}(x) = - \left [\frac{1}{n}
        \sum_{i=1}^{n} g_i(x_i) \right ]^{\beta}


    Where :math:`\beta` is a fixed non-linearity factor; in this exercise,
    :math:`\beta = 2`. The function :math:`g_i(x_i)` is given by:

    .. math::

    g_i(x_i) = \begin{cases}
    - \frac{x}{\alpha_i} + \frac{4}{5} &
    \textrm{if} \hspace{5pt} 0 \leq x_i \leq \frac{4}{5} \alpha_i \\
    \frac{5x}{\alpha_i} -4 &
    \textrm{if} \hspace{5pt} \frac{4}{5} \alpha_i \le x_i \leq \alpha_i \\
    \frac{5(x - \alpha_i)}{\alpha_i-1} &
    \textrm{if} \hspace{5pt} \alpha_i \le x_i \leq \frac{1 + 4\alpha_i}{5} \\
    \frac{x - 1}{1 - \alpha_i} &
    \textrm{if} \hspace{5pt} \frac{1 + 4\alpha_i}{5} \le x_i \leq 1
    \end{cases}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \in
    [0, 1]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = -1` for :math:`x_i = \alpha_i` for
    :math:`i = 1, ..., n`

    .. [1] Gavana, A. Global Optimization Benchmarks and AMPGO retrieved 2015

    TODO: this function was taken from the Gavana website. The following code
    is based on his code.  His code and the website don't match, the equations
    are wrong.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([0.0] * self.N, [1.0] * self.N))

        alpha = arange(1.0, self.N + 1.0) / (self.N + 1.0)

        self.global_optimum = [alpha]
        self.fglob = -1.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        alpha = arange(1.0, self.N + 1.0) / (self.N + 1.0)
        beta = 2.0

        g = zeros((self.N, ))

        for i in range(self.N):
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

        return -((1.0 / self.N) * sum(g)) ** beta


class DeckkersAarts(Benchmark):
    r"""
    Deckkers-Aarts objective function.

    This class defines the Deckkers-Aarts [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{DeckkersAarts}}(x) = 10^5x_1^2 + x_2^2 - (x_1^2 + x_2^2)^2
        + 10^{-5}(x_1^2 + x_2^2)^4


    with :math:`x_i \in [-20, 20]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = -24776.518242168` for
    :math:`x = [0, \pm 14.9451209]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    TODO: jamil solution and global minimum are slightly wrong.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-20.0] * self.N, [20.0] * self.N))
        self.custom_bounds = ([-1, 1], [14, 16])

        self.global_optimum = [[0.0, 14.9451209]]
        self.fglob = -24776.518342168

    def fun(self, x, *args):
        self.nfev += 1
        return (1.e5 * x[0] ** 2 + x[1] ** 2 - (x[0] ** 2 + x[1] ** 2) ** 2
                + 1.e-5 * (x[0] ** 2 + x[1] ** 2) ** 4)


class DeflectedCorrugatedSpring(Benchmark):
    r"""
    DeflectedCorrugatedSpring objective function.

    This class defines the Deflected Corrugated Spring [1]_ function global
    optimization problem. This is a multimodal minimization problem defined as
    follows:

    .. math::

       f_{\text{DeflectedCorrugatedSpring}}(x) = 0.1\sum_{i=1}^n \left[ (x_i -
       \alpha)^2 - \cos \left( K \sqrt {\sum_{i=1}^n (x_i - \alpha)^2}
       \right ) \right ]


    Where, in this exercise, :math:`K = 5` and :math:`\alpha = 5`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \in
    [0, 2\alpha]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = -1` for :math:`x_i = \alpha` for
    :math:`i = 1, ..., n`

    .. [1] Gavana, A. Global Optimization Benchmarks and AMPGO retrieved 2015

    TODO: website has a different equation to the gavana codebase. The function
    below is different to the equation above.  Also, the global minimum is
    wrong.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        alpha = 5.0
        self._bounds = list(zip([0] * self.N, [2 * alpha] * self.N))

        self.global_optimum = [[alpha for _ in range(self.N)]]
        self.fglob = -1.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1
        K, alpha = 5.0, 5.0

        return (-cos(K * sqrt(sum((x - alpha) ** 2)))
                + 0.1 * sum((x - alpha) ** 2))


class DeVilliersGlasser01(Benchmark):
    r"""
    DeVilliers-Glasser 1 objective function.

    This class defines the DeVilliers-Glasser 1 [1]_ function global optimization
    problem. This is a multimodal minimization problem defined as follows:

    .. math::

       f_{\text{DeVilliersGlasser01}}(x) = \sum_{i=1}^{24} \left[ x_1x_2^{t_i}
       \sin(x_3t_i + x_4) - y_i \right ]^2


    Where, in this exercise, :math:`t_i = 0.1(i - 1)` and
    :math:`y_i = 60.137(1.371^{t_i}) \sin(3.112t_i + 1.761)`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \in
    [1, 100]` for :math:`i = 1, ..., 4`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x_i = 0` for
    :math:`x = [60.137, 1.371, 3.112, 1.761]`.

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    """

    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([1.0] * self.N, [100.0] * self.N))

        self.global_optimum = [[60.137, 1.371, 3.112, 1.761]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        t = 0.1 * arange(24)
        y = 60.137 * (1.371 ** t) * sin(3.112 * t + 1.761)

        return sum((x[0] * (x[1] ** t) * sin(x[2] * t + x[3]) - y) ** 2.0)


class DeVilliersGlasser02(Benchmark):
    r"""
    DeVilliers-Glasser 2 objective function.

    This class defines the DeVilliers-Glasser 2 [1]_ function global optimization
    problem. This is a multimodal minimization problem defined as follows:

    .. math::

       f_{\text{DeVilliersGlasser01}}(x) = \sum_{i=1}^{24} \left[ x_1x_2^{t_i}
       \tanh \left [x_3t_i + \sin(x_4t_i) \right] \cos(t_ie^{x_5}) -
       y_i \right ]^2


    Where, in this exercise, :math:`t_i = 0.1(i - 1)` and
    :math:`y_i = 53.81(1.27^{t_i}) \tanh (3.012t_i + \sin(2.13t_i))
    \cos(e^{0.507}t_i)`.

    with :math:`x_i \in [1, 60]` for :math:`i = 1, ..., 5`.

    *Global optimum*: :math:`f(x) = 0` for
    :math:`x = [53.81, 1.27, 3.012, 2.13, 0.507]`.

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=5):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([1.0] * self.N, [60.0] * self.N))

        self.global_optimum = [[53.81, 1.27, 3.012, 2.13, 0.507]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        t = 0.1 * arange(16)
        y = (53.81 * 1.27 ** t * tanh(3.012 * t + sin(2.13 * t))
             * cos(exp(0.507) * t))

        return sum((x[0] * (x[1] ** t) * tanh(x[2] * t + sin(x[3] * t))
                   * cos(t * exp(x[4])) - y) ** 2.0)


class DixonPrice(Benchmark):
    r"""
    Dixon and Price objective function.

    This class defines the Dixon and Price global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{DixonPrice}}(x) = (x_i - 1)^2
        + \sum_{i=2}^n i(2x_i^2 - x_{i-1})^2


    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-10, 10]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x_i) = 0` for
    :math:`x_i = 2^{- \frac{(2^i - 2)}{2^i}}` for :math:`i = 1, ..., n`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    TODO: Gavana code not correct.  i array should start from 2.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-10.0] * self.N, [10.0] * self.N))
        self.custom_bounds = [(-2, 3), (-2, 3)]

        self.global_optimum = [[2.0 ** (-(2.0 ** i - 2.0) / 2.0 ** i)
                               for i in range(1, self.N + 1)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        i = arange(2, self.N + 1)
        s = i * (2.0 * x[1:] ** 2.0 - x[:-1]) ** 2.0
        return sum(s) + (x[0] - 1.0) ** 2.0


class Dolan(Benchmark):
    r"""
    Dolan objective function.

    This class defines the Dolan [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Dolan}}(x) = \lvert (x_1 + 1.7 x_2)\sin(x_1) - 1.5 x_3
        - 0.1 x_4\cos(x_5 + x_5 - x_1) + 0.2 x_5^2 - x_2 - 1 \rvert


    with :math:`x_i \in [-100, 100]` for :math:`i = 1, ..., 5`.

    *Global optimum*: :math:`f(x_i) = 10^{-5}` for
    :math:`x = [8.39045925, 4.81424707, 7.34574133, 68.88246895, 3.85470806]`

    .. [1] Gavana, A. Global Optimization Benchmarks and AMPGO retrieved 2015

    TODO Jamil equation is missing the absolute brackets around the entire
    expression.
    """

    def __init__(self, dimensions=5):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-100.0] * self.N,
                                [100.0] * self.N))

        self.global_optimum = [[-74.10522498, 44.33511286, 6.21069214,
                               18.42772233, -16.5839403]]
        self.fglob = 0

    def fun(self, x, *args):
        self.nfev += 1

        return (abs((x[0] + 1.7 * x[1]) * sin(x[0]) - 1.5 * x[2]
                - 0.1 * x[3] * cos(x[3] + x[4] - x[0]) + 0.2 * x[4] ** 2
                - x[1] - 1))


class DropWave(Benchmark):
    r"""
    DropWave objective function.

    This class defines the DropWave [1]_ global optimization problem. This is a
    multimodal minimization problem defined as follows:

    .. math::

        f_{\text{DropWave}}(x) = - \frac{1 + \cos\left(12 \sqrt{\sum_{i=1}^{n}
        x_i^{2}}\right)}{2 + 0.5 \sum_{i=1}^{n} x_i^{2}}


    with :math:`x_i \in [-5.12, 5.12]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = -1` for :math:`x = [0, 0]`

    .. [1] Gavana, A. Global Optimization Benchmarks and AMPGO retrieved 2015
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-5.12] * self.N, [5.12] * self.N))

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = -1.0

    def fun(self, x, *args):
        self.nfev += 1

        norm_x = sum(x ** 2)
        return -(1 + cos(12 * sqrt(norm_x))) / (0.5 * norm_x + 2)
