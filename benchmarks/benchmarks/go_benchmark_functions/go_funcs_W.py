from numpy import atleast_2d, arange, sum, cos, exp, pi
from .go_benchmark import Benchmark


class Watson(Benchmark):

    r"""
    Watson objective function.

    This class defines the Watson [1]_ global optimization problem. This is a
    unimodal minimization problem defined as follows:

    .. math::

        f_{\text{Watson}}(x) = \sum_{i=0}^{29} \left\{
                               \sum_{j=0}^4 ((j + 1)a_i^j x_{j+1})
                               - \left[ \sum_{j=0}^5 a_i^j
                               x_{j+1} \right ]^2 - 1 \right\}^2
                               + x_1^2


    Where, in this exercise, :math:`a_i = i/29`.

    with :math:`x_i \in [-5, 5]` for :math:`i = 1, ..., 6`.

    *Global optimum*: :math:`f(x) = 0.002288` for
    :math:`x = [-0.0158, 1.012, -0.2329, 1.260, -1.513, 0.9928]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    TODO Jamil #161 writes equation using (j - 1).  According to code in Adorio
    and Gavana it should be (j+1). However the equations in those papers
    contain (j - 1) as well.  However, I've got the right global minimum!!!
    """

    def __init__(self, dimensions=6):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-5.0] * self.N, [5.0] * self.N))

        self.global_optimum = [[-0.0158, 1.012, -0.2329, 1.260, -1.513,
                                0.9928]]
        self.fglob = 0.002288

    def fun(self, x, *args):
        self.nfev += 1

        i = atleast_2d(arange(30.)).T
        a = i / 29.
        j = arange(5.)
        k = arange(6.)

        t1 = sum((j + 1) * a ** j * x[1:], axis=1)
        t2 = sum(a ** k * x, axis=1)

        inner = (t1 - t2 ** 2 - 1) ** 2

        return sum(inner) + x[0] ** 2


class Wavy(Benchmark):

    r"""
    Wavy objective function.

    This class defines the W / Wavy [1]_ global optimization problem. This is a
    multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Wavy}}(x) = 1 - \frac{1}{n} \sum_{i=1}^{n}
                             \cos(kx_i)e^{-\frac{x_i^2}{2}}


    Where, in this exercise, :math:`k = 10`. The number of local minima is
    :math:`kn` and :math:`(k + 1)n` for odd and even :math:`k` respectively.

    Here, :math:`x_i \in [-\pi, \pi]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x = [0, 0]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-pi] * self.N, [pi] * self.N))

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return 1.0 - (1.0 / self.N) * sum(cos(10 * x) * exp(-x ** 2.0 / 2.0))


class WayburnSeader01(Benchmark):

    r"""
    Wayburn and Seader 1 objective function.

    This class defines the Wayburn and Seader 1 [1]_ global optimization
    problem. This is a unimodal minimization problem defined as follows:

    .. math::

        f_{\text{WayburnSeader01}}(x) = (x_1^6 + x_2^4 - 17)^2
                                        + (2x_1 + x_2 - 4)^2


    with :math:`x_i \in [-5, 5]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x = [1, 2]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-5.0] * self.N, [5.0] * self.N))
        self.custom_bounds = ([-2, 2], [-2, 2])

        self.global_optimum = [[1.0, 2.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return (x[0] ** 6 + x[1] ** 4 - 17) ** 2 + (2 * x[0] + x[1] - 4) ** 2


class WayburnSeader02(Benchmark):

    r"""
    Wayburn and Seader 2 objective function.

    This class defines the Wayburn and Seader 2 [1]_ global optimization
    problem. This is a unimodal minimization problem defined as follows:

    .. math::

        f_{\text{WayburnSeader02}}(x) = \left[ 1.613 - 4(x_1 - 0.3125)^2
                                        - 4(x_2 - 1.625)^2 \right]^2
                                        + (x_2 - 1)^2


    with :math:`x_i \in [-500, 500]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x = [0.2, 1]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-500.0] * self.N,
                           [500.0] * self.N))
        self.custom_bounds = ([-1, 2], [-1, 2])

        self.global_optimum = [[0.2, 1.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        u = (1.613 - 4 * (x[0] - 0.3125) ** 2 - 4 * (x[1] - 1.625) ** 2) ** 2
        v = (x[1] - 1) ** 2
        return u + v


class Weierstrass(Benchmark):

    r"""
    Weierstrass objective function.

    This class defines the Weierstrass [1]_ global optimization problem.
    This is a multimodal minimization problem defined as follows:

    .. math::

       f_{\text{Weierstrass}}(x) = \sum_{i=1}^{n} \left [
                                   \sum_{k=0}^{kmax} a^k \cos 
                                   \left( 2 \pi b^k (x_i + 0.5) \right) - n
                                   \sum_{k=0}^{kmax} a^k \cos(\pi b^k) \right ]


    Where, in this exercise, :math:`kmax = 20`, :math:`a = 0.5` and
    :math:`b = 3`.

    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-0.5, 0.5]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 4` for :math:`x_i = 0` for
    :math:`i = 1, ..., n`

    .. [1] Mishra, S. Global Optimization by Differential Evolution and
    Particle Swarm Methods: Evaluation on Some Benchmark Functions.
    Munich Personal RePEc Archive, 2006, 1005

    TODO line 1591.
    TODO Jamil, Gavana have got it wrong.  The second term is not supposed to
    be included in the outer sum. Mishra code has it right as does the
    reference referred to in Jamil#166.
    """
    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-0.5] * self.N, [0.5] * self.N))

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        kmax = 20
        a, b = 0.5, 3.0

        k = atleast_2d(arange(kmax + 1.)).T
        t1 = a ** k * cos(2 * pi * b ** k * (x + 0.5))
        t2 = self.N * sum(a ** k.T * cos(pi * b ** k.T))

        return sum(sum(t1, axis=0)) - t2


class Whitley(Benchmark):

    r"""
    Whitley objective function.

    This class defines the Whitley [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Whitley}}(x) = \sum_{i=1}^n \sum_{j=1}^n
                                \left[\frac{(100(x_i^2-x_j)^2
                                + (1-x_j)^2)^2}{4000} - \cos(100(x_i^2-x_j)^2
                                + (1-x_j)^2)+1 \right]


    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-10.24, 10.24]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x_i = 1` for
    :math:`i = 1, ..., n`

    .. [1] Gavana, A. Global Optimization Benchmarks and AMPGO retrieved 2015

    TODO Jamil#167 has '+ 1' inside the cos term, when it should be outside it.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-10.24] * self.N,
                           [10.24] * self.N))
        self.custom_bounds = ([-1, 2], [-1, 2])

        self.global_optimum = [[1.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        XI = x
        XJ = atleast_2d(x).T

        temp = 100.0 * ((XI ** 2.0) - XJ) + (1.0 - XJ) ** 2.0
        inner = (temp ** 2.0 / 4000.0) - cos(temp) + 1.0
        return sum(sum(inner, axis=0))


class Wolfe(Benchmark):

    r"""
    Wolfe objective function.

    This class defines the Wolfe [1]_ global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\text{Wolfe}}(x) = \frac{4}{3}(x_1^2 + x_2^2 - x_1x_2)^{0.75} + x_3


    with :math:`x_i \in [0, 2]` for :math:`i = 1, 2, 3`.

    *Global optimum*: :math:`f(x) = 0` for :math:`x = [0, 0, 0]`

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.
    """

    def __init__(self, dimensions=3):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([0.0] * self.N, [2.0] * self.N))

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return 4 / 3 * (x[0] ** 2 + x[1] ** 2 - x[0] * x[1]) ** 0.75 + x[2]
