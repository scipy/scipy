from __future__ import division
"""
==============================================================================
`go_benchmark_functions` --  Problems for testing global optimization routines
==============================================================================

This module provides a comprehensive set of problems for benchmarking global
optimization routines, such as scipy.optimize.basinhopping, or
scipy.optimize.differential_evolution.  The purpose is to see whether a given
optimization routine can find the global minimum, and how many function
evaluations it requires to do so.
The range of problems is extensive, with a range of difficulty. The problems are multivariate, with N=2 to N=17 provided.

References
----------
.. [1] Momin Jamil and Xin-She Yang, A literature survey of benchmark
    functions for global optimization problems, Int. Journal of Mathematical
    Modelling and Numerical Optimisation, Vol. 4, No. 2, pp. 150--194 (2013).
    http://arxiv.org/pdf/1308.4008v1.pdf
    (and references contained within)
.. [2] http://infinity77.net/global_optimization/index.html
.. [3] S. K. Mishra, Global Optimization By Differential Evolution and
    Particle Swarm Methods: Evaluation On Some Benchmark Functions, Munich
    Research Papers in Economics
.. [4] E. P. Adorio, U. P. Dilman, MVF - Multivariate Test Function Library
    in C for Unconstrained Global Optimization Methods, [Available Online]:
    http://www.geocities.ws/eadorio/mvf.pdf
.. [5] S. K. Mishra, Some New Test Functions For Global Optimization And
    Performance of Repulsive Particle Swarm Method, [Available Online]:
    http://mpra.ub.uni-muenchen.de/2718/
.. [6] NIST StRD Nonlinear Regression Problems, retrieved on 1 Oct, 2014
    http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml

"""

"""
Copyright 2013 Andrea Gavana
Author: <andrea.gavana@gmail.com>

Modifications 2014 Andrew Nelson
<andyfaff@gmail.com>
"""

import numpy as np
from numpy import abs, arange, arctan2, asarray, atleast_1d, cos, exp, floor, inf, log, ones, log10, arange
from numpy import pi, prod, roll, seterr, sign, sin, sqrt, sum, where, zeros, zeros_like, tan, tanh, dot

from scipy.misc import factorial

seterr(all='ignore')


class Benchmark(object):

    """
    Defines a global optimization benchmark problem.

    This abstract class defines the basic structure of a global
    optimization problem. Subclasses should implement the ``fun`` method
    for a particular optimization problem.

    Attributes
    ----------
    N
    bounds
    xmin
    xmax
    fglob : float
        The global minimum of the evaluated function.
    global_optimum : sequence
        A list of vectors that provide the locations of the global minimum.
        Note that some problems have multiple global minima, not all of which
        may be listed.
    nfev : int
        the number of function evaluations that the object has been asked to
        calculate.
    change_dimensionality : bool
        Whether we can change the benchmark function `x` variable length (i.e.,
        the dimensionality of the problem)
    custom_bounds : sequence
        a list of tuples that contain lower/upper bounds for use in plotting.
    """

    def __init__(self, dimensions):
        self.dimensions = dimensions
        self.nfev = 0
        self.fglob = np.nan
        self.global_optimum = None
        self.change_dimensionality = False
        self.custom_bounds = None

    def __str__(self):
        return '{0} ({1} dimensions)'.format(self.__class__.__name__, self.N)

    def __repr__(self):
        return self.__class__.__name__

    def initial_vector(self):
        """
        Random initialisation for the benchmark problem.

        Returns
        -------
        x : sequence
            a vector of length ``N`` that contains random floating point
            numbers that lie between the lower and upper bounds for a given
            parameter.
        """

        return asarray([np.random.uniform(l, u) for l, u in self._bounds])

    def success(self, x, tol=1.e-5):
        """
        Tests if a candidate solution at the global minimum.
        The default test is

        Parameters
        ----------
        x : sequence
            The candidate vector for testing if the global minimum has been
            reached. Must have ``len(x) == self.N``
        tol : float
            The evaluated function and known global minimum must differ by less
            than this amount to be at a global minimum.

        Returns
        -------
        bool : is the candidate vector at the global minimum?
        """
        val = self.fun(asarray(x))
        if abs(val - self.fglob) < tol:
            return True

        return False

    def fun(self, x):
        """
        Evaluation of the benchmark function.

        Parameters
        ----------
        x : sequence
            The candidate vector for evaluating the benchmark problem. Must
            have ``len(x) == self.N``.

        Returns
        -------
        val : float
              the evaluated benchmark function
        """

        raise NotImplementedError

    def change_dimensions(self, ndim):
        """
        Changes the dimensionality of the benchmark problem

        The dimensionality will only be changed if the problem is suitable

        Parameters
        ----------
        ndim - int
               The new dimensionality for the problem.
        """

        if self.change_dimensionality:
            self.dimensions = ndim

    @property
    def bounds(self):
        """
        The lower/upper bounds to be used for minimizing the problem.
        This a list of (lower, upper) tuples that contain the lower and upper
        bounds for the problem.  The problem should not be asked for evaluation
        outside these bounds. ``len(bounds) == N``.
        """
        if self.change_dimensionality:
            return [self._bounds[0]] * self.N
        else:
            return self._bounds

    @property
    def N(self):
        """
        The dimensionality of the problem.
        """
        return self.dimensions

    @property
    def xmin(self):
        """
        The lower bounds for the problem

        Returns
        -------
        xmin - sequence
            The lower bounds for the problem
        """
        return asarray([b[0] for b in self.bounds])

    @property
    def xmax(self):
        """
        The upper bounds for the problem

        Returns
        -------
        xmax - sequence
            The upper bounds for the problem
        """
        return asarray([b[1] for b in self.bounds])

#-----------------------------------------------------------------------
#                     SINGLE-OBJECTIVE PROBLEMS
#-----------------------------------------------------------------------

class Ackley01(Benchmark):
    """
    Ackley01 objective function.

    The Ackley01 [1]_ global optimization problem is a multimodal minimization
    problem defined as follows:

    .. math::

        f_{{Ackley01}}({x}) = -20e^{-0.2\sqrt{\frac{1}{n} \sum_{i=1}^n
         x_i^2}} - e^{ \frac{1}{n} \sum_{i=1}^n \cos(2 \pi x_i)} + 20 + e

    Here, :math:`n` represents the number of dimensions and :math:`x_i \in [-35,
    35]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for
    :math:`i=1,...,n`

    .. [1] Momin Jamil and Xin-She Yang, A literature survey of benchmark
     functions for global optimization problems, Int. Journal of Mathematical
     Modelling and Numerical Optimisation, Vol. 4, No. 2, pp. 150--194 (2013)

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-35.0] * self.N, [35.0] * self.N)
        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1
        u = sum(x ** 2)
        v = sum(cos(2 * pi * x))
        return (-20. * exp(-0.2 * sqrt(u / self.N))
                - exp(v / self.N) + 20. + exp(1.))


class Ackley02(object):
    """
    Ackley02 objective function.

    The Ackley02 [1]_ global optimization problem is a multimodal minimization
    problem defined as follows:

    .. math::

        f_{Ackley02}(\mathbf{x}) = -200e^{-0.02\sqrt{x_1^2 + x_2^2}}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \in [-32,
     32]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -200` for :math:`\mathbf{x} = [0, 0]`

    .. [1] Momin Jamil and Xin-She Yang, A literature survey of benchmark
     functions for global optimization problems, Int. Journal of Mathematical
     Modelling and Numerical Optimisation, Vol. 4, No. 2, pp. 150--194 (2013)

    """
    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-32.0] * self.N, [32.0] * self.N)
        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = -200.

    def fun(self, x):
        self.nfev += 1
        return -200 * exp(-0.02 * sqrt(x[0] ** 2 + x[1] ** 2))


class Ackley03(object):
    """
    Ackley03 [1]_ objective function.

    The Ackley03 global optimization problem is a multimodal minimization
    problem defined as follows:

    .. math::

        f_{Ackley03}(\mathbf{x}) = -200e^{-0.02\sqrt{x_1^2 + x_2^2}} +
            5e^{\cos(3x_1) + \sin(3x_2)}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \in [-32,
    32]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -195.62902825923879` for :math:`\mathbf{x}
    =[-0.68255758, -0.36070859]`

    .. [1] Momin Jamil and Xin-She Yang, A literature survey of benchmark
     functions for global optimization problems, Int. Journal of Mathematical
     Modelling and Numerical Optimisation, Vol. 4, No. 2, pp. 150--194 (2013)

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-32.0] * self.N, [32.0] * self.N)
        self.global_optimum = [[-0.68255758, -0.36070859]]
        self.fglob = -195.62902825923879

    def fun(self, x):
        a = -200 * exp(-0.02 * sqrt(x[0] ** 2 + x[1] ** 2))
        a += 5 * exp(cos(3 * x[0]) + sin(3 * x[1]))
        return a


class Adjiman(Benchmark):

    """
    Adjiman objective function.

    The Adjiman [1]_ global optimization problem is a multimodal minimization
    problem defined as follows:

    .. math::

        f_{{Adjiman}}(\mathbf{x}) = \cos(x_1)\sin(x_2) - \frac{x_1}{(x_2^2 +
        1)}

    Here, :math:`x_1 \in [-1, 2]` and :math:`x_2 \in [-1, 1]`.

    *Global optimum*: :math:`f(x_i) = -2.02181` for :math:`\mathbf{x} = [2.0,
    0.10578]`

    .. [1] Momin Jamil and Xin-She Yang, A literature survey of benchmark
     functions for global optimization problems, Int. Journal of Mathematical
     Modelling and Numerical Optimisation, Vol. 4, No. 2, pp. 150--194 (2013)

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = ([-1.0, 2.0], [-1.0, 1.0])
        self.global_optimum = [[2.0, 0.10578]]
        self.fglob = -2.02180678

    def fun(self, x, *args):
        self.nfev += 1
        return cos(x[0]) * sin(x[1]) - x[0] / (x[1] ** 2 + 1)


class Alpine01(Benchmark):

    """
    Alpine01 objective function.

    The Alpine01 [1]_ global optimization problem is a multimodal minimization
    problem defined as follows:

    .. math::

        f_{Alpine01}(\mathbf{x}) = \sum_{i=1}^{n} \lvert {x_i \sin \left( x_i
        \right) + 0.1 x_i} \rvert

    Here, :math:`n` represents the number of dimensions and :math:`x_i \in [-10,
    10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,
    n`

    .. [1] Momin Jamil and Xin-She Yang, A literature survey of benchmark
     functions for global optimization problems, Int. Journal of Mathematical
     Modelling and Numerical Optimisation, Vol. 4, No. 2, pp. 150--194 (2013)

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum(abs(x * sin(x) + 0.1 * x))


class Alpine02(Benchmark):

    """
    Alpine02 objective function.

    The Alpine02 [1]_ global optimization problem is a multimodal minimization
    problem defined as follows:

    .. math::

        f_{Alpine02}(\mathbf{x}) = \prod_{i=1}^{n} \sqrt{x_i} \sin(x_i)

    Here, :math:`n` represents the number of dimensions and :math:`x_i \in [0,
    10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = -6.1295` for :math:`\mathbf{x}=
    [7.91705268, 4.81584232]` for :math:`i=1,2`

    .. [1] Momin Jamil and Xin-She Yang, A literature survey of benchmark
     functions for global optimization problems, Int. Journal of Mathematical
     Modelling and Numerical Optimisation, Vol. 4, No. 2, pp. 150--194 (2013)

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N, [10.0] * self.N)
        #TODO check minima as a function of dimensionality
        self.global_optimum = [[7.91705268, 4.81584232]]
        self.fglob = -6.12950
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return prod(sqrt(x) * sin(x))


class AMGM(Benchmark):

    """
    AMGM objective function.

    The AMGM (Arithmetic Mean - Geometric Mean Equality) global optimization
    problem is a multimodal minimization problem defined as follows

    .. math::

        f_{{AMGM}}(\mathbf{x}) = \left ( \frac{1}{n} \sum_{i=1}^{n} x_i -
         \sqrt[n]{ \prod_{i=1}^{n} x_i} \right )^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \in [0,
    10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_1 = x_2 = ... = x_n` for
    :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N, [10.0] * self.N)
        self.global_optimum = [[1, 1]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        f1 = sum(x)
        f2 = prod(x)
        f1 = f1 / self.N
        f2 = f2 ** (1.0 / self.N)
        f = (f1 - f2) ** 2

        return f


class BartelsConn(Benchmark):

    """
    Bartels-Conn objective function.

    The BartelsConn [1]_ global optimization problem is a multimodal
    minimization problem defined as follows:

    .. math::

        f_{{BartelsConn}}(\mathbf{x}) = \lvert {x_1^2 + x_2^2 + x_1x_2} \rvert +
         \lvert {\sin(x_1)} \rvert + \lvert {\cos(x_2)} \rvert

    Here, :math:`n` represents the number of dimensions and :math:`x_i \in [-5,
    5]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 1` for :math:`\mathbf{x} = [0, 0]`

    .. [1] Momin Jamil and Xin-She Yang, A literature survey of benchmark
     functions for global optimization problems, Int. Journal of Mathematical
     Modelling and Numerical Optimisation, Vol. 4, No. 2, pp. 150--194 (2013)

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.0] * self.N, [5.0] * self.N)
        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 1.0

    def fun(self, x, *args):
        self.nfev += 1

        return (abs(x[0] ** 2.0 + x[1] ** 2.0 + x[0] * x[1]) + abs(sin(x[1]))
                + abs(cos(x[1])))


class Beale(Benchmark):

    """
    Beale objective function.

    The Beale [1]_ global optimization problem is a multimodal
    minimization problem defined as follows:

    .. math::

        f_{\text{Beale}}(\mathbf{x}) = \left(x_1 x_2 - x_1 + 1.5\right)^{2} +
        \left(x_1 x_2^{2} - x_1 + 2.25\right)^{2} + \left(x_1 x_2^{3} - x_1 +
        2.625\right)^{2}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \in [-4.5
    , 4.5]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\mathbf{x}=[3, 0.5]`

    .. [1] Momin Jamil and Xin-She Yang, A literature survey of benchmark
     functions for global optimization problems, Int. Journal of Mathematical
     Modelling and Numerical Optimisation, Vol. 4, No. 2, pp. 150--194 (2013)
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


class BiggsExp02(object):
    """
    BiggsExp02 objective function.

    The BiggsExp02 [1]_ global optimization problem is a multimodal minimization
    problem defined as follows

    .. math::

        \begin{array}\\ f_{{BiggsExp02}}(\mathbf{x}) = \sum_{i=1}^{10}
        (e^{-t_ix_1} - 5e^{-t_ix_2} - y_i)^2\\
        t_i = 0.1i\\
        y_i = e^{-t_i} - 5e^{-10t_i}
        \end{array}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \in [0,
     20]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\mathbf{x}=[1, 10]`

    .. [1] Momin Jamil and Xin-She Yang, A literature survey of benchmark
     functions for global optimization problems, Int. Journal of Mathematical
     Modelling and Numerical Optimisation, Vol. 4, No. 2, pp. 150--194 (2013)

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0] * 2,
                           [20] * 2)
        self.global_optimum = [[1., 10.]]
        self.fglob = 0

    def fun(self, x):
        self.nfev += 1

        t = arange(1, 11.) * 0.1
        y = exp(-t) - 5 * exp(-10 * t)
        vec = (exp(-t * x[0]) - 5 * exp(-t * x[1]) - y) ** 2

        return sum(vec)


class BiggsExp03(object):

    """
    BiggsExp03 objective function.

    The BiggsExp03 [1]_ global optimization problem is a multimodal minimization
    problem defined as follows

    .. math::

        \begin{array}\\ f_{BiggsExp03}(\mathbf{x}) = \sum_{i=1}^{10}
        (e^{-t_ix_1} - x_3e^{-t_ix_2} - y_i)^2\\
        t_i = 0.1i\\
        y_i = e^{-t_i} - 5e^{-10t_i}
        \end{array}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \in [0,
    20]` for :math:`i=1,2,3`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\mathbf{x}=[1, 10, 5]`

    .. [1] Momin Jamil and Xin-She Yang, A literature survey of benchmark
     functions for global optimization problems, Int. Journal of Mathematical
     Modelling and Numerical Optimisation, Vol. 4, No. 2, pp. 150--194 (2013)

    """

    def __init__(self, dimensions=3):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0] * 3,
                           [20] * 3)
        self.global_optimum = [[1., 10., 5.]]
        self.fglob = 0

    def fun(self, x):
        self.nfev += 1

        t = arange(1., 11.) * 0.1
        y = exp(-t) - 5 * exp(-10 * t)
        vec = (exp(-t * x[0]) - x[2] * exp(-t * x[1]) - y) ** 2

        return sum(vec)


class BiggsExp04(object):

    """
    BiggsExp04 objective function.

    The BiggsExp04 [1]_ global optimization problem is a multimodal minimization
    problem defined as follows

    .. math::

        \begin{array}\\ f_{BiggsExp04}(\mathbf{x}) = \sum_{i=1}^{10}
        (x_3e^{-t_ix_1} - x_4e^{-t_ix_2} - y_i)^2\\
        t_i = 0.1i\\
        y_i = e^{-t_i} - 5e^{-10t_i}
        \end{array}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \in [0,
    20]` for :math:`i=1,2,3,4`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\mathbf{x}=[1, 10, 1, 5]`

    .. [1] Momin Jamil and Xin-She Yang, A literature survey of benchmark
     functions for global optimization problems, Int. Journal of Mathematical
     Modelling and Numerical Optimisation, Vol. 4, No. 2, pp. 150--194 (2013)

    """

    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.] * 4,
                           [20.] * 4)
        self.global_optimum = [[1., 10., 1., 5.]]
        self.fglob = 0

    def fun(self, x):
        self.nfev += 1

        t = arange(1, 11.) * 0.1
        y = exp(-t) - 5 * exp(-10 * t)
        vec = (x[2] * exp(-t * x[0]) - x[3] * exp(-t * x[1]) - y) ** 2

        return sum(vec)


class BiggsExp05(object):

    """
    BiggsExp05 objective function.

    The BiggsExp05 [1]_ global optimization problem is a multimodal minimization
    problem defined as follows

    .. math::

        \begin{array}\\ f_{BiggsExp04}(\mathbf{x}) = \sum_{i=1}^{11}
        (x_3e^{-t_ix_1} - x_4e^{-t_ix_2} + 3e^{-t_ix_5} - y_i)^2\\
        t_i = 0.1i\\
        y_i = e^{-t_i} - 5e^{-10t_i} + 3e^{-4t_i}
        \end{array}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \in [0,
     20]` for :math:`i=1,...,5`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\mathbf{x}=[1, 10, 1, 5, 4]`

    .. [1] Momin Jamil and Xin-She Yang, A literature survey of benchmark
     functions for global optimization problems, Int. Journal of Mathematical
     Modelling and Numerical Optimisation, Vol. 4, No. 2, pp. 150--194 (2013)

    """

    def __init__(self, dimensions=5):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.] * 5,
                           [20.] * 5)
        self.global_optimum = [[1., 10., 1., 5., 4.]]
        self.fglob = 0

    target_E = 0
    solution = [1., 10., 1., 5., 4.]
    xmin = np.array([0., 0., 0., 0., 0.])
    xmax = np.array([20., 20., 20., 20., 20.])

    def fun(self, x):
        self.nfev += 1
        t = arange(1, 12.) * 0.1
        y = exp(-t) - 5 * exp(-10 * t) + 3 * exp(-4 * t)
        vec = (x[2] * exp(-t * x[0]) - x[3] * exp(-t * x[1])
               + 3 * exp(-t * x[4]) - y) ** 2

        return sum(vec)


class Bird(Benchmark):

    """
    Bird objective function.

    The Bird global optimization problem is a multimodal minimization
    problem defined as follows

    .. math::

        f_{Bird}(\mathbf{x}) = \left(x_1 - x_2\right)^{2} + e^{\left[1 -
         \sin\left(x_1\right) \right]^{2}} \cos\left(x_2\right) + e^{\left[1 -
          \cos\left(x_2\right)\right]^{2}} \sin\left(x_1\right)

    for :math:`x_i \in [-2\pi, 2\pi]`

    *Global optimum*: :math:`f(x_i) = -106.7645367198034` for :math:`\mathbf{x}
    = [4.701055751981055 , 3.152946019601391]` or :math:`\mathbf{x} =
    [-1.582142172055011, -3.130246799635430]`

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


class Bohachevsky(Benchmark):

    """
    Bohachevsky objective function.

    The Bohachevsky [1]_ global optimization problem is a multimodal
    minimization problem defined as follows

        .. math::

        f_{Bohachevsky}(\mathbf{x}) = \sum_{i=1}^{n-1}\left[x_i^2 + 2x_{i+1}^2 -
        0.3\cos(3\pi x_i) - 0.4\cos(4\pi x_{i+1}) + 0.7\right]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \in [-15,
    15]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for
    :math:`i=1,...,n`

    .. [1] Momin Jamil and Xin-She Yang, A literature survey of benchmark
     functions for global optimization problems, Int. Journal of Mathematical
     Modelling and Numerical Optimisation, Vol. 4, No. 2, pp. 150--194 (2013)

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-15.0] * self.N, [15.0] * self.N)
        self.custom_bounds = [(-2, 2), (-2, 2)]
        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        x0 = x[:-1]
        x1 = roll(x, -1)[:-1]

        return sum(x0 ** 2 + 2 * x1 ** 2 - 0.3 * cos(3 * pi * x0)
                   - 0.4 * cos(4 * pi * x1) + 0.7)


class BoxBetts(Benchmark):

    """
    BoxBetts objective function.

    The BoxBetts global optimization problem is a multimodal
    minimization problem defined as follows

    .. math::

        f_{BoxBetts}(\mathbf{x}) = \sum_{i=1}^k g(x_i)^2

    Where, in this exercise:

    .. math::
        g(\mathbf{x}) = e^{-0.1ix_1} - e^{-0.1ix_2} - x_3\left[e^{-0.1i}
        - e^{-i}\right]


    And :math:`k = 10`.

    Here, :math:`n` represents the number of dimensions and :math:`x_1 \in [0.9,
    1.2], x_2 \in [9, 11.2], x_3 \in [0.9, 1.2]`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\mathbf{x} = [1, 10, 1]`

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
             - (exp(-0.1 * i) - exp(-1.0 * i)) * x[2])
        return sum(g**2)


class Branin01(Benchmark):

    """
    Branin01  objective function.

    The Branin01 global optimization problem is a multimodal minimization
    problem defined as follows

    .. math::

        f_{Branin01}(\mathbf{x}) = \left(- 1.275 \frac{x_1^{2}}{\pi^{2}} + 5
        \frac{x_1}{\pi} + x_2 -6\right)^{2} + \left(10 - \frac{5}{4 \pi} \right)
        \cos\left(x_1\right) + 10

    Here, :math:`n` represents the number of dimensions and :math:`x_1 \\in [-5,
    10], x_2 \\in [0, 15]`

    *Global optimum*: :math:`f(x_i) = 0.39788735772973816` for
    :math:`\mathbf{x} = [-\pi, 12.275]` or :math:`\mathbf{x} = [\pi, 2.275]`
    or :math:`\mathbf{x} = [9.42478, 2.475]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = [(-5., 10.), (0., 15.)]

        self.global_optimum = [[-pi, 12.275], [pi, 2.275], [9.42478, 2.475]]
        self.fglob = 0.39788735772973816

    def fun(self, x, *args):
        self.nfev += 1

        return ((x[1] - (5.1 / (4 * pi ** 2)) * x[0] ** 2
                + 5 * x[0] / pi - 6) ** 2
                + 10 * (1 - 1 / (8 * pi)) * cos(x[0]) + 10)


class Branin02(Benchmark):

    """
    Branin02 objective function.

    The Branin02 global optimization problem is a multimodal minimization
    problem defined as follows


    .. math::

        f_{\text{Branin02}}(\mathbf{x}) = \left(- 1.275 \frac{x_1^{2}}{\pi^{2}}
        + 5 \frac{x_1}{\pi} + x_2 -6\right)^{2} + \left(10 - \frac{5}{4 \pi}
        \right) \cos\left(x_1\right) \cos\left(x_2\right) + \log(x_1^2+x_2^2 +1)
        + 10

    Here, :math:`n` represents the number of dimensions and :math:`x_i \in [-5,
    15]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 5.559037` for :math:`\mathbf{x} = [-3.2,
    12.53]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = [(-5.0, 15.0), (-5.0, 15.0)]

        self.global_optimum = [[-3.2, 12.53]]
        self.fglob = 5.559037

    def fun(self, x, *args):
        self.nfev += 1

        return ((x[1] - (5.1 / (4 * pi ** 2)) * x[0] ** 2
                + 5 * x[0] / pi - 6) ** 2
                + 10 * (1 - 1 / (8 * pi)) * cos(x[0]) * cos(x[1])
                + log(x[0] ** 2.0 + x[1] ** 2.0 + 1.0) + 10)


class Brent(Benchmark):

    """
    Brent objective function.

    The Brent [1]_ global optimization problem is a multimodal minimization
    problem defined as follows:

    .. math::

        f_{\text{Brent}}(\mathbf{x}) = (x_1 + 10)^2 + (x_2 + 10)^2 +
        e^{(-x_1^2-x_2^2)}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \in [-10,
    10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\mathbf{x} = [-10, -10]`

    .. [1] Momin Jamil and Xin-She Yang, A literature survey of benchmark
     functions for global optimization problems, Int. Journal of Mathematical
     Modelling and Numerical Optimisation, Vol. 4, No. 2, pp. 150--194 (2013)

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

    """
    Brown objective function.

    The Brown [1]_ global optimization problem is a multimodal minimization
    problem defined as follows:

    .. math::

        f_{\text{Brown}}(\mathbf{x}) = \sum_{i=1}^{n-1}\left[
        \left(x_i^2\right)^{x_{i+1}^2+1} + \left(x_{i+1}^2\right)^{x_i^2+1}
        \right]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \in [-1,
    4]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for
    :math:`i=1,...,n`

    .. [1] Momin Jamil and Xin-She Yang, A literature survey of benchmark
     functions for global optimization problems, Int. Journal of Mathematical
     Modelling and Numerical Optimisation, Vol. 4, No. 2, pp. 150--194 (2013)

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

    """
    Bukin02 objective function.

    The Bukin02 [1]_ global optimization problem is a multimodal minimization
    problem defined as follows:

    .. math::

        f_{\text{Bukin02}}(\mathbf{x}) = 100 (x_2 - 0.01x_1^2 + 1)
        + 0.01(x_1 + 10)^2

    Here, :math:`n` represents the number of dimensions and :math:`x_1 \in [-15,
    -5], x_2 \in [-3, 3]`

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\mathbf{x} = [-10, 0]`

    .. [1] Momin Jamil and Xin-She Yang, A literature survey of benchmark
     functions for global optimization problems, Int. Journal of Mathematical
     Modelling and Numerical Optimisation, Vol. 4, No. 2, pp. 150--194 (2013)

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = [(-15.0, -5.0), (-3.0, 3.0)]

        self.global_optimum = [[-10.0, 0.0]]
        self.fglob = 0.0

    def fun(self, x, *args):

        self.nfev += 1
        return (100 * (x[1] ** 2 - 0.01 * x[0] ** 2 + 1.0)
                + 0.01 * (x[0] + 10.0) ** 2.0)


class Bukin04(Benchmark):

    """
    Bukin04 objective function.

    The Bukin04 [1]_ global optimization problem is a multimodal minimization
    problem defined as follows:

    .. math::

        f_{\text{Bukin04}}(\mathbf{x}) = 100 x_2^{2} + 0.01 \lvert{x_1 + 10}
        \rvert

    Here, :math:`n` represents the number of dimensions and :math:`x_1 \in [-15,
    -5], x_2 \in [-3, 3]`

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\mathbf{x} = [-10, 0]`

    .. [1] Momin Jamil and Xin-She Yang, A literature survey of benchmark
     functions for global optimization problems, Int. Journal of Mathematical
     Modelling and Numerical Optimisation, Vol. 4, No. 2, pp. 150--194 (2013)

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

    """
    Bukin06 objective function.

    The Bukin06 [1]_ global optimization problem is a multimodal minimization
    problem defined as follows:

    .. math::

        f_{\text{Bukin06}}(\mathbf{x}) = 100 \sqrt{ \lvert{x_2 - 0.01 x_1^{2}}
        \rvert} + 0.01 \lvert{x_1 + 10} \rvert

    Here, :math:`n` represents the number of dimensions and :math:`x_1 \in [-15,
    -5], x_2 \in [-3, 3]`

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\mathbf{x} = [-10, 1]`

    .. [1] Momin Jamil and Xin-She Yang, A literature survey of benchmark
     functions for global optimization problems, Int. Journal of Mathematical
     Modelling and Numerical Optimisation, Vol. 4, No. 2, pp. 150--194 (2013)

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = [(-15.0, -5.0), (-3.0, 3.0)]
        self.global_optimum = [[-10.0, 1.0]]
        self.fglob = 0.0

    def fun(self, x, *args):

        self.nfev += 1
        return 100 * sqrt(abs(x[1] - 0.01 * x[0] ** 2)) + 0.01 * abs(x[0] + 10)


class CarromTable(Benchmark):

    """
    CarromTable objective function.

    The CarromTable [1]_ global optimization problem is a multimodal
    minimization problem defined as follows:

    .. math::

        f_{\text{CarromTable}}(\mathbf{x}) = - \frac{1}{30}\left(\cos(x_1)
        cos(x_2) e^{\left|1 - \frac{\sqrt{x_1^2 + x_2^2}}{\pi}\right|}\right)^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \in [-10,
    10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -24.15681551650653` for :math:`x_i = \pm
    9.646157266348881` for :math:`i=1,...,n`

    ..[1] S. K. Mishra, Global Optimization By Differential Evolution and
     Particle Swarm Methods: Evaluation On Some Benchmark Functions, Munich
     Research Papers in Economics

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.global_optimum = [(9.646157266348881, 9.646134286497169),
                               (-9.646157266348881, 9.646134286497169),
                               (9.646157266348881, -9.646134286497169),
                               (-9.646157266348881, -9.646134286497169)]
        self.fglob = -24.15681551650653

    def fun(self, x, *args):
        self.nfev += 1

        u = cos(x[0]) * cos(x[1])
        v = sqrt(x[0] ** 2 + x[1] ** 2)
        return -((u * exp(abs(1 - v / pi))) ** 2) / 30


class Chichinadze(Benchmark):

    """
    Chichinadze objective function.

    This class defines the Chichinadze global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Chichinadze}}(\\mathbf{x}) = x_{1}^{2} - 12 x_{1} + 8 \\sin\\left(\\frac{5}{2} \\pi x_{1}\\right) + 10 \\cos\\left(\\frac{1}{2} \\pi x_{1}\\right) + 11 - 0.2 \\frac{\\sqrt{5}}{e^{\\frac{1}{2} \\left(x_{2} -0.5\\right)^{2}}}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-30, 30]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -42.94438701899098` for :math:`\\mathbf{x} = [6.189866586965680, 0.5]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-30.0] * self.N, [30.0] * self.N)
        self.custom_bounds = [(-10, 10), (-10, 10)]

        self.global_optimum = [[6.189866586965680, 0.5]]
        self.fglob = -42.94438701899098

    def fun(self, x, *args):
        self.nfev += 1

        return (x[0] ** 2 - 12 * x[0] + 11 + 10 * cos(pi * x[0] / 2)
                + 8 * sin(5 * pi * x[0] / 2)
                - 1.0 / sqrt(5) * exp(-((x[1] - 0.5) ** 2) / 2))


class Cigar(Benchmark):

    """
    Cigar objective function.

    This class defines the Cigar global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Cigar}}(\\mathbf{x}) = x_1^2 + 10^6\\sum_{i=2}^{n} x_i^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = [(-5, 5), (-5, 5)]

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return x[0] ** 2 + 1e6 * sum(x[1:] ** 2)


class Cola(Benchmark):

    """
    Cola objective function.

    This class defines the Cola global optimization problem. The 17-dimensional function computes
    indirectly the formula :math:`f(n, u)` by setting :math:`x_0 = y_0, x_1 = u_0, x_i = u_{2(i2)}, y_i = u_{2(i2)+1}` :

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

        self._bounds = [[0.0, 4.0]] + \
            list(zip([-4.0] * (self.N - 1),
                 [4.0] * (self.N - 1)))

        self.global_optimum = [
            [0.651906, 1.30194, 0.099242, -0.883791, -0.8796,
             0.204651, -3.28414, 0.851188, -3.46245, 2.53245,
             -0.895246, 1.40992, -3.07367, 1.96257, -2.97872,
             -0.807849, -1.68978]]
        self.fglob = 11.7464

    def fun(self, x, *args):

        self.nfev += 1

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


class Colville(Benchmark):

    """
    Colville objective function.

    This class defines the Colville global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Colville}}(\\mathbf{x}) = \\left(x_{1} -1\\right)^{2} + 100 \\left(x_{1}^{2} - x_{2}\\right)^{2} + 10.1 \\left(x_{2} -1\\right)^{2} + \\left(x_{3} -1\\right)^{2} + 90 \\left(x_{3}^{2} - x_{4}\\right)^{2} + 10.1 \\left(x_{4} -1\\right)^{2} + 19.8 \\frac{x_{4} -1}{x_{2}}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,4`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 1` for :math:`i=1,...,4`

    """

    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[1 for _ in range(self.N)]]
        self.fglob = 0.0

    def fun(self, x, *args):

        self.nfev += 1
        return (100 * (x[0] ** 2 - x[1]) ** 2
                + (x[0] - 1) ** 2 + (x[2] - 1) ** 2
                + 90 * (x[2] ** 2 - x[3]) ** 2
                + 10.1 * ((x[1] - 1) ** 2 + (x[3] - 1) ** 2)
                + 19.8 * (1 / x[1]) * (x[3] - 1))


class Corana(Benchmark):

    """
    Corana objective function.

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

        self._bounds = zip([-5.0] * self.N, [5.0] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0

    def fun(self, x, *args):

        self.nfev += 1

        d = [1., 1000., 10., 100.]
        r = 0
        for j in range(4):
            zj = floor(abs(x[j] / 0.2) + 0.49999) * sign(x[j]) * 0.2
            if abs(x[j] - zj) < 0.05:
                r += 0.15 * ((zj - 0.05 * sign(zj)) ** 2) * d[j]
            else:
                r += d[j] * x[j] * x[j]
        return r


class CosineMixture(Benchmark):

    """
    Cosine Mixture objective function.

    This class defines the Cosine Mixture global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{CosineMixture}}(\\mathbf{x}) = -0.1 \\sum_{i=1}^n \\cos(5 \\pi x_i) - \\sum_{i=1}^n x_i^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,...,N`.

    *Global optimum*: :math:`f(x_i) = -0.1N` for :math:`x_i = 0` for :math:`i=1,...,N`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.change_dimensionality = True
        self._bounds = zip([-1.0] * self.N, [1.0] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = -0.1 * self.N

    def fun(self, x, *args):
        self.nfev += 1

        return -0.1 * sum(cos(5.0 * pi * x)) - sum(x ** 2.0)


class CrossInTray(Benchmark):

    """
    Cross-in-Tray objective function.

    This class defines the Cross-in-Tray global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{CrossInTray}}(\\mathbf{x}) = - 0.0001 \\left(\\left|{e^{\\left|{100 - \\frac{\\sqrt{x_{1}^{2} + x_{2}^{2}}}{\\pi}}\\right|} \\sin\\left(x_{1}\\right) \\sin\\left(x_{2}\\right)}\\right| + 1\\right)^{0.1}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-15, 15]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -2.062611870822739` for :math:`x_i = \\pm 1.349406608602084` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [(1.349406685353340, 1.349406608602084),
                               (-1.349406685353340, 1.349406608602084),
                               (1.349406685353340, -1.349406608602084),
                               (-1.349406685353340, -1.349406608602084)]
        self.fglob = -2.062611870822739

    def fun(self, x, *args):

        self.nfev += 1
        return (-0.0001 * (abs(sin(x[0]) * sin(x[1])
                           * exp(abs(100 - sqrt(x[0] ** 2 + x[1] ** 2) / pi)))
                           + 1) ** (0.1))


class CrossLegTable(Benchmark):

    """
    Cross-Leg-Table objective function.

    This class defines the Cross-Leg-Table global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{CrossLegTable}}(\\mathbf{x}) = - \\frac{1}{\\left(\\left|{e^{\\left|{100 - \\frac{\\sqrt{x_{1}^{2} + x_{2}^{2}}}{\\pi}}\\right|} \\sin\\left(x_{1}\\right) \\sin\\left(x_{2}\\right)}\\right| + 1\\right)^{0.1}}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -1`. The global minimum is found on the planes :math:`x_1 = 0` and :math:`x_2 = 0`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[0., 0.]]
        self.fglob = -1.0

    def fun(self, x, *args):

        self.nfev += 1
        u = 100 - sqrt(x[0] ** 2 + x[1] ** 2) / pi
        v = sin(x[0]) * sin(x[1])
        return -(abs(v * exp(abs(u))) + 1) ** (-0.1)


class CrownedCross(Benchmark):

    """
    Crowned Cross objective function.

    This class defines the Crowned Cross global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{CrownedCross}}(\\mathbf{x}) = 0.0001 \\left(\\left|{e^{\\left|{100- \\frac{\\sqrt{x_{1}^{2} + x_{2}^{2}}}{\\pi}}\\right|} \\sin\\left(x_{1}\\right) \\sin\\left(x_{2}\\right)}\\right| + 1\\right)^{0.1}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0.0001`. The global minimum is found on the planes :math:`x_1 = 0` and :math:`x_2 = 0`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[0, 0]]
        self.fglob = 0.0001

    def fun(self, x, *args):

        self.nfev += 1
        u = 100 - sqrt(x[0] ** 2 + x[1] ** 2) / pi
        v = sin(x[0]) * sin(x[1])
        return 0.0001 * (abs(v * exp(abs(u))) + 1) ** (0.1)


class Csendes(Benchmark):

    """
    Csendes objective function.

    This class defines the Csendes global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Csendes}}(\\mathbf{x}) = \\sum_{i=1}^n x_i^6 \\left[ 2 + \\sin \\left( \\frac{1}{x_i} \\right ) \\right]


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,...,N`.

    *Global optimum*: :math:`f(x_i) = 0.0` for :math:`x_i = 0` for :math:`i=1,...,N`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.change_dimensionality = True
        self._bounds = zip([-1.0] * self.N, [1.0] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        try:
            return sum((x ** 6.0) * (2.0 + sin(1.0 / x)))
        except ZeroDivisionError, FloatingPointError:
            return np.nan

    def success(self, x):
        """Is a candidate solution at the global minimum"""
        val = self.fun(asarray(x))
        if np.isnan(val):
            return True

        return False


class Cube(Benchmark):

    """
    Cube objective function.

    This class defines the Cube global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Cube}}(\\mathbf{x}) = 100(x_2 - x_1^3)^2 + (1 - x1)^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,N`.

    *Global optimum*: :math:`f(x_i) = 0.0` for :math:`\\mathbf{x} = [1, 1]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.custom_bounds = ([0, 2], [0, 2])

        self.global_optimum = [[1.0, 1.0]]
        self.fglob = 0.0

    def fun(self, x, *args):

        self.nfev += 1
        return 100.0 * (x[1] - x[0] ** 3.0) ** 2.0 + (1.0 - x[0]) ** 2.0


class Damavandi(Benchmark):

    """
    Damavandi objective function.

    This class defines the Damavandi global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Damavandi}}(\\mathbf{x}) = \\left[ 1 - \\lvert{\\frac{\\sin[\\pi(x_1-2)]\\sin[\\pi(x2-2)]}{\\pi^2(x_1-2)(x_2-2)}} \\rvert^5 \\right] \\left[2 + (x_1-7)^2 + 2(x_2-7)^2 \\right]


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 14]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0.0` for :math:`x_i = 2` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.change_dimensionality = True
        self._bounds = zip([0.0] * self.N, [14.0] * self.N)

        self.global_optimum = [[2 for _ in range(self.N)]]
        self.fglob = 0.0

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
            np.testing.assert_almost_equal(val, 0., 4)
            return True
        except AssertionError:
            return False

        return False


class Deb01(Benchmark):

    """
    Deb 1 objective function.

    This class defines the Deb 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Deb01}}(\\mathbf{x}) = - \\frac{1}{N} \\sum_{i=1}^n \\sin^6(5 \\pi x_i)


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0.0`. The number of global minima is :math:`5^n` that are evenly spaced
    in the function landscape, where :math:`n` represents the dimension of the problem.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.change_dimensionality = True

        self._bounds = zip([-1.0] * self.N, [1.0] * self.N)

        self.global_optimum = [[0.3, -0.3]]
        self.fglob = -1.0

    def fun(self, x, *args):
        self.nfev += 1
        return -(1.0 / self.N) * sum(sin(5 * pi * x) ** 6.0)


class Deb02(Benchmark):

    """
    Deb 2 objective function.

    This class defines the Deb 2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Deb02}}(\\mathbf{x}) = - \\frac{1}{N} \\sum_{i=1}^n \\sin^6 \\left[ 5 \\pi \\left ( x_i^{3/4} - 0.05 \\right) \\right ]


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 1]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0.0`. The number of global minima is :math:`5^n` that are evenly spaced
    in the function landscape, where :math:`n` represents the dimension of the problem.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self.change_dimensionality = True

        self._bounds = zip([0.0] * self.N, [1.0] * self.N)

        self.global_optimum = [[0.93388314, 0.68141781]]
        self.fglob = -1.0

    def fun(self, x, *args):
        self.nfev += 1

        return -(1.0 / self.N) * sum(sin(5 * pi * (x ** 0.75 - 0.05)) ** 6.0)


class Decanomial(Benchmark):

    """
    Decanomial objective function.

    This class defines the Decanomial function global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Decanomial}}(\\mathbf{x}) = 0.001 \\left(\\lvert{x_{2}^{4} + 12 x_{2}^{3} + 54 x_{2}^{2} + 108 x_{2} + 81.0}\\rvert + \\lvert{x_{1}^{10} - 20 x_{1}^{9} + 180 x_{1}^{8} - 960 x_{1}^{7} + 3360 x_{1}^{6} - 8064 x_{1}^{5} + 13340 x_{1}^{4} - 15360 x_{1}^{3} + 11520 x_{1}^{2} - 5120 x_{1} + 2624.0}\\rvert\\right)^{2}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [2, -3]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
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

    """
    Deceptive objective function.

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

    *Global optimum*: :math:`f(x_i) = -1` for :math:`x_i = \\alpha_i` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N, [1.0] * self.N)

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

    """
    Deckkers-Aarts objective function.

    This class defines the Deckkers-Aarts global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{DeckkersAarts}}(\\mathbf{x}) = 10^5x_1^2 + x_2^2 - (x_1^2 + x_2^2)^2 + 10^{-5}(x_1^2 + x_2^2)^4


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-20, 20]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -24776.518242168` for :math:`\\mathbf{x} = [0, \\pm 14.9451209]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-20.0] * self.N, [20.0] * self.N)
        self.custom_bounds = ([-1, 1], [14, 16])

        self.global_optimum = [[0.0, 14.9451209]]
        self.fglob = -24776.518342168

    def fun(self, x, *args):
        self.nfev += 1
        return (1.e5 * x[0] ** 2 + x[1] ** 2 - (x[0] ** 2 + x[1] ** 2) ** 2
                + 1.e-5 * (x[0] ** 2 + x[1] ** 2) ** 4)


class DeflectedCorrugatedSpring(Benchmark):

    """
    DeflectedCorrugatedSpring objective function.

    This class defines the Deflected Corrugated Spring function global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{DeflectedCorrugatedSpring}}(\\mathbf{x}) = 0.1\\sum_{i=1}^n \\left[ (x_i - \\alpha)^2 - \\cos \\left( K \\sqrt {\\sum_{i=1}^n (x_i - \\alpha)^2} \\right ) \\right ]


    Where, in this exercise, :math:`K = 5` and :math:`\\alpha = 5`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 2\\alpha]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = \\alpha` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        alpha = 5.0
        self._bounds = zip([0] * self.N, [2 * alpha] * self.N)

        self.global_optimum = [[alpha for _ in range(self.N)]]
        self.fglob = -1.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1
        K, alpha = 5.0, 5.0

        return (-cos(K * sqrt(sum((x - alpha) ** 2)))
                + 0.1 * sum((x - alpha) ** 2))


class DeVilliersGlasser01(Benchmark):

    """
    DeVilliers-Glasser 1 objective function.

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

        self._bounds = zip([1.0] * self.N, [100.0] * self.N)

        self.global_optimum = [[60.137, 1.371, 3.112, 1.761]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        t = 0.1 * arange(24)
        y = 60.137 * (1.371 ** t) * sin(3.112 * t + 1.761)

        return sum((x[0] * (x[1] ** t) * sin(x[2] * t + x[3]) - y) ** 2.0)


class DeVilliersGlasser02(Benchmark):

    """
    DeVilliers-Glasser 2 objective function.

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

        self._bounds = zip([1.0] * self.N, [60.0] * self.N)

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

    """
    Dixon and Price objective function.

    This class defines the Dixon and Price global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{DixonPrice}}(\\mathbf{x}) = (x_i - 1)^2 + \\sum_{i=2}^n i(2x_i^2 - x_{i-1})^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 2^{- \\frac{(2^i-2)}{2^i}}` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.custom_bounds = [(-2, 3), (-2, 3)]

        self.global_optimum = [[2.0 ** (-(2.0 ** i - 2.0) / 2.0 ** i)
                               for i in range(1, self.N + 1)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        i = arange(1, self.N)
        s = i * (2.0 * x[1:] ** 2.0 - x[:-1]) ** 2.0
        return sum(s) + (x[0] - 1.0) ** 2.0


class Dolan(Benchmark):

    """
    Dolan objective function.

    This class defines the Dolan global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Dolan}}(\\mathbf{x}) = \\lvert (x_1 + 1.7x_2)\\sin(x_1) - 1.5x_3 - 0.1x_4\\cos(x_5 + x_5 - x_1) + 0.2x_5^2 - x_2 - 1 \\rvert

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 10^{-5}` for :math:`\\mathbf{x} = [8.39045925, 4.81424707, 7.34574133, 68.88246895, 3.85470806]`

    """

    def __init__(self, dimensions=5):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)

        self.global_optimum = [
            [8.39045925, 4.81424707, 7.34574133, 68.88246895,
             3.85470806]]
        self.fglob = 1e-5

    def fun(self, x, *args):
        self.nfev += 1

        return (abs((x[0] + 1.7 * x[1]) * sin(x[0]) - 1.5 * x[2]
                - 0.1 * x[3] * cos(x[3] + x[4] - x[0]) + 0.2 * x[4] ** 2
                - x[1] - 1))


class DropWave(Benchmark):

    """
    DropWave objective function.

    This class defines the DropWave global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{DropWave}}(\\mathbf{x}) = - \\frac{1 + \\cos\\left(12 \\sqrt{\\sum_{i=1}^{n} x_i^{2}}\\right)}{2 + 0.5 \\sum_{i=1}^{n} x_i^{2}}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5.12, 5.12]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -1` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.12] * self.N, [5.12] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = -1.0

    def fun(self, x, *args):
        self.nfev += 1

        norm_x = sum(x ** 2)
        return -(1 + cos(12 * sqrt(norm_x))) / (0.5 * norm_x + 2)


class Easom(Benchmark):

    """
    Easom objective function.

    This class defines the Easom global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Easom}}(\\mathbf{x}) = a - \\frac{a}{e^{b \\sqrt{\\frac{\\sum_{i=1}^{n} x_i^{2}}{n}}}} + e - e^{\\frac{\\sum_{i=1}^{n} \\cos\\left(c x_i\\right)}{n}}

    Where, in this exercise, :math:`a = 20, b = 0.2` and :math:`c = 2\\pi`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        a = 20.0
        b = 0.2
        c = 2 * pi

        return (-a * exp(-b * sqrt(sum(x ** 2) / self.N))
                - exp(sum(cos(c * x)) / self.N) + a + exp(1))


class Eckerle4(Benchmark):
    """
    Eckerle4 objective function.
    Eckerle, K., NIST (1979).
    Circular Interference Transmittance Study.
    """

    # TODO, this is a NIST regression standard dataset
    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0., 1., 0., 0.1],
                           [1000, 20., 3., 6.])
        self.global_optimum = [[1.5543827178, 4.0888321754, 4.5154121844e2]]
        self.fglob = 1.4635887487E-03

    def fun(self, x, *args):
        self.nfev += 1

        a = asarray([1.5750000E-04, 1.6990000E-04, 2.3500000E-04, 3.1020000E-04,
                     4.9170000E-04, 8.7100000E-04, 1.7418000E-03, 4.6400000E-03,
                     6.5895000E-03, 9.7302000E-03, 1.4900200E-02, 2.3731000E-02,
                     4.0168300E-02, 7.1255900E-02, 1.2644580E-01, 2.0734130E-01,
                     2.9023660E-01, 3.4456230E-01, 3.6980490E-01, 3.6685340E-01,
                     3.1067270E-01, 2.0781540E-01, 1.1643540E-01, 6.1676400E-02,
                     3.3720000E-02, 1.9402300E-02, 1.1783100E-02, 7.4357000E-03,
                     2.2732000E-03, 8.8000000E-04, 4.5790000E-04, 2.3450000E-04,
                     1.5860000E-04, 1.1430000E-04, 7.1000000E-05])
        b = asarray([4.0000000E+02, 4.0500000E+02, 4.1000000E+02, 4.1500000E+02,
                     4.2000000E+02, 4.2500000E+02, 4.3000000E+02, 4.3500000E+02,
                     4.3650000E+02, 4.3800000E+02, 4.3950000E+02, 4.4100000E+02,
                     4.4250000E+02, 4.4400000E+02, 4.4550000E+02, 4.4700000E+02,
                     4.4850000E+02, 4.5000000E+02, 4.5150000E+02, 4.5300000E+02,
                     4.5450000E+02, 4.5600000E+02, 4.5750000E+02, 4.5900000E+02,
                     4.6050000E+02, 4.6200000E+02, 4.6350000E+02, 4.6500000E+02,
                     4.7000000E+02, 4.7500000E+02, 4.8000000E+02, 4.8500000E+02,
                     4.9000000E+02, 4.9500000E+02, 5.0000000E+02])

        vec = x[0] / x[1] * exp(-(b - x[2]) ** 2 / (2 * x[1] ** 2))
        return sum((a - vec) ** 2)


class EggCrate(Benchmark):

    """
    Egg Crate objective function.

    This class defines the Egg Crate global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{EggCrate}}(\\mathbf{x}) = x_1^2 + x_2^2 + 25 \\left[ \\sin^2(x_1) + \\sin^2(x_2) \\right]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.0] * self.N, [5.0] * self.N)

        self.global_optimum = [[0.0, 0.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return x[0] ** 2 + x[1] ** 2 + 25 * (sin(x[0]) ** 2 + sin(x[1]) ** 2)


class EggHolder(Benchmark):

    """
    Egg Holder objective function.

    This class defines the Egg Holder global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{EggHolder}}(\\mathbf{x}) = - x_{1} \\sin\\left(\\sqrt{\\lvert{x_{1} - x_{2} -47}\\rvert}\\right) - \\left(x_{2} + 47\\right) \\sin\\left(\\sqrt{\\left|{\\frac{1}{2} x_{1} + x_{2} + 47}\\right|}\\right)

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-512, 512]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -959.640662711` for :math:`\\mathbf{x} = [512, 404.2319]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-512.1] * self.N,
                           [512.0] * self.N)

        self.global_optimum = [[512.0, 404.2319]]
        self.fglob = -959.640662711
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        vec = (-(x[1:] + 47) * sin(sqrt(abs(x[1:] + x[:-1] / 2. + 47)))
               - x[:-1] * sin(sqrt(abs(x[:-1] - (x[1:] + 47)))))
        return sum(vec)


class ElAttarVidyasagarDutta(Benchmark):

    """
    El-Attar-Vidyasagar-Dutta objective function.

    This class defines the El-Attar-Vidyasagar-Dutta function global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{ElAttarVidyasagarDutta}}(\\mathbf{x}) = (x_1^2 + x_2 - 10)^2 + (x_1 + x_2^2 - 7)^2 + (x_1^2 + x_2^3 - 1)^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 1.712780354` for :math:`\\mathbf{x} = [3.40918683, -2.17143304]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = [(-4, 4), (-4, 4)]

        self.global_optimum = [[3.40918683, -2.17143304]]
        self.fglob = 1.712780354

    def fun(self, x, *args):
        self.nfev += 1

        return ((x[0] ** 2 + x[1] - 10) ** 2 + (x[0] + x[1] ** 2 - 7) ** 2
                + (x[0] ** 2 + x[1] ** 3 - 1) ** 2)


class Exp2(Benchmark):

    """
    Exp2 objective function.

    This class defines the Exp2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Exp2}}(\\mathbf{x}) = \\sum_{i=0}^9 \\left ( e^{-ix_1/10} - 5e^{-ix_2/10} -e^{-i/10} + 5e^{-i} \\right )^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 20]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = [1, 10.]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N, [20.0] * self.N)
        self.custom_bounds = [(0, 2), (0, 20)]

        self.global_optimum = [[1.0, 10.]]
        self.fglob = 0

    def fun(self, x, *args):
        self.nfev += 1

        i = arange(10.)
        vec = (exp(-i * x[0] / 10.) - 5 * exp(-i * x[1] / 10.) - exp(-i / 10)
               + 5 * exp(-i)) ** 2

        return sum(vec)


class Exponential(Benchmark):

    """
    Exponential objective function.

    This class defines the Exponential global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Exponential}}(\\mathbf{x}) = -e^{-0.5 \\sum_{i=1}^n x_i^2}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = -1` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-1.0] * self.N, [1.0] * self.N)

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = -1.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return -exp(-0.5 * sum(x ** 2.0))


class FreudensteinRoth(Benchmark):

    """
    FreudensteinRoth objective function.

    This class defines the Freudenstein & Roth global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{FreudensteinRoth}}(\\mathbf{x}) =  \\left\{x_1 - 13 + \\left[(5 - x_2)x_2 - 2 \\right] x_2 \\right\}^2 + \\left \{x_1 - 29 + \\left[(x_2 + 1)x_2 - 14 \\right] x_2 \\right\}^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [5, 4]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.custom_bounds = [(-3, 3), (-5, 5)]

        self.global_optimum = [[5.0, 4.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        f1 = (-13.0 + x[0] + ((5.0 - x[1]) * x[1] - 2.0) * x[1]) ** 2
        f2 = (-29.0 + x[0] + ((x[1] + 1.0) * x[1] - 14.0) * x[1]) ** 2

        return f1 + f2


class Gear(Benchmark):

    """
    Gear objective function.

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

        self._bounds = zip([12.0] * self.N, [60.0] * self.N)
        self.global_optimum = [[16, 19, 43, 49]]
        self.fglob = 2.7e-12

    def fun(self, x, *args):
        self.nfev += 1

        return (1. / 6.931
                - floor(x[0]) * floor(x[1]) / floor(x[2]) / floor(x[3])) ** 2


class Giunta(Benchmark):

    """
    Giunta objective function.

    This class defines the Giunta global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Giunta}}(\\mathbf{x}) = 0.6 + \\sum_{i=1}^{n} \\left[\\sin^{2}\\left(1 - \\frac{16}{15} x_i\\right) - \\frac{1}{50} \\sin\\left(4 - \\frac{64}{15} x_i\\right) - \\sin\\left(1 - \\frac{16}{15} x_i\\right)\\right]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0.06447042053690566` for :math:`\\mathbf{x} = [0.4673200277395354, 0.4673200169591304]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-1.0] * self.N, [1.0] * self.N)

        self.global_optimum = [[0.4673200277395354, 0.4673200169591304]]
        self.fglob = 0.06447042053690566

    def fun(self, x, *args):
        self.nfev += 1

        arg = 16 * x / 15.0 - 1
        return 0.6 + sum(sin(arg) + sin(arg) ** 2 + sin(4 * arg) / 50)


class GoldsteinPrice(Benchmark):

    """
    Goldstein-Price objective function.

    This class defines the Goldstein-Price global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{GoldsteinPrice}}(\\mathbf{x}) = \\left[ 1+(x_1+x_2+1)^2(19-14x_1+3x_1^2-14x_2+6x_1x_2+3x_2^2) \\right] \\left[ 30+(2x_1-3x_2)^2(18-32x_1+12x_1^2+48x_2-36x_1x_2+27x_2^2) \\right]


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-2, 2]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 3` for :math:`\\mathbf{x} = [0, -1]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-2.0] * self.N, [2.0] * self.N)

        self.global_optimum = [[0., -1.]]
        self.fglob = 3.0

    def fun(self, x, *args):
        self.nfev += 1

        a = 1 + (x[0] + x[1] + 1) ** 2 * \
            (19 - 14 * x[0] + 3 * x[0] ** 2 -
             14 * x[1] + 6 * x[0] * x[1] + 3 * x[1] ** 2)
        b = 30 + (2 * x[0] - 3 * x[1]) ** 2 * \
            (18 - 32 * x[0] + 12 * x[0] ** 2 +
             48 * x[1] - 36 * x[0] * x[1] + 27 * x[1] ** 2)
        return a * b


class Griewank(Benchmark):

    """
    Griewank objective function.

    This class defines the Griewank global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Griewank}}(\\mathbf{x}) = \\frac{1}{4000}\\sum_{i=1}^n x_i^2 - \\prod_{i=1}^n\\cos\\left(\\frac{x_i}{\\sqrt{i}}\\right) + 1

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-600, 600]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-600.0] * self.N,
                           [600.0] * self.N)
        self.custom_bounds = [(-50, 50), (-50, 50)]

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        i = arange(1., np.size(x))
        return sum(x ** 2 / 4000) - prod(cos(x / sqrt(i))) + 1


class Gulf(Benchmark):

    """
    Gulf objective function.

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

        self._bounds = zip([0.0] * self.N, [50.0] * self.N)

        self.global_optimum = [[50.0, 25.0, 1.5]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        m = 99.
        i = arange(1., m + 1)
        y = 25 + (-50 * log(i / 100.)) ** (2 / 3.)
        vec = (exp(-((abs(y - x[1])) ** x[2] / x[0])) - i / 100.)
        return sum(vec ** 2)


class Hansen(Benchmark):

    """
    Hansen objective function.

    This class defines the Hansen global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Hansen}}(\\mathbf{x}) = \\left[ \\sum_{i=0}^4(i+1)\\cos(ix_1+i+1)\\right ] \\left[\\sum_{j=0}^4(j+1)\\cos[(j+2)x_2+j+1])\\right ]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -176.54179` for :math:`\\mathbf{x} = [-7.58989583, -7.70831466]`.

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

    """
    Hartmann3 objective function.

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

        self._bounds = zip([0.0] * self.N, [1.0] * self.N)

        self.global_optimum = [[0.11461292,  0.55564907,  0.85254697]]
        self.fglob = -3.8627821478

    def fun(self, x, *args):
        self.nfev += 1

        a = asarray([[3.0,  0.1,  3.0,  0.1],
                     [10.0, 10.0, 10.0, 10.0],
                     [30.0, 35.0, 30.0, 35.0]])
        p = asarray([[0.36890, 0.46990, 0.10910, 0.03815],
                     [0.11700, 0.43870, 0.87320, 0.57430],
                     [0.26730, 0.74700, 0.55470, 0.88280]])
        c = asarray([1.0, 1.2, 3.0, 3.2])

        XX = np.atleast_2d(x).T
        d = sum(a * (XX - p) ** 2, axis=0)
        return -sum(c * exp(-d))


class Hartmann6(Benchmark):

    """
    Hartmann6 objective function.

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

        self._bounds = zip([0.0] * self.N, [1.0] * self.N)

        self.global_optimum = [[0.20168952, 0.15001069,
                               0.47687398, 0.27533243, 0.31165162, 0.65730054]]
        self.fglob = -3.32236801141551

    def fun(self, x, *args):
        self.nfev += 1

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

        XX = np.atleast_2d(x).T
        d = sum(a * (XX - p) ** 2, axis=0)
        return -sum(c * exp(-d))

        return -sum(c * exp(-d))


class HelicalValley(Benchmark):

    """
    HelicalValley objective function.

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

        self._bounds = zip([-100.0] * self.N, [100] * self.N)

        self.global_optimum = [[1.0, 0.0, 0.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return (100 * ((x[2] - 10 * arctan2(x[1], x[0]) / 2 / pi) ** 2
                + (sqrt(x[0] ** 2 + x[1] ** 2) - 1) ** 2) + x[2] ** 2)


class HimmelBlau(Benchmark):

    """
    HimmelBlau objective function.

    This class defines the HimmelBlau global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{HimmelBlau}}(\\mathbf{x}) = (x_1^2 + x_2 - 11)^2 + (x_1 + x_2^2 -7)^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-6, 6]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [3, 2]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-6] * self.N, [6] * self.N)

        self.global_optimum = [[3.0, 2.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return (x[0] ** 2 + x[1] - 11) ** 2 + (x[0] + x[1] ** 2 - 7) ** 2


class HolderTable(Benchmark):

    """
    HolderTable objective function.

    This class defines the HolderTable global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{HolderTable}}(\\mathbf{x}) = - \\left|{e^{\\left|{1 - \\frac{\\sqrt{x_{1}^{2} + x_{2}^{2}}}{\\pi} }\\right|} \\sin\\left(x_{1}\\right) \\cos\\left(x_{2}\\right)}\\right|

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -19.20850256788675` for :math:`x_i = \\pm 9.664590028909654` for :math:`i=1,2`

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


class Holzman(Benchmark):

    """
    Holzman objective function.

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

        self._bounds = ([0.0, 100.0], [0.0, 25.6], [0.0, 5.0])

        self.global_optimum = [[50.0, 25.0, 1.5]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        val = 0
        i = arange(1, 101)
        t = 2 / 3.
        u = 25 + (-50 * log(0.01 * i)) ** t
        v = (u - x[1]) ** x[2]
        w = exp(-v / x[0])
        return sum(-0.01 * i + w)


class Hosaki(Benchmark):

    """
    Hosaki objective function.

    This class defines the Hosaki global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Hosaki}}(\\mathbf{x}) = \\left ( 1 - 8x_1 + 7x_1^2 - \\frac{7}{3}x_1^3 + \\frac{1}{4}x_1^4 \\right )x_2^2e^{-x_1}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0,
    10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -2.3458115` for :math:`\\mathbf{x} = [4,
    2]`.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N, [10.0] * self.N)
        self.custom_bounds = [(0, 5), (0, 5)]

        self.global_optimum = [[4, 2]]
        self.fglob = -2.3458115

    def fun(self, x, *args):
        self.nfev += 1

        val = (1 - 8 * x[0] + 7 * x[0] ** 2 - 7 / 3. * x[0] ** 3
               + 0.25 * x[0] ** 4)
        return val * x[1] ** 2 * exp(-x[1])


class Infinity(Benchmark):

    """
    Infinity objective function.

    This class defines the Infinity global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Infinity}}(\\mathbf{x}) = \\sum_{i=1}^{n} x_i^{6} \\left [ \\sin\\left ( \\frac{1}{x_i} \\right )+2 \\right ]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-1.0] * self.N, [1.0] * self.N)

        self.global_optimum = [[1e-16 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum(x ** 6.0 * (sin(1.0 / x) + 2.0))


class JennrichSampson(Benchmark):

    """
    Jennrich-Sampson objective function.

    This class defines the Jennrich-Sampson global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{JennrichSampson}}(\\mathbf{x}) = \\sum_{i=1}^{10} \\left [2 + 2i - (e^{ix_1} + e^{ix_2}) \\right ]^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 124.3621824` for :math:`\\mathbf{x} = [0.257825, 0.257825]`.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-1.0] * self.N, [1.0] * self.N)

        self.global_optimum = [[0.257825, 0.257825]]
        self.custom_bounds = [(-1, 0.34), (-1, 0.34)]
        self.fglob = 124.3621824

    def fun(self, x, *args):
        self.nfev += 1

        i = arange(1, 11)
        return sum((2 + 2 * i - (exp(i * x[0]) + exp(i * x[1]))) ** 2)


class Judge(Benchmark):

    """
    Judge objective function.

    This class defines the Judge global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Judge}}(\\mathbf{x}) = \\sum_{i=1}^{20} \\left [ \\left (x_1 + A_i x_2 + B x_2^2 \\right ) - C_i \\right ]^2

    Where, in this exercise:

    .. math::

        \\begin{cases} C = [4.284, 4.149, 3.877, 0.533, 2.211, 2.389, 2.145,  3.231, 1.998, 1.379, 2.106, 1.428, 1.011, 2.179, 2.858, 1.388, 1.651, 1.593, 1.046, 2.152] \\\\
        A = [0.286, 0.973, 0.384, 0.276, 0.973, 0.543, 0.957, 0.948, 0.543, 0.797, 0.936, 0.889, 0.006, 0.828, 0.399, 0.617, 0.939, 0.784, 0.072, 0.889] \\\\
        B = [0.645, 0.585, 0.310, 0.058, 0.455, 0.779, 0.259, 0.202, 0.028, 0.099, 0.142, 0.296, 0.175, 0.180, 0.842, 0.039, 0.103, 0.620, 0.158, 0.704] \\end{cases}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 16.0817307` for :math:`\\mathbf{x} = [0.86479, 1.2357]`.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[0.86479, 1.2357]]
        self.custom_bounds = [(-2.0, 2.0), (-2.0, 2.0)]
        self.fglob = 16.0817307

    def fun(self, x, *args):
        self.nfev += 1

        C = asarray([4.284, 4.149, 3.877, 0.533, 2.211, 2.389, 2.145,  3.231,
                     1.998, 1.379, 2.106, 1.428, 1.011, 2.179, 2.858, 1.388,
                     1.651, 1.593, 1.046, 2.152])

        A = asarray([0.286, 0.973, 0.384, 0.276, 0.973, 0.543, 0.957, 0.948,
                     0.543, 0.797, 0.936, 0.889, 0.006, 0.828, 0.399, 0.617,
                     0.939, 0.784, 0.072, 0.889])

        B = asarray([0.645, 0.585, 0.310, 0.058, 0.455, 0.779, 0.259, 0.202,
                     0.028, 0.099, 0.142, 0.296, 0.175, 0.180, 0.842, 0.039,
                     0.103, 0.620, 0.158, 0.704])

        return sum(((x[0] + x[1] * A + (x[1] ** 2.0) * B) - C) ** 2.0)


class Katsuura(Benchmark):

    """
    Katsuura objective function.

    This class defines the Katsuura global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Katsuura}}(\\mathbf{x}) = \\prod_{i=0}^{n-1} \\left [ 1 + (i+1) \\sum_{k=1}^{d} \\lfloor (2^k x_i) \\rfloor 2^{-k} \\right ]

    Where, in this exercise, :math:`d = 32`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 100]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 1` for :math:`x_i = 0` for :math:`i=1,...,n`.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N, [100.0] * self.N)

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.custom_bounds = [(0, 1), (0, 1)]
        self.fglob = 1.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        d = 32
        k = np.atleast_2d(arange(1, d + 1)).T
        i = arange(0., self.N * 1.)
        inner = floor(2 ** k * x) * (2. ** (-k))
        return prod(sum(inner, axis=0) * (i + 1) + 1)


class Keane(Benchmark):

    """
    Keane objective function.

    This class defines the Keane global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Keane}}(\\mathbf{x}) = \\frac{\\sin^2(x_1 - x_2)\\sin^2(x_1 + x_2)}{\\sqrt{x_1^2 + x_2^2}}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0.0` for :math:`\\mathbf{x} = [7.85396153, 7.85396135]`.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[7.85396153, 7.85396135]]
        self.custom_bounds = [(-1, 0.34), (-1, 0.34)]
        self.fglob = 0.

    def fun(self, x, *args):
        self.nfev += 1

        val = sin(x[0] - x[1]) ** 2 * sin(x[0] + x[1]) ** 2
        return val / sqrt(x[0] ** 2 + x[1] ** 2)


class Kowalik(Benchmark):

    """
    Kowalik objective function.

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
    # TODO, this is a NIST regression standard dataset

    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.0] * self.N, [5.0] * self.N)
        self.global_optimum = [[0.192833, 0.190836, 0.123117, 0.135766]]
        self.fglob = 0.00030748610

    def fun(self, x, *args):
        self.nfev += 1

        b = asarray([4.0, 2.0, 1.0, 1 / 2.0, 1 / 4.0, 1 / 6.0, 1 / 8.0,
                     1 / 10.0, 1 / 12.0, 1 / 14.0, 1 / 16.0])
        a = asarray([0.1957, 0.1947, 0.1735, 0.1600, 0.0844,
                     0.0627, 0.0456, 0.0342, 0.0323, 0.0235,
                     0.0246])

        vec = a - (x[0] * (b ** 2 + b * x[1])
                   / (b ** 2 + b * x[2] + x[3]))
        return sum(vec ** 2)


class Langermann(Benchmark):

    """
    Langermann objective function.

    This class defines the Langermann global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Langermann}}(\\mathbf{x}) = - \\sum_{i=1}^{5} \\frac{c_i \\cos\\left\{\\pi \\left[\\left(x_{1}- a_i\\right)^{2} + \\left(x_{2} - b_i \\right)^{2}\\right]\\right\}}{e^{\\frac{\\left( x_{1} - a_i\\right)^{2} + \\left( x_{2} - b_i\\right)^{2}}{\\pi}}}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -5.1621259` for :math:`\\mathbf{x} = [2.00299219, 1.006096]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[2.00299219, 1.006096]]
        self.fglob = -5.1621259

    def fun(self, x, *args):
        self.nfev += 1

        a = [3, 5, 2, 1, 7]
        b = [5, 2, 1, 4, 9]
        c = [1, 2, 5, 2, 3]

        return (-sum(c * exp(-(1 / pi) * ((x[0] - a) ** 2 +
                    (x[1] - b) ** 2)) * cos(pi * ((x[0] - a) ** 2
                                            + (x[1] - b) ** 2))))


class LennardJones(Benchmark):

    """
    LennardJones objective function.

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

        self._bounds = zip([-4.0] * self.N, [4.0] * self.N)

        self.global_optimum = [[]]

        minima = [
            -1.0, -3.0, -6.0, -9.103852, -12.712062, -
            16.505384, -19.821489, -24.113360, -28.422532,
            -32.765970, -37.967600, -44.326801, -
            47.845157, -52.322627, -56.815742, -61.317995,
            -66.530949, -72.659782, -77.1777043]

        k = int(dimensions / 3)
        self.fglob = minima[k - 2]
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        k = int(self.N / 3)
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


class Leon(Benchmark):

    """
    Leon objective function.

    This class defines the Leon global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Leon}}(\\mathbf{x}) = \\left(1 - x_{1}\\right)^{2} + 100 \\left(x_{2} - x_{1}^{2} \\right)^{2}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1.2, 1.2]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 1` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-1.2] * self.N, [1.2] * self.N)

        self.global_optimum = [[1 for _ in range(self.N)]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return 100 * (x[1] - x[0] ** 2.0) ** 2.0 + (1 - x[0]) ** 2.0


class Levy03(Benchmark):

    """
    Levy 3 objective function.

    This class defines the Levy 3 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Levy03}}(\\mathbf{x}) = \\sin^2(\\pi y_1)+\\sum_{i=1}^{n-1}(y_i-1)^2[1+10\\sin^2(\\pi y_{i+1})]+(y_n-1)^2

    Where, in this exercise:

    .. math::

        y_i=1+\\frac{x_i-1}{4}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 1` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.custom_bounds = [(-5, 5), (-5, 5)]

        self.global_optimum = [[1 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        y = 1 + (x - 1) / 4
        v = sum((y[: -1] - 1) ** 2 * (1 + 10 * sin(pi * y[1:]) ** 2))
        return sin(pi * y[0]) ** 2 + v + (y[-1] - 1) ** 2


class Levy05(Benchmark):

    """
    Levy 5 objective function.

    This class defines the Levy 5 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Levy05}}(\\mathbf{x}) = \\sum_{i=1}^{5} i \\cos \\left[(i-1)x_1 + i \\right] \\times \\sum_{j=1}^{5} j \\cos \\left[(j+1)x_2 + j \\right] + (x_1 + 1.42513)^2 + (x_2 + 0.80032)^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = -176.1375779` for :math:`\\mathbf{x} = [-1.30685, -1.42485]`.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.custom_bounds = ([-2.0, 2.0], [-2.0, 2.0])

        self.global_optimum = [[-1.30685, -1.42485]]
        self.fglob = -176.1375779

    def fun(self, x, *args):
        self.nfev += 1

        i = arange(1, 6)
        a = i * cos((i - 1) * x[0] + i)
        b = i * cos((i + 1) * x[1] + i)

        return sum(a) * sum(b) + (x[0] + 1.42513) ** 2 + (x[1] + 0.80032) ** 2


class Levy13(Benchmark):

    """
    Levy13 objective function.

    This class defines the Levy13 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Levy13}}(\\mathbf{x}) = \\left(x_{1} -1\\right)^{2} \\left[\sin^{2}\\left(3 \\pi x_{2}\\right) + 1\\right] + \\left(x_{2} -1\\right)^{2} \\left[\\sin^{2}\\left(2 \\pi x_{2}\\right) + 1\\right] + \\sin^{2}\\left(3 \\pi x_{1}\\right)


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 1` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.custom_bounds = [(-5, 5), (-5, 5)]

        self.global_optimum = [[1 for _ in range(self.N)]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        u = sin(3 * pi * x[0]) ** 2
        v = (x[0] - 1) ** 2 * (1 + (sin(3 * pi * x[1])) ** 2)
        w = (x[1] - 1) ** 2 * (1 + (sin(2 * pi * x[1])) ** 2)
        return (u + v + w)


class Matyas(Benchmark):

    """
    Matyas objective function.

    This class defines the Matyas global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Matyas}}(\\mathbf{x}) = 0.26(x_1^2 + x_2^2) - 0.48x_1x_2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return 0.26 * (x[0] ** 2 + x[1] ** 2) - 0.48 * x[0] * x[1]


class McCormick(Benchmark):

    """
    McCormick objective function.

    This class defines the McCormick global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{McCormick}}(\\mathbf{x}) = - x_{1} + 2 x_{2} + \\left(x_{1} - x_{2}\\right)^{2} + \\sin\\left(x_{1} + x_{2}\\right) + 1


    Here, :math:`n` represents the number of dimensions and :math:`x_1 \\in [-1.5, 4]`, :math:`x_2 \\in [-3, 4]`.

    *Global optimum*: :math:`f(x_i) = -1.913222954981037` for :math:`\\mathbf{x} = [-0.5471975602214493, -1.547197559268372]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = [(-1.5, 4.0), (-3.0, 4.0)]

        self.global_optimum = [[-0.5471975602214493, -1.547197559268372]]
        self.fglob = -1.913222954981037

    def fun(self, x, *args):
        self.nfev += 1

        return (sin(x[0] + x[1]) + (x[0] - x[1]) ** 2 - 1.5 * x[0]
                + 2.5 * x[1] + 1)


class Meyer(Benchmark):

    """
    Meyer objective function.

    """

    # TODO, this is a NIST regression standard dataset
    def __init__(self, dimensions=3):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0., 100., 100.],
                           [1, 1000., 500.])
        self.global_optimum = [[5.6096364710e-3, 6.1813463463e3,
                                3.4522363462e2]]
        self.fglob = 8.7945855171e1

    def fun(self, x, *args):
        self.nfev += 1

        a = asarray([3.478E+04, 2.861E+04, 2.365E+04, 1.963E+04, 1.637E+04,
                     1.372E+04, 1.154E+04, 9.744E+03, 8.261E+03, 7.030E+03,
                     6.005E+03, 5.147E+03, 4.427E+03, 3.820E+03, 3.307E+03,
                     2.872E+03])
        b = asarray([5.000E+01, 5.500E+01, 6.000E+01, 6.500E+01, 7.000E+01,
                     7.500E+01, 8.000E+01, 8.500E+01, 9.000E+01, 9.500E+01,
                     1.000E+02, 1.050E+02, 1.100E+02, 1.150E+02, 1.200E+02,
                     1.250E+02])

        vec = x[0] * exp(x[1] / (b + x[2]))

        return sum((a - vec) ** 2)


class Michalewicz(Benchmark):

    """
    Michalewicz objective function.

    This class defines the Michalewicz global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Michalewicz}}(\\mathbf{x}) = - \\sum_{i=1}^{2} \\sin\\left(x_i\\right) \\sin^{2 m}\\left(\\frac{i x_i^{2}}{\\pi}\\right)


    Where, in this exercise, :math:`m = 10`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, \\pi]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -1.8013` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N, [pi] * self.N)

        self.global_optimum = [[2.20290555, 1.570796]]
        self.fglob = -1.8013

    def fun(self, x, *args):
        self.nfev += 1

        m = 10.0
        i = arange(1, self.N + 1)
        return -sum(sin(x) * sin(i * x ** 2 / pi) ** (2 * m))


class MieleCantrell(Benchmark):

    """
    Miele-Cantrell objective function.

    This class defines the Miele-Cantrell global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{MieleCantrell}}(\\mathbf{x}) = (e^{-x_1} - x_2)^4 + 100(x_2 - x_3)^6 + \\tan^4(x_3 - x_4) + x_1^8


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,...,4`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [0, 1, 1, 1]`

    """

    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-1.0] * self.N, [1.0] * self.N)

        self.global_optimum = [[0.0, 1.0, 1.0, 1.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return ((exp(-x[0]) - x[1]) ** 4 + 100 * (x[1] - x[2]) ** 6
                + tan(x[2] - x[3]) ** 4 + x[0] ** 8)


class Mishra01(Benchmark):

    """
    Mishra 1 objective function.

    This class defines the Mishra 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra01}}(\\mathbf{x}) = (1 + x_n)^{x_n} \\hspace{10pt} ; \\hspace{10pt} x_n = n - \\sum_{i=1}^{n-1} x_i


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 1]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 2` for :math:`x_i = 1` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N,
                           [1.0 + 1e-9] * self.N)

        self.global_optimum = [[1.0 for _ in range(self.N)]]
        self.fglob = 2.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        xn = self.N - sum(x[0:-1])
        return (1 + xn) ** xn


class Mishra02(Benchmark):

    """
    Mishra 2 objective function.

    This class defines the Mishra 2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra02}}(\\mathbf{x}) = (1 + x_n)^{x_n} \\hspace{10pt} ; \\hspace{10pt} x_n = n - \\sum_{i=1}^{n-1} \\frac{(x_i + x_{i+1})}{2}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 1]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 2` for :math:`x_i = 1` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N,
                           [1.0 + 1e-9] * self.N)

        self.global_optimum = [[1.0 for _ in range(self.N)]]
        self.fglob = 2.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        xn = self.N - sum((x[:-1] + x[1:]) / 2.0)
        return (1 + xn) ** xn


class Mishra03(Benchmark):

    """
    Mishra 3 objective function.

    This class defines the Mishra 3 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra03}}(\\mathbf{x}) = \\sqrt{\\lvert \\cos{\\sqrt{\\lvert x_1^2 + x_2^2 \\rvert}} \\rvert} + 0.01(x_1 + x_2)


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -0.1999` for :math:`x_i = {-9.99378322, -9.99918927}`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[-9.99378322, -9.99918927]]
        self.fglob = -0.19990562

    def fun(self, x, *args):
        self.nfev += 1

        return ((0.01 * (x[0] + x[1])
                + sqrt(abs(cos(sqrt(abs(x[0] ** 2 + x[1] ** 2)))))))


class Mishra04(Benchmark):

    """
    Mishra 4 objective function.

    This class defines the Mishra 4 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra04}}(\\mathbf{x}) = \\sqrt{\\lvert \\sin{\\sqrt{\\lvert x_1^2 + x_2^2 \\rvert}} \\rvert} + 0.01(x_1 + x_2)


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -0.17767` for :math:`x_i = {-8.71499636, -9.0533148}`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[-8.71499636, -9.0533148]]
        self.fglob = -0.17767

    def fun(self, x, *args):
        self.nfev += 1

        return ((0.01 * (x[0] + x[1])
                + sqrt(abs(sin(sqrt(abs(x[0] ** 2 + x[1] ** 2)))))))


class Mishra05(Benchmark):

    """
    Mishra 5 objective function.

    This class defines the Mishra 5 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra05}}(\\mathbf{x}) = \\left [ \\sin^2 ((\\cos(x_1) + \\cos(x_2))^2) + \\cos^2 ((\\sin(x_1) + \\sin(x_2))^2) + x_1 \\right ]^2 + 0.01(x_1 + x_2)


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -0.119829` for :math:`\\mathbf{x} = [-1.98682, -10]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[-1.98682, -10.0]]
        self.fglob = -0.119829

    def fun(self, x, *args):
        self.nfev += 1

        return (0.01 * (x[0] + x[1])
                + (sin((cos(x[0]) + cos(x[1])) ** 2) ** 2
                   + cos((sin(x[0]) + sin(x[1])) ** 2) ** 2 + x[0]) ** 2)


class Mishra06(Benchmark):

    """
    Mishra 6 objective function.

    This class defines the Mishra 6 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra06}}(\\mathbf{x}) = -\\log{\\left [ \\sin^2 ((\\cos(x_1) + \\cos(x_2))^2) - \\cos^2 ((\\sin(x_1) + \\sin(x_2))^2) + x_1 \\right ]^2} + 0.01 \\left[(x_1 -1)^2 + (x_2 - 1)^2 \\right]


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -2.28395` for :math:`\\mathbf{x} = [2.88631, 1.82326]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[2.88631, 1.82326]]
        self.fglob = -2.28395

    def fun(self, x, *args):
        self.nfev += 1

        a = 0.1 * ((x[0] - 1) ** 2 + (x[1] - 1) ** 2)
        u = (cos(x[0]) + cos(x[1])) ** 2
        v = (sin(x[0]) + sin(x[1])) ** 2
        return a - log((sin(u) ** 2 - cos(v) ** 2 + x[0]) ** 2)


class Mishra07(Benchmark):

    """
    Mishra 7 objective function.

    This class defines the Mishra 7 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra07}}(\\mathbf{x}) = \\left [\\prod_{i=1}^{n} x_i - n! \\right]^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = \\sqrt{n}` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.custom_bounds = [(-2, 2), (-2, 2)]
        self.global_optimum = [[sqrt(self.N)
                               for i in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return (prod(x) - factorial(self.N)) ** 2.0


class Mishra08(Benchmark):

    """
    Mishra 8 objective function.

    This class defines the Mishra 8 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra08}}(\\mathbf{x}) = 0.001 \\left[\\lvert x_1^{10} - 20x_1^9 + 180x_1^8 - 960 x_1^7 + 3360x_1^6 - 8064x_1^5 + 13340x_1^4 - 15360x_1^3 + 11520x_1^2 - 5120x_1 + 2624 \\rvert \\lvert x_2^4 + 12x_2^3 + 54x_2^2 + 108x_2 + 81 \\rvert \\right]^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [2, -3]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.custom_bounds = [(1.0, 2.0), (-4.0, 1.0)]
        self.global_optimum = [[2.0, 3.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        val = abs(x[0] ** 10 - 20 * x[0] ** 9 + 180 * x[0] ** 8
                  - 960 * x[0] ** 7 + 3360 * x[0] ** 6 - 8064 * x[0] ** 5
                  + 13340 * x[0] ** 4 - 15360 * x[0] ** 3 + 11520 * x[0] ** 2
                  - 5120 * x[0] + 2624)
        val *= abs(x[1] ** 2 + 12 * x[1] ** 2 +
                   54 * x[1] ** 2 + 108 * x[1] + 81)
        return 0.001 * val ** 2


class Mishra09(Benchmark):

    """
    Mishra 9 objective function.

    This class defines the Mishra 9 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra09}}(\\mathbf{x}) = \\left[ ab^2c + abc^2 + b^2 + (x_1 + x_2 - x_3)^2 \\right]^2


    Where, in this exercise:

    .. math::

        \\begin{cases} a = 2x_1^3 + 5x_1x_2 + 4x_3 - 2x_1^2x_3 - 18 \\\\
        b = x_1 + x_2^3 + x_1x_2^2 + x_1x_3^2 - 22 \\\\
        c = 8x_1^2 + 2x_2x_3 + 2x_2^2 + 3x_2^3 - 52 \\end{cases}


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2,3`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [1, 2, 3]`

    """

    def __init__(self, dimensions=3):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.global_optimum = [[1.0, 2.0, 3.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        a = (2 * x[0] ** 3 + 5 * x[0] * x[1]
             + 4 * x[2] - 2 * x[0] ** 2 * x[2] - 18)
        b = x[0] + x[1] ** 3 + x[0] * x[1] ** 2 + x[0] * x[2] ** 2 - 22.0
        c = 8 * x[0] ** 2 + 2 * x[1] * x[2] + \
            2 * x[1] ** 2 + 3 * x[1] ** 3 - 52

        return (a * c * b ** 2 + a * b * c ** 2 + b ** 2
                + (x[0] + x[1] - x[2]) ** 2) ** 2


class Mishra10(Benchmark):

    """
    Mishra 10 objective function.

    This class defines the Mishra 10 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::
    TODO - int(x) should be used instead of floor(x)!!!!!
       f_{\\text{Mishra10}}(\\mathbf{x}) = \\left[ \\lfloor x_1 \\perp x_2 \\rfloor - \\lfloor x_1 \\rfloor - \\lfloor x_2 \\rfloor \\right]^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [2, 2]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.global_optimum = [[2.0, 2.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        x1, x2 = int(x[0]), int(x[1])
#         TODO rewrite equation above with nint(x)
        f1 = x1 + x2
        f2 = x1 * x2
        return (f1 - f2) ** 2.0


class Mishra11(Benchmark):

    """
    Mishra 11 objective function.

    This class defines the Mishra 11 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Mishra11}}(\\mathbf{x}) = \\left [ \\frac{1}{n} \\sum_{i=1}^{n} \\lvert x_i \\rvert - \\left(\\prod_{i=1}^{n} \\lvert x_i \\rvert \\right )^{\\frac{1}{n}} \\right]^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.custom_bounds = [(-3, 3), (-3, 3)]

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        N = self.N
        return ((1.0 / N) * sum(abs(x)) - (prod(abs(x))) ** 1.0 / N) ** 2.0


class MultiModal(Benchmark):

    """
    MultiModal objective function.

    This class defines the MultiModal global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{MultiModal}}(\\mathbf{x}) = \\left( \\sum_{i=1}^n \\lvert x_i \\rvert \\right) \\left( \\prod_{i=1}^n \\lvert x_i \\rvert \\right)


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.custom_bounds = [(-5, 5), (-5, 5)]

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum(abs(x)) * prod(abs(x))


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
#         self.global_optimum = [self.a[0: self.N]]
        self.global_optimum = [[1.09263477,  1.39263477]]

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
#     TODO: is equation correct? If it's not then a factor of 10 needs to appear
#     before first sin**2 term.
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


class Qing(Benchmark):

    """
    Qing objective function.

    This class defines the Qing global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Qing}}(\\mathbf{x}) = \\sum_{i=1}^{n} (x_i^2 - i)^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-500, 500]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = \\pm \\sqrt(i)` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-500.0] * self.N,
                           [500.0] * self.N)
        self.custom_bounds = [(-2, 2), (-2, 2)]
        self.global_optimum = [[sqrt(_) for _ in range(1, self.N + 1)]]
        self.fglob = 0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        i = arange(1, self.N + 1)
        return sum((x ** 2.0 - i) ** 2.0)


class Quadratic(Benchmark):

    """
    Quadratic objective function.

    This class defines the Quadratic global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Quadratic}}(\\mathbf{x}) = -3803.84 - 138.08x_1 - 232.92x_2 + 128.08x_1^2 + 203.64x_2^2 + 182.25x_1x_2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -3873.72418` for :math:`\\mathbf{x} = [0.19388, 0.48513]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.custom_bounds = [(0, 1), (0, 1)]
        self.global_optimum = [[0.19388, 0.48513]]
        self.fglob = -3873.72418
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return (-3803.84 - 138.08 * x[0] - 232.92 * x[1] + 128.08 * x[0] ** 2.0
                + 203.64 * x[1] ** 2.0 + 182.25 * x[0] * x[1])


class Quintic(Benchmark):

    """
    Quintic objective function.

    This class defines the Quintic global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Quintic}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\left|{x_{i}^{5} - 3 x_{i}^{4} + 4 x_{i}^{3} + 2 x_{i}^{2} - 10 x_{i} -4}\\right|


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = -1` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.custom_bounds = [(-2, 2), (-2, 2)]

        self.global_optimum = [[-1.0 for _ in range(self.N)]]
        self.fglob = 0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum(abs(x ** 5 - 3 * x ** 4 + 4 * x ** 3 + 2 * x ** 2
                       - 10 * x - 4))


class Rana(Benchmark):

    """
    Rana objective function.

    This class defines the Rana global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Rana}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\left[x_{i} \\sin\\left(\\sqrt{\\lvert{x_{1} - x_{i} + 1}\\rvert}\\right) \\cos\\left(\\sqrt{\\lvert{x_{1} + x_{i} + 1}\\rvert}\\right) + \\left(x_{1} + 1\\right) \\sin\\left(\\sqrt{\\lvert{x_{1} + x_{i} + 1}\\rvert}\\right) \\cos\\left(\\sqrt{\\lvert{x_{1} - x_{i} + 1}\\rvert}\\right)\\right]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-500.000001, 500.000001]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = -928.5478` for :math:`x_i = -500` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-500.000001] * self.N,
                           [500.000001] * self.N)

        self.global_optimum = [[-300.3376, 500.]]
        self.fglob = -500.8021602966615
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        t1 = sqrt(abs(x[1:] + x[: -1] + 1))
        t2 = sqrt(abs(x[1:] - x[: -1] + 1))
        return sum((x[1:] + 1) * cos(t2) * sin(t1) + x[:-1] * cos(t1) * sin(t2))


class Rastrigin(Benchmark):

    """
    Rastrigin objective function.

    This class defines the Rastrigin global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Rastrigin}}(\\mathbf{x}) = 10n \\sum_{i=1}^n \\left[ x_i^2 - 10 \\cos(2\\pi x_i) \\right]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5.12, 5.12]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([-5.12] * self.N, [5.12] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return 10.0 * self.N + sum(x ** 2.0 - 10.0 * cos(2.0 * pi * x))


class Ratkowsky01(Benchmark):

    """
    Ratkowsky objective function.

    """

    # TODO, this is a NIST regression standard dataset
    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0., 1., 0., 0.1],
                           [1000, 20., 3., 6.])
        self.global_optimum = [[6.996415127e2, 5.2771253025, 7.5962938329e-1,
                                1.2792483859]]
        self.fglob = 8.786404908e3

    def fun(self, x, *args):
        self.nfev += 1

        a = asarray(
            [16.08, 33.83, 65.80, 97.20, 191.55, 326.20, 386.87, 520.53,
             590.03, 651.92, 724.93, 699.56, 689.96, 637.56, 717.41])
        b = arange(1, 16.)

        vec = x[0] / ((1 + exp(x[1] - x[2] * b)) ** (1 / x[3]))

        return sum((a - vec) ** 2)


class Ratkowsky02(Benchmark):

    """
    Ratkowsky02 objective function.

    This class defines the Ratkowsky 2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::
#TODO fix equation
       f_{\\text{Ratkowsky02}}(\\mathbf{x}) = \\sum_{i=1}^2 -e^{-2 \\log 2 (\\frac{x_i-0.1}{0.8})^2} \\left[\\sin^6(5 \\pi x_i) + 0.1\\cos^2(500 \\pi x_i) \\right]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 1]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 8.0565229338` for :math:`x_i = [7.2462237576e1, 2.6180768402, 6.7359200066e-2]`

    """

    # TODO, this is a NIST regression standard dataset
    def __init__(self, dimensions=4):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([10, 0.5,  0.01],
                           [200, 5., 0.5])
        self.global_optimum = [[7.2462237576e1, 2.6180768402, 6.7359200066e-2]]
        self.fglob = 8.0565229338

    def fun(self, x, *args):
        self.nfev += 1

        a = asarray([8.93, 10.8, 18.59, 22.33, 39.35, 56.11, 61.73, 64.62,
                     67.08])
        b = asarray([9., 14., 21., 28., 42., 57., 63., 70., 79.])
        vec = x[0] / (1 + exp(x[1] - x[2] * b))
        return sum((a - vec) ** 2)


class Ripple01(Benchmark):

    """
    Ripple 1 objective function.

    This class defines the Ripple 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Ripple01}}(\\mathbf{x}) = \\sum_{i=1}^2 -e^{-2 \\log 2 (\\frac{x_i-0.1}{0.8})^2} \\left[\\sin^6(5 \\pi x_i) + 0.1\\cos^2(500 \\pi x_i) \\right]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 1]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = -2.2` for :math:`x_i = 0.1` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([0.0] * self.N, [1.0] * self.N)

        self.global_optimum = [[0.1 for _ in range(self.N)]]
        self.fglob = -2.2

    def fun(self, x, *args):
        self.nfev += 1

        u = -2.0 * log(2.0) * ((x - 0.1) / 0.8) ** 2.0
        v = sin(5.0 * pi * x) ** 6.0 + 0.1 * cos(500.0 * pi * x) ** 2.0
        return sum(-exp(u) * v)


class Ripple25(Benchmark):

    """
    Ripple 25 objective function.

    This class defines the Ripple 25 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Ripple25}}(\\mathbf{x}) = \\sum_{i=1}^2 -e^{-2 \\log 2 (\\frac{x_i-0.1}{0.8})^2} \\left[\\sin^6(5 \\pi x_i) \\right]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 1]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = -2` for :math:`x_i = 0.1` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([0.0] * self.N, [1.0] * self.N)

        self.global_optimum = [[0.1 for _ in range(self.N)]]
        self.fglob = -2.0

    def fun(self, x, *args):
        self.nfev += 1

        u = -2.0 * log(2.0) * ((x - 0.1) / 0.8) ** 2.0
        v = sin(5.0 * pi * x) ** 6.0
        return sum(-exp(u) * v)


class Rosenbrock(Benchmark):

    """
    Rosenbrock objective function.

    This class defines the Rosenbrock global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Rosenbrock}}(\\mathbf{x}) = \\sum_{i=1}^{n-1} [100(x_i^2 - x_{i+1})^2 + (x_i - 1)^2]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 1` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.0] * self.N, [10.0] * self.N)
        self.custom_bounds = [(-2, 2), (-2, 2)]

        self.global_optimum = [[1 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum(100.0 * (x[1:] - x[:-1] ** 2.0) ** 2.0 + (1 - x[:-1]) ** 2.0)


class RosenbrockModified(Benchmark):

    """
    Modified Rosenbrock objective function.

    This class defines the Modified Rosenbrock global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{RosenbrockModified}}(\\mathbf{x}) = 74 + 100(x_2 - x_1^2)^2 + (1 - x_1)^2 - 400 e^{-\\frac{(x_1+1)^2 + (x_2 + 1)^2}{0.1}}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-2, 2]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 34.04024310` for :math:`\\mathbf{x} = [-0.90955374, -0.95057172]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-2.0] * self.N, [2.0] * self.N)
        self.custom_bounds = ([-1.0, 0.5], [-1.0, 1.0])

        self.global_optimum = [[-0.90955374, -0.95057172]]
        self.fglob = 34.040243106640844

    def fun(self, x, *args):
        self.nfev += 1

        a = 74 + 100. * (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2
        a -= 400 * exp(-((x[0] + 1.) ** 2 + (x[1] + 1.) ** 2) / 0.1)
        return a


class RotatedEllipse01(Benchmark):

    """
    Rotated Ellipse 1 objective function.

    This class defines the Rotated Ellipse 1 global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{RotatedEllipse01}}(\\mathbf{x}) = 7x_1^2 - 6 \\sqrt{3} x_1x_2 + 13x_2^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-500, 500]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [0, 0]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-500.0] * self.N,
                           [500.0] * self.N)
        self.custom_bounds = ([-2.0, 2.0], [-2.0, 2.0])

        self.global_optimum = [[0.0, 0.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return (7.0 * x[0] ** 2.0 - 6.0 * sqrt(3) * x[0] * x[1]
                + 13 * x[1] ** 2.0)


class RotatedEllipse02(Benchmark):

    """
    Rotated Ellipse 2 objective function.

    This class defines the Rotated Ellipse 2 global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{RotatedEllipse02}}(\\mathbf{x}) = x_1^2 - x_1x_2 + x_2^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-500, 500]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [0, 0]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-500.0] * self.N,
                           [500.0] * self.N)
        self.custom_bounds = ([-2.0, 2.0], [-2.0, 2.0])

        self.global_optimum = [[0.0, 0.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return x[0] ** 2.0 - x[0] * x[1] + x[1] ** 2.0


class Salomon(Benchmark):

    """
    Salomon objective function.

    This class defines the Salomon global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Salomon}}(\\mathbf{x}) = 1 - \\cos \\left (2 \\pi \\sqrt{\\sum_{i=1}^{n} x_i^2} \\right) + 0.1 \\sqrt{\\sum_{i=1}^n x_i^2}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = [(-50, 50), (-50, 50)]

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        u = sum(x ** 2)
        return 1 - cos(2 * pi * sqrt(u)) + 0.1 * sqrt(u)


class Sargan(Benchmark):

    """
    Sargan objective function.

    This class defines the Sargan global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Sargan}}(\\mathbf{x}) = \\sum_{i=1}^{n} n \\left (x_i^2 + 0.4 \\sum_{i \\neq j}^{n} x_ix_j \\right)

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
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

    """
    Schaffer 1 objective function.

    This class defines the Schaffer 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schaffer01}}(\\mathbf{x}) = 0.5 + \\frac{\\sin^2 (x_1^2 + x_2^2)^2 - 0.5}{1 + 0.001(x_1^2 + x_2^2)^2}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = [(-10, 10), (-10, 10)]

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        num = sin((x[0] ** 2 + x[1] ** 2)) ** 2 - 0.5
        den = (1 + 0.001 * (x[0] ** 2 + x[1] ** 2)) ** 2
        return 0.5 + num / den


class Schaffer02(Benchmark):

    """
    Schaffer 2 objective function.

    This class defines the Schaffer 2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schaffer02}}(\\mathbf{x}) = 0.5 + \\frac{\\sin^2 (x_1^2 - x_2^2)^2 - 0.5}{1 + 0.001(x_1^2 + x_2^2)^2}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = [(-10, 10), (-10, 10)]

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        num = sin((x[0] ** 2 - x[1] ** 2)) ** 2 - 0.5
        den = (1 + 0.001 * (x[0] ** 2 + x[1] ** 2)) ** 2
        return 0.5 + num / den


class Schaffer03(Benchmark):

    """
    Schaffer 3 objective function.

    This class defines the Schaffer 3 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schaffer03}}(\\mathbf{x}) = 0.5 + \\frac{\\sin^2 \\left( \\cos \\lvert x_1^2 - x_2^2 \\rvert \\right ) - 0.5}{1 + 0.001(x_1^2 + x_2^2)^2}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0.00156685` for :math:`\\mathbf{x} = [0, 1.253115]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = [(-10, 10), (-10, 10)]

        self.global_optimum = [[0.0, 1.253115]]
        self.fglob = 0.00156685

    def fun(self, x, *args):
        self.nfev += 1

        num = sin(cos(abs(x[0] ** 2 - x[1] ** 2))) ** 2 - 0.5
        den = (1 + 0.001 * (x[0] ** 2 + x[1] ** 2)) ** 2
        return 0.5 + num / den


class Schaffer04(Benchmark):

    """
    Schaffer 4 objective function.

    This class defines the Schaffer 4 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schaffer04}}(\\mathbf{x}) = 0.5 + \\frac{\\cos^2 \\left( \\sin(x_1^2 - x_2^2) \\right ) - 0.5}{1 + 0.001(x_1^2 + x_2^2)^2}^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0.292579` for :math:`\\mathbf{x} = [0, 1.253115]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = [(-10, 10), (-10, 10)]

        self.global_optimum = [[0.0, 1.253115]]
        self.fglob = 0.292579

    def fun(self, x, *args):
        self.nfev += 1

        num = cos(sin(abs(x[0] ** 2 - x[1] ** 2))) ** 2 - 0.5
        den = (1 + 0.001 * (x[0] ** 2 + x[1] ** 2)) ** 2
        return 0.5 + num / den


class SchmidtVetters(Benchmark):

    """
    Schmidt-Vetters objective function.

    This class defines the Schmidt-Vetters global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{SchmidtVetters}}(\\mathbf{x}) = \\frac{1}{1 + (x_1 - x_2)^2} + \\sin \\left(\\frac{\\pi x_2 + x_3}{2} \\right) + e^{\\left(\\frac{x_1+x_2}{x_2} - 2\\right)^2}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 10]` for :math:`i=1,2,3`.

    *Global optimum*: :math:`f(x_i) = 2.99643266` for :math:`x_i = [0.79876108,  0.79962581,  0.79848824]`

    """

    def __init__(self, dimensions=3):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([0.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[0.79876108,  0.79962581,  0.79848824]]
        self.fglob = 2.99643266

    def fun(self, x, *args):
        self.nfev += 1

        return (1 / (1 + (x[0] - x[1]) ** 2) + sin((pi * x[1] + x[2]) / 2)
                + exp(((x[0] + x[1]) / x[1] - 2) ** 2))


class Schwefel01(Benchmark):

    """
    Schwefel 1 objective function.

    This class defines the Schwefel 1 global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel01}}(\\mathbf{x}) = \\left(\\sum_{i=1}^n x_i^2 \\right)^{\\alpha}

    Where, in this exercise, :math:`\\alpha = \\sqrt{\\pi}`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = ([-4.0, 4.0], [-4.0, 4.0])

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        alpha = sqrt(pi)
        return (sum(x ** 2.0)) ** alpha


class Schwefel02(Benchmark):

    """
    Schwefel 2 objective function.

    This class defines the Schwefel 2 global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel02}}(\\mathbf{x}) = \\sum_{i=1}^n \\left(\\sum_{j=1}^i x_i \\right)^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = ([-4.0, 4.0], [-4.0, 4.0])

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        mat = np.repeat(np.atleast_2d(x), self.N, axis=0)
        inner = sum(np.tril(mat), axis=1)
        return sum(inner ** 2)


class Schwefel04(Benchmark):

    """
    Schwefel 4 objective function.

    This class defines the Schwefel 4 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel04}}(\\mathbf{x}) = \\sum_{i=1}^n \\left[(x_i - 1)^2 + (x_1 - x_i^2)^2 \\right]


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 1` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([0.0] * self.N, [10.0] * self.N)
        self.custom_bounds = ([0.0, 2.0], [0.0, 2.0])

        self.global_optimum = [[1.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum((x - 1.0) ** 2.0 + (x[0] - x ** 2.0) ** 2.0)


class Schwefel06(Benchmark):

    """
    Schwefel 6 objective function.

    This class defines the Schwefel 6 global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel06}}(\\mathbf{x}) = \\max(\\lvert x_1 + 2x_2 - 7 \\rvert, \\lvert 2x_1 + x_2 - 5 \\rvert)


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [1, 3]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = ([-10.0, 10.0], [-10.0, 10.0])

        self.global_optimum = [[1.0, 3.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return max(abs(x[0] + 2 * x[1] - 7), abs(2 * x[0] + x[1] - 5))


class Schwefel20(Benchmark):

    """
    Schwefel 20 objective function.

    This class defines the Schwefel 20 global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel20}}(\\mathbf{x}) = \\sum_{i=1}^n \\lvert x_i \\rvert


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum(abs(x))


class Schwefel21(Benchmark):

    """
    Schwefel 21 objective function.

    This class defines the Schwefel 21 global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel21}}(\\mathbf{x}) = \\smash{\\displaystyle\\max_{1 \leq i \leq n}} \\lvert x_i \\rvert


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return max(abs(x))


class Schwefel22(Benchmark):

    """
    Schwefel 22 objective function.

    This class defines the Schwefel 22 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel22}}(\\mathbf{x}) = \\sum_{i=1}^n \\lvert x_i \\rvert + \\prod_{i=1}^n \\lvert x_i \\rvert


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = ([-10.0, 10.0], [-10.0, 10.0])

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum(abs(x)) + prod(abs(x))


class Schwefel26(Benchmark):

    """
    Schwefel 26 objective function.

    This class defines the Schwefel 26 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel26}}(\\mathbf{x}) = 418.9829n - \\sum_{i=1}^n x_i \\sin(\\sqrt{|x_i|})

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-500, 500]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 420.968746` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([-500.0] * self.N,
                           [500.0] * self.N)

        self.global_optimum = [[420.968746 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return 418.982887 * self.N - sum(x * sin(sqrt(abs(x))))


class Schwefel36(Benchmark):

    """
    Schwefel 36 objective function.

    This class defines the Schwefel 36 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Schwefel36}}(\\mathbf{x}) = -x_1x_2(72 - 2x_1 - 2x_2)


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 500]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -3456` for :math:`\\mathbf{x} = [12, 12]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([0.0] * self.N, [500.0] * self.N)
        self.custom_bounds = ([0.0, 20.0], [0.0, 20.0])

        self.global_optimum = [[12.0, 12.0]]
        self.fglob = -3456.0

    def fun(self, x, *args):
        self.nfev += 1

        return -x[0] * x[1] * (72 - 2 * x[0] - 2 * x[1])


class Shekel05(Benchmark):

    """
    Shekel 5 objective function.

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

        self._bounds = zip([0.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[4.0 for _ in range(self.N)]]
        self.fglob = -10.15319585

    def fun(self, x, *args):
        self.nfev += 1

        m = 5
        A = asarray([[4.0, 4.0, 4.0, 4.0],
                     [1.0, 1.0, 1.0, 1.0],
                     [8.0, 8.0, 8.0, 8.0],
                     [6.0, 6.0, 6.0, 6.0],
                     [3.0, 7.0, 3.0, 7.0]])

        C = asarray([0.1, 0.2, 0.2, 0.4, 0.4])

        return -sum(1.0 / (dot(x - a, x - a) + c) for a, c in zip(A, C))


class Shekel07(Benchmark):

    """
    Shekel 7 objective function.

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

        self._bounds = zip([0.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[4.0 for _ in range(self.N)]]
        self.fglob = -10.4028188

    def fun(self, x, *args):
        self.nfev += 1

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


class Shekel10(Benchmark):

    """
    Shekel 10 objective function.

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

        self._bounds = zip([0.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[4.0 for _ in range(self.N)]]
        self.fglob = -10.5362837262

    def fun(self, x, *args):
        self.nfev += 1

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


class Shubert01(Benchmark):

    """
    Shubert 1 objective function.

    This class defines the Shubert 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Shubert01}}(\\mathbf{x}) = \\left( \\sum\\limits_{i=1}^{5} i\\cos[(i+1)x_1 + i] \\right) \\left( \\sum\\limits_{i=1}^{5} i\\cos[(i+1)x_2 + i] \\right)

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = -186.7309` for :math:`\\mathbf{x} = [-7.0835, 4.8580]` (and many others).

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.global_optimum = [[-7.0835, 4.8580]]

        self.fglob = -186.7309

        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        j = np.atleast_2d(arange(1, 6)).T
        y = j * cos((j + 1) * x + j)
        return prod(sum(y, axis=0))
        # TODO change equation to reflect higher dimensions are possible


class Shubert03(Benchmark):

    """
    Shubert 3 objective function.

    This class defines the Shubert 3 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Shubert03}}(\\mathbf{x}) = \\sum_{i=1}^n \\sum_{j=1}^5 j \\sin \\left[(j+1)x_i \\right] + j

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = -24.062499` for :math:`\\mathbf{x} = [5.791794, 5.791794]` (and many others).

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[5.791794, 5.791794]]
        self.fglob = -24.062499

        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        j = np.atleast_2d(arange(1, 6)).T
        y = -j * sin((j + 1) * x + j)
        # TODO change equation to reflect higher dimensions are possible
        return sum(sum(y))


class Shubert04(Benchmark):

    """
    Shubert 4 objective function.

    This class defines the Shubert 4 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Shubert04}}(\\mathbf{x}) = \\sum_{i=1}^n \\sum_{j=1}^5 j \\cos \\left[(j+1)x_i \\right] + j

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = -29.016015` for :math:`\\mathbf{x} = [-0.80032121, -7.08350592]` (and many others).

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[-0.80032121, -7.08350592]]
        self.fglob = -29.016015

        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        j = np.atleast_2d(arange(1, 6)).T
        y = -j * cos((j + 1) * x + j)
        # TODO change equation to reflect higher dimensions are possible
        return sum(sum(y))


class SineEnvelope(Benchmark):

    """
    SineEnvelope objective function.

    This class defines the SineEnvelope global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{SineEnvelope}}(\\mathbf{x}) = -\\sum_{i=1}^{n-1}\\left[\\frac{\\sin^2(\\sqrt{x_{i+1}^2+x_{i}^2}-0.5)}{(0.001(x_{i+1}^2+x_{i}^2)+1)^2}+0.5\\right]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
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

    """
    Six Hump Camel objective function.

    This class defines the Six Hump Camel global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{SixHumpCamel}}(\\mathbf{x}) = 4x_1^2+x_1x_2-4x_2^2-2.1x_1^4+4x_2^4+\\frac{1}{3}x_1^6

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -1.031628453489877` for :math:`\\mathbf{x} = [0.08984201368301331 , -0.7126564032704135]`
    or :math:`\\mathbf{x} = [-0.08984201368301331, 0.7126564032704135]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.0] * self.N, [5.0] * self.N)
        self.custom_bounds = [(-2, 2), (-1.5, 1.5)]

        self.global_optimum = [(0.08984201368301331, -0.7126564032704135),
                               (-0.08984201368301331, 0.7126564032704135)]
        self.fglob = -1.031628

    def fun(self, x, *args):
        self.nfev += 1
        return ((4 - 2.1 * x[0] ** 2 + x[0] ** 4 / 3) * x[0] ** 2 + x[0] * x[1]
                + (4 * x[1] ** 2 - 4) * x[1] ** 2)


class Sodp(Benchmark):

    """
    Sodp objective function.

    This class defines the Sum Of Different Powers global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Sodp}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\lvert{x_{i}}\\rvert^{i + 1}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-1.0] * self.N, [1.0] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        i = arange(1, self.N + 1)
        return sum(abs(x) ** (i + 1))


class Sphere(Benchmark):

    """
    Sphere objective function.

    This class defines the Sphere global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Sphere}}(\\mathbf{x}) = \\sum_{i=1}^{n} x_i^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 1]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([-5.12] * self.N, [5.12] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum(x ** 2)


class Step(Benchmark):

    """
    Step objective function.

    This class defines the Step global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Step}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\left ( \\lfloor x_i  + 0.5 \\rfloor \\right )^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0.5` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)
        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)
        self.custom_bounds = ([-5, 5], [-5, 5])

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum((floor(x + 0.5)) ** 2.0)


class Stochastic(Benchmark):

    """
    Stochastic objective function.

    This class defines a Stochastic global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

        f_{\\text{Stochastic}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\epsilon_i \\left | {x_i - \\frac{1}{i}} \\right |

    The variable :math:`\\epsilon_i, (i=1,...,n)` is a random variable uniformly distributed in :math:`[0, 1]`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = [1/n]` for :math:`i=1,...,n`
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.0] * self.N, [5.0] * self.N)

        self.global_optimum = [[1.0 / _ for _ in range(1, self.N + 1)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        rnd = np.random.uniform(0.0, 1.0, size=(self.N, ))
        i = arange(1, self.N + 1)

        return sum(rnd * abs(x - 1.0 / i))


class StretchedV(Benchmark):

    """
    StretchedV objective function.

    This class defines the Stretched V global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{StretchedV}}(\\mathbf{x}) = \sum_{i=1}^{n-1} t^{1/4} [\sin (50t^{0.1}) + 1]^2

    Where, in this exercise:

    .. math::

       t = x_{i+1}^2 + x_i^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [0., 0.]` when :math:`n = 2`.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10] * self.N, [10] * self.N)

        self.global_optimum = [[0, 0]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        t = x[1:] ** 2 + x[: -1] ** 2
#         TODO: fix equation in docs
        return sum(t ** 0.25 * (sin(50.0 * t ** 0.1) + 1) ** 2)


class StyblinskiTang(Benchmark):

    """
    StyblinskiTang objective function.

    This class defines the Styblinski-Tang global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{StyblinskiTang}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\left(x_i^4 - 16x_i^2 + 5x_i \\right)

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = -39.16616570377142n` for :math:`x_i = -2.903534018185960` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.0] * self.N, [5.0] * self.N)

        self.global_optimum = [[-2.903534018185960 for _ in range(self.N)]]
        self.fglob = -39.16616570377142 * self.N
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum(x ** 4 - 16 * x ** 2 + 5 * x) / 2


class TestTubeHolder(Benchmark):

    """
    TestTubeHolder objective function.

    This class defines the TestTubeHolder global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{TestTubeHolder}}(\\mathbf{x}) = - 4 \\left | {e^{\\left|{\\cos\\left(\\frac{1}{200} x_{1}^{2} + \\frac{1}{200} x_{2}^{2}\\right)}\\right|} \\sin\\left(x_{1}\\right) \\cos\\left(x_{2}\\right)}\\right |

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -10.872299901558` for :math:`\\mathbf{x} = [-\\pi/2, 0]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[-pi / 2, 0.0]]
        self.fglob = -10.87229990155800

    def fun(self, x, *args):
        self.nfev += 1

        u = sin(x[0]) * cos(x[1])
        v = (x[0] ** 2 + x[1] ** 2) / 200
        return -4 * abs(u * exp(abs(cos(v))))


class Thurber(Benchmark):

    """
    Thurber objective function.

    """

    # TODO, this is a NIST regression standard dataset
    def __init__(self, dimensions=7):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip(
            [500., 500., 100., 10., 0.1, 0.1, 0.],
            [2000., 2000., 1000., 150., 2., 1., 0.2])
        self.global_optimum = [[1.288139680e3, 1.4910792535e3, 5.8323836877e2,
                                75.416644291, 0.96629502864, 0.39797285797,
                                4.9727297349e-2]]
        self.fglob = 5642.7082397

    def fun(self, x, *args):
        self.nfev += 1

        a = asarray([80.574, 84.248, 87.264, 87.195, 89.076, 89.608, 89.868,
                     90.101, 92.405, 95.854, 100.696, 101.06, 401.672, 390.724,
                     567.534, 635.316, 733.054, 759.087, 894.206, 990.785,
                     1090.109, 1080.914, 1122.643, 1178.351, 1260.531, 1273.514,
                     1288.339, 1327.543, 1353.863, 1414.509, 1425.208, 1421.384,
                     1442.962, 1464.350, 1468.705, 1447.894, 1457.628])
        b = asarray([-3.067, -2.981, -2.921, -2.912, -2.840, -2.797, -2.702,
                     -2.699, -2.633, -2.481, -2.363, -2.322, -1.501, -1.460,
                     -1.274, -1.212, -1.100, -1.046, -0.915, -0.714, -0.566,
                     -0.545, -0.400, -0.309, -0.109, -0.103, 0.010, 0.119,
                     0.377, 0.790, 0.963, 1.006, 1.115, 1.572, 1.841, 2.047,
                     2.200])

        vec = x[0] + x[1] * b + x[2] * b ** 2 + x[3] * b ** 3
        vec /= 1 + x[4] * b + x[5] * b ** 2 + x[6] * b ** 3

        return sum((a - vec) ** 2)


class Treccani(Benchmark):

    """
    Treccani objective function.

    This class defines the Treccani global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Treccani}}(\\mathbf{x}) = x_1^4 + 4x_1^3 + 4x_1^2 + x_2^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [-2, 0]` or :math:`\\mathbf{x} = [0, 0]`.

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.0] * self.N, [5.0] * self.N)
        self.custom_bounds = [(-2, 2), (-2, 2)]

        self.global_optimum = [[-2.0, 0.0]]
        self.fglob = 0

    def fun(self, x, *args):
        self.nfev += 1

        return x[0] ** 4 + 4.0 * x[0] ** 3 + 4.0 * x[0] ** 2 + x[1] ** 2


class Trefethen(Benchmark):

    """
    Trefethen objective function.

    This class defines the Trefethen global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Trefethen}}(\\mathbf{x}) = 0.25 x_{1}^{2} + 0.25 x_{2}^{2} + e^{\\sin\\left(50 x_{1}\\right)} - \\sin\\left(10 x_{1} + 10 x_{2}\\right) + \\sin\\left(60 e^{x_{2}}\\right) + \\sin\\left[70 \\sin\\left(x_{1}\\right)\\right] + \\sin\\left[\\sin\\left(80 x_{2}\\right)\\right]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -3.3068686474` for :math:`\\mathbf{x} = [-0.02440307923, 0.2106124261]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.custom_bounds = [(-5, 5), (-5, 5)]

        self.global_optimum = [[-0.02440307923, 0.2106124261]]
        self.fglob = -3.3068686474

    def fun(self, x, *args):
        self.nfev += 1

        val = 0.25 * x[0] ** 2 + 0.25 * x[1] ** 2
        val += exp(sin(50. * x[0])) - sin(10 * x[0] + 10 * x[1])
        val += sin(60 * exp(x[1]))
        val += sin(70 * sin(x[0]))
        val += sin(sin(80 * x[1]))
        return val


class ThreeHumpCamel(Benchmark):

    """
    Three Hump Camel objective function.

    This class defines the Three Hump Camel global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{ThreeHumpCamel}}(\\mathbf{x}) = 2x_1^2 - 1.05x_1^4 + \\frac{x_1^6}{6} + x_1x_2 + x_2^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [0, 0]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.0] * self.N, [5.0] * self.N)
        self.custom_bounds = [(-2, 2), (-1.5, 1.5)]

        self.global_optimum = [[0.0, 0.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return (2.0 * x[0] ** 2.0 - 1.05 * x[0] ** 4.0 + x[0] ** 6 / 6.0
                + x[0] * x[1] + x[1] ** 2.0)


class Trid(Benchmark):

    """
    Trid objective function.

    This class defines the Trid global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Trid}}(\\mathbf{x}) = \\sum_{i=1}^{n}(x_i - 1)^2 - \\sum_{i=2}^{n} x_ix_{i-1}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-20, 20]` for :math:`i=1,...,6`.

    *Global optimum*: :math:`f(x_i) = -50` for :math:`\\mathbf{x} = [6, 10, 12, 12, 10, 6]`

    """

    def __init__(self, dimensions=6):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N, [20.0] * self.N)

        self.global_optimum = [[6, 10, 12, 12, 10, 6]]
        self.fglob = -50.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum((x - 1.0) ** 2.0) - sum(x[1:] * x[:-1])


class Trigonometric01(Benchmark):

    """
    Trigonometric 1 objective function.

    This class defines the Trigonometric 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Trigonometric01}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\left [n - \\sum_{j=1}^{n} \\cos(x_j) + i \\left(1 - cos(x_i) - sin(x_i) \\right ) \\right]^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, \\pi]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N, [pi] * self.N)

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        i = arange(1.0, self.N + 1)
        return sum((self.N - sum(cos(x) + i * (1 - cos(x) - sin(x)))) ** 2.0)


class Trigonometric02(Benchmark):

    """
    Trigonometric 2 objective function.

    This class defines the Trigonometric 2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Trigonometric2}}(\\mathbf{x}) = 1 + \\sum_{i=1}^{n} 8 \\sin^2 \\left[7(x_i - 0.9)^2 \\right] + 6 \\sin^2 \\left[14(x_i - 0.9)^2 \\right] + (x_i - 0.9)^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-500, 500]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 1` for :math:`x_i = 0.9` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-500.0] * self.N,
                           [500.0] * self.N)
        self.custom_bounds = [(0, 2), (0, 2)]

        self.global_optimum = [[0.9 for _ in range(self.N)]]
        self.fglob = 1.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        vec = (8 * sin(7 * (x - 0.9) ** 2) ** 2
               + 6 * sin(14 * (x - 0.9) ** 2) ** 2
               + (x - 0.9) ** 2)
        return 1.0 + sum(vec)


class Tripod(Benchmark):

    """
    Tripod objective function.

    This class defines the Tripod global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Tripod}}(\\mathbf{x}) = p(x_2) \\left[1 + p(x_1) \\right] + \\lvert x_1 + 50p(x_2) \\left[1 - 2p(x_1) \\right] \\rvert + \\lvert x_2 + 50\\left[1 - 2p(x_2)\\right] \\rvert

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-100, 100]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [0, -50]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-100.0] * self.N,
                           [100.0] * self.N)

        self.global_optimum = [[0.0, -50.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        p1 = float(x[0] >= 0)
        p2 = float(x[1] >= 0)

        return (p2 * (1.0 + p1) + abs(x[0] + 50.0 * p2 * (1.0 - 2.0 * p1))
                + abs(x[1] + 50.0 * (1.0 - 2.0 * p2)))


class Ursem01(Benchmark):

    """
    Ursem 1 objective function.

    This class defines the Ursem 1 global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Ursem01}}(\\mathbf{x}) = - \\sin(2x_1 - 0.5 \\pi) - 3 \\cos(x_2) - 0.5x_1

    Here, :math:`n` represents the number of dimensions and :math:`x_1 \\in [-2.5, 3]`, :math:`x_2 \\in [-2, 2]`.

    *Global optimum*: :math:`f(x_i) = -4.81681406371` for :math:`\\mathbf{x} = [1.69714, 0.0]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = [(-2.5, 3.0), (-2.0, 2.0)]

        self.global_optimum = [[1.69714, 0.0]]
        self.fglob = -4.81681406371

    def fun(self, x, *args):
        self.nfev += 1

        return (-sin(2 * x[0] - 0.5 * pi) - 3.0 * cos(x[1]) - 0.5 * x[0])


class Ursem03(Benchmark):

    """
    Ursem 3 objective function.

    This class defines the Ursem 3 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Ursem03}}(\\mathbf{x}) = - \\sin(2.2 \\pi x_1 + 0.5 \\pi) \\frac{2 - \\lvert x_1 \\rvert}{2} \\frac{3 - \\lvert x_1 \\rvert}{2} - \\sin(2.2 \\pi x_2 + 0.5 \\pi) \\frac{2 - \\lvert x_2 \\rvert}{2} \\frac{3 - \\lvert x_2 \\rvert}{2}

    Here, :math:`n` represents the number of dimensions and :math:`x_1 \\in [-2, 2]`, :math:`x_2 \\in [-1.5, 1.5]`.

    *Global optimum*: :math:`f(x_i) = -3` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = [(-2, 2), (-1.5, 1.5)]

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = -3.0

    def fun(self, x, *args):
        self.nfev += 1

        u = -(sin(2.2 * pi * x[0] + 0.5 * pi)
              * ((2.0 - abs(x[0])) / 2.0) * ((3.0 - abs(x[0])) / 2))
        v = -(sin(2.2 * pi * x[1] + 0.5 * pi)
              * ((2.0 - abs(x[1])) / 2) * ((3.0 - abs(x[1])) / 2))
        return u + v


class Ursem04(Benchmark):

    """
    Ursem 4 objective function.

    This class defines the Ursem 4 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Ursem04}}(\\mathbf{x}) = -3 \\sin(0.5 \\pi x_1 + 0.5 \\pi) \\frac{2 - \\sqrt{x_1^2 + x_2 ^ 2}}{4}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-2, 2]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -1.5` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-2.0] * self.N, [2.0] * self.N)

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = -1.5

    def fun(self, x, *args):
        self.nfev += 1

        return (-3 * sin(0.5 * pi * x[0] + 0.5 * pi)
                * (2 - sqrt(x[0] ** 2 + x[1] ** 2)) / 4)


class UrsemWaves(Benchmark):

    """
    Ursem Waves objective function.

    This class defines the Ursem Waves global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{UrsemWaves}}(\\mathbf{x}) = -0.9x_1^2 + (x_2^2 - 4.5x_2^2)x_1x_2 + 4.7 \\cos \\left[ 2x_1 - x_2^2(2 + x_1) \\right ] \\sin(2.5 \\pi x_1)

    Here, :math:`n` represents the number of dimensions and :math:`x_1 \\in [-0.9, 1.2]`, :math:`x_2 \\in [-1.2, 1.2]`.

    *Global optimum*: :math:`f(x_i) = -8.5536` for :math:`x_i = 1.2` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = [(-0.9, 1.2), (-1.2, 1.2)]

        self.global_optimum = [[1.2 for _ in range(self.N)]]
        self.fglob = -8.5536

    def fun(self, x, *args):
        self.nfev += 1

        u = -0.9 * x[0] ** 2
        v = (x[1] ** 2 - 4.5 * x[1] ** 2) * x[0] * x[1]
        w = 4.7 * cos(2 * x[0] - x[1] ** 2 * (2 + x[0])) * sin(2.5 * pi * x[0])
        return (u + v + w)


class VenterSobiezcczanskiSobieski(Benchmark):

    """
    Venter Sobiezcczanski-Sobieski objective function.

    This class defines the Venter Sobiezcczanski-Sobieski global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{VenterSobiezcczanskiSobieski}}(\\mathbf{x}) = x_1^2 - 100 \\cos^2(x_1) - 100 \\cos(x_1^2/30) + x_2^2 - 100 \\cos^2(x_2) - 100 \\cos(x_2^2/30)

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-50, 50]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -400` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-50.0] * self.N, [50.0] * self.N)
        self.custom_bounds = ([-10, 10], [-10, 10])

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = -400

    def fun(self, x, *args):
        self.nfev += 1

        u = x[0] ** 2.0 - 100.0 * cos(x[0]) ** 2.0
        v = -100.0 * cos(x[0] ** 2.0 / 30.0) + x[1] ** 2.0
        w = - 100.0 * cos(x[1]) ** 2.0 - 100.0 * cos(x[1] ** 2.0 / 30.0)
        return u + v + w


class Vincent(Benchmark):

    """
    Vincent objective function.

    This class defines the Vincent global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Vincent}}(\\mathbf{x}) = - \\sum_{i=1}^{n} \\sin(10 \\log(x))

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0.25, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = -n` for :math:`x_i = 7.70628098` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.25] * self.N, [10.0] * self.N)

        self.global_optimum = [[7.70628098 for _ in range(self.N)]]
        self.fglob = -float(self.N)
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return -sum(sin(10.0 * log(x)))


class Watson(Benchmark):

    """
    Watson objective function.

    This class defines the Watson global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Watson}}(\\mathbf{x}) = \\sum_{i=0}^{29} \\left\\{ \\sum_{j=0}^4 ((j + 1)a_i^j x_{j+1}) - \\left[ \\sum_{j=0}^5 a_i^j x_{j+1} \\right ]^2 - 1 \\right\\}^2 + x_1^2


    Where, in this exercise, :math:`a_i = i/29`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,...,6`.

    *Global optimum*: :math:`f(x_i) = 0.002288` for :math:`\\mathbf{x} = [-0.0158, 1.012, -0.2329, 1.260, -1.513, 0.9928]`

    """

    def __init__(self, dimensions=6):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.0] * self.N, [5.0] * self.N)

        self.global_optimum = [
            [-0.0158, 1.012, -0.2329, 1.260, -1.513, 0.9928]]
        self.fglob = 0.002288

    def fun(self, x, *args):
        self.nfev += 1

        i = np.atleast_2d(arange(30.)).T
        a = i / 29.
        j = arange(5.)
        k = arange(6.)

        t1 = sum((j + 1) * a ** j * x[1:], axis=1)
        t2 = sum(a ** k * x, axis=1)

        inner = (t1 - t2 ** 2 - 1) ** 2

        return sum(inner) + x[0] ** 2


class Wavy(Benchmark):

    """
    W / Wavy objective function.

    This class defines the W / Wavy global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Wavy}}(\\mathbf{x}) = 1 - \\frac{1}{n} \\sum_{i=1}^{n} \\cos(kx_i)e^{-\\frac{x_i^2}{2}}


    Where, in this exercise, :math:`k = 10`. The number of local minima is :math:`kn` and :math:`(k + 1)n` for odd and even :math:`k` respectively.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-\\pi, \\pi]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,2`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-pi] * self.N, [pi] * self.N)

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return 1.0 - (1.0 / self.N) * sum(cos(10 * x) * exp(-x ** 2.0 / 2.0))


class WayburnSeader01(Benchmark):

    """
    Wayburn and Seader 1 objective function.

    This class defines the Wayburn and Seader 1 global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{WayburnSeader01}}(\\mathbf{x}) = (x_1^6 + x_2^4 - 17)^2 + (2x_1 + x_2 - 4)^2


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [1, 2]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.0] * self.N, [5.0] * self.N)
        self.custom_bounds = ([-2, 2], [-2, 2])

        self.global_optimum = [[1.0, 2.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return (x[0] ** 6 + x[1] ** 4 - 17) ** 2 + (2 * x[0] + x[1] - 4) ** 2


class WayburnSeader02(Benchmark):

    """
    Wayburn and Seader 2 objective function.

    This class defines the Wayburn and Seader 2 global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{WayburnSeader02}}(\\mathbf{x}) = \\left[ 1.613 - 4(x_1 - 0.3125)^2 - 4(x_2 - 1.625)^2 \\right]^2 + (x_2 - 1)^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-500, 500]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [0.2, 1]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-500.0] * self.N,
                           [500.0] * self.N)
        self.custom_bounds = ([-1, 2], [-1, 2])

        self.global_optimum = [[0.2, 1.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        u = (1.613 - 4 * (x[0] - 0.3125) ** 2 - 4 * (x[1] - 1.625) ** 2) ** 2
        v = (x[1] - 1) ** 2
        return u + v


class Weierstrass(Benchmark):

    """
    Weierstrass objective function.

    This class defines the Weierstrass global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Weierstrass}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\left [ \\sum_{k=0}^{kmax} a^k \\cos \\left( 2 \\pi b^k (x_i + 0.5) \\right) - n \\sum_{k=0}^{kmax} a^k \\cos(\\pi b^k) \\right ]


    Where, in this exercise, :math:`kmax = 20`, :math:`a = 0.5` and :math:`b = 3`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-0.5, 0.5]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 4` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-0.5] * self.N, [0.5] * self.N)

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 4.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        kmax = 20
        a, b = 0.5, 3.0

        i = arange(self.N)
        k = np.atleast_2d(arange(kmax + 1.)).T

        t1 = a ** k * cos(2 * pi * b ** k * (x + 0.5))
        t2 = self.N * sum(a ** k.T * cos(pi * b ** k.T))

        inner = sum(t1, axis=0) - t2
        return sum(inner)


class Whitley(Benchmark):

    """
    Whitley objective function.

    This class defines the Whitley global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Whitley}}(\\mathbf{x}) = \\sum_{i=1}^n \\sum_{j=1}^n \\left[\\frac{(100(x_i^2-x_j)^2 + (1-x_j)^2)^2}{4000} - \\cos(100(x_i^2-x_j)^2 + (1-x_j)^2)+1 \\right]


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10.24, 10.24]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 1` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.24] * self.N,
                           [10.24] * self.N)
        self.custom_bounds = ([-1, 2], [-1, 2])

        self.global_optimum = [[1.0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        XI = x
        XJ = np.atleast_2d(x).T

        temp = 100.0 * ((XI ** 2.0) - XJ) + (1.0 - XJ) ** 2.0
        inner = (temp ** 2.0 / 4000.0) - cos(temp) + 1.0
        return sum(sum(inner, axis=0))


class Wolfe(Benchmark):

    """
    Wolfe objective function.

    This class defines the Wolfe global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Wolfe}}(\\mathbf{x}) = \\frac{4}{3}(x_1^2 + x_2^2 - x_1x_2)^{0.75} + x_3


    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [0, 2]` for :math:`i=1,2,3`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,2,3`

    """

    def __init__(self, dimensions=3):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N, [2.0] * self.N)

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        return 4 / 3 * (x[0] ** 2 + x[1] ** 2 - x[0] * x[1]) ** 0.75 + x[2]


class XinSheYang01(Benchmark):

    """
    Xin-She Yang 1 objective function.

    This class defines the Xin-She Yang 1 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{XinSheYang01}}(\\mathbf{x}) = \\sum_{i=1}^{n} \\epsilon_i \\lvert x_i \\rvert^i


    The variable :math:`\\epsilon_i, (i=1,...,n)` is a random variable uniformly distributed in :math:`[0, 1]`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 5]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.0] * self.N, [5.0] * self.N)
        self.custom_bounds = ([-2, 2], [-2, 2])

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        i = arange(1.0, self.N + 1.0)
        return sum(np.random.random() * (abs(x) ** i))


class XinSheYang02(Benchmark):

    """
    Xin-She Yang 2 objective function.

    This class defines the Xin-She Yang 2 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{XinSheYang02}}(\\mathbf{x}) = \\frac{\\sum_{i=1}^{n} \\lvert{x_{i}}\\rvert}{e^{\\sum_{i=1}^{n} \\sin\\left(x_{i}^{2.0}\\right)}}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-2\\pi, 2\\pi]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-2 * pi] * self.N,
                           [2 * pi] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum(abs(x)) * exp(-sum(sin(x ** 2.0)))


class XinSheYang03(Benchmark):

    """
    Xin-She Yang 3 objective function.

    This class defines the Xin-She Yang 3 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{XinSheYang03}}(\\mathbf{x}) = e^{-\\sum_{i=1}^{n} (x_i/\\beta)^{2m}} - 2e^{-\\sum_{i=1}^{n} x_i^2} \\prod_{i=1}^{n} \\cos^2(x_i)


    Where, in this exercise, :math:`\\beta = 15` and :math:`m = 3`.

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-20, 20]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = -1` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-20.0] * self.N, [20.0] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = -1.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        beta, m = 15.0, 5.0
        u = sum((x / beta) ** (2 * m))
        v = sum(x ** 2)
        w = prod(cos(x) ** 2)

        return exp(-u) - 2 * exp(-v) * w


class XinSheYang04(Benchmark):

    """
    Xin-She Yang 4 objective function.

    This class defines the Xin-She Yang 4 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{XinSheYang04}}(\\mathbf{x}) = \\left[ \\sum_{i=1}^{n} \\sin^2(x_i) - e^{-\\sum_{i=1}^{n} x_i^2} \\right ] e^{-\\sum_{i=1}^{n} \\sin^2 \\sqrt{ \\lvert x_i \\rvert }}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = -1` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = -1.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        u = sum(sin(x) ** 2)
        v = sum(x ** 2)
        w = sum(sin(sqrt(abs(x))) ** 2)
        return (u - exp(-v)) * exp(-w)


class Xor(Benchmark):

    """
    Xor objective function.

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

        self._bounds = zip([-1.0] * self.N, [1.0] * self.N)

        self.global_optimum = [[1.0, -1.0, 1.0,
                               -1.0, -1.0, 1.0, 1.0, -1.0, 0.421134]]
        self.fglob = 0.9597588

    def fun(self, x, *args):
        self.nfev += 1

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


class YaoLiu04(Benchmark):

    """
    Yao-Liu 4 objective function.

    This class defines the Yao-Liu function 4 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{YaoLiu04}}(\\mathbf{x}) = {max}_i \\left\{ \\left | x_i \\right | , 1 \\leq i \\leq n \\right\}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return abs(x.max())


class YaoLiu09(Benchmark):

    """
    Yao-Liu 9 objective function.

    This class defines the Yao-Liu function 9 global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{YaoLiu09}}(\\mathbf{x}) = \\sum_{i=1}^n \\left [ x_i^2 - 10 \\cos(2 \\pi x_i ) + 10 \\right ]

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5.12, 5.12]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.12] * self.N, [5.12] * self.N)

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        return sum(x ** 2.0 - 10.0 * cos(2 * pi * x) + 10)


class Zacharov(Benchmark):

    """
    Zacharov objective function.

    This class defines the Zacharov global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

         f_{\\text{Zacharov}}(\\mathbf{x}) = \\sum_{i=1}^{n} x_i^2 + \\left ( \\frac{1}{2} \\sum_{i=1}^{n} i x_i \\right )^2 + \\left ( \\frac{1}{2} \\sum_{i=1}^{n} i x_i \\right )^4

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-5, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` for :math:`x_i = 0` for :math:`i=1,...,n`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-5.0] * self.N, [10.0] * self.N)
        self.custom_bounds = ([-1, 1], [-1, 1])

        self.global_optimum = [[0 for _ in range(self.N)]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        u = sum(x ** 2)
        v = sum(arange(1, self.N + 1) * x)
        return u + (0.5 * v) ** 2 + (0.5 * v) ** 4


class ZeroSum(Benchmark):

    """
    ZeroSum objective function.

    This class defines the ZeroSum global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

         f_{\\text{ZeroSum}}(\\mathbf{x}) = \\begin{cases}0 & \\textrm{if} \\sum_{i=1}^n x_i = 0 \\\\
                1 + \\left(10000 \\left |\\sum_{i=1}^n x_i\\right| \\right)^{0.5} & \\textrm{otherwise}\\end{cases}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,...,n`.

    *Global optimum*: :math:`f(x_i) = 0` where :math:`\\sum_{i=1}^n x_i = 0`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)

        self.global_optimum = [[]]
        self.fglob = 0.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        if abs(sum(x)) < 3e-16:
            return 0.0
#         TODO: the abs term doesn't appear to be in the equation.
        return 1.0 + (10000.0 * abs(sum(x))) ** 0.5


class Zettl(Benchmark):

    """
    Zettl objective function.

    This class defines the Zettl global optimization problem. This
    is a multimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Zettl}}(\\mathbf{x}) = \\frac{1}{4} x_{1} + \\left(x_{1}^{2} - 2 x_{1} + x_{2}^{2}\\right)^{2}

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-1, 5]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -0.0037912` for :math:`\\mathbf{x} = [-0.029896, 0.0]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-1.0] * self.N, [5.0] * self.N)

        self.global_optimum = [[-0.02989597760285287, 0.0]]
        self.fglob = -0.003791237220468656

    def fun(self, x, *args):
        self.nfev += 1

        return (x[0] ** 2 + x[1] ** 2 - 2 * x[0]) ** 2 + x[0] / 4


class Zimmerman(Benchmark):

    """
    Zimmerman objective function.

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

    *Global optimum*: :math:`f(x_i) = 0` for :math:`\\mathbf{x} = [7, 2]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([0.0] * self.N, [100.0] * self.N)
        self.custom_bounds = ([0.0, 8.0], [0.0, 8.0])

        self.global_optimum = [[7.0, 2.0]]
        self.fglob = 0.0

    def fun(self, x, *args):
        self.nfev += 1

        Zh1 = lambda x: 9.0 - x[0] - x[1]
        Zh2 = lambda x: (x[0] - 3.0) ** 2.0 + (x[1] - 2.0) ** 2.0 - 16.0
        Zh3 = lambda x: x[0] * x[1] - 14.0
        Zp = lambda x: 100.0 * (1.0 + x)

        return max(Zh1(x),
                   Zp(Zh2(x)) * sign(Zh2(x)),
                   Zp(Zh3(x)) * sign(Zh3(x)),
                   Zp(-x[0]) * sign(x[0]),
                   Zp(-x[1]) * sign(x[1]))


class Zirilli(Benchmark):

    """
    Zettl objective function.

    This class defines the Zirilli global optimization problem. This
    is a unimodal minimization problem defined as follows:

    .. math::

       f_{\\text{Zirilli}}(\\mathbf{x}) = 0.25x_1^4 - 0.5x_1^2 + 0.1x_1 + 0.5x_2^2

    Here, :math:`n` represents the number of dimensions and :math:`x_i \\in [-10, 10]` for :math:`i=1,2`.

    *Global optimum*: :math:`f(x_i) = -0.3523` for :math:`\\mathbf{x} = [-1.0465, 0]`

    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = zip([-10.0] * self.N, [10.0] * self.N)
        self.custom_bounds = ([-2.0, 2.0], [-2.0, 2.0])

        self.global_optimum = [[-1.0465, 0.0]]
        self.fglob = -0.35238603

    def fun(self, x, *args):
        self.nfev += 1

        return 0.25 * x[0] ** 4 - 0.5 * x[0] ** 2 + 0.1 * x[0] + 0.5 * x[1] ** 2


#-----------------------------------------------------------------------
#                 UNIVARIATE SINGLE-OBJECTIVE PROBLEMS
#-----------------------------------------------------------------------

class Problem02(Benchmark):

    """
    Univariate Problem02 objective function.

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

        self._bounds = [(2.7, 7.5)]

        self.global_optimum = 5.145735
        self.fglob = -1.899599

    def fun(self, x, *args):
        self.nfev += 1

        x = x[0]
        return sin(x) + sin(10.0 / 3.0 * x)


class Problem03(Benchmark):

    """
    Univariate Problem03 objective function.

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

        self._bounds = [(-10, 10)]

        self.global_optimum = -6.7745761
        self.fglob = -12.03124

    def fun(self, x, *args):
        self.nfev += 1

        x = x[0]
        y = 0.0
        for k in range(1, 6):
            y += k * sin((k + 1) * x + k)

        return -y


class Problem04(Benchmark):

    """
    Univariate Problem04 objective function.

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

        self._bounds = [(1.9, 3.9)]

        self.global_optimum = 2.868034
        self.fglob = -3.85045

    def fun(self, x, *args):
        self.nfev += 1

        x = x[0]
        return -(16 * x ** 2 - 24 * x + 5) * exp(-x)


class Problem05(Benchmark):

    """
    Univariate Problem05 objective function.

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

        self._bounds = [(0.0, 1.2)]

        self.global_optimum = 0.96609
        self.fglob = -1.48907

    def fun(self, x, *args):
        self.nfev += 1

        x = x[0]
        return -(1.4 - 3 * x) * sin(18.0 * x)


class Problem06(Benchmark):

    """
    Univariate Problem06 objective function.

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

        self._bounds = [(-10.0, 10.0)]

        self.global_optimum = 0.67956
        self.fglob = -0.824239

    def fun(self, x, *args):
        self.nfev += 1

        x = x[0]
        return -(x + sin(x)) * exp(-x ** 2.0)


class Problem07(Benchmark):

    """
    Univariate Problem07 objective function.

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

        self._bounds = [(2.7, 7.5)]

        self.global_optimum = 5.19978
        self.fglob = -1.6013

    def fun(self, x, *args):
        self.nfev += 1

        x = x[0]
        return sin(x) + sin(10.0 / 3.0 * x) + log(x) - 0.84 * x + 3


class Problem08(Benchmark):

    """
    Univariate Problem08 objective function.

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

        self._bounds = [(-10, 10)]

        self.global_optimum = -7.083506
        self.fglob = -14.508

    def fun(self, x, *args):
        self.nfev += 1

        x = x[0]

        y = 0.0
        for k in range(1, 6):
            y += k * cos((k + 1) * x + k)

        return -y


class Problem09(Benchmark):

    """
    Univariate Problem09 objective function.

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

        self._bounds = [(3.1, 20.4)]

        self.global_optimum = 17.039
        self.fglob = -1.90596

    def fun(self, x, *args):
        self.nfev += 1

        x = x[0]
        return sin(x) + sin(2.0 / 3.0 * x)


class Problem10(Benchmark):

    """
    Univariate Problem10 objective function.

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

        self._bounds = [(0, 10)]

        self.global_optimum = 7.9787
        self.fglob = -7.916727

    def fun(self, x, *args):
        self.nfev += 1

        x = x[0]
        return -x * sin(x)


class Problem11(Benchmark):

    """
    Univariate Problem11 objective function.

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

        self._bounds = [(-pi / 2, 2 * pi)]

        self.global_optimum = 2.09439
        self.fglob = -1.5

    def fun(self, x, *args):
        self.nfev += 1

        x = x[0]
        return 2 * cos(x) + cos(2 * x)


class Problem12(Benchmark):

    """
    Univariate Problem12 objective function.

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

        self._bounds = [(0, 2 * pi)]

        self.global_optimum = pi
        self.fglob = -1

    def fun(self, x, *args):
        self.nfev += 1

        x = x[0]
        return (sin(x)) ** 3.0 + (cos(x)) ** 3.0


class Problem13(Benchmark):

    """
    Univariate Problem13 objective function.

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

        self._bounds = [(0.001, 0.99)]

        self.global_optimum = 1.0 / sqrt(2)
        self.fglob = -1.5874

    def fun(self, x, *args):
        self.nfev += 1

        x = x[0]
        return -x ** (2.0 / 3.0) - (1.0 - x ** 2) ** (1.0 / 3.0)


class Problem14(Benchmark):

    """
    Univariate Problem14 objective function.

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

        self._bounds = [(0.0, 4.0)]

        self.global_optimum = 0.224885
        self.fglob = -0.788685

    def fun(self, x, *args):
        self.nfev += 1

        x = x[0]
        return -exp(-x) * sin(2.0 * pi * x)


class Problem15(Benchmark):

    """
    Univariate Problem15 objective function.

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

        self._bounds = [(-5.0, 5.0)]

        self.global_optimum = 2.41422
        self.fglob = -0.03553

    def fun(self, x, *args):
        self.nfev += 1

        x = x[0]
        return -(-x ** 2.0 + 5 * x - 6) / (x ** 2 + 1)


class Problem18(Benchmark):

    """
    Univariate Problem18 objective function.

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

        self._bounds = [(0.0, 6.0)]

        self.global_optimum = 2
        self.fglob = 0

    def fun(self, x, *args):
        self.nfev += 1

        x = x[0]

        if x <= 3:
            return (x - 2.0) ** 2.0

        return 2 * log(x - 2.0) + 1


class Problem20(Benchmark):

    """
    Univariate Problem20 objective function.

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

        self._bounds = [(-10, 10)]

        self.global_optimum = 1.195137
        self.fglob = -0.0634905

    def fun(self, x, *args):
        self.nfev += 1

        x = x[0]
        return -(x - sin(x)) * exp(-x ** 2.0)


class Problem21(Benchmark):

    """
    Univariate Problem21 objective function.

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

        self._bounds = [(0, 10)]

        self.global_optimum = 4.79507
        self.fglob = -9.50835

    def fun(self, x, *args):
        self.nfev += 1

        x = x[0]
        return x * sin(x) + x * cos(2.0 * x)


class Problem22(Benchmark):

    """
    Univariate Problem22 objective function.

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

        self._bounds = [(0, 20)]

        self.global_optimum = 9.0 * pi / 2.0
        self.fglob = exp(-27.0 * pi / 2.0) - 1.0

    def fun(self, x, *args):
        self.nfev += 1

        x = x[0]
        return exp(-3.0 * x) - (sin(x)) ** 3.0
