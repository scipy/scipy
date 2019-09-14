# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import

from numpy import cos, sqrt, sin, abs
from .go_benchmark import Benchmark


class NeedleEye(Benchmark):

    r"""
    NeedleEye objective function.

    This class defines the Needle-Eye [1]_ global optimization problem. This is a
    a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{NeedleEye}}(x) =
            \begin{cases}
            1 & \textrm{if }\hspace{5pt} \lvert x_i \rvert  <  eye \hspace{5pt}
            \forall i \\
            \sum_{i=1}^n (100 + \lvert x_i \rvert) & \textrm{if } \hspace{5pt}
            \lvert x_i \rvert > eye \\
            0 & \textrm{otherwise}\\
            \end{cases}


    Where, in this exercise, :math:`eye = 0.0001`.

    Here, :math:`n` represents the number of dimensions and
    :math:`x_i \in [-10, 10]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 1` for :math:`x_i = 0` for
    :math:`i = 1, ..., n`

    .. [1] Gavana, A. Global Optimization Benchmarks and AMPGO retrieved 2015
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-10.0] * self.N, [10.0] * self.N))

        self.global_optimum = [[0.0 for _ in range(self.N)]]
        self.fglob = 1.0
        self.change_dimensionality = True

    def fun(self, x, *args):
        self.nfev += 1

        f = fp = 0.0
        eye = 0.0001

        for val in x:
            if abs(val) >= eye:
                fp = 1.0
                f += 100.0 + abs(val)
            else:
                f += 1.0

        if fp < 1e-6:
            f = f / self.N

        return f


class NewFunction01(Benchmark):

    r"""
    NewFunction01 objective function.

    This class defines the NewFunction01 [1]_ global optimization problem. This is a
    multimodal minimization problem defined as follows:

    .. math::

       f_{\text{NewFunction01}}(x) = \left | {\cos\left(\sqrt{\left|{x_{1}^{2}
       + x_{2}}\right|}\right)} \right |^{0.5} + (x_{1} + x_{2})/100


    with :math:`x_i \in [-10, 10]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = -0.18459899925` for
    :math:`x = [-8.46669057, -9.99982177]`

    .. [1] Mishra, S. Global Optimization by Differential Evolution and
    Particle Swarm Methods: Evaluation on Some Benchmark Functions.
    Munich Personal RePEc Archive, 2006, 1005

    TODO line 355
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-10.0] * self.N, [10.0] * self.N))

        self.global_optimum = [[-8.46668984648, -9.99980944557]]
        self.fglob = -0.184648852475

    def fun(self, x, *args):
        self.nfev += 1

        return ((abs(cos(sqrt(abs(x[0] ** 2 + x[1]))))) ** 0.5
                + 0.01 * (x[0] + x[1]))


class NewFunction02(Benchmark):

    r"""
    NewFunction02 objective function.

    This class defines the NewFunction02 global optimization problem. This is a
    multimodal minimization problem defined as follows:

    .. math::

       f_{\text{NewFunction02}}(x) = \left | {\sin\left(\sqrt{\lvert{x_{1}^{2}
       + x_{2}}\rvert}\right)} \right |^{0.5} + (x_{1} + x_{2})/100


    with :math:`x_i \in [-10, 10]` for :math:`i = 1, 2`.

    *Global optimum*: :math:`f(x) = -0.19933159253` for
    :math:`x = [-9.94103375, -9.99771235]`

    .. [1] Mishra, S. Global Optimization by Differential Evolution and
    Particle Swarm Methods: Evaluation on Some Benchmark Functions.
    Munich Personal RePEc Archive, 2006, 1005

    TODO Line 368
    TODO WARNING, minimum value is estimated from running many optimisations and
    choosing the best.
    """

    def __init__(self, dimensions=2):
        Benchmark.__init__(self, dimensions)

        self._bounds = list(zip([-10.0] * self.N, [10.0] * self.N))

        self.global_optimum = [[-9.94114736324, -9.99997128772]]
        self.fglob = -0.199409030092

    def fun(self, x, *args):
        self.nfev += 1

        return ((abs(sin(sqrt(abs(x[0] ** 2 + x[1]))))) ** 0.5
                + 0.01 * (x[0] + x[1]))


#Newfunction 3 from Gavana is entered as Mishra05.
