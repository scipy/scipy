# -*- coding: utf-8 -*-
from numpy import cos, exp, log, pi, sin, sqrt

from .go_benchmark import Benchmark, safe_import

with safe_import():
    try:
        from scipy.special import factorial  # new
    except ImportError:
        from scipy.misc import factorial  # old


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

    This class defines the Univariate Problem15 global optimization problem.
    This is a multimodal minimization problem defined as follows:

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

    This class defines the Univariate Problem18 global optimization problem.
    This is a multimodal minimization problem defined as follows:

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

    This class defines the Univariate Problem20 global optimization problem.
    This is a multimodal minimization problem defined as follows:

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

    This class defines the Univariate Problem21 global optimization problem.
    This is a multimodal minimization problem defined as follows:

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

    This class defines the Univariate Problem22 global optimization problem.
    This is a multimodal minimization problem defined as follows:

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
