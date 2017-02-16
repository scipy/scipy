
.. _continuous-frechet_l:

Fréchet (left-skewed, Extreme Value Type III, Weibull maximum) Distribution
============================================================================

Defined for :math:`x<0` and :math:`c>0` .

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & c\left(-x\right)^{c-1}\exp\left(-\left(-x\right)^{c}\right)\\ F\left(x;c\right) & = & \exp\left(-\left(-x\right)^{c}\right)\\ G\left(q;c\right) & = & -\left(-\log q\right)^{1/c}\end{eqnarray*}

The mean is the negative of the right-skewed Frechet distribution
given above, and the other statistical parameters can be computed from

.. math::

     \mu_{n}^{\prime}=\left(-1\right)^{n}\Gamma\left(1+\frac{n}{c}\right).

.. math::

     h\left[X\right]=-\frac{\gamma}{c}-\log\left(c\right)+\gamma+1

where :math:`\gamma` is Euler's constant and equal to

.. math::

     \gamma\approx0.57721566490153286061.

Implementation: `scipy.stats.frechet_l`
