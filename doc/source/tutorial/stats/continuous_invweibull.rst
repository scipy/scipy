
.. _continuous-invweibull:

Inverted Weibull Distribution
=============================

Shape parameter :math:`c>0` and :math:`x>0` . Then

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & cx^{-c-1}\exp\left(-x^{-c}\right)\\ F\left(x;c\right) & = & \exp\left(-x^{-c}\right)\\ G\left(q;c\right) & = & \left(-\log q\right)^{-1/c}\end{eqnarray*}

.. math::

     h\left[X\right]=1+\gamma+\frac{\gamma}{c}-\log\left(c\right)

where :math:`\gamma` is Euler's constant.

Implementation: `scipy.stats.invweibull`
