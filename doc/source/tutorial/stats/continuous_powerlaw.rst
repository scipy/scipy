
.. _continuous-powerlaw:

Power-function Distribution
===========================

A special case of the beta distribution with :math:`b=1` : defined for :math:`x\in\left[0,1\right]`

.. math::

     a>0

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;a\right) & = & ax^{a-1}\\ F\left(x;a\right) & = & x^{a}\\ G\left(q;a\right) & = & q^{1/a}\\ \mu & = & \frac{a}{a+1}\\ \mu_{2} & = & \frac{a\left(a+2\right)}{\left(a+1\right)^{2}}\\ \gamma_{1} & = & 2\left(1-a\right)\sqrt{\frac{a+2}{a\left(a+3\right)}}\\ \gamma_{2} & = & \frac{6\left(a^{3}-a^{2}-6a+2\right)}{a\left(a+3\right)\left(a+4\right)}\\ m_{d} & = & 1\end{eqnarray*}

.. math::

     h\left[X\right]=1-\frac{1}{a}-\log\left(a\right)

Implementation: `scipy.stats.powerlaw`
