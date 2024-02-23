
.. _continuous-lomax:

Pareto Second Kind (Lomax) Distribution
=======================================

This is Pareto of the first kind with :math:`L=-1.0` .  There is one shape parameter
:math:`c>0` and support :math:`x\geq0`.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \frac{c}{\left(1+x\right)^{c+1}}\\ F\left(x;c\right) & = & 1-\frac{1}{\left(1+x\right)^{c}}\\ G\left(q;c\right) & = & \left(1-q\right)^{-1/c}-1\end{eqnarray*}

.. math::

     h\left[X\right]=\frac{1}{c}+1-\log\left(c\right).

Implementation: `scipy.stats.lomax`
