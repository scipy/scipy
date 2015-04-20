
.. _continuous-gompertz:

Gompertz (Truncated Gumbel) Distribution
========================================

For :math:`x\geq0` and :math:`c>0` . In JKB the two shape parameters :math:`b,a` are reduced to the single shape-parameter :math:`c=b/a` . As :math:`a` is just a scale parameter when :math:`a\neq0` . If :math:`a=0,` the distribution reduces to the exponential distribution scaled by :math:`1/b.` Thus, the standard form is given as

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & ce^{x}\exp\left[-c\left(e^{x}-1\right)\right]\\ F\left(x;c\right) & = & 1-\exp\left[-c\left(e^{x}-1\right)\right]\\ G\left(q;c\right) & = & \log\left[1-\frac{1}{c}\log\left(1-q\right)\right]\end{eqnarray*}

.. math::

     h\left[X\right]=1-\log\left(c\right)-e^{c}\mathrm{Ei}\left(1,c\right),

where

.. math::

     \mathrm{Ei}\left(n,x\right)=\int_{1}^{\infty}t^{-n}\exp\left(-xt\right)dt

Implementation: `scipy.stats.gompertz`
