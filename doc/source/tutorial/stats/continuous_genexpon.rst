
.. _continuous-genexpon:

Generalized Exponential Distribution
====================================

Three positive shape parameters :math:`a,b,c>0` with support :math:`x\geq0`.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;a,b,c\right) & = & \left(a+b\left(1-e^{-cx}\right)\right)\exp\left(ax-bx+\frac{b}{c}\left(1-e^{-cx}\right)\right)\\
    F\left(x;a,b,c\right) & = & 1-\exp\left(ax-bx+\frac{b}{c}\left(1-e^{-cx}\right)\right)\\
    G\left(q;a,b,c\right) & = & F^{-1}\end{eqnarray*}

Implementation: `scipy.stats.genexpon`
