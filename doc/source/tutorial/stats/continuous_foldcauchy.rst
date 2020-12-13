
.. _continuous-foldcauchy:

Folded Cauchy Distribution
==========================

This formula can be expressed in terms of the standard formulas for
the Cauchy distribution (call the cdf :math:`C\left(x\right)` and
the pdf :math:`d\left(x\right)` ).
If :math:`Y` is cauchy then :math:`\left|Y\right|` is folded cauchy.
There is one shape parameter :math:`c` and the support is :math:`x\geq0.`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \frac{1}{\pi\left(1+\left(x-c\right)^{2}\right)}+\frac{1}{\pi\left(1+\left(x+c\right)^{2}\right)}\\
    F\left(x;c\right) & = & \frac{1}{\pi}\tan^{-1}\left(x-c\right)+\frac{1}{\pi}\tan^{-1}\left(x+c\right)\\
    G\left(q;c\right) & = & F^{-1}\left(q;c\right)\end{eqnarray*}

No moments

Implementation: `scipy.stats.foldcauchy`
