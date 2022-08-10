
.. _continuous-truncweibull_min:

Truncated Weibull Minimum Extreme Value Distribution
====================================================

A doubly truncated version of Weibull minimum extreme value distribution.
Defined for :math:`a<x<=b` and :math:`c>0`.

.. math::
   :nowrap:

    \begin{eqnarray*}
        f\left(x;c,a,b\right) & = & \frac{cx^{c-1}\exp\left(-x^{c}\right)}{\exp\left(-a^{c}\right) - \exp\left(-b^{c}\right)} \\
        F\left(x;c,a,b\right) & = & \frac{\exp\left(-a^{c}\right) - \exp\left(-x^{c}\right)}{\exp\left(-a^{c}\right) - \exp\left(-b^{c}\right)} \\
        G\left(q;c,a,b\right) & = & \left[-\log\left(\left(1-q\right)\exp\left(-a^{c}\right)+q\exp\left(-b^{c}\right)\right)\right]^{1/c}
    \end{eqnarray*}

.. math::

     \mu_{n}^{\prime}=\frac{\exp\left(a^{c}\right)}{1-\exp\left(-b^{c}\right)}\left[\gamma\left(\frac{n}{c}+1,b^{c}\right)-\gamma\left(\frac{n}{c}+1,a^{c}\right)\right]

where :math:`\gamma\left(\right)` is the lower incomplete gamma function.

Implementation: `scipy.stats.truncweibull_min`
