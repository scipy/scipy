
.. _continuous-exponweib:

Exponentiated Weibull Distribution
==================================

Two positive shape parameters :math:`a,c>0`, and the support is :math:`x\in\left[0,\infty\right)`.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;a,c\right) & = & ac\left[1-\exp\left(-x^{c}\right)\right]^{a-1}\exp\left(-x^{c}\right)x^{c-1}\\ F\left(x;a,c\right) & = & \left[1-\exp\left(-x^{c}\right)\right]^{a}\\ G\left(q;a,c\right) & = & \left[-\log\left(1-q^{1/a}\right)\right]^{1/c}\end{eqnarray*}

Implementation: `scipy.stats.exponweib`
