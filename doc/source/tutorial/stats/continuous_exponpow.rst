
.. _continuous-exponpow:

Exponential Power Distribution
==============================

One positive shape parameter :math:`b`. The support is :math:`x\geq0.`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;b\right) & = & ebx^{b-1}\exp\left(x^{b}-e^{x^{b}}\right)\\
    F\left(x;b\right) & = & 1-\exp\left(1-e^{x^{b}}\right)\\
    G\left(q;b\right) & = & \log\left(1-\log\left(1-q\right)\right)^{1/b}\end{eqnarray*}

Implementation: `scipy.stats.exponpow`
