
.. _continuous-truncpareto:

Truncated Pareto Distribution
=============================

Two shape parameters, exponent :math:`b>0` and upper truncation :math:`c>1`. The support is :math:`1\leq x \leq c`.

The standard form is

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;b, c\right) & = & \frac{b}{1 - c^{-b}}\frac{1}{x^{b+1}}\\ F\left(x;b, c\right) & = & \frac{1 - x^{-b}}{1 - c^{-b}} \\ G\left(q;b, c\right) & = & \left(1 - \left(1 - c^{-b}\right)q\right)^{-1/b}.\end{eqnarray*}

The non-central moments are given by

.. math::

    \mu_{n}^{\prime} = \begin{cases}\displaystyle\frac{b\log c}{1 - c^{-b}} & \mbox{if $n = b$,}\\ \displaystyle\frac{b}{b-n}\frac{c^b - c^n}{c^b - 1} & \mbox{otherwise.}\end{cases}

The entropy is

.. math::

     h\left(X\right)= -\left[\log\left(\frac{b}{1 - c^{-b}}\right) + (b+1)\left(\frac{\log c}{c^b - 1} - \frac{1}{b}\right)\right].


Implementation: `scipy.stats.truncpareto`
