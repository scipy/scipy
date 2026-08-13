.. _continuous-skew-cauchy:

Skewed Cauchy Distribution
==========================

This distribution is a generalization of the Cauchy distribution. It
has a single shape parameter :math:`-1 < a < 1` that skews the distribution.
The special case :math:`a=0` yields the Cauchy distribution.

Functions
---------

.. math::
   :nowrap:

    \begin{eqnarray*}
    f(x, a) & = & \frac{1}{\pi \left(\frac{x^2}{\left(a x + 1 \right)^2} + 1 \right)},\quad x\ge0; \\
                 & = & \frac{1}{\pi \left(\frac{x^2}{\left(-a x + 1 \right)^2} + 1 \right)},\quad x<0. \\
    F(x, a) & = & \frac{1 - a}{2} + \frac{1 + a}{\pi} \arctan\left(\frac{x}{1 + a} \right),\quad x\ge0; \\
                 & = & \frac{1 - a}{2} + \frac{1 - a}{\pi} \arctan\left(\frac{x}{1 - a} \right),\quad x<0.
    \end{eqnarray*}

The mean, variance, skewness, and kurtosis are all undefined.

References
----------

-  "Skewed generalized *t* distribution", Wikipedia
   https://en.wikipedia.org/wiki/Skewed_generalized_t_distribution#Skewed_Cauchy_distribution

Implementation: `scipy.stats.skewcauchy`
