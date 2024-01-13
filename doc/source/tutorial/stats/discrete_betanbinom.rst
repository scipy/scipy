
.. _discrete-betanbinom:

Beta-Negative Binomial Distribution
===================================

The beta-negative binomial distribution is a negative binomial distribution with a probability of success `p` that follows a beta distribution. The probability mass function for `betanbinom`, defined for :math:`k\geq 0`, is:

.. math::

    f(k; n, a, b) = \binom{n + k - 1}{k} \frac{B(a + n, b + k)}{B(a, b)}

for ``k`` in ``{0, 1,...}``, where :math:`B(a, b)` is the Beta function.

In the limiting case of :math:`n = 1`, the beta-negative binomial distribution reduces to a beta-geometric distribution with the probability mass function:

.. math::

    f(k; a, b) = \frac{B(a + 1, b + k)}{B(a, b)}

Implementation: `scipy.stats.betanbinom`
