
.. _discrete-betabinom:

Beta-Binomial Distribution
==========================

The beta-binomial distribution is a binomial distribution with a probability of success `p` that follows a beta distribution. The probability mass function for `betabinom`, defined for :math:`0 \leq k \leq n`, is:

.. math::

    f(k; n, a, b) = \binom{n}{k} \frac{B(k + a, n - k + b)}{B(a, b)}

for ``k`` in ``{0, 1,..., n}``, where :math:`B(a, b)` is the Beta function.

In the limiting case of :math:`a = b = 1`, the beta-binomial distribution reduces to a discrete uniform distribution:

.. math::

    f(k; n, 1, 1) = \frac{1}{n + 1}

In the limiting case of :math:`n = 1`, the beta-binomial distribution reduces to a Bernoulli distribution with the shape parameter :math:`p = a / (a + b)`:

.. math::

    f(k; 1, a, b) = \begin{cases}a / (a + b) & \text{if}\; k = 0 \\b / (a + b) & \text{if}\; k = 1\end{cases}

Implementation: `scipy.stats.betabinom`
