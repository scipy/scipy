
.. _discrete-wnch:

Wallenius' Noncentral Hypergeometric Distribution
===============================================

A random variable has Wallenius' Noncentral Hypergeometric distribution with
parameters

:math:`N \in {\mathbb N}`,
:math:`m_1 \in [0, N]`,
:math:`n \in [0, N]`,
:math:`\omega > 0`,

if its probability mass function is given by

.. math::

    p(x; n, m_1, N) = \binom{m_1}{x} \binom{m_2}{n-x}\int_0^1 \left(1-t^{\omega/D}\right)^x\left(1-t^{1/D}\right)^{n-x} dt

for
:math:`x \in [x_l, x_u]`,
where
:math:`m_2 = N - m_1`,
:math:`x_l = \max(0, n - m_2)`,
:math:`x_u = \min(n, m_1)`,

.. math::

    D = \omega(m_1 - x) + (m_2-(n-x)),

and the binomial coefficients are

.. math::

    \binom{n}{k} \equiv \frac{n!}{k! (n - k)!}.

References
----------
-  Agner Fog, "Biased Urn Theory", https://cran.r-project.org/web/packages/BiasedUrn/vignettes/UrnTheory.pdf
-  "Wallenius' noncentral hypergeometric distribution", Wikipedia, https://en.wikipedia.org/wiki/Wallenius'_noncentral_hypergeometric_distribution

Implementation: `scipy.stats.wnch`
