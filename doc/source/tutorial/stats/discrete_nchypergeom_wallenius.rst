
.. _discrete-nchypergeom-wallenius:

Wallenius' Noncentral Hypergeometric Distribution
=================================================

A random variable has Wallenius' Noncentral Hypergeometric distribution with
parameters

:math:`M \in {\mathbb N}`,
:math:`n \in [0, M]`,
:math:`N \in [0, M]`,
:math:`\omega > 0`,

if its probability mass function is given by

.. math::

    p(x; N, n, M) = \binom{n}{x} \binom{M - n}{N-x}\int_0^1 \left(1-t^{\omega/D}\right)^x\left(1-t^{1/D}\right)^{N-x} dt

for
:math:`x \in [x_l, x_u]`,
where
:math:`x_l = \max(0, N - (M - n))`,
:math:`x_u = \min(N, n)`,

.. math::

    D = \omega(n - x) + ((M - n)-(N-x)),

and the binomial coefficients are

.. math::

    \binom{n}{k} \equiv \frac{n!}{k! (n - k)!}.

References
----------
-  Agner Fog, "Biased Urn Theory", https://cran.r-project.org/web/packages/BiasedUrn/vignettes/UrnTheory.pdf
-  "Wallenius' noncentral hypergeometric distribution", Wikipedia, https://en.wikipedia.org/wiki/Wallenius'_noncentral_hypergeometric_distribution

Implementation: `scipy.stats.nchypergeom_wallenius`
