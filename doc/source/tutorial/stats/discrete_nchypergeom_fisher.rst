
.. _discrete-nchypergeom-fisher:

Fisher's Noncentral Hypergeometric Distribution
===============================================

A random variable has Fisher's Noncentral Hypergeometric distribution with
parameters

:math:`M \in {\mathbb N}`,
:math:`n \in [0, M]`,
:math:`N \in [0, M]`,
:math:`\omega > 0`,

if its probability mass function is given by

.. math::

    p(x; M, n, N, \omega) = \frac{\binom{n}{x}\binom{M - n}{N-x}\omega^x}{P_0},

for
:math:`x \in [x_l, x_u]`,
where
:math:`x_l = \max(0, N - (M - n))`,
:math:`x_u = \min(N, n)`,

.. math::

    P_k = \sum_{y=x_l}^{x_u} \binom{n}{y} \binom{M - n}{N-y} \omega^y y^k,

and the binomial coefficients are

.. math::

    \binom{n}{k} \equiv \frac{n!}{k! (n - k)!}.

Other functions of this distribution are

.. math::
   :nowrap:

    \begin{eqnarray*}
    \mu & = & \frac{P_0}{P_1},\\
    \mu_{2} & = & \frac{P_2}{P_0} - \left(\frac{P_1}{P_0}\right)^2,\\
    \end{eqnarray*}

References
----------
-  Agner Fog, "Biased Urn Theory", https://cran.r-project.org/web/packages/BiasedUrn/vignettes/UrnTheory.pdf
-  "Fisher's noncentral hypergeometric distribution", Wikipedia, https://en.wikipedia.org/wiki/Fisher's_noncentral_hypergeometric_distribution

Implementation: `scipy.stats.nchypergeom_fisher`
