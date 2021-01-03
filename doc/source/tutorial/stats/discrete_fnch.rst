
.. _discrete-fnch:

Fisher's Noncentral Hypergeometric Distribution
===============================================

A random variable has Fisher's Noncentral Hypergeometric distribution with
parameters

:math:`N \in {\mathbb N}`,
:math:`m_1 \in [0, N]`,
:math:`n \in [0, N]`,
:math:`\omega > 0`,

if its probability mass function is given by

.. math::

    p(x; N, m_1, n, \omega) = \frac{\binom{m_1}{x}\binom{m_2}{n-x}\omega^x}{P_0},

for
:math:`x \in [x_l, x_u]`,
where
:math:`m_2 = N - m_1`,
:math:`x_l = \max(0, n - m_2)`,
:math:`x_u = \min(n, m_1)`,

.. math::

    P_k = \sum_{y=x_l}^{x_u} \binom{m_1}{y} \binom{m_2}{n-y} \omega^y y^k,

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

Implementation: `scipy.stats.fnch`
