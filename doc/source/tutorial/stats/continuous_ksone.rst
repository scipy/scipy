
.. _continuous-ksone:

KSone Distribution
==================


This is the distribution of maximum positive differences between an
empirical distribution function, computed from :math:`n` samples or observations,
and a comparison (or target) cumulative distribution function.

Writing :math:`D_n^+ = \sup_t \left(F_{empirical,n}(t)-F_{target}(t)\right)`,
``ksone`` is the distribution of the :math:`D_n^+` values.
(The distribution of :math:`D_n^- = \sup_t \left(F_{target}(t)-F_{empirical,n}(t)\right)`
differences follows the same distribution, so ``ksone`` can be used for one-sided tests on either side.)


There is one shape parameter :math:`n`, a positive integer, and the support is :math:`x\in\left[0,1\right]`.

.. math::
   :nowrap:

    \begin{eqnarray*} F\left(n, x\right) & = & 1 - \sum_{j=0}^{\lfloor n(1-x)\rfloor} \dbinom{n}{j} x \left(x+\frac{j}{n}\right)^{j-1} \left(1-x-\frac{j}{n}\right)^{n-j}\\
    & = & 1 - \textrm{scipy.special.smirnov}(n, x) \\
    \lim_{n \rightarrow\infty} F\left(n, \frac{x}{\sqrt n}\right) & = & e^{-2 x^2} \end{eqnarray*}


References
----------

-  "Kolmogorov-Smirnov test", Wikipedia
   https://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test

-  Birnbaum, Z. W.; Tingey, Fred H. "One-Sided Confidence Contours for Probability Distribution Functions."
   *Ann. Math. Statist*. 22 (1951), no. 4, 592--596.


Implementation: `scipy.stats.ksone`
