
.. _continuous-ksone:

KSone Distribution
==================


This is the distribution of maximum positive "gaps", when comparing a theoretical distribution with an empirical one.
Writing :math:`D_n^+ = \sup_t [F_{theoretical}(t)-S_n(t)]`,  KSOne is the distribution of the :math:`D_n^+` values.

One shape parameter

.. math::
   n, \textrm{a positive integer}

.. math::
   :nowrap:

    \begin{eqnarray*} F\left(n, x\right) & = & 1 - \sum_{j=0}^{\lfloor n(1-x)\rfloor} \dbinom{n}{v} x (x+\frac{j}{n})^{j-1} (1-x-\frac{j}{n})^{n-j}\\ & = & 1 - \textrm{scipy.special.smirnov}(n, x) \\\lim_{n \rightarrow\infty} F\left(n, \frac{x}{\sqrt n}\right) & = & e^{-2 x^2} \end{eqnarray*}


References
----------

-  "Kolmogorov-Smirnov test", Wikipedia
   https://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test

-  Birnbaum, Z. W. and Fred H. Tingey (1951), One-sided confidence contours for probability distribution functions. *The Annals of Mathematical Statistics*, **22**/4, 592-596.



Implementation: `scipy.stats.ksone`
