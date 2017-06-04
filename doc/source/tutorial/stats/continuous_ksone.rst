
.. _continuous-ksone:

KSone Distribution
==================


This is the distribution of maximum positive differences between an
empirical distribution function, computed from n samples or observations,
and a comparison (or target) cumulative distribution function.

Writing :math:`D_n^+ = \sup_t \left(F_{target}(t)-F_{empirical,n}(t)\right)`,
KSone is the distribution of the :math:`D_n^+` values.
(The distribution of :math:`D_n^- = -\inf_t \left(F_{target}(t)-F_{empirical,n}(t)\right)`
 differences follows the same distribution, so ksone can be used for one-sided tests on either side.)


One shape parameter

.. math::
   n, \textrm{a positive integer}


.. math::
   :nowrap:

    \begin{eqnarray*} F\left(n, x\right) & = & 1 - \sum_{j=0}^{\lfloor n(1-x)\rfloor} \dbinom{n}{j} x (x+\frac{j}{n})^{j-1} (1-x-\frac{j}{n})^{n-j}\\ & = & 1 - \textrm{scipy.special.smirnov}(n, x) \\ \lim_{n \rightarrow\infty} F\left(n, \frac{x}{\sqrt n}\right) & = & e^{-2 x^2} \end{eqnarray*}


:math:`x\in\left[0,1\right]`

References
----------

-  "Kolmogorov-Smirnov test", Wikipedia
   https://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test

-  Birnbaum, Z. W. and Fred H. Tingey (1951), One-sided confidence contours for probability distribution functions. *The Annals of Mathematical Statistics*, **22**/4, 592-596.



Implementation: `scipy.stats.ksone`
