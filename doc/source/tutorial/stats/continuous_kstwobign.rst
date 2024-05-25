
.. _continuous-kstwobign:

KStwobign Distribution
======================

This is the limiting distribution of the normalized maximum absolute differences between an
empirical distribution function, computed from :math:`n` samples or observations,
and a comparison (or target) cumulative distribution function.  (``ksone`` is the distribution
of the unnormalized positive differences, :math:`D_n^+`.)

Writing :math:`D_n = \sup_t \left|F_{empirical,n}(t) - F_{target}(t)\right|`,
the normalization factor is :math:`\sqrt{n}`, and ``kstwobign`` is the limiting distribution
of the :math:`\sqrt{n} D_n` values as :math:`n\rightarrow\infty`.

Note that :math:`D_n=\max(D_n^+, D_n^-)`, but :math:`D_n^+` and :math:`D_n^-` are not independent.

``kstwobign`` can also be used with the differences between two empirical distribution functions,
for sets of observations with :math:`m` and :math:`n` samples respectively,
where :math:`m` and :math:`n` are "big".
Writing :math:`D_{m,n} = \sup_t \left|F_{1,m}(t)-F_{2,n}(t)\right|`,  where
:math:`F_{1,m}` and :math:`F_{2,n}` are the two empirical distribution functions, then
``kstwobign`` is also the limiting distribution of the :math:`\sqrt{\frac{mn}{m+n}}D_{m,n}` values,
as :math:`m,n\rightarrow\infty` and :math:`m/n\rightarrow a \ne 0, \infty`.

There are no shape parameters, and the support is :math:`x\in\left[0,\infty\right)`.


.. math::
    :nowrap:

    \begin{eqnarray*}  F\left(x\right) & = & 1 - 2 \sum_{k=1}^{\infty} (-1)^{k-1} e^{-2k^2 x^2}\\
    & = & \frac{\sqrt{2\pi}}{x} \sum_{k=1}^{\infty} e^{-(2k-1)^2 \pi^2/(8x^2)}\\
    & = & 1 - \textrm{scipy.special.kolmogorov}(n, x) \\
    f\left(x\right) & = & 8x \sum_{k=1}^{\infty} (-1)^{k-1} k^2 e^{-2k^2 x^2} \end{eqnarray*}


References
----------

-  "Kolmogorov-Smirnov test", Wikipedia
   https://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test

-  Kolmogoroff, A. "Confidence Limits for an Unknown Distribution Function.""
   *Ann. Math. Statist.* 12 (1941), no. 4, 461--463.

-  Smirnov, N. "On the estimation of the discrepancy between empirical curves of distribution for two independent samples"
   *Bull. Math. Univ. Moscou.*, 2 (1039), 2-26.

-  Feller, W. "On the Kolmogorov-Smirnov Limit Theorems for Empirical Distributions."
   *Ann. Math. Statist.* 19 (1948), no. 2, 177--189. and "Errata"  *Ann. Math. Statist.* 21 (1950), no. 2, 301--302.


Implementation: `scipy.stats.kstwobign`
