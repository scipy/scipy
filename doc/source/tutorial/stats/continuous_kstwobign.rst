
.. _continuous-kstwobign:

KStwo Distribution
==================

This is the distribution of maximum absolute gaps, when comparing
two observations with m and n samples respectively, when m and n are "big".

Writing :math:`D_{m,n} = \sup_t (F_{1,m}(t)-F_{2,n}(t))`,  where
:math:`F_{1,m}` and :math:`F_{2,n}` are the empirical distribution functions, then
KStwo is the distribution of the :math:`\sqrt\frac{mn}{m+n}D_{m,n}` values.

.. math::
   :nowrap:

    \begin{eqnarray*}  F\left(x\right) & = & 1 - 2 \sum_{k=1}^{\infty} (-1)^{k-1} e^{-2k^2 x^2}\\  & = & \frac{\sqrt{2\pi}}{x} \sum_{k=1}^{\infty} e^{-2(k-1)^2 \pi^2/(8x^2)}\\  & = & 1 - \textrm{scipy.special.kolmogorov}(n, x) \\ f\left(x\right) & = & 8x \sum_{k=1}^{\infty} (-1)^{k-1} k^2 e^{-2k^2 x^2} \end{eqnarray*}


References
----------

-  "Kolmogorov-Smirnov test", Wikipedia
   https://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test


Implementation: `scipy.stats.kstwobign`
