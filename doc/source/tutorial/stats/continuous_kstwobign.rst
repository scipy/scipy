
.. _continuous-kstwobign:

KStwo Distribution
==================

This is the distribution of maximum absolute differences
between two empirical distribution functions for sets of
observations with m and n samples respectively, when m and n are "big".

Writing :math:`D_{m,n} = \sup_t \left|F_{1,m}(t)-F_{2,n}(t)\right|`,  where
:math:`F_{1,m}` and :math:`F_{2,n}` are the two empirical distribution functions, then
KStwo is the limiting distribution of the :math:`\sqrt{\left(\frac{mn}{m+n}\right)D_{m,n}}` values,
as :math:`m,n\rightarrow\infty`.

.. math::
   :nowrap:

    \begin{eqnarray*}  F\left(x\right) & = & 1 - 2 \sum_{k=1}^{\infty} (-1)^{k-1} e^{-2k^2 x^2}\\  & = & \frac{\sqrt{2\pi}}{x} \sum_{k=1}^{\infty} e^{-(2k-1)^2 \pi^2/(8x^2)}\\  & = & 1 - \textrm{scipy.special.kolmogorov}(n, x) \\ f\left(x\right) & = & 8x \sum_{k=1}^{\infty} (-1)^{k-1} k^2 e^{-2k^2 x^2} \end{eqnarray*}


:math:`x\in\left[0,1\right]`

References
----------

-  "Kolmogorov-Smirnov test", Wikipedia
   https://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test


Implementation: `scipy.stats.kstwobign`
