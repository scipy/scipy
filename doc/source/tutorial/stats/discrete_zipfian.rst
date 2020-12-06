
.. _discrete-zipfian:

Zipfian Distribution
========================

A random variable has the Zipfian distribution with parameters
:math:`s \ge 0` and :math:`N \in {1, 2, 3, \dots}` if its probability
mass function is given by

.. math::
   :nowrap:

    \begin{eqnarray*} p\left(k; s, N \right) & = & \frac{1}{H_{N, s}k^{s}}\quad k\geq1\end{eqnarray*}

where

.. math::

    H_{N, s}=\sum_{n=1}^{N}\frac{1}{n^{s}}=

is the :math:`N^\mbox{th}` generalized harmonic number of order :math:`s`.
Other functions of this distribution are

.. math::
   :nowrap:

    \begin{eqnarray*}
     F\left(x; s, N\right) & = & \frac{H_{k, s}}{H_{N, s}}, \\
    \mu & = & \frac{H_{N, s-1}}{H_{N, s}}, \mbox{and}\\
    \mu_{2} & = & \frac{H_{N, s-2}}{H_{N, s}} - \frac{H^2_{N, s-1}}{H^2_{N, s}}.\\
    \end{eqnarray*}

References
----------
-  "Zipf's Law", Wikipedia, https://en.wikipedia.org/wiki/Zipf%27s_law

Implementation: `scipy.stats.zipfian`
