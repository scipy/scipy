
.. _discrete-zipfian:

Zipfian Distribution
========================

A random variable has the Zipfian distribution with parameters
:math:`s \ge 0` and :math:`N \in \{1, 2, 3, \dots\}` if its probability
mass function is given by

.. math::
   :nowrap:

    \begin{eqnarray*} p\left(k; s, N \right) & = & \frac{1}{H_{N, s}k^{s}}\quad k \in \{1, 2, \dots, n-1, n\} \end{eqnarray*}

where

.. math::

    H_{N, s}=\sum_{n=1}^{N}\frac{1}{n^{s}}

is the :math:`N`:sup:`th` generalized harmonic number of order
:math:`s`. Other functions of this distribution are

.. math::
   :nowrap:

    \begin{eqnarray*}
     F\left(x; s, N\right) & = & \frac{H_{k, s}}{H_{N, s}}, \\
    \mu & = & \frac{H_{N, s-1}}{H_{N, s}},\\
    \mu_{2} & = & \frac{H_{N, s-2}}{H_{N, s}} - \frac{H^2_{N, s-1}}{H^2_{N, s}},\\
    \gamma_1 & = & \frac{\frac{H_{N, s-3}}{H_{N, s}} - 3 \frac{H_{N, s-1}H_{N, s-2}}{H_{N, s}^2} + 2\frac{H_{N, s-1}^3}{H_{N, s}^3}}{\left(\frac{H_{N, s-2}H_{N, s}- H_{N, s-1}^2}{H_{N, s}^2}\right)^{\frac{3}{2}}}, \mbox{and}\\
    \gamma_2 & = & \frac{H_{N, s}^3 H_{N, s-4} - 4 H_{N, s}^2 H_{N, s-1} H_{N, s-3} + 6 H_{N, s} H_{N, s-1}^2 H_{N, s-2} - 3 H_{N, s-1}^4}{\left(H_{N, s-2} H_{N, s} - H_{N, s-1}^2 \right)^2}.
    \end{eqnarray*}

References
----------
-  "Zipf's Law", Wikipedia, https://en.wikipedia.org/wiki/Zipf%27s_law
-  Larry Leemis, "Zipf Distribution", Univariate Distribution Relationships. http://www.math.wm.edu/~leemis/chart/UDR/PDFs/Zipf.pdf

Implementation: `scipy.stats.zipfian`
