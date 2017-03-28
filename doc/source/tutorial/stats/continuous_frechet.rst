
.. _continuous-frechet:

Fréchet (Type II Extreme Value) Distribution
============================================

A type of extreme-value distribution.  Defined for :math:`x > 0` and :math:`\alpha > 0`

.. math::
   :nowrap:

    \begin{eqnarray*}
        f\left(x; \alpha\right) & = & \alpha x^{-\alpha - 1}\exp\left(-x^{-\alpha}\right) \\
        F\left(x; \alpha\right) & = & \exp\left(-x^{-\alpha}\right) \\
        G\left(q; \alpha\right) & = & \left[\frac{-1}{\log\left(q\right)}\right]^{1/\alpha}
    \end{eqnarray*}

.. math::

     \mu_{n}^{\prime} = \begin{cases}
                            \Gamma\left(1 - \frac{n}{\alpha}\right) & \text{if}\; \alpha > n \\
                            \infty & \text{if}\; 0 < \alpha \le n
                        \end{cases}

.. math::
   :nowrap:

    \begin{eqnarray*}
        \mu & = & \begin{cases}
                      \Gamma\left(1-\frac{1}{\alpha}\right) & \text{if}\; \alpha > 1 \\
                      \infty & \text{if}\; 0 < \alpha \le 1
                  \end{cases} \\
        \mu_{2} & = & \begin{cases}
                          \Gamma\left(1-\frac{2}{\alpha}\right) -
                              \Gamma^{2}\left(1-\frac{1}{\alpha}\right) & \text{if}\; \alpha > 2 \\
                          \infty & \text{if}\; 0 < \alpha \le 2
                      \end{cases} \\
        \gamma_{1} & = & \begin{cases}
                             \frac{\Gamma\left(1-\frac{3}{\alpha}\right) -
                                   3\Gamma\left(1-\frac{2}{\alpha}\right)\Gamma\left(1-\frac{1}{\alpha}\right) +
                                   2\Gamma^{3}\left(1-\frac{1}{\alpha}\right)}
                                  {\mu_{2}^{3/2}} & \text{if}\; \alpha > 3 \\
                             \infty & \text{if}\; 0 < \alpha \le 3
                         \end{cases} \\
        \gamma_{2} & = & \begin{cases}
                             \frac{\Gamma\left(1-\frac{4}{\alpha}\right) -
                                   4\Gamma\left(1-\frac{1}{\alpha}\right)\Gamma\left(1-\frac{3}{\alpha}\right) +
                                   6\Gamma^{2}\left(1-\frac{1}{\alpha}\right)\Gamma\left(1-\frac{2}{\alpha}\right) -
                                   3\Gamma^{4}\left(1-\frac{1}{\alpha}\right)}
                                  {\mu_{2}^{2}} - 3 & \text{if}\; \alpha > 4 \\
                             \infty & \text{if}\; 0 < \alpha \le 4
                         \end{cases} \\
        m_{d} & = & \left(\frac{\alpha}{1+\alpha}\right)^{\frac{1}{\alpha}} \\
        m_{n} & = & \left(\frac{1}{\ln\left(2\right)}\right)^{\frac{1}{\alpha}}
    \end{eqnarray*}

.. math::

     h\left[X\right] = \frac{\gamma}{\alpha}-\log\left(\alpha\right)+\gamma+1

where :math:`\gamma` is Euler's constant and equal to

.. math::

     \gamma\approx0.57721566490153286061.

Implementation: `scipy.stats.frechet`
